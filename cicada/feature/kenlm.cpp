//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/kenlm.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"
#include "cicada/cluster.hpp"

#include "kenlm/lm/model.hh"
#include "kenlm/lm/left.hh"

#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/atomicop.hpp"
#include "utils/array_power2.hpp"
#include "utils/unordered_map.hpp"

//
// ngram language model feature via kenlm
//

//
// Here, we use KenLMNGram to wrap kenlm ngram language model, especially to manage
// dynamically accessing vocabulary mapping.

namespace cicada
{
  namespace feature
  {
    template <typename Model>
    struct KenLMNGram
    {
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol symbol_type;
      typedef cicada::Symbol word_type;

      typedef boost::filesystem::path path_type;
      
      typedef Model model_type;
      
      struct Cache
      {
	typedef int64_t value_type;
	
	Cache() : value(value_type(-1)) {}
	
	volatile value_type value;
      };
      typedef Cache cache_type;
      typedef utils::array_power2<cache_type, 1024 * 32, std::allocator<cache_type> > cache_set_type;
      
      lm::WordIndex vocabulary(const word_type& word) const
      {
	cache_set_type& caches = const_cast<cache_set_type&>(caches_);
	
	cache_type cache;
	cache_type cache_new;
	
	cache.value = utils::atomicop::fetch_and_add(caches[word.id() & (caches.size() - 1)].value, int64_t(0));
	
	uint32_t __word = (cache.value >> 32) & 0xffffffff;
	uint32_t __id   = (cache.value) & 0xffffffff;
	
	if (__word == word.id())
	  return lm::WordIndex(__id);

        __id = model_.GetVocabulary().Index(static_cast<const std::string&>(word));
	
	cache_new.value = (uint64_t(word.id()) << 32) | (uint64_t(__id) & 0xffffffff);
	
	utils::atomicop::compare_and_swap(caches[word.id() & (caches.size() - 1)].value, cache.value, cache_new.value);
	
	return lm::WordIndex(__id);
      }
      
      KenLMNGram(const path_type& file,
		 const lm::ngram::Config& config) : caches_(), model_(file.string().c_str(), config) {}

      static
      KenLMNGram<Model>& create(const path_type& path, const bool populate)
      {
	typedef boost::mutex            mutex_type;
	typedef mutex_type::scoped_lock lock_type;
	
	typedef boost::shared_ptr<KenLMNGram<Model> > ngram_ptr_type;
	typedef typename utils::unordered_map<std::string, ngram_ptr_type, boost::hash<utils::piece>, std::equal_to<std::string>,
					      std::allocator<std::pair<const std::string, ngram_ptr_type> > >::type ngram_map_type;
	
	static mutex_type     ngram_mutex;
	static ngram_map_type ngram_map;
	
	lock_type lock(ngram_mutex);
	
	typename ngram_map_type::iterator iter = ngram_map.find(path.string());
	if (iter == ngram_map.end()) {
	  lm::ngram::Config config;
	  
	  config.load_method = populate ? util::POPULATE_OR_READ : util::LAZY;
	  
	  iter = ngram_map.insert(std::make_pair(path.string(), ngram_ptr_type(new KenLMNGram<Model>(path, config)))).first;
	}
	
	return *(iter->second);
      }

      cache_set_type caches_;
      model_type     model_;
    };


    template <typename Model>
    class KenLMImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::Cluster cluster_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;

      typedef boost::filesystem::path path_type;

      typedef KenLMNGram<Model> ngram_type;

      typedef lm::ngram::ChartState ngram_state_type;
      
    public:
      KenLMImpl(const path_type& __path, const bool populate)
	: ngram(&ngram_type::create(__path, populate)),
	  cluster(0), no_bos_eos(false), skip_sgml_tag(false),
	  log10(std::log(10))
      {
	id_oov = 0;
	id_bos = ngram->vocabulary(vocab_type::BOS);
	id_eos = ngram->vocabulary(vocab_type::EOS);
      }
      
      KenLMImpl(const KenLMImpl& x)
	: ngram(x.ngram),
	  cluster(x.cluster ? &cluster_type::create(x.cluster->path()) : 0),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  feature_name(x.feature_name),
	  feature_name_oov(x.feature_name_oov),
	  id_oov(x.id_oov),
	  id_bos(x.id_bos),
	  id_eos(x.id_eos),
	  log10(std::log(10))
      { }
      
      KenLMImpl& operator=(const KenLMImpl& x)
      {
	ngram = x.ngram;
	cluster = (x.cluster ? &cluster_type::create(x.cluster->path()) : 0);
	no_bos_eos = x.no_bos_eos;
	skip_sgml_tag = x.skip_sgml_tag;
	
	feature_name     = x.feature_name;
	feature_name_oov = x.feature_name_oov;
	id_oov           = x.id_oov;
	id_bos           = x.id_bos;
	id_eos           = x.id_eos;
	
	return *this;
      }
            
      
      struct extract_cluster
      {
	extract_cluster(const cluster_type* __cluster): cluster(__cluster) {}
	
	const cluster_type* cluster;

	symbol_type operator()(const symbol_type& word) const
	{
	  return cluster->operator[](word);
	}
      };

      struct extract_word
      {
	const symbol_type& operator()(const symbol_type& word) const
	{
	  return word;
	}
      };

      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON;
	}
      };
      
      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || (word != vocab_type::BOS && word != vocab_type::EOS && word.is_sgml_tag());
	}
      };
      
      double ngram_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 int& oov)
      {
	if (cluster) {
	  if (skip_sgml_tag)
	    return ngram_score(state, states, edge, oov, extract_cluster(cluster), skipper_sgml());
	  else
	    return ngram_score(state, states, edge, oov, extract_cluster(cluster), skipper_epsilon());
	} else {
	  if (skip_sgml_tag)
	    return ngram_score(state, states, edge, oov, extract_word(), skipper_sgml());
	  else
	    return ngram_score(state, states, edge, oov, extract_word(), skipper_epsilon());
	}
      }
      
      template <typename Extract, typename Skipper>
      double ngram_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 int& oov,
			 Extract extract,
			 Skipper skipper)
      {
	const rule_type& rule = *(edge.rule);
        const phrase_type& target = rule.rhs;
	
	lm::ngram::RuleScore<Model> ruleScore(ngram->model_, *reinterpret_cast<ngram_state_type*>(state));

	phrase_type::const_iterator titer_begin = target.begin();
        phrase_type::const_iterator titer_end   = target.end();
	
	int non_terminal_pos = 0;
	bool initial = true;
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    const int __non_terminal_index = titer->non_terminal_index();
            const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
            ++ non_terminal_pos;
	    
	    if (initial) // special case for left-edge NT
	      ruleScore.BeginNonTerminal(*reinterpret_cast<const ngram_state_type*>(states[antecedent_index]));
	    else
	      ruleScore.NonTerminal(*reinterpret_cast<const ngram_state_type*>(states[antecedent_index]));
	    
	    initial = false;
	  } else if (! skipper(*titer)) {
	    
	    if (initial && *titer == vocab_type::BOS)
	      ruleScore.BeginSentence();
	    else {
	      const lm::WordIndex id = ngram->vocabulary(extract(*titer));
	      
	      oov += (id == id_oov);
	      
	      ruleScore.Terminal(id);
	    }
	    
	    initial = false;
	  }
	}
	
	const double ret = ruleScore.Finish();
	
	reinterpret_cast<ngram_state_type*>(state)->ZeroRemaining();
	
	return ret * log10;
      }
      
      double ngram_final_score(const state_ptr_type& state)
      {
	if (no_bos_eos) {
	  // we do not need scoring by bos and eos...
	  return 0.0;
	} else {
	  const ngram_state_type& ngram_state = *reinterpret_cast<const ngram_state_type*>(state);
	  
	  ngram_state_type ngram_state_next;
	  lm::ngram::RuleScore<Model> ruleScore(ngram->model_, ngram_state_next);
	  
	  ruleScore.BeginSentence();
	  ruleScore.NonTerminal(ngram_state);
	  ruleScore.Terminal(id_eos);
	  
	  return ruleScore.Finish() * log10;
	}
      }

      size_type reserve_state_size() const
      {
	return sizeof(ngram_state_type);
      }      
      
      ngram_type*     ngram;
      
      // cluster...
      cluster_type* cluster;
      
      bool no_bos_eos;
      bool skip_sgml_tag;
      
      // names...
      feature_type feature_name;
      feature_type feature_name_oov;
      
      lm::WordIndex id_oov;
      lm::WordIndex id_bos;
      lm::WordIndex id_eos;
      
      double log10;
    };
    
    template <typename Model>
    KenLM<Model>::KenLM(const std::string& parameter)
      : pimpl(0), pimpl_coarse(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "kenlm")
	throw std::runtime_error("is this really kenlm feature function? " + parameter);

      path_type   path;
      bool        populate = false;
      path_type   cluster_path;
      bool        skip_sgml_tag = false;
      bool        no_bos_eos = false;
      
      path_type   coarse_path;
      bool        coarse_populate = false;
      path_type   coarse_cluster_path;
      
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "populate")
	  populate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "cluster")
	  cluster_path = piter->second;
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "coarse-file")
	  coarse_path = piter->second;
	else if (utils::ipiece(piter->first) == "coarse-populate")
	  coarse_populate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "coarse-cluster")
	  coarse_cluster_path = piter->second;
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (path.empty())
	throw std::runtime_error("no ngram file? " + path.string());
      
      if (! coarse_path.empty() && ! boost::filesystem::exists(coarse_path))
	throw std::runtime_error("no coarse ngram language model? " + coarse_path.string());
      
      std::auto_ptr<impl_type> ngram_impl(new impl_type(path, populate));

      ngram_impl->no_bos_eos = no_bos_eos;
      ngram_impl->skip_sgml_tag = skip_sgml_tag;
      
      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.string());
	
	ngram_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      // two contexts (order - 1) for each edge, with two separator..
      // state_size..
      base_type::__state_size   = ngram_impl->reserve_state_size();
      base_type::__feature_name = (name.empty() ? std::string("kenlm") : name);
      
      ngram_impl->feature_name     = base_type::__feature_name;
      ngram_impl->feature_name_oov = static_cast<const std::string&>(base_type::__feature_name) + ":oov-penalty";
      
      pimpl = ngram_impl.release();

      // coarse ngram
      if (! coarse_path.empty()) {
	std::auto_ptr<impl_type> ngram_impl(new impl_type(coarse_path, coarse_populate));
	
	ngram_impl->no_bos_eos = no_bos_eos;
	ngram_impl->skip_sgml_tag = skip_sgml_tag;
	
	if (! coarse_cluster_path.empty()) {
	  if (! boost::filesystem::exists(coarse_cluster_path))
	    throw std::runtime_error("no cluster file: " + coarse_cluster_path.string());
	  
	  ngram_impl->cluster = &cicada::Cluster::create(coarse_cluster_path);
	}
	
	pimpl_coarse = ngram_impl.release();
      }
    }
    
    template <typename Model>
    KenLM<Model>::~KenLM()
    {
      std::auto_ptr<impl_type> tmp(pimpl);
      if (pimpl_coarse)
	std::auto_ptr<impl_type> tmp_coarse(pimpl_coarse);
    }
    
    template <typename Model>
    KenLM<Model>::KenLM(const KenLM& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl)),
	pimpl_coarse(x.pimpl_coarse ? new impl_type(*x.pimpl_coarse) : 0)
    {}
    
    template <typename Model>
    KenLM<Model>& KenLM<Model>::operator=(const KenLM<Model>& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      if (x.pimpl_coarse) {
	if (pimpl_coarse)
	  *pimpl_coarse = *x.pimpl_coarse;
	else
	  pimpl_coarse = new impl_type(*x.pimpl_coarse);
      } else {
	if (pimpl_coarse)
	  delete pimpl_coarse;
	pimpl_coarse = 0;
      }
      
      return *this;
    }
    
    
    template <typename Model>
    void KenLM<Model>::apply(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     const bool final) const
    {
      int oov = 0;
      double score = pimpl->ngram_score(state, states, edge, oov);
      
      if (final)
	score += pimpl->ngram_final_score(state);
      
      if (score != 0.0)
	features[pimpl->feature_name] = score;
      else
	features.erase(pimpl->feature_name);

      if (oov)
	features[pimpl->feature_name_oov] = - oov;
      else
	features.erase(pimpl->feature_name_oov);
    }

    template <typename Model>
    void KenLM<Model>::apply_coarse(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    const bool final) const
    {
      // use of coarse, but it is mandatory...
      if (pimpl_coarse) {
	int oov = 0;
	double score = pimpl_coarse->ngram_score(state, states, edge, oov);
	
	if (final)
	  score += pimpl_coarse->ngram_final_score(state);
	
      if (score != 0.0)
	features[pimpl->feature_name] = score;
      else
	features.erase(pimpl->feature_name);

      if (oov)
	features[pimpl->feature_name_oov] = - oov;
      else
	features.erase(pimpl->feature_name_oov);	
      } else
	apply(state, states, edge, features, final);
    }

    // temporarily assigned feature function...
    
    template <typename Model>
    void KenLM<Model>::apply_predict(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     const bool final) const
    {
      
    }
    
    template <typename Model>
    void KenLM<Model>::apply_scan(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  const int dot,
				  feature_set_type& features,
				  const bool final) const
    {
      
    }
    
    template <typename Model>
    void KenLM<Model>::apply_complete(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    FeatureFunction::feature_function_ptr_type KenLMFactory::create(const std::string& parameter) const
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      typedef FeatureFunction::feature_function_ptr_type feature_function_ptr_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "kenlm")
	throw std::runtime_error("is this really kenlm feature function? " + parameter);
      
      path_type path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	if (utils::ipiece(piter->first) == "file") {
	  path = piter->second;
	  break;
	}
      
      if (path.empty() || ! boost::filesystem::exists(path))
	throw std::runtime_error("no filename? " + path.string());
      
      lm::ngram::ModelType m;
      if (! lm::ngram::RecognizeBinary(path.string().c_str(), m))
	m = lm::ngram::PROBING;
      
      switch (m) { 
      case lm::ngram::PROBING:
	return feature_function_ptr_type(new KenLM<lm::ngram::ProbingModel>(parameter));
      case lm::ngram::REST_PROBING:
	return feature_function_ptr_type(new KenLM<lm::ngram::RestProbingModel>(parameter));
      case lm::ngram::TRIE:
	return feature_function_ptr_type(new KenLM<lm::ngram::TrieModel>(parameter));
      case lm::ngram::ARRAY_TRIE:
	return feature_function_ptr_type(new KenLM<lm::ngram::ArrayTrieModel>(parameter));
      case lm::ngram::QUANT_TRIE:
	return feature_function_ptr_type(new KenLM<lm::ngram::QuantTrieModel>(parameter));
      case lm::ngram::QUANT_ARRAY_TRIE:
	return feature_function_ptr_type(new KenLM<lm::ngram::QuantArrayTrieModel>(parameter));
      default:
	throw std::runtime_error("Unrecognized kenlm binary file type:" + boost::lexical_cast<std::string>(m));
      }
    }
  };
};
