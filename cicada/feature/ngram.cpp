//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/ngram.hpp"
#include "cicada/ngram_scorer.hpp"

#include "cicada/feature/ngram.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"
#include "cicada/cluster.hpp"

#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"

// faster ngram state representation inspired by
//
// Unary Data Structures for Language Models
// Jeffrey Sorensen, Cyril Allauzen
// INTERSPEECH 2011
//
// and
// 
// @InProceedings{pauls-klein:2011:ACL-HLT2011,
//  author    = {Pauls, Adam  and  Klein, Dan},
//  title     = {Faster and Smaller N-Gram Language Models},
//  booktitle = {Proceedings of the 49th Annual Meeting of the Association for Computational Linguistics: Human Language Technologies},
//  month     = {June},
//  year      = {2011},
//  address   = {Portland, Oregon, USA},
//  publisher = {Association for Computational Linguistics},
//  pages     = {258--267},
//  url       = {http://www.aclweb.org/anthology/P11-1027}
// }

// and employes the idea from
//
// @InProceedings{heafield-koehn-lavie:2012:EMNLP-CoNLL,
//   author    = {Heafield, Kenneth  and  Koehn, Philipp  and  Lavie, Alon},
//   title     = {Language Model Rest Costs and Space-Efficient Storage},
//   booktitle = {Proceedings of the 2012 Joint Conference on Empirical Methods in Natural Language Processing and Computational Natural Language Learning},
//   month     = {July},
//   year      = {2012},
//   address   = {Jeju Island, Korea},
//   publisher = {Association for Computational Linguistics},
//   pages     = {1169--1178},
//   url       = {http://www.aclweb.org/anthology/D12-1107}
// }

namespace cicada
{
  namespace feature
  {
    class NGramImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Symbol word_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::NGram  ngram_type;
      
      typedef cicada::Cluster cluster_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      typedef std::vector<char, std::allocator<char> > buffer_type;
      
      typedef boost::filesystem::path path_type;
      
    public:
      
      NGramImpl(const path_type& __path, const bool populate)
	: ngram(&ngram_type::create(__path)),
	  order(0), cluster(0), approximate(false), no_bos_eos(false), skip_sgml_tag(false)
      {
	if (populate)
	  ngram->populate();
	
	// set up correct ordering...
	order = ngram->index.order();
	
	id_oov = ngram->index.vocab()[vocab_type::UNK];
	id_bos = ngram->index.vocab()[vocab_type::BOS];
	id_eos = ngram->index.vocab()[vocab_type::EOS];
	
	scorer.assign(*ngram);

	buffer_bos.reserve(scorer.buffer_size());
	buffer_tmp.reserve(scorer.buffer_size());
	
	buffer_bos.resize(scorer.buffer_size());
	buffer_tmp.resize(scorer.buffer_size());

	scorer.ngram_state_.clear(&(*buffer_bos.begin()));
	
	ngram->lookup_context(&id_bos, (&id_bos) + 1, scorer.ngram_state_.suffix(&(*buffer_bos.begin())));
      }

      NGramImpl(const NGramImpl& x)
	: ngram(&ngram_type::create(x.ngram->path())),
	  order(x.order),
	  cluster(x.cluster ? &cluster_type::create(x.cluster->path()) : 0),
	  approximate(x.approximate),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  feature_name(x.feature_name),
	  feature_name_oov(x.feature_name_oov),
	  id_oov(x.id_oov),
	  id_bos(x.id_bos),
	  id_eos(x.id_eos)
      {
	scorer.assign(*ngram);

	buffer_bos.reserve(scorer.buffer_size());
	buffer_tmp.reserve(scorer.buffer_size());
	
	buffer_bos.resize(scorer.buffer_size());
	buffer_tmp.resize(scorer.buffer_size());

	scorer.ngram_state_.clear(&(*buffer_bos.begin()));
	
	ngram->lookup_context(&id_bos, (&id_bos) + 1, scorer.ngram_state_.suffix(&(*buffer_bos.begin())));
      }

      NGramImpl& operator=(const NGramImpl& x)
      {
	ngram = &ngram_type::create(x.ngram->path());
	order = x.order;
	cluster = (x.cluster ? &cluster_type::create(x.cluster->path()) : 0);
	approximate = x.approximate;
	no_bos_eos = x.no_bos_eos;
	skip_sgml_tag = x.skip_sgml_tag;
	
	feature_name     = x.feature_name;
	feature_name_oov = x.feature_name_oov;
	id_oov           = x.id_oov;
	id_bos           = x.id_bos;
	id_eos           = x.id_eos;
	
	scorer.assign(*ngram);
	
	buffer_bos.reserve(scorer.buffer_size());
	buffer_tmp.reserve(scorer.buffer_size());
	
	buffer_bos.resize(scorer.buffer_size());
	buffer_tmp.resize(scorer.buffer_size());
	
	scorer.ngram_state_.clear(&(*buffer_bos.begin()));
	
	ngram->lookup_context(&id_bos, (&id_bos) + 1, scorer.ngram_state_.suffix(&(*buffer_bos.begin())));
	
	return *this;
      }

      size_type buffer_size() const
      {
	return scorer.buffer_size();
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
	
	scorer.assign(state);
	
	phrase_type::const_iterator titer_begin = target.begin();
        phrase_type::const_iterator titer_end   = target.end();
	
	int non_terminal_pos = 0;
	bool initial = true;
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  if (titer->is_non_terminal()) {
	    const int __non_terminal_index = titer->non_terminal_index();
            const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
            ++ non_terminal_pos;

	    if (initial)
	      scorer.initial_non_terminal(states[antecedent_index]);
	    else
	      scorer.non_terminal(states[antecedent_index]);
	    
	    initial = false;
	  } else if (! skipper(*titer)) {
	    
	    if (no_bos_eos && extract(*titer) == vocab_type::BOS && scorer.ngram_state_.empty(state))
	      scorer.initial_bos(&(*buffer_bos.begin()));
	    else {
	      const word_type::id_type id = ngram->index.vocab()[extract(*titer)];
	      
	      oov += (id == id_oov);
	      
	      scorer.terminal(id);
	    }
	    
	    initial = false;
	  }
	}

	return scorer.complete();
      }
      
      
      double ngram_final_score(const state_ptr_type& state)
      {
	if (no_bos_eos)
	  return 0.0;
	else {
	  scorer.assign(&(*buffer_tmp.begin()));
	  
	  scorer.initial_bos(&(*buffer_bos.begin()));
	  scorer.non_terminal(state);
	  scorer.terminal(id_eos);
	  
	  return scorer.complete();
	}
      }

      double ngram_coarse_score(const edge_type& edge,
				int& oov)
      {
	if (cluster) {
          if (skip_sgml_tag)
            return ngram_coarse_score(edge, oov, extract_cluster(cluster), skipper_sgml());
          else
            return ngram_coarse_score(edge, oov, extract_cluster(cluster), skipper_epsilon());
        } else {
          if (skip_sgml_tag)
            return ngram_coarse_score(edge, oov, extract_word(), skipper_sgml());
          else
            return ngram_coarse_score(edge, oov, extract_word(), skipper_epsilon());
        }
      }
      
      template <typename Extract, typename Skipper>
      double ngram_coarse_score(const edge_type& edge,
				int& oov,
				Extract extract,
				Skipper skipper)
      {
	const rule_type& rule = *(edge.rule);
	const phrase_type& phrase = rule.rhs;
	
	scorer.assign(&(*buffer_tmp.begin()));
	
	double score = 0.0;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) {
	  if (piter->is_non_terminal()) {
	    score += scorer.complete();
	    
	    scorer.assign(&(*buffer_tmp.begin()));
	  } else if (! skipper(*piter)) {
	    
	    if (no_bos_eos && extract(*piter) == vocab_type::BOS && scorer.ngram_state_.empty(&(*buffer_tmp.begin())))
	      scorer.initial_bos(&(*buffer_bos.begin()));
	    else {
	      const word_type::id_type id = ngram->index.vocab()[extract(*piter)];
	      
	      oov += (id == id_oov);
	      scorer.terminal(id);
	    }
	  }
	}
	
	return score + scorer.complete();
      }
      
      double ngram_predict_score(state_ptr_type& state)
      {
	scorer.assign(state);
	
	if (! no_bos_eos)
	  scorer.ngram_state_.suffix_.copy(scorer.ngram_state_.suffix(&(*buffer_bos.begin())),
					   scorer.ngram_state_.suffix(state));
	
	return 0.0;
      }

      double ngram_scan_score(state_ptr_type& state,
			      const edge_type& edge,
			      const int dot,
			      int& oov)
      {
	if (cluster) {
          if (skip_sgml_tag)
            return ngram_scan_score(state, edge, dot, oov, extract_cluster(cluster), skipper_sgml());
          else
            return ngram_scan_score(state, edge, dot, oov, extract_cluster(cluster), skipper_epsilon());
        } else {
          if (skip_sgml_tag)
            return ngram_scan_score(state, edge, dot, oov, extract_word(), skipper_sgml());
          else
            return ngram_scan_score(state, edge, dot, oov, extract_word(), skipper_epsilon());
        }
      }
      
      template <typename Extract, typename Skipper>
      double ngram_scan_score(state_ptr_type& state,
                              const edge_type& edge,
                              const int dot,
                              int& oov,
                              Extract extract,
                              Skipper skipper)
      {
	const rule_type& rule = *(edge.rule);
	const phrase_type& phrase = rule.rhs;

	double score = 0.0;
	
	void* state_curr = state;
	void* state_next = &(*buffer_tmp.begin());

	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin() + dot; piter != piter_end && ! piter->is_non_terminal(); ++ piter)
	  if (! skipper(*piter)) {
	    
	    if (no_bos_eos && extract(*piter) == vocab_type::BOS && scorer.ngram_state_.size_suffix(state_curr) == 0)
	      scorer.ngram_state_.suffix_.copy(scorer.ngram_state_.suffix(&(*buffer_bos.begin())),
					       scorer.ngram_state_.suffix(state_next));
	    else {
	      const word_type::id_type id = ngram->index.vocab()[extract(*piter)];
	      
	      oov += (id == id_oov);
	      
	      score += ngram->ngram_score(scorer.ngram_state_.suffix(state_curr), 
					  id_eos,
					  scorer.ngram_state_.suffix(state_next)).prob;
	    }
	    
	    std::swap(state_curr, state_next);
	  }
	
	// copy state...
	if (state_curr != state)
	  scorer.ngram_state_.suffix_.copy(scorer.ngram_state_.suffix(state_next),
					   scorer.ngram_state_.suffix(state));
	
	scorer.ngram_state_.suffix_.fill(scorer.ngram_state_.suffix(state));
	
	return score;
      }
      
      double ngram_complete_score(state_ptr_type& state)
      {
	if (no_bos_eos)
	  return 0.0;
	
	void* state_tmp = &(*buffer_tmp.begin());
	
	const double prob = ngram->ngram_score(scorer.ngram_state_.suffix(state), 
					       id_eos,
					       scorer.ngram_state_.suffix(state_tmp)).prob;
	
	scorer.ngram_state_.suffix_.copy(scorer.ngram_state_.suffix(state_tmp),
					 scorer.ngram_state_.suffix(state));
	scorer.ngram_state_.suffix_.fill(scorer.ngram_state_.suffix(state));
	
	return prob;
      }
      
      // ngrams
      ngram_type* ngram;
      int         order;
      
      // cluster...
      cluster_type* cluster;
      
      bool approximate;
      bool no_bos_eos;
      bool skip_sgml_tag;
      
      // names...
      feature_type feature_name;
      feature_type feature_name_oov;
      
      // bos/eos etc...
      symbol_type::id_type id_oov;
      symbol_type::id_type id_bos;
      symbol_type::id_type id_eos;
      
      cicada::NGramScorer scorer;
      
      // bos state
      buffer_type buffer_bos;
      buffer_type buffer_tmp;
    };
    
    NGram::NGram(const std::string& parameter)
      : pimpl(0), pimpl_coarse(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "ngram")
	throw std::runtime_error("is this really ngram feature function? " + parameter);

      path_type   path;
      bool        populate = false;
      path_type   cluster_path;
      bool        approximate = false;
      bool        skip_sgml_tag = false;
      bool        no_bos_eos = false;
      
      path_type   coarse_path;
      bool        coarse_populate = false;
      path_type   coarse_cluster_path;
      bool        coarse_approximate = false;
      
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "populate")
	  populate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "cluster")
	  cluster_path = piter->second;
	else if (utils::ipiece(piter->first) == "approximate")
	  approximate = utils::lexical_cast<bool>(piter->second);
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
	else if (utils::ipiece(piter->first) == "coarse-approximate")
	  coarse_approximate = utils::lexical_cast<bool>(piter->second);
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

      if (ngram_impl->order <= 0)
	throw std::runtime_error("invalid ngram order: " + utils::lexical_cast<std::string>(ngram_impl->order));
            
      ngram_impl->approximate = approximate;
      ngram_impl->no_bos_eos = no_bos_eos;
      ngram_impl->skip_sgml_tag = skip_sgml_tag;
      
      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.string());
	
	ngram_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      // two contexts (order - 1) for each edge, with two separator..
      base_type::__state_size = ngram_impl->buffer_size();
      base_type::__feature_name = (name.empty() ? std::string("ngram") : name);
      
      ngram_impl->feature_name     = base_type::__feature_name;
      ngram_impl->feature_name_oov = static_cast<const std::string&>(base_type::__feature_name) + ":oov-penalty";
      
      pimpl = ngram_impl.release();

      // ...
      if (! coarse_path.empty()) {
	std::auto_ptr<impl_type> ngram_impl(new impl_type(coarse_path, coarse_populate));
	
	if (ngram_impl->order <= 0)
	  throw std::runtime_error("invalid coarse ngram order: " + utils::lexical_cast<std::string>(ngram_impl->order));
	
	ngram_impl->approximate = coarse_approximate;
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
    
    NGram::~NGram()
    {
      std::auto_ptr<impl_type> tmp(pimpl);
      if (pimpl_coarse)
	std::auto_ptr<impl_type> tmp_coarse(pimpl_coarse);
    }
    
    NGram::NGram(const NGram& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl)),
	pimpl_coarse(x.pimpl_coarse ? new impl_type(*x.pimpl_coarse) : 0)
    {}
    
    NGram& NGram::operator=(const NGram& x)
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
    
    
    void NGram::apply(state_ptr_type& state,
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

    void NGram::apply_coarse(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     const bool final) const
    {
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
      } else {
	// state-less.... here, we ignored final flag...do we add this...?
	int oov = 0;
	const double score = pimpl->ngram_coarse_score(edge, oov);
	
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);

	if (oov)
	  features[pimpl->feature_name_oov] = - oov;
	else
	  features.erase(pimpl->feature_name_oov);
      }
    }

    // temporarily assigned feature function...
    
    void NGram::apply_predict(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      const bool final) const
    {
      // add <s>
      if (final)
	pimpl->ngram_predict_score(state);
    }
    
    void NGram::apply_scan(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   const int dot,
			   feature_set_type& features,
			   const bool final) const
    {
      int oov = 0;
      const double score = pimpl->ngram_scan_score(state, edge, dot, oov);
      
      if (score != 0.0)
	features[pimpl->feature_name] = score;
      else
	features.erase(pimpl->feature_name);
      
      if (oov)
	features[pimpl->feature_name_oov] = - oov;
      else
	features.erase(pimpl->feature_name_oov);
    }
    
    void NGram::apply_complete(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       const bool final) const
    {
      // if final, add scoring for </s>
      
      if (final) {
	const double score = pimpl->ngram_complete_score(state);
	
	if (score != 0.0)
	  features[pimpl->feature_name] = score;
	else
	  features.erase(pimpl->feature_name);
      }
    }


  };
};
