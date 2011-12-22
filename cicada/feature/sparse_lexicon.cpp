//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <set>

#include "sparse_lexicon.hpp"

#include "cicada/lexicon.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/alloc_vector.hpp"

#include <google/dense_hash_set>
#include <google/dense_hash_map>

namespace cicada
{
  namespace feature
  {

    class SparseLexiconImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;
      
      typedef cicada::ClusterStemmer normalizer_type;
      typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef symbol_type word_type;
      typedef std::pair<word_type, word_type> word_pair_type;
      typedef std::vector<word_pair_type, std::allocator<word_pair_type> > word_pair_set_type;

      typedef std::vector<word_type, std::allocator<word_type> > word_set_type;
      typedef std::vector<word_set_type, std::allocator<word_set_type> > word_map_type;
      typedef std::set<size_type, std::less<size_type>, std::allocator<size_type> > pos_set_type;
      
      typedef google::dense_hash_set<word_pair_type, utils::hashmurmur<size_t>, std::equal_to<word_pair_type> > word_pair_unique_type;
      
      typedef utils::alloc_vector<feature_set_type, std::allocator<feature_set_type> > cache_set_type;
      
      
      typedef cicada::Lexicon lexicon_type;
      
      SparseLexiconImpl()
	: uniques_prev(), uniques_next(), caches(),
	  lexicon_prev(0), lexicon_next(0),
	  skip_sgml_tag(false), unique_source(false), prefix("sparse-lexicon"), forced_feature(false)
      {
	uniques_prev.set_empty_key(word_pair_type()); 
	uniques_next.set_empty_key(word_pair_type()); 
      }
      
      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS;
	}
      };

      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS || word.is_sgml_tag();
	}
      };
      
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features)
      {
	if (skip_sgml_tag)
	  lexicon_score(edge, features, skipper_sgml());
	else
	  lexicon_score(edge, features, skipper_epsilon());
      }
      
      
      template <typename Skipper>
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features,
			 Skipper skipper)
      {
	const phrase_type& phrase = edge.rule->rhs;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) 
	  if (piter->is_terminal() && ! skipper(*piter)) {
	    if (! caches.exists(piter->id())) {
	      feature_set_type& features = caches[piter->id()];
	      features.clear();
	      
	      {
		word_pair_set_type::const_iterator witer_end = sources_prev.end();
		for (word_pair_set_type::const_iterator witer = sources_prev.begin(); witer != witer_end; ++ witer) 
		  if (! lexicon_prev || exists(*lexicon_prev, *witer, *piter))
		    assign("prev", *witer, *piter, features);
	      }
	      
	      {
		word_pair_set_type::const_iterator witer_end = sources_next.end();
		for (word_pair_set_type::const_iterator witer = sources_next.begin(); witer != witer_end; ++ witer)
		  if (! lexicon_next || exists(*lexicon_next, *witer, *piter))
		    assign("next", *witer, *piter, features);
	      }
	    }
	    
	    features += caches[piter->id()];
	  }
      }

      word_set_type words_source_prev;
      word_set_type words_source_next;
      word_set_type words_target;
  
      void assign(const char* id, const word_pair_type& prev, const word_type& next, feature_set_type& features)
      {
	if (normalizers_source.empty() && normalizers_target.empty()) {
	  const std::string name = (prefix + ':' + id + ':'
				    + static_cast<const std::string&>(prev.first)
				    + ':' + static_cast<const std::string&>(prev.second)
				    + '_' + static_cast<const std::string&>(next));
	  
	  if (forced_feature || feature_type::exists(name))
	    features[name] += 1.0;
	} else {
	  words_source_prev.clear();
	  words_source_next.clear();
	  words_target.clear();
	  
	  words_source_prev.push_back(prev.first);
	  words_source_next.push_back(prev.second);
	  for (size_t i = 0; i != normalizers_source.size(); ++ i) {
	    {
	      const word_type normalized = normalizers_source[i](prev.first);
	      if (normalized != prev.first)
		words_source_prev.push_back(normalized);
	    }
	    {
	      const word_type normalized = normalizers_source[i](prev.second);
	      if (normalized != prev.second)
		words_source_next.push_back(normalized);
	    }
	  }
	  
	  words_target.push_back(next);
	  for (size_t i = 0; i != normalizers_target.size(); ++ i) {
	    const word_type normalized = normalizers_target[i](next);
	    if (normalized != next)
	      words_target.push_back(normalized);
	  }
	  
	  word_set_type::const_iterator piter_begin = words_source_prev.begin();
	  word_set_type::const_iterator piter_end   = words_source_prev.end();
	  word_set_type::const_iterator niter_begin = words_source_next.begin();
	  word_set_type::const_iterator niter_end   = words_source_next.end();
	  word_set_type::const_iterator titer_begin = words_target.begin();
	  word_set_type::const_iterator titer_end   = words_target.end();
	  
	  for (word_set_type::const_iterator piter = piter_begin; piter != piter_end; ++ piter)
	    for (word_set_type::const_iterator niter = niter_begin; niter != niter_end; ++ niter)
	      for (word_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
		const std::string name = (prefix + ':' + id + ':'
					  + static_cast<const std::string&>(*piter)
					  + ':' + static_cast<const std::string&>(*niter)
					  + '_' + static_cast<const std::string&>(*titer));
		
		if (forced_feature || feature_type::exists(name))
		  features[name] += 1.0;
	      }
	}
      }

      void assign(const lattice_type& lattice)
      {
	if (skip_sgml_tag)
	  assign(lattice, skipper_sgml());
	else
	  assign(lattice, skipper_epsilon());
      }
      
      template <typename Skipper>
      void assign(const lattice_type& lattice, Skipper skipper)
      {
	word_set_type words;
	pos_set_type  positions;
	
	clear();
	
	lattice_prev.clear();
	lattice_prev.resize(lattice.size() + 1);
	lattice_prev.front().push_back(vocab_type::BOS);
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	  positions.clear();
	  
	  lattice_type::arc_set_type::const_iterator aiter_end = lattice[pos].end();
	  for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter) {
	    if (! skipper(aiter->label)) {
	      words.clear();
	      words.push_back(aiter->label);
	      
	      lattice_prev[pos + aiter->distance].insert(lattice_prev[pos + aiter->distance].end(), words.begin(), words.end());
	      
	      // we will compute pair of lattice_prev[pos] and words
	      word_set_type::const_iterator piter_end = lattice_prev[pos].end();
	      for (word_set_type::const_iterator piter = lattice_prev[pos].begin(); piter != piter_end; ++ piter) {
		
		word_set_type::const_iterator niter_end = words.end();
		for (word_set_type::const_iterator niter = words.begin(); niter != niter_end; ++ niter) {
		  
		  if (! lexicon_prev || exists(*lexicon_prev, *piter, *niter))
		    sources_prev.push_back(std::make_pair(*piter, *niter));
		  
		  if (*piter != vocab_type::BOS && (! lexicon_next || exists(*lexicon_next, *piter, *niter)))
		    sources_next.push_back(std::make_pair(*piter, *niter));
		}
	      }
	    } else
	      positions.insert(pos + aiter->distance);
	  }
	  
	  // copy lattice_prev[pos] into  positons.
	  pos_set_type::const_iterator piter_end = positions.end();
	  for (pos_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter)
	    lattice_prev[*piter].insert(lattice_prev[*piter].end(), lattice_prev[pos].begin(), lattice_prev[pos].end());
	}
	
	// we will compute pair of lattice_prev[lattice.size()] and EOS
	word_set_type::const_iterator piter_end = lattice_prev[lattice.size()].end();
	for (word_set_type::const_iterator piter = lattice_prev[lattice.size()].begin(); piter != piter_end; ++ piter)
	  if (! lexicon_next || exists(*lexicon_next, *piter, vocab_type::EOS))
	    sources_next.push_back(std::make_pair(*piter, vocab_type::EOS));
	
	if (unique_source) {
	  uniques_prev.clear();
	  uniques_next.clear();
	  
	  uniques_prev.insert(sources_prev.begin(), sources_prev.end());
	  uniques_next.insert(sources_next.begin(), sources_next.end());
	  
	  sources_prev.clear();
	  sources_next.clear();
	  
	  sources_prev.insert(sources_prev.end(), uniques_prev.begin(), uniques_prev.end());
	  sources_next.insert(sources_next.end(), uniques_next.begin(), uniques_next.end());
	}
	
	std::sort(sources_prev.begin(), sources_prev.end());
	std::sort(sources_next.begin(), sources_next.end());
      }
      
      void clear()
      {
	sources_prev.clear();
	sources_next.clear();
	caches.clear();
      }

      bool exists(const lexicon_type& lexicon, const word_pair_type& prev, const word_type& next)
      {
	word_type codes[2];
	codes[0] = prev.first;
	codes[1] = prev.second;
	return lexicon.exists(codes, codes + 2, next);
      }

      bool exists(const lexicon_type& lexicon, const word_type& prev, const word_type& next)
      {
	word_type codes[2];
	codes[0] = prev;
	codes[1] = next;
	return lexicon.exists(codes, codes + 2);
      }
      
      normalizer_set_type normalizers_source;
      normalizer_set_type normalizers_target;
      
      word_pair_unique_type  uniques_prev;
      word_pair_unique_type  uniques_next;
      word_map_type          lattice_prev;
      
      word_pair_set_type sources_prev;
      word_pair_set_type sources_next;
      cache_set_type caches;

      lexicon_type* lexicon_prev;
      lexicon_type* lexicon_next;
      
      bool skip_sgml_tag;
      bool unique_source;
      
      std::string prefix;
      bool forced_feature;
    };
    
    SparseLexicon::SparseLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "sparse-lexicon")
	throw std::runtime_error("this is not sparse lexicon feature: " + parameter);
      
      impl_type::normalizer_set_type normalizers_source;
      impl_type::normalizer_set_type normalizers_target;
      
      bool skip_sgml_tag = false;
      bool unique_source = false;
      
      std::string name;

      std::string lexicon_prev;
      std::string lexicon_next;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "cluster-source") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers_source.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "cluster-target") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers_target.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "stemmer-source")
	  normalizers_source.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "stemmer-target")
	  normalizers_target.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "unique-source")
	  unique_source = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "lexicon-next")
	  lexicon_next = piter->second;
	else if (utils::ipiece(piter->first) == "lexicon-prev")
	  lexicon_prev = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for sparse lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> lexicon_impl(new impl_type());

      lexicon_impl->normalizers_source.swap(normalizers_source);
      lexicon_impl->normalizers_target.swap(normalizers_target);
      
      lexicon_impl->skip_sgml_tag = skip_sgml_tag;
      lexicon_impl->unique_source = unique_source;
      lexicon_impl->prefix = (name.empty() ? std::string("sparse-lexicon") : name);
      
      if (! lexicon_prev.empty())
	lexicon_impl->lexicon_prev = &cicada::Lexicon::create(lexicon_prev);
      if (! lexicon_next.empty())
	lexicon_impl->lexicon_next = &cicada::Lexicon::create(lexicon_next);

      
      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("sparse-lexicon") : name);
      base_type::__sparse_feature = true;
      
      pimpl = lexicon_impl.release();
    }
    
    SparseLexicon::~SparseLexicon() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    SparseLexicon::SparseLexicon(const SparseLexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {
      pimpl->lexicon_prev = 0;
      pimpl->lexicon_next = 0;
      if (x.pimpl->lexicon_prev)
	pimpl->lexicon_prev = &cicada::Lexicon::create(x.pimpl->lexicon_prev->path().string());
      if (x.pimpl->lexicon_next)
	pimpl->lexicon_next = &cicada::Lexicon::create(x.pimpl->lexicon_next->path().string());
    }

    SparseLexicon& SparseLexicon::operator=(const SparseLexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      pimpl->lexicon_prev = 0;
      pimpl->lexicon_next = 0;
      if (x.pimpl->lexicon_prev)
	pimpl->lexicon_prev = &cicada::Lexicon::create(x.pimpl->lexicon_prev->path().string());
      if (x.pimpl->lexicon_next)
	pimpl->lexicon_next = &cicada::Lexicon::create(x.pimpl->lexicon_next->path().string());
      
      return *this;
    }
    
    void SparseLexicon::apply(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));

      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->lexicon_score(edge, features);
    }

    void SparseLexicon::apply_coarse(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void SparseLexicon::apply_predict(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void SparseLexicon::apply_scan(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   const int dot,
				   feature_set_type& features,
				   const bool final) const
    {}
    void SparseLexicon::apply_complete(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       const bool final) const
    {}
    
    
    void SparseLexicon::assign(const size_type& id,
			       const hypergraph_type& hypergraph,
			       const lattice_type& lattice,
			       const span_set_type& spans,
			       const sentence_set_type& targets,
			       const ngram_count_set_type& ngram_counts)
    {
      //
      // how do we assign lexion feature from hypergraph...???
      // we assume that the lattice is always filled with the source-word...!
      //
      
      pimpl->clear();
      
      if (! lattice.empty())
	pimpl->assign(lattice);
    }
    
  };
};
