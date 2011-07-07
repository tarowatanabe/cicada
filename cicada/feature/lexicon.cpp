//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "lexicon.hpp"
#include "cicada/parameter.hpp"

#include "utils/piece.hpp"

#include <google/dense_hash_set>
#include <google/dense_hash_map>

namespace cicada
{
  namespace feature
  {

    class LexiconImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

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
      
      typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;
      typedef google::dense_hash_map<word_pair_type, feature_type, utils::hashmurmur<size_t>, std::equal_to<word_pair_type> > cache_set_type;
      
      
      LexiconImpl() : uniques(), words(), caches(), forced_feature(false) { uniques.set_empty_key(word_type()); caches.set_empty_key(word_pair_type()); }
      
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features)
      {
	const phrase_type& phrase = edge.rule->rhs;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) {
	  
	  sentence_type::const_iterator witer_end = words.end();
	  for (sentence_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) {
	    
	    std::pair<cache_set_type::iterator, bool> result = caches.insert(std::make_pair(word_pair_type(*witer, *piter), feature_type()));
	    if (result.second) {
	      if (forced_feature)
		result.first->second = "lexicon:" + static_cast<const std::string&>(*witer) + ":" + static_cast<const std::string&>(*piter);
	      else {
		const std::string name = "lexicon:" + static_cast<const std::string&>(*witer) + ":" + static_cast<const std::string&>(*piter);
		if (feature_type::exists(name))
		  result.first->second = name;
	      }
	    }
	    
	    if (! result.first->second.empty())
	      features[result.first->second] += 1.0;
	  }
	}
      }

      void assign(const lattice_type& lattice)
      {
	uniques.clear();
	caches.clear();
	
	lattice_type::const_iterator liter_end = lattice.end();
	for (lattice_type::const_iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	  for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	    if (aiter->label != vocab_type::EPSILON)
	      uniques.insert(aiter->label);
	}

	words.clear();
	words.insert(words.end(), uniques.begin(), uniques.end());
      }
      
      void assign(const hypergraph_type& forest)
      {
	uniques.clear();
	caches.clear();
	
	hypergraph_type::edge_set_type::const_iterator eiter_end = forest.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = forest.edges.begin(); eiter != eiter_end; ++ eiter)
	  if (eiter->rule) {
	    const rule_type& rule = *(eiter->rule);
	    
	    rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	    for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	      if (siter->is_terminal() && *siter != vocab_type::EPSILON && *siter != vocab_type::BOS && *siter != vocab_type::EOS)
		uniques.insert(*siter);
	  }
	
	words.clear();
	words.insert(words.end(), uniques.begin(), uniques.end());
      }
      
      void clear()
      {
	words.clear();
	caches.clear();
      }
      
      word_set_type  uniques;
      sentence_type  words;
      cache_set_type caches;
      
      bool forced_feature;
    };
    
    Lexicon::Lexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "lexicon")
	throw std::runtime_error("this is not Lexicon feature: " + parameter);
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for lexicon: " << piter->first << "=" << piter->second << std::endl;
      
      std::auto_ptr<impl_type> lexicon_impl(new impl_type());
      
      // two-side context + length + counts-id 
      base_type::__state_size = 0;
      base_type::__feature_name = "lexicon";
      base_type::__sparse_feature = true;
      
      pimpl = lexicon_impl.release();
    }
    
    Lexicon::~Lexicon() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    Lexicon::Lexicon(const Lexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    Lexicon& Lexicon::operator=(const Lexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Lexicon::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));

      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->lexicon_score(edge, features);
    }

    void Lexicon::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {
      
    }
    
    void Lexicon::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    void Lexicon::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Lexicon::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    
    void Lexicon::assign(const size_type& id,
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
      else if (hypergraph.is_valid())
	pimpl->assign(hypergraph);

    }
    
  };
};
