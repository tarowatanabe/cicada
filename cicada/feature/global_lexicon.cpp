//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "cicada/global_lexicon.hpp"
#include "cicada/feature/global_lexicon.hpp"
#include "cicada/parameter.hpp"

#include "utils/piece.hpp"

#include <google/dense_hash_set>

namespace cicada
{

  namespace feature
  {
  
    class GlobalLexiconImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::GlobalLexicon lexicon_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;

      typedef lexicon_type::path_type path_type;
      
      typedef symbol_type word_type;
      
      typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;
      
      typedef std::vector<double, std::allocator<double> > cache_set_type;
      
      GlobalLexiconImpl(const path_type& __path)
	: lexicon(__path)
      {
	words.set_empty_key(symbol_type());
      }
      
      double global_lexicon_score(const edge_type& edge)
      {
	const phrase_type& phrase = edge.rule->rhs;
	
	double score = 0.0;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	  if (*piter != vocab_type::EPSILON && *piter != vocab_type::BOS && *piter != vocab_type::EOS && piter->is_terminal())
	    score += score_lexicon(*piter);
	
	return score;
      }

      void assign(const lattice_type& lattice)
      {
	caches.clear();
	words.clear();
	
	lattice_type::const_iterator liter_end = lattice.end();
	for (lattice_type::const_iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	  for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	    if (aiter->label != vocab_type::EPSILON)
	      words.insert(aiter->label);
	}
      }
      
      void assign(const hypergraph_type& forest)
      {
	caches.clear();
	words.clear();
	
	hypergraph_type::edge_set_type::const_iterator eiter_end = forest.edges.end();
	for (hypergraph_type::edge_set_type::const_iterator eiter = forest.edges.begin(); eiter != eiter_end; ++ eiter)
	  if (eiter->rule) {
	    const rule_type& rule = *(eiter->rule);
	    
	    rule_type::symbol_set_type::const_iterator siter_end = rule.rhs.end();
	    for (rule_type::symbol_set_type::const_iterator siter = rule.rhs.begin(); siter != siter_end; ++ siter)
	      if (siter->is_terminal() && *siter != vocab_type::EPSILON && *siter != vocab_type::BOS && *siter != vocab_type::EOS)
		words.insert(*siter);
	  }
      }

      double score_lexicon(const symbol_type& word)
      {
	if (word.id() >= caches.size())
	  caches.resize(word.id() + 1, 0.0);
	
	if (caches[word.id()] == 0.0)
	  caches[word.id()] = lexicon(word, words.begin(), words.end());
	
	return caches[word.id()];
      }
      
      word_set_type  words;
      lexicon_type   lexicon;
      cache_set_type caches;
    };
  
  
    GlobalLexicon::GlobalLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (utils::ipiece(param.name()) != "global-lexicon")
	throw std::runtime_error("is this really a global-lexicon feature function? " + parameter);

      boost::filesystem::path path;
       
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for global lexicon: " << piter->first << "=" << piter->second << std::endl;
      }

      if (path.empty())
	throw std::runtime_error("no global lexicon file? " + path.string());
      
      std::auto_ptr<impl_type> global_lexicon_impl(new impl_type(path));
      
      // no-context...
      base_type::__state_size = 0;
      base_type::__feature_name = std::string("global-lexicon");
      
      pimpl = global_lexicon_impl.release();
    }
    
    GlobalLexicon::~GlobalLexicon() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    GlobalLexicon::GlobalLexicon(const GlobalLexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    GlobalLexicon& GlobalLexicon::operator=(const GlobalLexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void GlobalLexicon::apply(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      const double score = pimpl->global_lexicon_score(edge);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
    }

    void GlobalLexicon::apply_coarse(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     feature_set_type& estimates,
				     const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    void GlobalLexicon::apply_predict(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      feature_set_type& estimates,
				      const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void GlobalLexicon::apply_scan(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   const int dot,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
    void GlobalLexicon::apply_complete(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       feature_set_type& estimates,
				       const bool final) const
    {}


    void GlobalLexicon::assign(const size_type& id,
			       const hypergraph_type& hypergraph,
			       const lattice_type& lattice,
			       const span_set_type& spans,
			       const sentence_set_type& targets,
			       const ngram_count_set_type& ngram_counts)
    {
      pimpl->words.clear();
      pimpl->caches.clear();
      
      if (! lattice.empty())
	pimpl->assign(lattice);
      else if (hypergraph.is_valid())
	pimpl->assign(hypergraph);
    }

  };
  
};
