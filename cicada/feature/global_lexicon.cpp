
#include "cicada/global_lexicon.hpp"
#include "cicada/feature/global_lexicon.hpp"
#include "cicada/parameter.hpp"


#include "utils/sgi_hash_set.hpp"

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
      
#ifdef HAVE_TR1_UNORDERED_SET
      typedef std::tr1::unordered_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type > > word_set_type;
#else
      typedef sgi::hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type > > word_set_type;
#endif

      
      GlobalLexiconImpl(const path_type& __path, const bool __yield_source)
	: lexicon(__path), yield_source(__yield_source) {}
      
      double global_lexicon_score(const edge_type& edge)
      {
	const phrase_type& phrase = (yield_source ? edge.rule->source : edge.rule->target);

	double score = 0.0;

	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter)
	  if (*piter != vocab_type::EPSILON && piter->is_terminal())
	    score += lexicon(*piter, words.begin(), words.end());

	return score;
      }

      void assign(const lattice_type& lattice,
		  const hypergraph_type& hypergraph)
      {
	words.clear();
	
	if (! lattice.empty()) {
	  lattice_type::const_iterator liter_end = lattice.end();
	  for (lattice_type::const_iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	    lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	    for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	      words.insert(aiter->label);
	  }
	} else {
	  hypergraph_type::edge_set_type::const_iterator eiter_end = hypergraph.edges.end();
	  for (hypergraph_type::edge_set_type::const_iterator eiter = hypergraph.edges.begin(); eiter != eiter_end; ++ eiter) {
	    phrase_type::const_iterator piter_end = eiter->rule->source.end();
	    for (phrase_type::const_iterator piter = eiter->rule->source.begin(); piter != piter_end; ++ piter)
	      if (*piter != vocab_type::EPSILON && piter->is_terminal())
		words.insert(*piter);
	  }
	}
      }
      
      word_set_type words;
      lexicon_type  lexicon;

      bool yield_source;
    };
  
  
    GlobalLexicon::GlobalLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "global-lexicon")
	throw std::runtime_error("is this really a global-lexicon feature function? " + parameter);

      bool source = false;
      bool target = false;
      
      boost::filesystem::path path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	} else if (strcasecmp(piter->first.c_str(), "file") == 0)
	  path = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for global lexicon: " << piter->first << "=" << piter->second << std::endl;
      }

      if (path.empty())
	throw std::runtime_error("no global lexicon file? " + path.file_string());
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      
      // default to target
      if (! source && ! target)
	target = true;
      
      std::auto_ptr<impl_type> global_lexicon_impl(new impl_type(path, source));
      
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
      pimpl->assign(lattice, hypergraph);
    }

  };
  
};
