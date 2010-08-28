// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__GLOBAL_LEXICON__HPP__
#define __CICADA__FEATURE__GLOBAL_LEXICON__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

namespace cicada
{
  namespace feature
  {
    
    class GlobalLexiconImpl;
    
    class GlobalLexicon : public FeatureFunction
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
    private:
      typedef FeatureFunction   base_type;
      typedef GlobalLeixconImpl impl_type;
      
    public:
      GlobalLexicon(const std::string& parameter);
      GlobalLexicon(const GlobalLexicon&);
      
      ~GlobalLexicon();

      GlobalLexicon& operator=(const GlobalLexicon&);
      
    private:
      GlobalLexicon() {}
      
    public:
      virtual void apply(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const;
      virtual void apply_coarse(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const;
      virtual void initialize();
      
      virtual void assign(const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans);
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new GlobalLexicon(*this)); }
      
    private:
      impl_type* pimpl;
    };

    
  };
};


#endif


