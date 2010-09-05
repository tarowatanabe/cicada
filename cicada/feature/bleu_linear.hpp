// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__BLEU_LINEAR__HPP__
#define __CICADA__FEATURE__BLEU_LINEAR__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <cicada/eval.hpp>

namespace cicada
{
  namespace feature
  {
    
    class BleuLinearImpl;
    
    class BleuLinear : public FeatureFunction
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      typedef cicada::eval::Score score_type;
      typedef score_type::score_ptr_type score_ptr_type;

    private:
      typedef FeatureFunction base_type;
      typedef BleuLinearImpl        impl_type;
      
    public:
      BleuLinear(const std::string& parameter);
      BleuLinear(const BleuLinear&);

      ~BleuLinear();

      BleuLinear& operator=(const BleuLinear&);
      
    private:
      BleuLinear() {}
      
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
      virtual void assign(const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts);
      
      virtual void initialize();

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new BleuLinear(*this)); }
      
      void clear();
      
      void assign(const score_ptr_type& score);
      
    private:
      impl_type* pimpl;
    };
    
  };
};


#endif

