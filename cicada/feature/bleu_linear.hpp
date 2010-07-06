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
      ~BleuLinear();
      
    private:
      BleuLinear() {}
      BleuLinear(const BleuLinear&) {}
      BleuLinear& operator=(const BleuLinear&) { return *this; }
      
    public:
      virtual void operator()(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates) const;
      virtual void operator()(const state_ptr_type& state,
			      feature_set_type& features,
			      feature_set_type& estimates) const {}
      virtual void initialize();
      
      void clear();
      void insert(const int source_size, const sentence_type& sentence);
      void insert(const score_ptr_type& score);
      
    private:
      impl_type* pimpl;
    };
    
  };
};


#endif

