// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__BLEU__HPP__
#define __CICADA__FEATURE__BLEU__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

namespace cicada
{
  namespace feature
  {
    
    class BleuImpl;
    
    class Bleu : public FeatureFunction
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      
    private:
      typedef FeatureFunction base_type;
      typedef BleuImpl        impl_type;
      
    public:
      Bleu(const std::string& parameter);
      ~Bleu();
      
    private:
      Bleu() {}
      Bleu(const Bleu&) {}
      Bleu& operator=(const Bleu&) { return *this; }
      
    public:
      virtual void operator()(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates) const;
      virtual void operator()(state_ptr_type& state,
			      feature_set_type& features) const;
      
      void clear();
      void insert(const int source_size, const sentence_type& sentence);
      
    private:
      impl_type* pimpl;
    };
    
  };
};


#endif

