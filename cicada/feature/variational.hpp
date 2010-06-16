// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__VARIATIONAL__HPP__
#define __CICADA__FEATURE__VARIATIONAL__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

namespace cicada
{
  namespace feature
  {
    class VariationalImpl;

    class Variational : public FeatureFunction
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;

      typedef sentence_type ngram_type;
      
    private:
      typedef FeatureFunction base_type;
      typedef VariationalImpl impl_type;
      
    public:
      Variational(const std::string& parameter);
      ~Variational();
      
    private:
      Variational() {}
      Variational(const Variational&) {}
      Variational& operator=(const Variational&) { return *this; }
      
    public:
      virtual void operator()(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates) const;
      virtual void operator()(const state_ptr_type& state,
			      feature_set_type& features) const;
      
      void clear();

      template <typename Iterator>
      void insert(Iterator first, Iterator last, const double& logprob)
      {
	insert(ngram_type(first, last), logprob);
      }
      
      void insert(const ngram_type& ngram, const double& logprob);
      
    private:
      impl_type* pimpl;
    };
    
  };
};


#endif

