// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__VARIATIONAL__HPP__
#define __CICADA__FEATURE__VARIATIONAL__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/weight_vector.hpp>

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

      typedef cicada::WeightVector<double> weight_set_type;
      
    private:
      typedef FeatureFunction base_type;
      typedef VariationalImpl impl_type;
      
    public:
      Variational(const std::string& parameter);
      Variational(const Variational&);
      
      ~Variational();

      Variational& operator=(const Variational&);
      
    private:
      Variational() {}
      
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
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new Variational(*this)); }
      
      void clear();
      
      void insert(const hypergraph_type& graph, const weight_set_type& weights);
      
    private:
      impl_type* pimpl;
    };
    
  };
};


#endif

