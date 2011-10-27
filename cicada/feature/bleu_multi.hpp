// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__BLEU_MULTI__HPP__
#define __CICADA__FEATURE__BLEU_MULTI__HPP__ 1

#include <string>
#include <vector>

#include <cicada/feature_function.hpp>
#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <cicada/eval.hpp>

namespace cicada
{
  namespace feature
  {
    class BleuMulti : public FeatureFunction
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
      typedef FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_function_ptr_type feature_function_ptr_type;
      typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type> > feature_function_ptr_set_type;
      
    public:
      BleuMulti(const std::string& parameter);
      BleuMulti(const BleuMulti&);
      BleuMulti& operator=(const BleuMulti&);
      
    private:
      BleuMulti() {}
      
    public:
      virtual void apply(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 const bool final) const;
      virtual void apply_coarse(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const;
      virtual void apply_predict(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const;
      virtual void apply_scan(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      const int dot,
			      feature_set_type& features,
			      const bool final) const;
      virtual void apply_complete(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  const bool final) const;
      virtual void assign(const size_type& id,
			  const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts);
      
      virtual void initialize();
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new BleuMulti(*this)); }
      
    private:
      feature_function_ptr_set_type bleus;
    };
    
  };
};


#endif

