// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__ALIGNMENT__HPP__
#define __CICADA__FEATURE__ALIGNMENT__HPP__ 1

#include <cicada/cluster_stemmer.hpp>
#include <cicada/feature_function.hpp>

namespace cicada
{
  namespace feature
  {
    
    // alignment related feature with states...
    template <typename F>
    class AlignmentBaseState : public FeatureFunction
    {
    public:
      AlignmentBaseState(const std::string& parameter)
	: f(parameter, FeatureFunction::__state_size, FeatureFunction::__feature_name)
      {
	if (__state_size == 0)
	  throw std::runtime_error("Are we stateless feature...?");
      }
      
      virtual void apply(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	f(*this, state, states, edge, features, estimates, final);
      }
      virtual void apply_coarse(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const {}
      virtual void apply_predict(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const {}
      virtual void apply_scan(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      const int dot,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const {}
      virtual void apply_complete(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
      {
	f(*this, state, states, edge, features, estimates, final);
      }

      virtual void assign(const size_type& id,
			  const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts)
      {
	f(*this, id, hypergraph, lattice, spans, targets, ngram_counts);
      }
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new AlignmentBaseState<F>(*this)); }
      
    private:
      F f;
    };
    
    // alignment related feature without states...
    
    template <typename F>
    class AlignmentBase : public FeatureFunction
    {
    public:
      AlignmentBase(const std::string& parameter)
	: f(parameter, FeatureFunction::__state_size, FeatureFunction::__feature_name)
      {
	if (__state_size != 0)
	  throw std::runtime_error("we are not stateless feature...?");
      }
      
      virtual void apply(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	f(*this, state, states, edge, features, estimates, final);
      }
      virtual void apply_coarse(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
      {
	f(*this, state, states, edge, features, estimates, final);
      }
      virtual void apply_predict(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
      {
	f(*this, state, states, edge, features, estimates, final);
      }
      virtual void apply_scan(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      const int dot,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const {}
      virtual void apply_complete(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const {  }
      
      virtual void assign(const size_type& id,
			  const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts)
      {
	f(*this, id, hypergraph, lattice, spans, targets, ngram_counts);
      }
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new AlignmentBase<F>(*this)); }
      
    private:
      F f;
    };

    // we will put into align namespace...
    namespace align
    {
      struct Base
      {
	typedef size_t    size_type;
	typedef ptrdiff_t difference_type;
	
	typedef cicada::Symbol     symbol_type;
	typedef cicada::Vocab      vocab_type;
	typedef cicada::Sentence   sentence_type;
	typedef cicada::HyperGraph hypergraph_type;
	typedef cicada::SpanVector span_set_type;
	typedef cicada::Lattice    lattice_type;
	typedef cicada::Rule       rule_type;
	
	typedef cicada::SentenceVector sentence_set_type;
	typedef cicada::NGramCountSet  ngram_count_set_type;
	
	typedef hypergraph_type::feature_set_type   feature_set_type;
	typedef hypergraph_type::attribute_set_type attribute_set_type;

	typedef feature_set_type::feature_type     feature_type;
	typedef attribute_set_type::attribute_type attribute_type;
	
	typedef FeatureFunction feature_function_type;
	
	typedef feature_function_type::state_ptr_type     state_ptr_type;
	typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
	
	typedef feature_function_type::edge_type edge_type;
      };

      class RelativePosition : public Base
      {
      public:
	typedef ClusterStemmer normalizer_type;
	typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;
	
	RelativePosition(const std::string& parameter, size_type& __state_size, feature_type& __feature_name);
	
	void operator()(const feature_function_type& feature_function,
			state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const;
	
	void operator()(const feature_function_type& feature_function,
			const size_type& id,
			const hypergraph_type& hypergraph,
			const lattice_type& lattice,
			const span_set_type& spans,
			const sentence_set_type& targets,
			const ngram_count_set_type& ngram_counts);
	
      private:
	normalizer_set_type normalizers;
	const sentence_type* sentence;
	int source_size;
	int target_size;
      };
      
      struct NullJump : public Base
      {
	NullJump(const std::string& parameter, size_type& __state_size, feature_type& __feature_name);
	
	void operator()(const feature_function_type& feature_function,
			state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const;
	
	void operator()(const feature_function_type& feature_function,
			const size_type& id,
			const hypergraph_type& hypergraph,
			const lattice_type& lattice,
			const span_set_type& spans,
			const sentence_set_type& targets,
			const ngram_count_set_type& ngram_counts) {}

	feature_type feature_none_none;
	feature_type feature_none_word;
	feature_type feature_word_none;
	feature_type feature_word_word;
      };
        
      class TargetBigram : public Base
      {
      public:
	typedef ClusterStemmer normalizer_type;
	typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;

	TargetBigram(const std::string& parameter, size_type& __state_size, feature_type& __feature_name);
	
	void operator()(const feature_function_type& feature_function,
			state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const;
	
	void operator()(const feature_function_type& feature_function,
			const size_type& id,
			const hypergraph_type& hypergraph,
			const lattice_type& lattice,
			const span_set_type& spans,
			const sentence_set_type& targets,
			const ngram_count_set_type& ngram_counts);
	
      private:
	normalizer_set_type normalizers;
	const sentence_type* sentence;
      };
    };

    // actual classes...
    typedef AlignmentBase<align::RelativePosition>  RelativePosition;
    typedef AlignmentBaseState<align::NullJump>     NullJump;
    typedef AlignmentBaseState<align::TargetBigram> TargetBigram;
    
  };
};


#endif
