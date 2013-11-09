// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__NGRAM_NN__HPP__
#define __CICADA__FEATURE__NGRAM_NN__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>

#include <boost/filesystem.hpp>

namespace cicada
{
  namespace feature
  {
    
    class NGramNNImpl;

    class NGramNN : public FeatureFunction
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
    public:      
      typedef boost::filesystem::path path_type;
      
    private:
      typedef FeatureFunction base_type;
      typedef NGramNNImpl     impl_type;
      
    public:
      // parameter = key:[key=value (delimited by ',')]*
      
      // ngram parameter = ngram:file=file-name,name=feature-name,order=5
      // "ngram" is the key for this ngram-feature
      // file: file name
      // name: name of this feature function. default to ngram
      // order: ngram's order
      
      NGramNN(const std::string& parameter);
      NGramNN(const NGramNN&);
      ~NGramNN();
      
      NGramNN& operator=(const NGramNN&);

    private:
      NGramNN() {}
      
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
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new NGramNN(*this)); }
      
    private:
      
      impl_type* pimpl;
      impl_type* pimpl_coarse;
    };
    
  };
};


#endif
