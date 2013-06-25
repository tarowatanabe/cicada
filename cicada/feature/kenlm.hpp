// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__KENLM__HPP__
#define __CICADA__FEATURE__KENLM__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>

#include <boost/filesystem.hpp>

namespace cicada
{
  namespace feature
  {
    template <typename Model>
    class KenLMImpl;

    template <typename Model>
    class KenLM : public FeatureFunction
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
    public:      
      typedef boost::filesystem::path path_type;
      
    private:
      typedef FeatureFunction  base_type;
      typedef KenLMImpl<Model> impl_type;
      
    public:
      // parameter = key:[key=value (delimited by ',')]*
      
      // kenlm parameter = kenlm:file=file-name,name=feature-name,order=5
      // "kenlm" is the key for this kenlm-feature
      // file: file name
      // name: name of this feature function. default to kenlm
      // order: kenlm's order
      
      KenLM(const std::string& parameter);
      KenLM(const KenLM&);
      ~KenLM();
      
      KenLM& operator=(const KenLM&);

    private:
      KenLM() {}
      
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

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new KenLM(*this)); }
      
    private:
      
      impl_type* pimpl;
      impl_type* pimpl_coarse;
    };
    
    struct KenLMFactory
    {
      FeatureFunction::feature_function_ptr_type create(const std::string& parameter) const;
    };
  };
};


#endif
