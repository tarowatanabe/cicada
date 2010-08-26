// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__NGRAM__HPP__
#define __CICADA__FEATURE__NGRAM__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/ngram.hpp>

#include <boost/filesystem.hpp>

namespace cicada
{
  namespace feature
  {
    
    class NGramImpl;

    class NGram : public FeatureFunction
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
      typedef NGramImpl       impl_type;
      
    public:
      // parameter = key:[key=value (delimited by ',')]*
      
      // ngram parameter = ngram:file=file-name,name=feature-name,order=5
      // "ngram" is the key for this ngram-feature
      // file: file name
      // name: name of this feature function. default to ngram
      // order: ngram's order
      
      NGram(const std::string& parameter);
      NGram(const NGram&);
      ~NGram();
      
      NGram& operator=(const NGram&);

    private:
      NGram() {}
      
    public:
      virtual void operator()(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const;
      virtual void operator()(const edge_type& edge,
			      feature_set_type& features) const;
      

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new NGram(*this)); }
      
    private:
      
      impl_type* pimpl;
    };
    
  };
};


#endif
