// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__NUMBER__HPP__
#define __CICADA__FORMAT__NUMBER__HPP__ 1

#include <cicada/format.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{
  namespace format
  {
    class NumberImpl;
    
    class Number : public cicada::Format
    {
    private:
      typedef NumberImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> > pimpl_set_type;

    public:
      typedef boost::filesystem::path path_type;
      
    public:
      Number(const std::string& locale_str_source,
	     const std::string& locale_str_target)
      { initialize(path_type(), path_type(), locale_str_source, locale_str_target); }
      Number(const path_type& path_source,
	     const path_type& path_target,
	     const std::string& locale_str_source,
	     const std::string& locale_str_target)
      { initialize(path_source, path_target, locale_str_source, locale_str_target); }
      
      Number(const Number& x);
      Number& operator=(const Number& x);
      ~Number();
      
    public:
       virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const;

    private:
      void initialize(const path_type& path_source,
		      const path_type& path_target,
		      const std::string& locale_str_source,
		      const std::string& locale_str_target);
      
    private:
      pimpl_set_type pimpls;
    };
  };
};

#endif
