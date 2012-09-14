// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__DATE__HPP__
#define __CICADA__FORMAT__DATE__HPP__ 1

#include <cicada/format.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{
  namespace format
  {
    class DateImpl;
    
    class Date : public cicada::Format
    {
    private:
      typedef DateImpl impl_type;
      typedef std::vector<impl_type*, std::allocator<impl_type*> > pimpl_set_type;
      
    public:
      typedef boost::filesystem::path path_type;

    public:
      Date(const std::string& locale_str_source,
	   const std::string& locale_str_target)
      { initialize(path_type(), path_type(), locale_str_source, locale_str_target); }
      Date(const path_type& path_source,
	   const path_type& path_target,
	   const std::string& locale_str_source,
	   const std::string& locale_str_target)
      { initialize(path_source, path_target, locale_str_source, locale_str_target); }
      ~Date();

    private:
      Date(const Date& x) {}
      Date& operator=(const Date& x) { return *this; }
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
