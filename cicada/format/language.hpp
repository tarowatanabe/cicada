// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__LANGUAGE__HPP__
#define __CICADA__FORMAT__LANGUAGE__HPP__ 1

#include <cicada/format.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{
  namespace format
  {
    class LanguageImpl;
    
    class Language : public cicada::Format
    {
    private:
      typedef LanguageImpl impl_type;
      
    public:
      Language(const std::string& locale_str_source,
	       const std::string& locale_str_target);
      ~Language();

    private:      
      Language(const Language& x) {}
      Language& operator=(const Language& x) { return *this; }
      
    public:
      virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const;
      
    private:
      impl_type* pimpl;
    };
  };
};

#endif
