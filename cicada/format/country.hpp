// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__COUNTRY__HPP__
#define __CICADA__FORMAT__COUNTRY__HPP__ 1

#include <cicada/format.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{
  namespace format
  {
    class CountryImpl;
    
    class Country : public cicada::Format
    {
    private:
      typedef CountryImpl impl_type;
      
    public:
      Country(const std::string& locale_str_source,
	      const std::string& locale_str_target);
      ~Country();

    private:      
      Country(const Country& x) {}
      Country& operator=(const Country& x) { return *this; }
      
    public:
      virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const;
      
    private:
      impl_type* pimpl;
    };
  };
};

#endif
