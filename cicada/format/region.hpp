// -*- mode: c++ -*-
//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__REGION__HPP__
#define __CICADA__FORMAT__REGION__HPP__ 1

#include <cicada/format.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{
  namespace format
  {
    class RegionImpl;
    
    class Region : public cicada::Format
    {
    private:
      typedef RegionImpl impl_type;
      
    public:
      Region(const std::string& locale_str_source,
	     const std::string& locale_str_target);
      ~Region();

    private:      
      Region(const Region& x) {}
      Region& operator=(const Region& x) { return *this; }
      
    public:
      virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const;
      
    private:
      impl_type* pimpl;
    };
  };
};

#endif
