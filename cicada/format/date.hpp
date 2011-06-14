// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__DATE__HPP__
#define __CICADA__FORMAT__DATE__HPP__ 1

#include <cicada/format.hpp>

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
      Date(const std::string& locale_str_source,
	   const std::string& locale_str_target);
      Date(const Date& x);
      Date& operator=(const Date& x);
      ~Date();
      
    public:
       virtual void operator()(const phrase_type& phrase, phrase_set_type& phrases) const;
      
    private:
      pimpl_set_type pimpls;
    };
  };
};

#endif
