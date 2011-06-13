// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FORMAT__NUMBER__HPP__
#define __CICADA__FORMAT__NUMBER__HPP__ 1

#include <cicada/format.hpp>

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
      Number(const std::string& locale_str_source,
	     const std::string& locale_str_target);
      Number(const Number& x);
      Number& operator=(const Number& x);
      ~Number();
      
    public:
      virtual const phrase_set_type& operator()(const phrase_type& phrase) const;
      
    private:
      pimpl_set_type pimpls;
    };
  };
};

#endif
