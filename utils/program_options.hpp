// -*- mode: c++ -*-

#ifndef __UTILS__PROGRAM_OPTIONS__HPP__
#define __UTILS__PROGRAM_OPTIONS__HPP__ 1

// additional definition for program_options...

#include <memory>

#include <boost/program_options.hpp>

namespace utils
{
  
  boost::program_options::typed_value<bool>* true_false_switch(bool* value)
  {
    typedef boost::program_options::typed_value<bool> value_type;
    
    std::auto_ptr<value_type> ret(new value_type(value));
    if (value)
      ret->default_value(*value, *value ? "true" : "false");
    else
      ret->default_value(false, "false");
    ret->implicit_value(true, "true");
    
    return ret.release();
  }

  boost::program_options::typed_value<bool>* true_false_switch()
  {
    return true_false_switch(0);
  }
    
};


#endif
