// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PERMUTE__HPP__
#define __CICADA__OPERATION__PERMUTE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

#include <utils/dense_hash_set.hpp>

namespace cicada
{
  namespace operation
  {

    class Permute : public Operation
    {
      typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > exclude_set_type;
      typedef exclude_set_type deterministic_set_type;

      struct Filter
      {
	Filter(const exclude_set_type& __excludes)
	  : excludes(__excludes) {}
    
	const exclude_set_type& excludes;
    
	template <typename Cat>
	bool operator()(const Cat& x) const
	{
	  return ! excludes.empty() && excludes.find(x) != excludes.end();
	}
      };

      struct FilterDeterministic
      {
	FilterDeterministic(const deterministic_set_type& __deterministics)
	  : deterministics(__deterministics) {}
	
	const deterministic_set_type& deterministics;
	
	template <typename Cat>
	bool operator()(const Cat& x) const
	{
	  return deterministics.find(x) != deterministics.end();
	}	
      };

    public:
      Permute(const std::string& parameter, const int __debug);
      
      void operator()(data_type& data) const;
      
      exclude_set_type       excludes;
      deterministic_set_type deterministics;
      int size;
      
      int debug;
    };

  };
};


#endif
