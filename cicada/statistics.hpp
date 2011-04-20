// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__STATISTICS__HPP__
#define __CICADA__STATISTICS__HPP__ 1

#include <cstddef>
#include <stdint.h>

#include <iostream>

#include <cicada/attribute.hpp>

#include <boost/functional/hash/hash.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  class Statistics
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef int64_t count_type;
    typedef double  second_type;
    
    typedef Attribute attribute_type;
    
    struct Stat
    {
      count_type  count;
      
      count_type node;
      count_type edge;
      
      second_type user_time;
      second_type cpu_time;
      
      Stat() : count(0), node(0), edge(0), user_time(0), cpu_time(0) {}
      Stat(const count_type& __count,
	   const count_type& __node,
	   const count_type& __edge,
	   const second_type& __user_time,
	   const second_type& __cpu_time)
	: count(__count), node(__node), edge(__edge),
	  user_time(__user_time), cpu_time(__cpu_time) {}
      
      void clear()
      {
	count = 0;
	node = 0;
	edge = 0;
	user_time = 0;
	cpu_time  = 0;
      }
      
      Stat operator+() const
      {
	return *this;
      }
      
      Stat operator-() const
      {
	Stat stat;
	stat -= *this;
	return stat;
      }
      
      Stat& operator+=(const Stat& x)
      {
	count += x.count;
	node += x.node;
	edge += x.edge;
	user_time += x.user_time;
	cpu_time  += x.cpu_time;
	return *this;
      }
      
      Stat& operator-=(const Stat& x)
      {
	count -= x.count;
	node -= x.node;
	edge -= x.edge;
	user_time -= x.user_time;
	cpu_time  -= x.cpu_time;
	return *this;
      }

    public:
      friend
      std::ostream& operator<<(std::ostream& os, const Stat& stat);
    };
    
    typedef Stat stat_type;
    typedef Stat statistic_type;
    
  private:
    typedef google::dense_hash_map<attribute_type, stat_type, boost::hash<attribute_type>, std::equal_to<attribute_type> > stat_set_type;
    
  public:
    typedef stat_set_type::value_type     value_type;

    typedef stat_set_type::const_iterator const_iterator;
    typedef stat_set_type::iterator       iterator;
    
  public:
    Statistics() : stats() { stats.set_empty_key(attribute_type(attribute_type::id_type(-1))); }
    
    inline const_iterator begin() const { return stats.begin(); }
    inline       iterator begin()       { return stats.begin(); }
    
    inline const_iterator end() const { return stats.end(); }
    inline       iterator end()       { return stats.end(); }
    
    inline const_iterator find(const attribute_type& x) const { return stats.find(x); }
    inline       iterator find(const attribute_type& x)       { return stats.find(x); }
    
    std::pair<iterator, bool> insert(const value_type& x) { return stats.insert(x); }
    
    stat_type& operator[](const attribute_type& x) { return stats[x]; }

    void clear() { stats.clear(); }
    
    bool empty() const { return stats.empty(); }
    size_type size() const { return stats.size(); }
    
    Statistics& operator+=(const Statistics& x)
    {
      const_iterator iter_end = x.end();
      for (const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	operator[](iter->first) += iter->second;
      
      return *this;
    }
    
    Statistics& operator-=(const Statistics& x)
    {
      const_iterator iter_end = x.end();
      for (const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	operator[](iter->first) -= iter->second;
      
      return *this;
    }

  public:
    friend
    std::ostream& operator<<(std::ostream& os, const Statistics& stats);
    
  private:
    stat_set_type stats;
  };

  
  inline
  Statistics::stat_type operator+(const Statistics::stat_type& x, const Statistics::stat_type& y)
  {
    Statistics::stat_type stat(x);
    stat += y;
    return stat;
  }
  
  inline
  Statistics::stat_type operator-(const Statistics::stat_type& x, const Statistics::stat_type& y)
  {
    Statistics::stat_type stat(x);
    stat -= y;
    return stat;
  }

  inline
  Statistics operator+(const Statistics& x, const Statistics& y)
  {
    Statistics stats(x);
    stats += y;
    return stats;
  }

  inline
  Statistics operator-(const Statistics& x, const Statistics& y)
  {
    Statistics stats(x);
    stats -= y;
    return stats;
  }
  
};

#endif
