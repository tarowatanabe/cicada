// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__TABLE_COUNT__HPP__
#define __UTILS__TABLE_COUNT__HPP__ 1

#include <numeric>
#include <limits>
#include <cmath>
#include <stdexcept>

#include <boost/swap.hpp>
#include <boost/functional/hash/hash.hpp>

#include <utils/compact_map.hpp>

namespace utils
{

  namespace detail
  {
    template <typename Count, typename Alloc=std::allocator<Count> >
    class table_count_histogram
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      typedef Count     count_type;
      
      struct unassigned { count_type operator()() const { return count_type(-1); } };
      struct erased { count_type operator()() const { return count_type(-2); } };
      
      typedef typename Alloc::template rebind<std::pair<const count_type, count_type> >::other allocator_type;
      
      typedef utils::compact_map<count_type, count_type,
				 unassigned, erased,
				 boost::hash<count_type>,
				 std::equal_to<count_type>,
				 allocator_type> count_set_type;
      
    public:
      typedef typename count_set_type::const_iterator const_iterator;
      typedef typename count_set_type::iterator       iterator;
      
      void increment(const count_type& count)
      {
	++ counts_[count];
      }
      
      void decrement(const count_type& count)
      {
	typename count_set_type::iterator citer = counts_.find(count);
	if (citer == counts_.end())
	  throw std::runtime_error("invalid decrement");
	
	-- citer->second;
	
	if (! citer->second)
	  counts_.erase(citer);
      }
      
      bool empty() const { return counts_.empty(); }
      size_type size() const { return counts_.size(); }
      
      inline const_iterator begin() const { return counts_.begin(); }
      inline       iterator begin()       { return counts_.begin(); }
      inline const_iterator end() const { return counts_.end(); }
      inline       iterator end()       { return counts_.end(); }

      void clear()
      {
	counts_.clear();
      }

      void swap(table_count_histogram& x)
      {
	counts_.swap(x.counts_);
      }
      
    private:
      count_set_type counts_;
    };    
  };
};

namespace std
{
  template <typename C, typename A>
  inline
  void swap(utils::detail::table_count_histogram<C,A>& x,
	    utils::detail::table_count_histogram<C,A>& y)
  {
    x.swap(y);
  }
};

namespace utils
{
  template <typename Count, size_t Floors=1, typename Alloc=std::allocator<Count> >
  class table_count
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef Count     count_type;
    
    typedef detail::table_count_histogram<Count,Alloc> histogram_type;

    typedef const histogram_type* const_iterator;

  public:
    table_count() : customers_(0), tables_(0) {}
    
    const histogram_type& operator[](size_type floor) const { return floors_[floor]; }

    const_iterator begin() const { return floors_; }
    const_iterator end() const { return floors_ + Floors; }
    
    count_type customers() const { return customers_; }
    count_type tables() const { return tables_; }

    void clear()
    {
      customers_ = 0;
      tables_ = 0;
      
      for (size_type floor = 0; floor != Floors; ++ floor)
	floors_[floor].clear();
    }

    std::pair<size_type, bool> increment_new(const size_type floor=0)
    {
      floors_[floor].increment(1);
      ++ customers_;
      ++ tables_;
      
      return std::make_pair(floor, true);
    }
    
    template <typename Sampler>
    std::pair<size_type, bool> increment_existing(const double discount, Sampler& sampler)
    {
      double r = sampler.uniform() * (customers_ - discount * tables_);
      
      for (size_type floor = 0; floor != Floors; ++ floor) {
	typename histogram_type::iterator hiter_end = floors_[floor].end();
	typename histogram_type::iterator hiter     = floors_[floor].begin();
	
	for (/**/; hiter != hiter_end; ++ hiter) {
	  const double step = (hiter->first - discount) * hiter->second;
	  
	  if (step > r) {
	    const count_type count_customer = hiter->first;
	    
	    floors_[floor].decrement(count_customer);
	    floors_[floor].increment(count_customer + 1);
	    
	    ++ customers_;
	    
	    return std::make_pair(floor, true);
	  }
	  
	  r -= step;
	}
      }
      
      throw std::runtime_error("invalid increment for existing table");
    }
    
    template <typename Sampler>
    std::pair<size_type, bool> decrement(Sampler& sampler)
    {
      difference_type r = sampler.uniform() * customers_;
      
      for (size_type floor = 0; floor != Floors; ++ floor) {
	typename histogram_type::iterator hiter_end = floors_[floor].end();
	typename histogram_type::iterator hiter     = floors_[floor].begin();
	
	for (/**/; hiter != hiter_end; ++ hiter) {
	  const difference_type step = hiter->first * hiter->second;

	  if (step > r) {
	    -- customers_;
	    
	    const count_type count_customer = hiter->first;
	    
	    if (count_customer == 1) {
	      floors_[floor].decrement(count_customer);
	      -- tables_;
	    } else {
	      floors_[floor].decrement(count_customer);
	      floors_[floor].increment(count_customer - 1);
	    }
	    
	    return std::make_pair(floor, count_customer == 1);
	  }
	  
	  r -= step;
	}
      }
      
      throw std::runtime_error("invalid decrement");
    }
    

    void swap(table_count& x, table_count& y)
    {
      std::swap(customers_, x.customers_);
      std::swap(tables_, x.tables_);
      
      for (size_type floor = 0; floor != Floors; ++ floor)
	floors_[floor].swap(x.floors_[floor]);
    }
    
  private:
    count_type     customers_;
    count_type     tables_;
    histogram_type floors_[Floors];
  };
};

namespace std
{
  template <typename C, size_t F, typename A>
  inline
  void swap(utils::table_count<C,F,A>& x, utils::table_count<C,F,A>& y)
  {
    x.swap(y);
  }
};


#endif
