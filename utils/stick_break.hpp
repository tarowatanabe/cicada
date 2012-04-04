// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// a stick breaking representation of PYP
//

#ifndef __UTILS__STICK_BREAK__HPP__
#define __UTILS__STICK_BREAK__HPP__ 1

#include <stdexcept>
#include <vector>
#include <numeric>

#include <boost/lexical_cast.hpp>

namespace utils
{
  
  template <typename Alloc=std::allocator<double> >
  class stick_break
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef double stick_type;

  private:
    typedef typename Alloc::template rebind<stick_type >::other stick_alloc_type;
    typedef std::vector<stick_type, stick_alloc_type>           stick_set_type;
    typedef stick_break<Alloc> self_type;

  public:
    stick_break()
      : sticks(1, 1.0),
	m_discount(0.9),
	m_strength(1.0) {}
    
    stick_break(const double& __discount,
		const double& __strength)
      : sticks(1, 1.0),
	m_discount(__discount),
	m_strength(__strength)
    {
      verify_parameters();
    }
    
  public:    
    double& discount() { return m_discount; }
    double& strength() { return m_strength; }
    
    const double& discount() const { return m_discount; }
    const double& strength() const { return m_strength; }
    
    bool empty() const { return sticks.size() == 1; }
    size_type size() const { return sticks.size() - 1; }

    const double& operator[](size_type pos) const
    {
      return sticks[pos];
    }
    
    template <typename Sampler>
    void increment(Sampler& sampler)
    {
      const double stick = sticks.back();
      const double beta = sampler.beta(1.0 - m_discount, m_strength + m_discount * sticks.size());
      
      sticks.back() = stick * beta;
      sticks.push_back(stick * (1.0 - beta));
    }
    
    template <typename Sampler>
    void sample_parameters(size_type size, Sampler& sampler)
    {
      sticks.clear();
      sticks.push_back(1.0);
      
      for (size_t i = 0; i != size; ++ i) {
	const double stick = sticks.back();
	const double beta = sampler.beta(1.0 - m_discount, m_strength + m_discount * sticks.size());
	
	sticks.back() = stick * beta;
	sticks.push_back(stick * (1.0 - beta));
      }
    }

    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      sticks.clear();
      sticks.insert(sticks.end(), first, last);

      if (sticks.empty()) 
	sticks.push_back(1.0);
      else {
	const double sum = std::accumulate(sticks.begin(), sticks.end(), 0.0);
	
	if (sum >= 1.0 || sum < 0.0)
	  throw std::runtime_error("invalid sticks");
	
	sticks.push_back(1.0 - sum);
      }
    }

    template <typename Mapping>
    void permute(const Mapping& mapping)
    {
      stick_set_type sticks_new(sticks.size());
      
      for (size_type i = 0; i != sticks.size(); ++ i)
	sticks_new[i] = sticks[mapping[i]];
      
      sticks_new.swap(sticks);
    }

    void swap(stick_break& x)
    {
      sticks.swap(x.sticks);
      std::swap(m_discount, x.m_discount);
      std::swap(m_strength, x.m_strength);
    }

    bool verify_parameters()
    {
      if (m_discount < 0.0 || m_discount >= 1.0)
	throw std::runtime_error("invalid discount: " + boost::lexical_cast<std::string>(m_discount));
      
      if (m_strength <= - m_discount)
	throw std::runtime_error("invalid strength: " + boost::lexical_cast<std::string>(m_strength));
      
      return true;
    }
    
  private:
    stick_set_type sticks;

    double m_discount;
    double m_strength;
  };
};

namespace std
{
  template <typename A>
  inline
  void swap(utils::stick_break<A>& x, utils::stick_break<A>& y)
  {
    x.swap(y);
  }
  
};

#endif
