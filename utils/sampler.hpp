// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__SAMPLER__HPP__
#define __UTILS__SAMPLER__HPP__ 1

#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>

#include <boost/random.hpp>

#include <utils/random_seed.hpp>

namespace utils
{

  template <typename Generator>
  struct sampler
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    sampler()
      : dist(0,1), gen(), random(gen, dist)
    {
      gen.seed(utils::random_seed());
    }
    
    explicit sampler(const unsigned int seed)
      : dist(0,1), gen(), random(gen, dist)
    {
      gen.seed(seed);
    }
    
  public:
    Generator& generator() { return gen; }
    
    // draw from U(0,1)
    double operator()() { return random(); }
    
    // draw from U(0,1)
    double uniform() { return random(); }
    
    // draw from N(mean,var)
    double normal(const double& mean, const double& var) { return boost::normal_distribution<double>(mean, var)(random); }
    
    // draw from Poisson
    double poisson(const int lambda) { return boost::poisson_distribution<int>(lambda)(random); }
    
    bool accept_metropolis_hasting(const double& p_cur,
				   const double& p_prev,
				   const double& q_cur,
				   const double& q_prev)
    {
      const double a = (p_cur / p_prev) * (q_prev / q_cur);
      
      return (std::log(a) >= 0.0) || (uniform() < a);
    }
    
    size_type select(const double& a, const double& b)
    {
      return uniform() > (a / (a + b));
    }
    
  private:
    boost::uniform_real<> dist;
    Generator gen;
    boost::variate_generator<Generator&, boost::uniform_real<> > random;
  };
  
  
};

#endif
