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

    bool bernoulli_sample(const double& p)
    {
      return random() < p;
    }
    
    double expon_sample(const double& l)
    {
      return - std::log(1.0 - random()) / l;
    }

    double gamma_sample(const double& a, const double& scale)
    {
      double b, c, e, u, v, w, y, x, z;
      
      if (a > 1) { // Best's XG method
        b = a - 1;
        c = 3 * a - 0.75;
	
        bool accept = false;
        do {
	  u = random();
	  v = random();
	  w = u * (1 - u);
	  y = std::sqrt(c / w) * (u-  0.5);
	  x = b + y;
	  
	  if(x >= 0) {
	    z = 64 * w * w * w * v * v;
	    accept = (z <= 1 - 2 * y * y / x || std::log(z) <= 2 * (b * std::log(x / b) - y));
	  }
        } while (!accept);
      } else { // Johnk's method
        do {
	  x = std::pow(random(), 1 / a);
	  y = std::pow(random(), 1 / (1 - a));
        } while (x + y > 1);
	
        x = expon_sample(1.0) * x / (x + y);
      }
      return x * scale;
    }
    
    double beta_sample(const double& a, const double& b)
    {
      const double ga = gamma_sample(a, 1);
      const double gb = gamma_sample(b, 1);
      
      return ga / (ga + gb);
    }
    
  private:
    boost::uniform_real<> dist;
    Generator gen;
    boost::variate_generator<Generator&, boost::uniform_real<> > random;
  };
  
  
};

#endif
