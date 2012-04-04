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
    
    typedef Generator generator_type;

    sampler()
      : gen()
    {
      gen.seed(utils::random_seed());
    }
    
    explicit sampler(const unsigned int seed)
      : gen()
    {
      gen.seed(seed);
    }
    
  public:
    generator_type& generator() { return gen; }
    
    double operator()() { return uniform(); }
    
    double uniform()
    {
      return boost::uniform_01<double>()(gen);
    }

    double uniform(const double& min, const double& max)
    {
      return boost::random::uniform_real_distribution<double>(min, max)(gen);
    }
    
    double normal(const double& mean, const double& var)
    {
      return boost::normal_distribution<double>(mean, var)(gen);
    }
    
    int poisson(const double& lambda)
    {
      return boost::poisson_distribution<int, double>(lambda)(gen);
    }
    
    bool bernoulli(const double& p)
    {
      return boost::bernoulli_distribution<double>(p)(gen);
    }

    int binomial(const int t, const double& p)
    {
      return boost::binomial_distribution<double>(t, p)(gen);
    }
    
    int geometric(const double& p)
    {
      return boost::geometric_distribution<int, double>(p)(gen);
    }
    
    double exponential(const double& lambda)
    {
      // we will use PRML style
      
      return - std::log(1.0 - uniform()) / lambda;
      //return boost::exponential_distribution<double>(lambda)(gen);
    }

    double gamma(const double& alpha, const double& beta)
    {
      // in PRML, beta is inverse! I decided to use follow PRML style gamma distribution, not
      // boost::math style beta...
      return boost::gamma_distribution<double>(alpha, 1.0 / beta)(gen);
      
#if 0
      double b, c, u, v, w, y, x, z;
      
      if (a > 1) { // Best's rejection method. Devroye (1986) p.410
        b = a - 1;
        c = 3 * a - 0.75;
	
        bool accept = false;
        while (! accept) {
	  u = uniform();
	  v = uniform();
	  w = u * (1 - u);
	  y = std::sqrt(c / w) * (u - 0.5);
	  x = b + y;
	  
	  if (x >= 0) {
	    z = 64.0 * w * w * w * v * v;
	    accept = (z <= 1.0 - (2.0 * y * y) / x || std::log(z) <= 2.0 * (b * std::log(x / b) - y));
	  }
        } 
      } else { // Johnk's method. Devroye (1986) p.418
        do {
	  x = std::pow(uniform(), 1 / a);
	  y = std::pow(uniform(), 1 / (1 - a));
        } while (x + y > 1);
	
        x = exponential(1.0) * x / (x + y);
      }
      return x / scale;
#endif
    }
    
    double beta(const double& a, const double& b)
    {
      const double x = gamma(a, 1);
      const double y = gamma(b, 1);
      
      return x / (x + y);
    }
    
    bool draw(double a, double b, const double temperature=1.0)
    {
      using namespace std;

      if (temperature == 1.0)
	return bernoulli(a / (a + b));
      else {
	a = pow(a, 1.0 / temperature);
	b = pow(b, 1.0 / temperature);
	return bernoulli(a / (a + b));
      }
    }

    template <typename Iterator>
    Iterator draw(Iterator first, Iterator last, const double temperature=1.0)
    {
      using namespace std;

      if (std::distance(first, last) <= 1) return first;
      
      if (temperature == 1.0) {
	const double draw = uniform() * std::accumulate(first, last, 0.0);
	
	double sum = *first;
	++ first;
	for (/**/; first != last && sum < draw; ++ first)
	  sum += *first;
	
	return -- first;
      } else {
	const double anneal = 1.0 / temperature;
	
	double scale = 0.0;
	for (Iterator iter = first; iter != last; ++ iter)
	  scale += pow(*iter, anneal);
	
	const double draw = uniform() * scale;
	
	double sum = *first;
	++ first;
	for (/**/; first != last && sum < draw; ++ first)
	  sum += pow(*first, anneal);
	
	return -- first;
      }
    }
    
    
  private:
    Generator gen;
  };
  
  
};

#endif
