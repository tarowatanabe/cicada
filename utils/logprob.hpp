// -*- mode: c++ -*-

#ifndef __UTILS__LOGPROB__HPP__
#define __UTILS__LOGPROB__HPP__ 1

#include <limits>
#include <iostream>

#include <utils/mathop.hpp>
#include <utils/float_traits.hpp>

#include <boost/numeric/conversion/bounds.hpp>

namespace utils
{
  struct logprob
  {
  public:
    typedef float base_type;
    
  public:    
    struct proxy_type
    {
      friend struct logprob;
      
      proxy_type(const base_type& x) : __log(__round(x)) {}
      
      operator logprob() const { return logprob(*this); }
      
    private:
      base_type __log;
    };
    
    template <typename Tp>
    static inline proxy_type from_log(const Tp& x)
    {
      return proxy_type(x);
    }
    
  public:
    logprob() : __log(__zero()) {}
    explicit logprob(const proxy_type& x) : __log(x.__log) {}
    logprob(const double& x) : __log(x == 1.0 ? __one() : __round(mathop::log(x))) {}
    
  public:
    inline const base_type& log() const { return __log; }
    inline       base_type& log()       { return __log; }

    double prob() const { return mathop::exp(double(__log)); }
    operator double() const { return prob(); }
    operator float() const { return prob(); }
    
  public:
    logprob& operator*=(const logprob& x)
    {
      if (__log == __zero() || x.__log == __zero())
	__log = __zero();
      else
	__log = __log + x.__log;
      return *this; 
    }
    
    logprob& operator/=(const logprob& x)
    {
      if (__log != __zero())
	__log = __log - x.__log;
      return *this;
    }
    
    logprob& operator+=(const logprob& x)
    {
      __log = mathop::logsum(double(__log), double(x.__log));
      return *this;
    }
    
    logprob& operator-=(const logprob& x)
    {
      const double diff = prob() - x.prob();
      if (diff < 0) throw std::runtime_error("logprob do not assume minus value");
      __log = __round(mathop::log(diff));
      return *this;
    }
    
    template <typename Tp>
    logprob& operator^=(const Tp& n)
    {
      __log *= n;
      return *this;
    }
    
    logprob operator*(const logprob& x) const { logprob tmp = *this; return tmp *= x; }
    logprob operator/(const logprob& x) const { logprob tmp = *this; return tmp /= x; }
    logprob operator+(const logprob& x) const { logprob tmp = *this; return tmp += x; }
    logprob operator-(const logprob& x) const { logprob tmp = *this; return tmp -= x; }
    logprob operator^(const logprob& x) const { logprob tmp = *this; return tmp ^= x.prob(); } // it is non-sense...
    
    logprob operator*(const double& x) const { logprob tmp = *this; return tmp *= logprob(x); }
    logprob operator/(const double& x) const { logprob tmp = *this; return tmp /= logprob(x); }
    logprob operator+(const double& x) const { logprob tmp = *this; return tmp += logprob(x); }
    logprob operator-(const double& x) const { logprob tmp = *this; return tmp -= logprob(x); }
    logprob operator^(const double& x) const { logprob tmp = *this; return tmp ^= x; }
    
    friend
    std::istream& operator>>(std::istream& is, logprob& x);
    friend
    std::ostream& operator<<(std::ostream& os, const logprob& x);

    friend
    bool operator==(const logprob& x, const logprob& y);
    friend
    bool operator!=(const logprob& x, const logprob& y);
    friend
    bool operator<(const logprob& x, const logprob& y);
    friend
    bool operator>(const logprob& x, const logprob& y);
    friend
    bool operator<=(const logprob& x, const logprob& y);
    friend
    bool operator>=(const logprob& x, const logprob& y);
    
  private:
    static inline base_type __one() { return base_type(0); }
    static inline base_type __zero() { return boost::numeric::bounds<base_type>::lowest(); }
    static inline base_type __infinity() { return boost::numeric::bounds<base_type>::highest(); }
    static inline base_type __round(const base_type& x)
    {
      return std::max(__zero(), std::min(__infinity(), x));
    }
    
  private:
    base_type __log;
  };

  inline
  std::istream& operator>>(std::istream& is, logprob& x)
  {
    is >> x.__log;
    return is;
  }
  inline
  std::ostream& operator<<(std::ostream& os, const logprob& x)
  {
    os << x.__log;
    return os;
  }
  
  inline
  double operator*(const double& x, const logprob& y)
  {
    return x * y.prob();
  }
  inline
  double operator/(const double& x, const logprob& y)
  {
    return x / y.prob();
  }
  inline
  double operator+(const double& x, const logprob& y)
  {
    return x + y.prob();
  }
  inline
  double operator-(const double& x, const logprob& y)
  {
    return x - y.prob();
  }
  
  inline
  bool operator==(const logprob& x, const logprob& y)
  {
    return x.__log == y.__log;
  }
  inline
  bool operator!=(const logprob& x, const logprob& y)
  {
    return x.__log != y.__log;
  }
  inline
  bool operator<(const logprob& x, const logprob& y)
  {
    return x.__log < y.__log;
  }
  inline
  bool operator>(const logprob& x, const logprob& y)
  {
    return x.__log > y.__log;
  }
  inline
  bool operator<=(const logprob& x, const logprob& y)
  {
    return x.__log <= y.__log;
  }
  inline 
  bool operator>=(const logprob& x, const logprob& y)
  {
    return x.__log >= y.__log;
  }
  
  template <>
  struct float_traits<utils::logprob>
  {
    static inline utils::logprob zero() { return utils::logprob(); } 
    static inline utils::logprob max() { return utils::logprob::from_log(boost::numeric::bounds<float>::highest()); }
    static inline utils::logprob min() { return utils::logprob(); }
  };

};



#endif
