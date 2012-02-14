// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SEMIRING__TUPLE__HPP__
#define __CICADA__SEMIRING__TUPLE__HPP__ 1

#include <limits>
#include <algorithm>

#include <cicada/semiring/traits.hpp>

#include <utils/simple_vector.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  namespace semiring
  {
    template <typename Tp, typename Alloc=std::allocator<Tp> >
    class Tuple
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef Tp value_type;
      typedef utils::simple_vector<Tp, Alloc>  tuple_type;
      
    private:
      typedef Tuple<Tp, Alloc> self_type;

    public:
      typedef typename tuple_type::const_reference const_reference;
      typedef typename tuple_type::reference       reference;
      
    public:
      static inline self_type zero() { return self_type(); }
      
    public:
      Tuple() : __tuple() {}
      Tuple(size_type size) : __tuple(size) {}
      Tuple(size_type size, const Tp& x) : __tuple(size, x) {}
      
    public:
      void clear() { __tuple.clear(); }
      void resize(size_type __size) { __tuple.resize(__size); }
      
      bool empty() const { return __tuple.empty(); }
      size_type size() const { return __tuple.size(); }
      
    public:
      inline const_reference operator[](size_type pos) const { return __tuple[pos]; }
      inline       reference operator[](size_type pos)       { return __tuple[pos]; }
      
    public:
      template <typename T, typename A>
      Tuple<Tp,Alloc>& operator+=(const Tuple<T,A>& x)
      {
	const size_type __size = utils::bithack::max(__tuple.size(), x.__tuple.size());
	
	__tuple.resize(__size);
	std::transform(x.__tuple.begin(), x.__tuple.end(), __tuple.begin(), __tuple.begin(), std::plus<Tp>());
	
	return *this;
      }
      
      template <typename T, typename A>
      Tuple<Tp,Alloc>& operator*=(const Tuple<T,A>& x)
      {
	const size_type __size = utils::bithack::min(__tuple.size(), x.__tuple.size());
	
	std::fill(std::transform(x.__tuple.begin(), x.__tuple.begin() + __size, __tuple.begin(), __tuple.begin(), std::multiplies<Tp>()), __tuple.end(), Tp());
	
	return *this;
      }
      
      template <typename T>
      Tuple<Tp,Alloc>& operator*=(const T& x)
      {
	std::transform(__tuple.begin(), __tuple.end(), __tuple.begin(), std::bind2nd(std::multiplies<Tp>(), x));
	return *this;
      }
      
    private:
      tuple_type __tuple;
    };
    
    
    template <typename Tp, typename Alloc, typename T, typename A>
    inline
    Tuple<Tp, Alloc> operator+(const Tuple<Tp, Alloc>& x, const Tuple<T, A>& y)
    {
      Tuple<Tp, Alloc> __value(x);
      __value += y;
      return __value;
    }
    
    template <typename Tp, typename Alloc, typename T, typename A>
    inline
    Tuple<Tp, Alloc> operator*(const Tuple<Tp, Alloc>& x, const Tuple<T, A>& y)
    {
      Tuple<Tp, Alloc> __value(x);
      __value *= y;
      return __value;
    }

    template <typename Tp, typename Alloc, typename T>
    inline
    Tuple<Tp, Alloc> operator*(const Tuple<Tp, Alloc>& x, const T& y)
    {
      Tuple<Tp, Alloc> __value(x);
      __value *= y;
      return __value;
    }


    template <typename Tp, typename Alloc>
    struct traits<Tuple<Tp, Alloc> >
    {
      static inline Tuple<Tp,Alloc> zero() { return Tuple<Tp,Alloc>::zero();  }
    };

  };
};

#endif
