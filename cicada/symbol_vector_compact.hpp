// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__SYMBOL_VECTOR_COMPACT__HPP__
#define __CICADA__SYMBOL_VECTOR_COMPACT__HPP__ 1

#include <cicada/symbol.hpp>
#include <cicada/symbol_vector.hpp>

#include <utils/hashmurmur.hpp>
#include <utils/bithack.hpp>
#include <utils/byte_aligned_code.hpp>
#include <utils/simple_vector.hpp>

namespace cicada
{
  
  class SymbolVectorCompact
  {
  public:
    typedef cicada::Symbol  symbol_type;
    
    typedef size_t      size_type;
    typedef ptrdiff_t   difference_type;
    typedef symbol_type value_type;
    
  private:
    typedef uint8_t         byte_type;
    typedef utils::simple_vector<byte_type, std::allocator<byte_type> > impl_type;
    
    struct codec_type
    {
      typedef uint8_t byte_type;
      
      static size_t encode(byte_type* buffer, const symbol_type::id_type& value)
      {
	return utils::byte_aligned_encode(value, reinterpret_cast<char*>(buffer));
      }
      
      static size_t encode(byte_type* buffer, const symbol_type& value)
      {
	return encode(buffer, value.id());
      }
      
      static size_t decode(const byte_type* buffer, symbol_type::id_type& value)
      {
	return utils::byte_aligned_decode(value, reinterpret_cast<const char*>(buffer));
      }
      
      static size_t decode(const byte_type* buffer, symbol_type& value)
      {
	symbol_type::id_type value_id = 0;
	const size_t ret = utils::byte_aligned_decode(value_id, reinterpret_cast<const char*>(buffer));
	value = symbol_type(value_id);
	return ret;
      }
    };

  public:
    struct iterator
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      typedef std::input_iterator_tag   iterator_category;
      typedef symbol_type value_type;
      typedef const value_type* pointer;
      typedef const value_type& reference;
      typedef const byte_type* ptr_type;
      
    public:
      iterator(ptr_type iter, ptr_type last) : __iter(iter), __last(last), __impl() { increment(); }
      iterator() : __iter(0), __last(0), __impl() {}
      
      const value_type& operator*() const { return __impl; }
      const value_type* operator->() const { return &__impl; }
      
      iterator& operator++()
      {
	increment();
	return *this;
      }
      
      iterator operator++(int)
      {
	iterator tmp = *this;
	increment();
	return tmp;
      }
      
      friend
      bool operator==(const iterator& x, const iterator& y)
      {
	return (x.__iter == y.__iter) && (x.__last == y.__last);
      }
      
      friend
      bool operator!=(const iterator& x, const iterator& y)
      {
	return (x.__iter != y.__iter) || (x.__last != y.__last);
      }

    private:
      void increment()
      {
	if (__iter == __last) {
	  __iter = 0;
	  __last = 0;
	} else
	  std::advance(__iter, codec_type::decode(&(*__iter), __impl));
      }
      
    private:
      ptr_type   __iter;
      ptr_type   __last;
      value_type __impl;
    };

    typedef iterator const_iterator;
    
    typedef const value_type& reference;
    typedef const value_type& const_reference;

  public:
    SymbolVectorCompact() {}
    
    template <typename Iterator>
    SymbolVectorCompact(Iterator first, Iterator last) { assign(first, last); }
    SymbolVectorCompact(const SymbolVector& x) { assign(x); }
    SymbolVectorCompact(const SymbolVectorCompact& x) : impl(x.impl) {}
    
    SymbolVectorCompact& operator=(const SymbolVectorCompact& x)
    {
      impl = x.impl;
      return *this;
    }
    
    SymbolVectorCompact& operator=(const SymbolVector& x)
    {
      assign(x);
      return *this;
    }
    
    void assign(const SymbolVectorCompact& x)
    {
      impl.assign(x.impl);
    }
    
    template <typename Iterator>
    void assign(Iterator first, Iterator last)
    {
      impl.resize(std::distance(first, last) * 8);
      
      impl_type::iterator citer = impl.begin();
      for (/**/; first != last; ++ first)
	std::advance(citer, codec_type::encode(&(*citer), symbol_type(*first)));
      
      impl.erase(citer, impl.end());
    }
    
    void assign(const SymbolVector& x)
    {
      impl.resize(x.size() * 8);
      
      impl_type::iterator citer = impl.begin();
      
      SymbolVector::const_iterator last = x.end();
      for (SymbolVector::const_iterator iter = x.begin(); iter != last; ++ iter)
	std::advance(citer, codec_type::encode(&(*citer), *iter));
      
      impl.erase(citer, impl.end());
    }
	  
  public:
    const_iterator begin() const { return const_iterator(&(*impl.begin()), &(*impl.end())); }
    const_iterator end() const { return const_iterator(); }
    
    bool empty() const { return impl.empty(); }
    size_type size_compressed() const  { return impl.size(); }

    void clear() { impl.clear(); }

    void swap(SymbolVectorCompact& x)
    {
      impl.swap(x.impl);
    }

  public:
    friend size_t hash_value(SymbolVectorCompact const& x) { return utils::hashmurmur<size_t>()(x.impl.begin(), x.impl.end(), 0); }
    friend bool operator==(const SymbolVectorCompact& x, const SymbolVectorCompact& y) { return x.impl == y.impl; }
    friend bool operator!=(const SymbolVectorCompact& x, const SymbolVectorCompact& y) { return x.impl != y.impl; }
    friend bool operator<(const SymbolVectorCompact& x, const SymbolVectorCompact& y) { return x.impl < y.impl; }
    friend bool operator>(const SymbolVectorCompact& x, const SymbolVectorCompact& y) { return x.impl > y.impl; }
    friend bool operator<=(const SymbolVectorCompact& x, const SymbolVectorCompact& y) { return x.impl <= y.impl; }
    friend bool operator>=(const SymbolVectorCompact& x, const SymbolVectorCompact& y) { return x.impl >= y.impl; }

    
  private:
    impl_type impl;
  };
};

namespace std
{
  inline
  void swap(cicada::SymbolVectorCompact& x, cicada::SymbolVectorCompact& y)
  {
    x.swap(y);
  }
};

#endif
