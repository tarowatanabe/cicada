// -*- mode: c++ -*-

#ifndef __UTILS__ALLOC_VECTOR__HPP__
#define __UTILS__ALLOC_VECTOR__HPP__ 1

#include <vector>

#include <utils/memory.hpp>

// we perform allocation everytime operator[] is called...

namespace utils
{
  template <typename _Tp, typename _Alloc=std::allocator<_Tp> >
  class alloc_vector
  {
  private:
    typedef _Tp* map_type;
    typedef typename _Alloc::template rebind<_Tp*>::other map_alloc_type;
    typedef std::vector<map_type, map_alloc_type> map_vector_type;
    typedef _Alloc alloc_type;
    
    struct impl_type : public alloc_type
    {
      inline       map_vector_type& map()       { return __map; }
      inline const map_vector_type& map() const { return __map; }
      inline       alloc_type& alloc()       { return static_cast<alloc_type&>(*this); }
      inline const alloc_type& alloc() const { return static_cast<const alloc_type&>(*this); }

      void swap(impl_type& x)
      {
	__map.swap(x.__map);
	std::swap(alloc(), x.alloc());
      }

      void resize(size_t __size)
      {
	if (__size == __map.size()) return;
	
	if (__size > __map.size())
	  __map.resize(__size, 0);
	else {
	  deallocate_range(__map.begin() + __size, __map.end());
	  __map.resize(__size, 0);
	}
      }
      
      void clear()
      {
	deallocate_range(__map.begin(), __map.end());
	__map.clear();
      }
      
      template <typename Iterator>
      void deallocate_range(Iterator first, Iterator last)
      {
	for (/**/; first != last; ++ first)
	  if (*first)
	    alloc().deallocate(*first, 1);
      }
      
      impl_type() : __map() {}
      impl_type(size_t __n) : __map(__n, 0) { }
      ~impl_type() { clear(); }
      
    private:
      map_vector_type __map;
    };
    
  public:
    typedef _Tp                                      value_type;
    typedef _Tp                                      mapped_type;
    typedef size_t                                   size_type;
    typedef ptrdiff_t                                difference_type;
    
    typedef typename map_vector_type::iterator               iterator;
    typedef typename map_vector_type::const_iterator         const_iterator;
    typedef typename map_vector_type::reverse_iterator       reverse_iterator;
    typedef typename map_vector_type::const_reverse_iterator const_reverse_iterator;
    typedef       _Tp&                                       reference;
    typedef const _Tp&                                       const_reference;
    
  public:
    alloc_vector() : impl() {}
    alloc_vector(size_type __n) : impl(__n) {}
    alloc_vector(const alloc_vector& x) : impl(x.size()) {  __copy_initialize(x); }
    ~alloc_vector() { clear(); }
    
    alloc_vector& operator=(const alloc_vector& x)
    {
      __copy_assign(x);
      return *this;
    }
    
  public:
    
    inline const_iterator begin() const { return impl.map().begin(); }
    inline       iterator begin()       { return impl.map().begin(); }
    
    inline const_iterator end() const { return impl.map().end(); }
    inline       iterator end()       { return impl.map().end(); }
    
    inline const_reverse_iterator rbegin() const { return impl.map().rbegin(); }
    inline       reverse_iterator rbegin()       { return impl.map().rbegin(); }
    
    inline const_reverse_iterator rend() const { return impl.map().rend(); }
    inline       reverse_iterator rend()       { return impl.map().rend(); }
    
    size_type size() const { return impl.map().size(); }
    bool empty() const { return impl.map().empty(); }
    bool exists(size_type __pos) const { return (__pos < impl.map().size() && impl.map()[__pos]); }

    void push_back(const value_type& __value)
    {
      impl.map().push_back(impl.alloc().allocate(1));
      utils::construct_object(impl.map().back(), __value);
    }
    
    void resize(size_type __size)
    {
      if (__size < impl.map().size())
	__destroy_range(impl.map().begin() + __size, impl.map().end());
      impl.resize(__size);
    }
    void resize(size_type __size, const value_type& __value) { resize(__size); }
    void reserve(size_type __size) { impl.map().reserve(__size); }
    
    void shrink() { map_vector_type(impl.map()).swap(impl.map()); }
    void clear()
    { 
      __destroy_range(impl.map().begin(), impl.map().end());
      impl.clear();
    }

    
    void swap(alloc_vector& x) { impl.swap(x.impl); }
    
    void erase(iterator __iter)
    {
      if (! *__iter) return;
      
      utils::destroy_object(*__iter);
      impl.alloc().deallocate(*__iter, 1);
      *__iter = 0;
    }
    
    void erase(iterator __first, iterator __last)
    {
      for (/**/; __first != __last; ++ __first)
	if (*__first) {
	  utils::destroy_object(*__first);
	  impl.alloc().deallocate(*__first, 1);
	  *__first = 0;
	}
    }
    
    reference operator[](size_type __n)
    {
      if (__n >= impl.map().size())
	impl.map().resize(__n + 1, 0);
      if (! impl.map()[__n]) {
	// perform allocation and construction...
	impl.map()[__n] = impl.alloc().allocate(1);
	utils::construct_object(impl.map()[__n]);
      }
      return *(impl.map()[__n]);
    }
    
    const_reference operator[](size_type __n) const
    {
      if (__n >= impl.map().size())
	const_cast<impl_type&>(impl).map().resize(__n + 1, 0);
      if (! impl.map()[__n]) {
	// perform allocation and construction...
	const_cast<impl_type&>(impl).map()[__n] = const_cast<impl_type&>(impl).alloc().allocate(1);
	utils::construct_object(const_cast<impl_type&>(impl).map()[__n]);
      }
      return *(impl.map()[__n]);
    }
    
  private:
    void __copy_initialize(const alloc_vector& x)
    {
      iterator iter2 = begin();
      for (const_iterator iter = x.begin(); iter != x.end(); ++ iter, ++ iter2)
	if (*iter) {
	  *iter2 = impl.alloc().allocate(1);
	  utils::construct_object(*iter2, *(*iter));
	}
    }
    
    void __copy_assign(const alloc_vector& x)
    {
      resize(x.size());
      
      iterator iter2 = begin();
      for (const_iterator iter = x.begin(); iter != x.end(); ++ iter, ++ iter2)
	if (*iter && *iter2)
	  *(*iter2) = *(*iter);
	else if (*iter) {
	  *iter2 = impl.alloc().allocate(1);
	  utils::construct_object(*iter2, *(*iter));
	} else if (*iter2) {
	  utils::destroy_object(*iter2);
	  impl.alloc().deallocate(*iter2, 1);
	  *iter2 = 0;
	}
    }
    
    void __destroy_range(iterator first, iterator last)
    {
      typedef typename boost::has_trivial_destructor<_Tp>::type __trivial;
      __destroy_range(first, last, __trivial());
    }
    
    void __destroy_range(iterator first, iterator last, const boost::true_type&)
    {
      
    }
    void __destroy_range(iterator first, iterator last, const boost::false_type&)
    {
      for (/**/; first != last; ++ first)
	if (*first)
	  utils::destroy_object(*first);
    }
    
  private:
    impl_type impl;
  };
  
};

namespace std
{
  template <typename Tp, typename Alloc>
  inline
  void swap(utils::alloc_vector<Tp,Alloc>& x,
	    utils::alloc_vector<Tp,Alloc>& y)
  {
    x.swap(y);
  }
};

#endif
