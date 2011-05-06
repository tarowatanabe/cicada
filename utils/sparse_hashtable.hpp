// -*- mode: c++ -*-

#ifndef __UTILS__SPARSE_HASHTABLE__H__
#define __UTILS__SPARSE_HASHTABLE__H__ 1

#include <string>
#include <memory>

#include <google/sparse_hash_set>

#include <utils/memory.hpp>

namespace utils
{
  
  template <typename Key,
	    typename Value,
	    typename ExtractKey,
	    typename Hash,
	    typename Equal,
	    typename Alloc>
  class sparse_hashtable : public Alloc
  {
  public:
    typedef Key   key_type;
    typedef Value value_type;
    typedef Alloc allocator_type;
    
  private:
    struct hasher : public Hash, public ExtractKey
    {
      size_t operator()(const value_type* x) const
      {
	return Hash::operator()(ExtractKey::operator()(*x));
      }
    };
    struct equal : public Equal, public ExtractKey
    {
      bool operator()(const value_type* x, const value_type* y) const
      {
	return Equal::operator()(ExtractKey::operator()(*x), ExtractKey::operator()(*y));
      }
    };

    typedef typename Alloc::template rebind<value_type*>::other map_alloc_type;
    typedef google::sparse_hash_set<value_type*, hasher, equal, map_alloc_type> hashtable_type;
    
    typedef typename hashtable_type::iterator       iterator_base_type;
    typedef typename hashtable_type::const_iterator const_iterator_base_type;
    
  public:
    typedef typename hashtable_type::size_type       size_type;
    typedef typename hashtable_type::difference_type difference_type;
    
    class iterator : public iterator_base_type
    {
    private:
      typedef iterator_base_type base_type;
      
    public:
      typedef typename base_type::iterator_category iterator_category;
      typedef typename base_type::difference_type   difference_type;
      typedef Value       value_type;
      typedef value_type* pointer;
      
      iterator(const base_type& x) : base_type(x) {}
      
      value_type& operator*() { return (*base_type::operator*()); }
      value_type* operator->() { return &(*base_type::operator*()); }
    };
    
    class const_iterator : public const_iterator_base_type
    {
    public:
      typedef const_iterator_base_type base_type;
      typedef typename base_type::iterator_category iterator_category;
      typedef typename base_type::difference_type   difference_type;
      typedef Value       value_type;
      typedef value_type* pointer;
      
      const_iterator(const base_type& x) : base_type(x) {}
      
      const value_type& operator*() { return (*base_type::operator*()); }
      const value_type* operator->() { return &(*base_type::operator*()); }
    };
    
  public:
    sparse_hashtable() : hashtable() {  hashtable.set_deleted_key(0); }
    sparse_hashtable(const sparse_hashtable& x) : hashtable() { hashtable.set_deleted_key(0); assign(x); }
    ~sparse_hashtable() { clear(); }
    sparse_hashtable& operator=(const sparse_hashtable& x)
    {
      assign(x);
      return *this;
    }
    
  public:
    void swap(sparse_hashtable& x)
    {
      using namespace std;
      swap(alloc(), x.alloc());
      swap(hashtable, x.hashtable);
    }
    
    bool empty() const { return hashtable.empty(); }
    size_type size() const { return hashtable.size(); }
    
  public:
    inline const_iterator begin() const { return hashtable.begin(); }
    inline       iterator begin()       { return hashtable.begin(); }
    inline const_iterator end() const { return hashtable.end(); }
    inline       iterator end()       { return hashtable.end(); }
    inline const_iterator find(const value_type& x) const { return hashtable.find(const_cast<value_type*>(&x)); }
    inline       iterator find(const value_type& x)       { return hashtable.find(const_cast<value_type*>(&x)); }

    std::pair<iterator, bool> insert(const value_type& x)
    {
      typename hashtable_type::iterator iter = hashtable.find(const_cast<value_type*>(&x));
      if (iter != hashtable.end())
	return std::make_pair(iter, false);
      else {
	value_type* p = alloc().allocate(1);
	utils::construct_object(p, x);
	return hashtable.insert(p);
      }
    }
    
    void erase(const value_type& x)
    {
      typename hashtable_type::iterator iter = hashtable.find(const_cast<value_type*>(&x));
      if (iter != hashtable.end()) {
	value_type* p = *iter;
	
	hashtable.erase(iter);

	utils::destroy_object(p);
	alloc().deallocate(p, 1);
      }
    }
    
    void erase(iterator x)
    {
      value_type* p = &(*x);
      
      hashtable.erase(x);
      
      utils::destroy_object(p);
      alloc().deallocate(p, 1);
    }
    
    
    void clear()
    {
      typename hashtable_type::iterator iter_end = hashtable.end();
      for (typename hashtable_type::iterator iter = hashtable.begin(); iter != iter_end; ++ iter) {
	utils::destroy_object(*iter);
	alloc().deallocate(*iter, 1);
      }
      hashtable.clear();
    }
    
    void assign(const sparse_hashtable& x)
    {
      clear();
      
      typename hashtable_type::const_iterator iter_end = x.hashtable.end();
      for (typename hashtable_type::const_iterator iter = x.hashtable.begin(); iter != iter_end; ++ iter) {
	value_type* p = alloc().allocate(1);
	utils::construct_object(p, *(*iter));
	hashtable.insert(p);
      }
    }
    
  private:
    allocator_type& alloc() { return static_cast<allocator_type&>(*this); }
    const allocator_type& alloc() const { return static_cast<const allocator_type&>(*this); }
    
  private:
    hashtable_type hashtable;
  };
};

namespace std
{
  template <typename K, typename V, typename X, typename H, typename E, typename A>
  inline
  void swap(utils::sparse_hashtable<K,V,X,H,E,A>& x,
	    utils::sparse_hashtable<K,V,X,H,E,A>& y)
  {
    x.swap(y);
  }
};


#endif

