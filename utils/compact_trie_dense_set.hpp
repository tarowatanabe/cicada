// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__COMPACT_TRIE_DENSE_SET__HPP__
#define __UTILS__COMPACT_TRIE_DENSE_SET__HPP__ 1

#include <google/dense_hash_map>

#include <utils/chunk_vector.hpp>

#include <boost/functional/hash.hpp>

namespace utils
{
  struct __compact_trie_dense_set_base
  {
    typedef uint32_t                   id_type;
    
    static const id_type& npos() {
      static const id_type __npos(-1);
      return __npos;
    }
  };

  template <typename Key,
	    typename Hash=boost::hash<Key>,
	    typename Equal=std::equal_to<Key>,
	    typename Alloc=std::allocator<Key > >
  class compact_trie_dense_set : public __compact_trie_dense_set_base
  {
  public:
    typedef Key                        key_type;
    typedef Key                        value_type;
    
    typedef Hash                       hash_type;
    typedef Equal                      equal_type;
    
    typedef size_t                     size_type;
    typedef ptrdiff_t                  difference_type;
    
  private:  
    typedef typename Alloc::template rebind<std::pair<const key_type, id_type> >::other id_map_alloc_type;
    
    typedef google::dense_hash_map<key_type, id_type, hash_type, equal_type > id_map_type;
    typedef google::dense_hash_map<key_type, id_type, hash_type, equal_type > id_map_root_type;
    
    struct Node
    {
      id_map_type __map;
      
      Node(const key_type& __empty) : __map() { __map.set_empty_key(__empty); }
    };
    typedef Node node_type;
    
    typedef typename Alloc::template rebind<node_type>::other node_alloc_type;
    typedef utils::chunk_vector<node_type, 4096 / sizeof(node_type), node_alloc_type> node_set_type;

  public:
    typedef typename id_map_type::const_iterator const_iterator;
    typedef typename id_map_type::const_iterator       iterator;

    typedef typename id_map_root_type::const_iterator const_root_iterator;
    typedef typename id_map_root_type::const_iterator       root_iterator;
    
  public:
    compact_trie_dense_set(const key_type& __empty) { __root.set_empty_key(__empty); }

  private:
    compact_trie_dense_set() {}
    
  public:
    const_root_iterator begin() const { return __root.begin(); }
    const_root_iterator end() const { return __root.end(); }

    const_iterator begin(id_type __id) const { return __nodes[__id].__map.begin(); }
    const_iterator end(id_type __id) const { return __nodes[__id].__map.end(); }
    
    void clear() { __root.clear(); __nodes.clear(); }
    
    bool empty() const { return __nodes.empty(); }
    bool empty(id_type __id) const {
      if (__id == npos())
	return __root.empty();
      else
	return __nodes[__id].__map.empty();
    }
    
    bool is_root(id_type __id) const { return __id == npos(); }

    void swap(compact_trie_dense_set& x)
    {
      __root.swap(x.__root);
      __nodes.swap(x.__nodes);
    }
    
    id_type root() const { return npos(); }
    
    id_type find(id_type __id, const key_type& key) const
    {
      if (__nodes.empty())
	return npos();

      if (__id == npos()) {
	typename id_map_root_type::const_iterator riter = __root.find(key);
	return (riter != __root.end() ? riter->second : npos());
      } else {
	typename id_map_type::const_iterator niter = __nodes[__id].__map.find(key);
	return (niter != __nodes[__id].__map.end() ? niter->second : npos());
      }
    }
    
    template <typename Iterator>
    id_type find(Iterator first, Iterator last) const
    {
      if (__nodes.empty())
	return npos();
      
      id_type __id = npos();
      for (/**/; first != last; ++ first) {
	__id = find(__id, *first);
	if (__id == npos())
	  return __id;
      }
      return __id;
    }
    
    id_type insert(id_type __id, const key_type& key)
    {
      return __insert_key(__id, key);
    }
    
    template <typename Iterator>
    id_type insert(Iterator first, Iterator last)
    {
      typedef typename boost::is_integral<Iterator>::type __integral;
      return __insert_dispatch(first, last, __integral());
    }
    
    
    
  private:
    template <typename Integer>
    id_type __insert_dispatch(Integer __id, Integer __key, boost::true_type)
    {
      return __insert_key(__id, __key);
    }
    
    template <typename Iterator>
    id_type __insert_dispatch(Iterator first, Iterator last, boost::false_type)
    {
      return __insert_range(first, last);
    }
    
    
    id_type __insert_key(id_type __id, const key_type& key)
    {
      if (__id == npos())  {
	typename id_map_root_type::iterator riter = __root.find(key);
	if (riter != __root.end())
	  return riter->second;
	else {
	  __root.insert(std::make_pair(key, __nodes.size()));
	  __id = __nodes.size();
	  __nodes.push_back(node_type(__root.empty_key()));
	  return __id;
	}
      } else {
	typename id_map_type::iterator niter = __nodes[__id].__map.find(key);
	if (niter != __nodes[__id].__map.end())
	  return niter->second;
	else {
	  __nodes[__id].__map.insert(std::make_pair(key, __nodes.size()));
	  __id = __nodes.size();
	  __nodes.push_back(node_type(__root.empty_key()));
	  return __id;
	}
      }
    }

    template <typename Iterator>
    id_type __insert_range(Iterator first, Iterator last)
    {
      id_type __id = npos();
      for (/**/; first != last; ++ first)
	__id = __insert_key(__id, *first);
      
      return __id;
    }
    
    
  private:
    id_map_root_type __root;
    node_set_type    __nodes;
  };
};

namespace std
{
  template <typename Key, typename Hash, typename Equal, typename Alloc>
  inline
  void swap(utils::compact_trie_dense_set<Key,Hash,Equal,Alloc>& x,
	    utils::compact_trie_dense_set<Key,Hash,Equal,Alloc>& y)
  {
    x.swap(y);
  }
};

#endif
