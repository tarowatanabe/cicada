// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__COMPACT_TRIE_SET__HPP__
#define __UTILS__COMPACT_TRIE_SET__HPP__ 1

#include <utils/chunk_vector.hpp>
#include <utils/unordered_map.hpp>

#include <boost/functional/hash.hpp>

namespace utils
{
  struct __compact_trie_set_base
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
	    typename Alloc=std::allocator<Key> >
  class compact_trie_set : public __compact_trie_set_base
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
    
    typedef utils::unordered_map<key_type, id_type, hash_type, equal_type, id_map_alloc_type>::type id_map_type;
    
    struct Node
    {
      id_map_type __map;
      
      Node() : __map() { }
    };
    typedef Node node_type;
    
    typedef typename Alloc::template rebind<node_type>::other node_alloc_type;
    typedef utils::chunk_vector<node_type, 4096 / sizeof(node_type), node_alloc_type> node_set_type;

  public:
    typedef typename id_map_type::const_iterator const_iterator;
    typedef typename id_map_type::const_iterator       iterator;
    
  public:
    compact_trie_set() {}
    
  public:
    const_iterator begin() const { return __root.begin(); }
    const_iterator end() const { return __root.end(); }

    const_iterator begin(id_type __id) const { return (__id == npos() ? __root.begin() : __nodes[__id].__map.begin()); }
    const_iterator end(id_type __id) const { return (__id == npos() ? __root.end() : __nodes[__id].__map.end()); }
    
    void clear() { __root.clear(); __nodes.clear(); }

    size_type size() const { return __nodes.size(); }
    
    bool empty() const { return __nodes.empty(); }
    bool empty(id_type __id) const
    {
      return (__id == npos() ? __root.empty() : __nodes[__id].__map.empty());
    }
    
    bool is_root(id_type __id) const { return __id == npos(); }

    void swap(compact_trie_set& x)
    {
      __root.swap(x.__root);
      __nodes.swap(x.__nodes);
    }
    
    id_type root() const { return npos(); }
    
    id_type find(id_type __id, const key_type& key) const
    {
      return __find_key(__id, key);
    }
    
    template <typename Iterator>
    id_type find(Iterator first, Iterator last) const
    {
      typedef typename boost::is_integral<Iterator>::type __integral;
      return __find_dispatch(first, last, __integral());
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
    // find...
    template <typename Integer>
    id_type __find_dispatch(Integer __id, Integer __key, boost::true_type) const
    {
      return __find_key(__id, __key);
    }
    
    template <typename Iterator>
    id_type __find_dispatch(Iterator first, Iterator last, boost::false_type) const
    {
      return __find_range(first, last);
    }
    
    id_type __find_key(const id_type& __id, const key_type& key) const
    {
      if (__nodes.empty())
	return npos();
      
      const id_map_type& mapping = (__id == npos() ? __root : __nodes[__id].__map);
      
      typename id_map_type::const_iterator riter = mapping.find(key);
      return (riter != mapping.end() ? riter->second : npos());
    }
    
    template <typename Iterator>
    id_type __find_range(Iterator first, Iterator last) const
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
    
    // insertions...
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
    
    id_type __insert_key(const id_type& __id, const key_type& key)
    {
      id_map_type& mapping = (__id == npos() ? __root : __nodes[__id].__map);
      
      typename id_map_type::iterator riter = mapping.find(key);
      if (riter != mapping.end())
	return riter->second;
      else {
	const id_type __id_new = __nodes.size();
	
	__nodes.push_back(node_type());
	mapping.insert(std::make_pair(key, __id_new));
	
	return __id_new;
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
    id_map_type   __root;
    node_set_type __nodes;
  };
};

namespace std
{
  template <typename Key, typename Hash, typename Equal, typename Alloc>
  inline
  void swap(utils::compact_trie_set<Key,Hash,Equal,Alloc>& x,
	    utils::compact_trie_set<Key,Hash,Equal,Alloc>& y)
  {
    x.swap(y);
  }
};

#endif
