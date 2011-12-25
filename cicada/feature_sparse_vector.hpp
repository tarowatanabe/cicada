// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE_SPARSE_VECTOR__HPP__
#define __CICADA__FEATURE_SPARSE_VECTOR__HPP__ 1

#include <map>
#include <memory>
#include <utility>
#include <algorithm>
#include <iterator>
#include <iostream>

#include <cicada/feature.hpp>

#include <google/dense_hash_map>

//
// feature vector, for use as "temporary" sparse vector for faster access, w/o sorting
//

namespace cicada
{
   template <typename Tp, typename Alloc=std::allocator<Tp> >
   class FeatureSparseVector
   {
  public:
     typedef cicada::Feature feature_type;
     
     typedef size_t    size_type;
     typedef ptrdiff_t difference_type;

   private:
     typedef google::dense_hash_map<feature_type, data_type, boost::hash<feature_type>, std::equal_to<feature_type> > map_type;
     
   public:
     typedef typename map_type::key_type    key_type;
     typedef typename map_type::data_type   data_type;
     typedef typename map_type::mapped_type mapped_type;
     typedef typename map_type::value_type  value_type;
     
     typedef typename map_type::const_iterator  const_iterator;
     typedef typename map_type::iterator        iterator;
     
     typedef typename map_type::const_reference  const_reference;
     typedef typename map_type::reference        reference;
     
   public:
     typedef FeatureSparseVector<Tp, Alloc> self_type;
     
   public:
     FeatureSparseVector() : __map() { initialize(); }
     FeatureSparseVector(const self_type& x) : __map(x.__map) { }
     template <typename T, typename A>
     FeatureSparseVector(const FeatureSparseVector<T,A>& x) : __map() { initialize(); __map.insert(x.begin(), x.end()); }
     template <typename Iterator>
     FeatureSparseVector(Iterator first, Iterator last) : __map() { initialize(); __map.insert(first, last); } 
     
     FeatureSparseVector& operator=(const self_type& x)
     {
       __map = x.__map;
       return *this;
     }
     
     template <typename T, typename A>
     FeatureSparseVector& operator=(const FeatureSparseVector<T,A>& x)
     {
       __map.clear();
       __map.insert(x.begin(), x.end());
       return *this;
     }

   public:
     size_type size() const { return __map.size(); }
     bool empty() const { return __map.empty(); }
     
     void assign(const self_type& x)
     {
       __map = x.__map;
     }
     
     template <typename T, typename A>
     void assign(const FeatureSparseVector<T,A>& x)
     {
       __map.clear();
       __map.insert(x.begin(), x.end());
     }
     
     template <typename Iterator>
     void assign(Iterator first, Iterator last)
     {
       __map.clear();
       __map.insert(first, last);
     }
     
     Tp operator[](const key_type& x) const
     {
       const_iterator iter = find(x);
       return (iter == end() ? Tp() : iter->second);
     }
     
     Tp& operator[](const key_type& x)
     {
       return __map[x];
     }

     inline const_iterator begin() const { return __map.begin(); }
     inline       iterator begin()       { return __map.begin(); }
     
     inline const_iterator end() const { return __map.end(); }
     inline       iterator end()       { return __map.end(); }
     
     inline const_iterator find(const key_type& x) const { return __map.find(x); }
     inline       iterator find(const key_type& x)       { return __map.find(x); }

     void erase(const key_type& x) { __map.erase(x); }
     
     void swap(FeatureSparseVector& x) { __map.swap(x.__map); }
     
     void clear() { __map.clear(); }
     
   private:
     void initialize()
     {
       __map.set_empty_key(feature_type(feature_type::id_type(-1)));
       __map.set_deleted_key(feature_type(feature_type::id_type(-2)));
     }
     
   private:
     map_type __map;
   };
};

#endif
