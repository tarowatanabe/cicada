// -*- mode: c++ -*-

#ifndef __SUCCINCTDB__SUCCINCT_DB__HPP__
#define __SUCCINCTDB__SUCCINCT_DB__HPP__ 1


#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <vector>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include <succinct_db/succinct_trie.hpp>

#include <utils/repository.hpp>
#include <utils/tempfile.hpp>

namespace succinctdb
{
  template <typename Key, typename Data, typename Alloc=std::allocator<std::pair<Key, Data> > >
  struct __succinct_db_writer
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef int64_t   off_type;
    typedef uint32_t  pos_type;
    
    typedef Key       key_type;
    typedef Data      data_type;
    typedef Data      mapped_type;
    
    typedef boost::filesystem::path path_type;

    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<data_type>::other data_alloc_type;
    typedef typename Alloc::template rebind<char>::other      byte_alloc_type;
    
    __succinct_db_writer(const path_type& path) { open(path); }
    ~__succinct_db_writer() { close(); }
    
    size_type insert(const key_type* key, size_type key_size, const data_type* data, size_type data_size)
    {
      __os_key->write((char*) key, sizeof(key_type) * key_size);
      __os_data->write((char*) data, sizeof(data_type) * data_size);
      
      __key_offset += key_size;
      __os_key_offset->write((char*) &__key_offset, sizeof(__key_offset));
      
      __data_offset += data_size;
      __os_data_offset->write((char*) &__data_offset, sizeof(__data_offset));
      
      return __size ++;
    }
    
    size_type size() const { return __size; }
    
    void open(const path_type& path)
    {
      
    }
    
    
    struct __value_type
    {
      const key_type* first;
      const key_type* last;
      pos_type        pos;
      
      size_type size() const { return last - first; }
      const key_type& operator[](size_type __pos) const { return *(first + __pos); }

      __value_type() {}
    };
    typedef typename Alloc::template rebind<__value_type>::other __value_alloc_type;
    typedef std::vector<__value_type, __value_alloc_type> __value_set_type;
    
    struct __extract_key
    {
      const __value_type& operator()(const __value_type& x) const { return x; }
    };

    struct __extract_data
    {
      const pos_type& operator()(const __value_type& x) const { return x.pos; }
    };
    
    struct __less_value
    {
      bool operator()(const __value_type& x, const __value_type& y) const {
	return std::lexicographical_compare(x.first, x.last, y.first, x.last);
      }
    };
    
    void close()
    {
      
    }
    
  private:
    path_type path_output;
    path_type path_key;
    path_type path_key_offset;
    path_type path_data;
    path_type path_data_offset;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_key;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_key_offset;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_data;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_data_offset;
    off_type __key_offset;
    off_type __data_offset;
    size_type __size;
  };
  

  template <typename Key, typename Data, typename Alloc=std::allocator<std::pair<Key, Data> > >
  class succinct_db
  {
    
    
  };
  
};

#endif
