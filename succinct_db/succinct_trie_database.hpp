// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __SUCCINCT_DB__SUCCINCT_TRIE_DATABASE__HPP__
#define __SUCCINCT_DB__SUCCINCT_TRIE_DATABASE__HPP__ 1


#include <stdint.h>

#include <iostream>
#include <algorithm>
#include <vector>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include <succinct_db/succinct_trie.hpp>

#include <utils/map_file.hpp>
#include <utils/vertical_coded_vector.hpp>
#include <utils/vertical_coded_device.hpp>
#include <utils/repository.hpp>
#include <utils/tempfile.hpp>
#include <utils/bithack.hpp>

namespace succinctdb
{
  // writer
  template <typename Key, typename Data, typename Alloc=std::allocator<std::pair<Key, Data> > >
  struct __succinct_trie_database_writer
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint64_t  off_type;
    typedef uint64_t  pos_type;
    
    typedef Key       key_type;
    typedef Data      data_type;
    typedef Data      mapped_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<data_type>::other data_alloc_type;
    typedef typename Alloc::template rebind<off_type>::other  off_alloc_type;
    typedef typename Alloc::template rebind<char>::other      byte_alloc_type;
    typedef typename Alloc::template rebind<std::pair<key_type, pos_type> >::other trie_alloc_type;

    // we assume that pointer size is multiple of two!
    static const size_type pointer_size = sizeof(void*);
    static const size_type pointer_mask = ~(pointer_size - 1);
    
    __succinct_trie_database_writer(const path_type& path) 
      : path_output(), path_key_data(), path_size(), __offset(0), __size(0) { open(path); }
    ~__succinct_trie_database_writer() { close(); }

    path_type path() const { return path_output; }
    
    size_type insert(const key_type* buf, size_type buf_size, const data_type* data, size_type data_size)
    {
      const size_type buf_size_bytes   = buf_size * sizeof(key_type);
      const size_type buf_size_aligned = (buf_size_bytes + pointer_size - 1) & pointer_mask;
      
      const size_type pos_size_bytes   = sizeof(pos_type);
      const size_type pos_size_aligned = (pos_size_bytes + pointer_size - 1) & pointer_mask;
      

      // key value
      __os_key_data->write((char*) buf, buf_size * sizeof(key_type));
      if (buf_size_aligned > buf_size_bytes) {
	char __buf[pointer_size];
	__os_key_data->write((char*) __buf, buf_size_aligned - buf_size_bytes);
      }
      
      __os_key_data->write((char*) &__size, sizeof(pos_type));
      if (pos_size_aligned > pos_size_bytes) {
	char __buf[pointer_size];
	__os_key_data->write((char*) __buf, pos_size_aligned - pos_size_bytes);
      }
      
      __os_key_size->write((char*) &buf_size, sizeof(size_type));
      
      // mapped value
      __os_data->write((char*) data, sizeof(data_type) * data_size);
      __offset += data_size;
      __os_data_off->write((char*) &__offset, sizeof(off_type));
      
      return __size ++;
    }
    size_type size() const { return __size; }
    
    
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      close();
      
      repository_type rep(path, repository_type::write);
      rep["type"] = "succinct-trie-database";

      const path_type tmp_dir = utils::tempfile::tmp_dir();

      path_output = path;
      
      path_key_data = utils::tempfile::file_name(tmp_dir / "succinct-db.key-data.XXXXXX");
      path_size     = utils::tempfile::file_name(tmp_dir / "succinct-db.size.XXXXXX");
      
      utils::tempfile::insert(path_key_data);
      utils::tempfile::insert(path_size);
      
      __offset = 0;
      __size = 0;
      
      __os_key_data.reset(new boost::iostreams::filtering_ostream());
      __os_key_size.reset(new boost::iostreams::filtering_ostream());
      
      __os_data.reset(new boost::iostreams::filtering_ostream());
      __os_data_off.reset(new boost::iostreams::filtering_ostream());
      
      __os_key_data->push(boost::iostreams::file_sink(path_key_data.file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      __os_key_size->push(boost::iostreams::zlib_compressor());
      __os_key_size->push(boost::iostreams::file_sink(path_size.file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      
      __os_data->push(boost::iostreams::file_sink(rep.path("mapped").file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      __os_data_off->push(utils::vertical_coded_sink<off_type, off_alloc_type>(rep.path("offset"), 1024 * 1024));
      
      // initial data offset...
      __os_data_off->write((char*) &__offset, sizeof(off_type));
    }


    struct __value_type
    {
      size_type size() const { return last - first; }
      const key_type& operator[](size_type pos) const { return *(first + pos); }
      
      __value_type() : first(0), last(0) {}
      
      const key_type* first;
      const key_type* last;
    };
    typedef typename Alloc::template rebind<__value_type>::other __value_alloc_type;
    typedef std::vector<__value_type, __value_alloc_type> __value_set_type;
    
    struct __extract_key
    {
      const __value_type& operator()(const __value_type& x) const { return x; }
    };
    
    struct __extract_data
    {
      const pos_type& operator()(const __value_type& x) const
      {
	const size_type key_size_bytes = sizeof(key_type) * (x.last - x.first);
	const size_type key_size_aligned = (key_size_bytes + pointer_size - 1) & pointer_mask;
	
	return *reinterpret_cast<const pos_type*>(((char*) x.first) + key_size_aligned);
      }
    };

    struct __less_value
    {
      bool operator()(const __value_type& x, const __value_type& y) const {
	return std::lexicographical_compare(x.first, x.last, y.first, x.last);
      }
    };

    void clear() { close(); }
    
    void close()
    {
      typedef utils::map_file<char, byte_alloc_type> map_file_type;
      typedef succinct_trie<key_type, pos_type, trie_alloc_type> succinct_trie_type;
      
      if (__os_key_data)
	__os_key_data->flush();
      if (__os_key_size)
	__os_key_size->flush();
      
      if (__os_data)
	__os_data->flush();
      if (__os_data_off)
	__os_data_off->flush();

      __os_key_data.reset();
      __os_key_size.reset();
      
      __os_data.reset();
      __os_data_off.reset();
      
      if (boost::filesystem::exists(path_key_data) && boost::filesystem::exists(path_size) && ! path_output.empty()) {
	
	::sync();

	typedef utils::repository repository_type;
	
	repository_type rep(path_output, repository_type::read);
	
	map_file_type map_key_data(path_key_data);
	__value_set_type values(__size);
	
	if (__size > 0) {
	  const char* iter = reinterpret_cast<const char*>(&(*map_key_data.begin()));
	  boost::iostreams::filtering_istream is;
	  is.push(boost::iostreams::zlib_decompressor());
	  is.push(boost::iostreams::file_source(path_size.file_string()));
	  for (size_type i = 0; i < __size; ++ i) {
	    size_type key_size = 0;
	    is.read((char*) &key_size, sizeof(size_type));
	    
	    values[i].first = reinterpret_cast<const key_type*>(iter);
	    values[i].last  = reinterpret_cast<const key_type*>(iter) + key_size;
	    
	    const size_type key_size_bytes = key_size * sizeof(key_type);
	    const size_type key_size_aligned = (key_size_bytes + pointer_size - 1) & pointer_mask;
	    
	    const size_type pos_size_bytes   = sizeof(pos_type);
	    const size_type pos_size_aligned = (pos_size_bytes + pointer_size - 1) & pointer_mask;
	    
	    iter += key_size_aligned + pos_size_aligned;
	  }
	}
	boost::filesystem::remove(path_size);
	utils::tempfile::erase(path_size);
	
	// sorting
	std::sort(values.begin(), values.end(), __less_value());
	
	{
	  succinct_trie_type succinct_trie;
	  succinct_trie.build(rep.path("index"), values.begin(), values.end(), __extract_key(), __extract_data());
	}
	
	map_key_data.clear();
	boost::filesystem::remove(path_key_data);
	utils::tempfile::erase(path_key_data);
      }
      
      path_output = path_type();
      path_key_data = path_type();
      path_size = path_type();
            
      __offset = 0;
      __size = 0;
    }
    
  private:
    path_type path_output;
    
    path_type path_key_data;
    path_type path_size;
        
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_key_data;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_key_size;
    
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_data;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __os_data_off;
    
    off_type __offset;
    pos_type __size;
  };
  
  
  
  // mapped value access...
  template <typename Iterator>
  struct __succinct_trie_database_mapped
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef Iterator iterator;
    
  public:
    __succinct_trie_database_mapped()
      : __first(), __last() {}
    __succinct_trie_database_mapped(iterator first, iterator last)
      : __first(first), __last(last) {}
    
  public:
    iterator begin() const { return __first; }
    iterator end() const { return __last; }
    
    size_type size() const { return __last - __first; }
    bool empty() const { return __first == __last; }
    
  private:
    iterator __first;
    iterator __last;
  };
  
  template <typename Key, typename Data, typename Alloc=std::allocator<std::pair<Key, Data> > >
  class succinct_trie_database
  {
  public:
    typedef enum {
      READ,
      WRITE,
    } mode_type;
      
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint64_t  off_type;
    typedef uint64_t  pos_type;
      
    typedef Key  key_type;
    typedef Data data_type;
    typedef Data mapped_type;
      
    typedef boost::filesystem::path path_type;

  private:
    typedef typename Alloc::template rebind<data_type>::other data_alloc_type;
    typedef typename Alloc::template rebind<off_type>::other  off_alloc_type;
    typedef typename Alloc::template rebind<std::pair<key_type, pos_type> >::other trie_alloc_type;
    
  private:
    typedef succinct_trie_mapped<Key, pos_type, trie_alloc_type>          succinct_trie_type;
    typedef utils::map_file<data_type, data_alloc_type>                   data_set_type;
    typedef utils::vertical_coded_vector_mapped<off_type, off_alloc_type> off_set_type;
    
    typedef __succinct_trie_database_writer<Key,Data,Alloc>                     succinct_writer_type;

  public:
    typedef __succinct_trie_database_mapped<typename data_set_type::const_iterator> value_type;
      
  public:
    succinct_trie_database() {}
    succinct_trie_database(const path_type& path, const mode_type mode=READ) { open(path, mode); }
    ~succinct_trie_database() { close(); }
      
  public:
    path_type path() const
    {
      if (__succinct_trie)
	return __succinct_trie->path().parent_path();
      else if (__succinct_writer)
	return __succinct_writer->path();
      else
	return path_type();
    }

    // methods supported by both read/write mode
    void open(const path_type& path, const mode_type mode=READ)
    {
      clear();
	
      if (mode == READ) {
	typedef utils::repository repository_type;
	  
	repository_type rep(path, repository_type::read);
	  
	__succinct_trie.reset(new succinct_trie_type(rep.path("index")));
	__mapped.open(rep.path("mapped"));
	__offsets.open(rep.path("offset"));
	  
      } else
	__succinct_writer.reset(new succinct_writer_type(path));
    }

    void write(const path_type& file) const
    {
      if (file == path()) return;
      
      if (__succinct_trie) {
	typedef utils::repository repository_type;
	
	repository_type rep(file, repository_type::write);
	rep["type"] = "succinct-trie-database";
	
	__succinct_trie->write(rep.path("index"));
	__mapped.write(rep.path("mapped"));
	__offsets.write(rep.path("offset"));
      }
    }
    
    void close() { __succinct_trie.reset(); __succinct_writer.reset(); __mapped.clear(); __offsets.clear(); }
    void clear() { close(); }
      
    bool is_open() const { return __succinct_trie || __succinct_writer; }
    bool is_writer() const { return __succinct_writer; }
    bool is_reader() const { return __succinct_trie; }
      
    size_type size() const { return (is_open() ? (is_reader() ? __succinct_trie->size() : __succinct_writer->size()) : size_type(0)); }
    bool empty() const { return size() == 0; }
      
  public:
    size_type insert(const key_type* key, size_type key_size, const data_type* data, size_type data_size)
    {
      return __succinct_writer->insert(key, key_size, data, data_size);
    }
    
  public:    
    // structures supported by read-mode
    
    struct cursor : public succinct_trie_type::cursor
    {
      typedef typename succinct_trie_type::cursor base_type;
      typedef succinct_trie_database<Key,Data,Alloc> impl_type;
      
      cursor() : base_type(), __impl() {}
      cursor(const base_type& x, const impl_type* impl) : base_type(x), __impl(impl) {}
      
      value_type data() const { return __impl->operator[](base_type::node()); }
      value_type operator*() const { return data(); }

      const impl_type* __impl;
    };
    typedef cursor const_cursor;
    
    struct iterator : public succinct_trie_type::iterator
    {
      typedef typename succinct_trie_type::iterator base_type;
      typedef succinct_trie_database<Key,Data,Alloc> impl_type;
     
      iterator() : base_type(), __impl() {}
      iterator(const base_type& x, const impl_type* impl) : base_type(x), __impl(impl) {}

      cursor begin() const { return cursor(base_type::begin(), __impl); }
      cursor end() const { return cursor(base_type::end(), __impl); }
      
      value_type data() const { return __impl->operator[](base_type::node()); }
      value_type operator*() const { return data(); }
      
      const impl_type* __impl;
    };
    typedef iterator const_iterator;
    
    struct reverse_iterator : public succinct_trie_type::reverse_iterator
    {
      typedef typename succinct_trie_type::reverse_iterator base_type;
      typedef succinct_trie_database<Key,Data,Alloc> impl_type;
      
      reverse_iterator() : base_type(), __impl() {}
      reverse_iterator(const base_type& x, const impl_type* impl) : base_type(x), __impl(impl) {}
      
      cursor begin() const { return cursor(base_type::begin(), __impl); }
      cursor end() const { return cursor(base_type::end(), __impl); }
      
      value_type data() const { return __impl->operator[](base_type::node()); }
      value_type operator*() const { return data(); }

      const impl_type* __impl;
    };
    typedef reverse_iterator const_reverse_iterator;
    
    
    
  public:
    // operations supported by read-mode
    
    const_iterator begin(size_type node_pos) const { return const_iterator(__succinct_trie->begin(node_pos), this); }
    const_iterator end(size_type node_pos)   const { return const_iterator(__succinct_trie->end(node_pos), this); }
    const_iterator begin() const { return const_iterator(__succinct_trie->begin(), this); }
    const_iterator end()   const { return const_iterator(__succinct_trie->end(), this); }
    
    const_reverse_iterator rbegin(size_type node_pos) const { return const_reverse_iterator(__succinct_trie->rbegin(node_pos), this); }
    const_reverse_iterator rend(size_type node_pos)   const { return const_reverse_iterator(__succinct_trie->rend(node_pos), this); }
    const_reverse_iterator rbegin() const { return const_reverse_iterator(__succinct_trie->rbegin(), this); }
    const_reverse_iterator rend()   const { return const_reverse_iterator(__succinct_trie->rend(), this); }
    
    const_cursor cbegin(size_type node_pos) const { return (is_valid(node_pos) && exists(node_pos)
							    ? const_cursor(__succinct_trie->cbegin(node_pos), this)
							    : const_cursor(__succinct_trie->cend(), this)); }
    const_cursor cend(size_type node_pos)   const { return (is_valid(node_pos) && exists(node_pos)
							    ? const_cursor(__succinct_trie->cend(node_pos), this)
							    : const_cursor(__succinct_trie->cend(), this)); }
    const_cursor cbegin() const { return const_cursor(__succinct_trie->cbegin(), this); }
    const_cursor cend()   const { return const_cursor(__succinct_trie->cend(), this); }
     
    size_type find(const key_type* key_buf, size_type key_size, size_type node_pos) const
    {
      size_type key_pos = 0;
      return traverse(key_buf, node_pos, key_pos, key_size);
    }
 
    size_type find(const key_type* key_buf, size_type key_size) const
    {
      size_type node_pos = 0;
      size_type key_pos = 0;
      return traverse(key_buf, node_pos, key_pos, key_size);
    }
    
    size_type traverse(const key_type* key_buf, size_type& node_pos, size_type& key_pos, size_type key_size) const
    {
      return __succinct_trie->traverse(key_buf, node_pos, key_pos, key_size);
    }
    
    std::pair<size_type, size_type> range(const size_type node_pos) const { return __succinct_trie->range(node_pos); }
    size_type parent(size_type node_pos) const { return __succinct_trie->parent(node_pos); }
    
    bool is_next_sibling(size_type node_pos) const { return __succinct_trie->is_next_sibling(node_pos); }
    bool exists(size_type node_pos) const { return __succinct_trie->exists(node_pos); }
    bool has_children(size_type node_pos) const { return __succinct_trie->has_children(node_pos); }
    bool is_valid(size_type node_pos) const { return node_pos != succinct_trie_type::out_of_range(); }
    
    value_type operator[](size_type node_pos) const
    {
      const pos_type pos = __succinct_trie->operator[](node_pos);
      
      typename data_set_type::const_iterator first = __mapped.begin() + __offsets[pos];
      typename data_set_type::const_iterator last  = __mapped.begin() + __offsets[pos + 1];
      
      return value_type(first, last);
    }

    uint64_t size_bytes() const { return __succinct_trie->size_bytes() + __mapped.size_bytes() + __offsets.size_bytes(); }
    uint64_t size_compressed() const { return __succinct_trie->size_compressed() + __mapped.size_compressed() + __offsets.size_compressed(); }
    uint64_t size_cache() const { return __succinct_trie->size_cache() + __mapped.size_cache() + __offsets.size_cache(); }
    
  public:
    boost::shared_ptr<succinct_trie_type>   __succinct_trie;
    boost::shared_ptr<succinct_writer_type> __succinct_writer;
    
    data_set_type __mapped;
    off_set_type  __offsets;
  };  
};

#endif
