// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __SUCCINCT_DB__SUCCINCT_HASH__HPP__
#define __SUCCINCT_DB__SUCCINCT_HASH__HPP__ 1

#include <stdint.h>

#include <vector>
#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <memory>
#include <utility>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>

#include <utils/bithack.hpp>
#include <utils/map_file.hpp>
#include <utils/vertical_coded_vector.hpp>
#include <utils/vertical_coded_device.hpp>
#include <utils/packed_vector.hpp>
#include <utils/packed_device.hpp>
#include <utils/repository.hpp>
#include <utils/filesystem.hpp>

namespace succinctdb
{

  // we will use keys and offs field to access the hash db contents..
  template <typename Key, typename KeyIterator, typename OffIterator>
  struct __succinct_hash_iterator
  {
    typedef Key       key_type;
    typedef uint32_t  pos_type;
    typedef uint64_t  off_type;
    typedef uint64_t  hash_value_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef std::random_access_iterator_tag iterator_category;
    
    typedef __succinct_hash_iterator<Key,KeyIterator,OffIterator> __self_type;
    
    __succinct_hash_iterator(KeyIterator key_iter, OffIterator off_iter)
      : __key_iter(key_iter), __off_iter(off_iter) {}
    __succinct_hash_iterator(const __succinct_hash_iterator& x)
      : __key_iter(x.__key_iter), __off_iter(x.__off_iter) {}
    __succinct_hash_iterator()
      : __key_iter(), __off_iter() {}
    
    __self_type& operator++()
    { 
      __key_iter += *(__off_iter + 1) - *(__off_iter);
      ++ __off_iter;
      return *this;
    }
    
    __self_type& operator--()
    {
      __key_iter -= *(__off_iter) - *(__off_iter - 1);
      -- __off_iter;
      return *this;
    }

    __self_type& operator+=(difference_type __n)
    {
      __key_iter += *(__off_iter + __n) - *(__off_iter);
      __off_iter += __n;
      return *this;
    }
    
    __self_type& operator-=(difference_type __n)
    {
      __key_iter -= *(__off_iter) - *(__off_iter - __n);
      __off_iter -= __n;
      return *this;
    }
    
    __self_type operator++(int) { __self_type __tmp = *this; ++ *this; return __tmp; }
    __self_type operator--(int) { __self_type __tmp = *this; -- *this; return __tmp; }
    
    __self_type operator+(difference_type __n) const { __self_type __tmp = *this; return __tmp += __n; }
    __self_type operator-(difference_type __n) const { __self_type __tmp = *this; return __tmp -= __n; }
    
    KeyIterator begin() const { return __key_iter; }
    KeyIterator end() const { return __key_iter + size(); }
    size_type size() const { return *(__off_iter + 1) - *(__off_iter); }
    
    KeyIterator __key_iter;
    OffIterator __off_iter;
  };
  
  template <typename Key, typename KeyIterator, typename OffIterator>
  inline
  bool operator==(const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& x,
		  const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& y)
  {
    return x.__key_iter == y.__key_iter && x.__off_iter == y.__off_iter;
  }
  
  template <typename Key, typename KeyIterator, typename OffIterator>
  inline
  bool operator!=(const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& x,
		  const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& y)
  {
    return x.__key_iter != y.__key_iter || x.__off_iter != y.__off_iter;
  }
  
  template <typename Key, typename KeyIterator, typename OffIterator>
  inline
  bool operator<(const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& x,
		 const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& y)
  {
    return x.__key_iter < y.__key_iter || (! y.__key_iter < x.__key_iter &&  x.__off_iter < y.__off_iter);
  }
  
  template <typename Key, typename KeyIterator, typename OffIterator>
  inline
  bool operator>(const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& x,
		 const __succinct_hash_iterator<Key,KeyIterator,OffIterator>& y)
  {
    return y < x;
  }
  
  struct __succinct_hash_base
  {
    typedef uint32_t  pos_type;
    typedef uint64_t  off_type;
    typedef uint64_t  hash_value_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef boost::filesystem::path path_type;
    
    static const pos_type& npos() {
      static const pos_type __npos(-1);
      return __npos;
    }
  };

  template <typename Key, typename Alloc=std::allocator<Key> >
  class succinct_hash_mapped : public __succinct_hash_base
  {
  public:
    typedef Key       key_type;

  private:
    typedef __succinct_hash_base base_type;
    
  public:
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::pos_type        pos_type;
    typedef base_type::off_type        off_type;
    typedef base_type::hash_value_type hash_value_type;
    typedef base_type::path_type       path_type;
    
    using base_type::npos;
    
  private:
    typedef typename Alloc::template rebind<pos_type>::other  pos_alloc_type;
    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<off_type>::other  off_alloc_type;
    
    typedef utils::packed_vector_mapped<pos_type, pos_alloc_type> bin_set_type;
    typedef utils::packed_vector_mapped<pos_type, pos_alloc_type> next_set_type;
    typedef utils::map_file<key_type, key_alloc_type> key_set_type;
    typedef utils::vertical_coded_vector_mapped<off_type, off_alloc_type> off_set_type;

  public:
    typedef __succinct_hash_iterator<Key,typename key_set_type::const_iterator, typename off_set_type::const_iterator> iterator;
    typedef __succinct_hash_iterator<Key,typename key_set_type::const_iterator, typename off_set_type::const_iterator> const_iterator;
    
  public:
    succinct_hash_mapped(size_type __bucket_size=0)
      : bins(), nexts(), keys(), offs() {}
    succinct_hash_mapped(const path_type& path)
      : bins(), nexts(), keys(), offs() { open(path); }
    
    const_iterator begin() const
    { 
      return const_iterator(keys.begin(), offs.begin());
    }
    const_iterator end() const
    {
      return const_iterator(keys.end(), offs.end() - 1);
    }
    const_iterator operator[](size_type pos) const
    { 
      return const_iterator(keys.begin() + offs[pos], offs.begin() + pos);
    }
    size_type size() const
    {
      return nexts.size();
    }
    bool empty() const
    {
      return nexts.empty();
    }
    size_type bucket_size() const
    {
      return bins.size();
    }
    path_type path() const { return bins.path().parent_path(); }
    
    uint64_t size_bytes() const { return bins.size_bytes() + nexts.size_bytes() + keys.size_bytes() + offs.size_bytes(); }
    uint64_t size_compressed() const { return bins.size_compressed() + nexts.size_compressed() + keys.size_compressed() + offs.size_compressed(); }
    uint64_t size_cache() const { return bins.size_cache() + nexts.size_cache() + keys.size_cache() + offs.size_cache(); }
    
    pos_type find(const key_type* buf, size_type size, hash_value_type hash) const
    {
      const size_type hash_mask = bins.size() - 1;
      const size_type key = hash & hash_mask;
      pos_type i = bins[key];
      for (/**/; i && ! equal_to(i - 1, buf, size); i = nexts[i - 1]);
      return i - 1;
    }

    void close() { clear(); }
    void clear()
    {
      bins.clear();
      nexts.clear();
      keys.clear();
      offs.clear();
    }

    static bool exists(const path_type& path)
    {
      if (! utils::repository::exists(path)) return false;
      if (! bin_set_type::exists(path / "bins")) return false;
      if (! next_set_type::exists(path / "nexts")) return false;
      if (! key_set_type::exists(path / "keys")) return false;
      if (! off_set_type::exists(path / "offs")) return false;
      return true;
    }
    
    void open(const path_type& path)
    {
      typedef utils::repository repository_type;

      clear();
      
      repository_type rep(path, repository_type::read);
      bins.open(rep.path("bins"));
      nexts.open(rep.path("nexts"));
      keys.open(rep.path("keys"));
      offs.open(rep.path("offs"));
    }

    void write(const path_type& file) const
    {
      if (path() == file) return;
      
      // remove first...
      if (boost::filesystem::exists(file) && ! boost::filesystem::is_directory(file))
	boost::filesystem::remove_all(file);
      
      // create directory
      if (! boost::filesystem::exists(file))
	boost::filesystem::create_directories(file);
      
      // remove all the files...
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(file); iter != iter_end; ++ iter)
	boost::filesystem::remove_all(*iter);
      
      // copy all...
      for (boost::filesystem::directory_iterator iter(path()); iter != iter_end; ++ iter)
	utils::filesystem::copy_files(*iter, file);
    }
  
  private:
    bool equal_to(pos_type pos, const key_type* buf, size_type size) const
    {
      typename key_set_type::const_iterator first = keys.begin() + offs[pos];
      typename key_set_type::const_iterator last  = keys.begin() + offs[pos + 1];
      
      return (last - first == difference_type(size)) && std::equal(first, last, buf);
    }
    
  private:
    bin_set_type  bins;
    next_set_type nexts;
    key_set_type  keys;
    off_set_type  offs;
  };

  template <typename Key, typename Alloc=std::allocator<Key> >
  class succinct_hash_stream : public __succinct_hash_base
  {
  public:
    typedef Key       key_type;

  private:
    typedef __succinct_hash_base base_type;
    
  public:
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::pos_type        pos_type;
    typedef base_type::off_type        off_type;
    typedef base_type::hash_value_type hash_value_type;
    typedef base_type::path_type       path_type;
    
    using base_type::npos;
    
  private:
    typedef typename Alloc::template rebind<pos_type>::other  pos_alloc_type;
    typedef typename Alloc::template rebind<off_type>::other  off_alloc_type;
    
    typedef std::vector<pos_type, pos_alloc_type> bin_set_type;
    
  public:
    succinct_hash_stream() { }
    succinct_hash_stream(const path_type& path, size_type bin_size = 1024 * 1024 * 4) { open(path, bin_size); }
    ~succinct_hash_stream() { close(); }
    
    size_type size() const
    {
      return __size;
    }
    bool empty() const
    {
      return __size == 0;
    }
    size_type bucket_size() const
    {
      return bins.size();
    }
    path_type path() const { return __path; }
    
    pos_type insert(const key_type* buf, size_type size, hash_value_type hash)
    {
      const size_type hash_mask = bins.size() - 1;
      const size_type pos = hash & hash_mask;

      os_nexts->write((char*) &bins[pos], sizeof(pos_type));
      bins[pos] = __size + 1;
      ++ __size;
      
      os_keys->write((char*) buf, sizeof(key_type) * size);
      
      __offset += size;
      os_offs->write((char*) &__offset, sizeof(__offset));
      
      return __size - 1;
    }
    
    void clear() { close(); }
    void close() {
      os_nexts.reset();
      os_keys.reset();
      os_offs.reset();
      
      __size = 0;
      __offset = 0;
      
      if (! bins.empty()) {
	typedef utils::repository repository_type;
	
	repository_type rep(__path, repository_type::read);
	dump_file(rep.path("bins"), bins, true);
      }
      
      bins.clear();
      __path = path_type();
    }
    
    void open(const path_type& path, size_type bin_size)
    {
      typedef utils::repository repository_type;

      clear();
      
      
      repository_type rep(path, repository_type::write);      
      
      __path = path;
      __size = 0;
      __offset = 0;
      
      const size_type bin_size_power2 = (utils::bithack::is_power2(bin_size)
					 ? bin_size
					 : size_type(utils::bithack::next_largest_power2(bin_size)));
      
      bins.clear();
      bins.reserve(bin_size_power2);
      bins.resize(bin_size_power2, 0);
      
      os_nexts.reset(new boost::iostreams::filtering_ostream());
      os_keys.reset(new boost::iostreams::filtering_ostream());
      os_offs.reset(new boost::iostreams::filtering_ostream());
      
      //os_nexts->push(boost::iostreams::file_sink(rep.path("nexts").file_string()), 1024 * 1024);
      os_nexts->push(utils::packed_sink<pos_type, pos_alloc_type>(rep.path("nexts")));
#if BOOST_FILESYSTEM_VERSION == 2
      os_keys->push(boost::iostreams::file_sink(rep.path("keys").file_string()), 1024 * 1024);
#else
      os_keys->push(boost::iostreams::file_sink(rep.path("keys").string()), 1024 * 1024);
#endif
      os_offs->push(utils::vertical_coded_sink<off_type, off_alloc_type>(rep.path("offs")));
      
      os_keys->exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
      
      // initial offset...
      os_offs->write((char*) &__offset, sizeof(__offset));
    }

  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data, const bool packed=false) const
    {
      typedef typename _Data::value_type value_type;
      typedef typename Alloc::template rebind<value_type>::other value_alloc_type;
      
      std::auto_ptr<boost::iostreams::filtering_ostream> os(new boost::iostreams::filtering_ostream());
#if BOOST_FILESYSTEM_VERSION == 2
      if (packed)
	os->push(utils::packed_sink<value_type, value_alloc_type>(file));
      else
	os->push(boost::iostreams::file_sink(file.file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
#else
      if (packed)
	os->push(utils::packed_sink<value_type, value_alloc_type>(file));
      else
	os->push(boost::iostreams::file_sink(file.string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);	
#endif
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os->write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("error: succinct hash: write");
    }
    
  private:
    path_type __path;
    pos_type  __size;
    off_type  __offset;
    
    bin_set_type  bins;
    boost::shared_ptr<boost::iostreams::filtering_ostream> os_nexts;
    boost::shared_ptr<boost::iostreams::filtering_ostream> os_keys;
    boost::shared_ptr<boost::iostreams::filtering_ostream> os_offs;
  };
  
  template <typename Key, typename Alloc=std::allocator<Key> >
  class succinct_hash : public __succinct_hash_base
  {
  public:
    typedef Key       key_type;

  private:
    typedef __succinct_hash_base base_type;
    
  public:
    typedef base_type::size_type       size_type;
    typedef base_type::difference_type difference_type;
    typedef base_type::pos_type        pos_type;
    typedef base_type::off_type        off_type;
    typedef base_type::hash_value_type hash_value_type;
    typedef base_type::path_type       path_type;
    
    using base_type::npos;
    
  private:
    typedef typename Alloc::template rebind<pos_type>::other  pos_alloc_type;
    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<off_type>::other  off_alloc_type;
    
    typedef std::vector<pos_type, pos_alloc_type> bin_set_type;
    typedef std::vector<pos_type, pos_alloc_type> next_set_type;
    typedef std::vector<key_type, key_alloc_type> key_set_type;
    typedef utils::vertical_coded_vector<off_type, off_alloc_type> off_set_type;

  public:
    typedef __succinct_hash_iterator<Key,typename key_set_type::const_iterator, typename off_set_type::const_iterator> iterator;
    typedef __succinct_hash_iterator<Key,typename key_set_type::const_iterator, typename off_set_type::const_iterator> const_iterator;
    
  public:
    succinct_hash(size_type __bucket_size = 1024 * 1024 * 4)
      : bins(utils::bithack::is_power2(__bucket_size) ? __bucket_size : utils::bithack::next_largest_power2(__bucket_size), 0),
	nexts(), keys(), offs() { clear(); }
    

    void close() { clear(); }
    void clear()
    {
      std::fill(bins.begin(), bins.end(), 0);
      nexts.clear();
      keys.clear();
      offs.clear();
      offs.push_back(off_type(0));
    }
    
    const_iterator begin() const
    { 
      return const_iterator(keys.begin(), offs.begin());
    }
    const_iterator end() const
    {
      return const_iterator(keys.end(), offs.end() - 1);
    }
    const_iterator operator[](size_type pos) const
    { 
      return const_iterator(keys.begin() + offs[pos], offs.begin() + pos);
    }
    size_type size() const
    {
      return nexts.size();
    }
    bool empty() const
    {
      return nexts.empty();
    }
    size_type bucket_size() const
    {
      return bins.size();
    }

    uint64_t size_bytes() const
    {
      return bins.size() * sizeof(pos_type) + nexts.size() * sizeof(pos_type) + keys.size() * sizeof(key_type) + offs.size_bytes();
    }
    uint64_t size_compressed() const
    {
      return bins.size() * sizeof(pos_type) + nexts.size() * sizeof(pos_type) + keys.size() * sizeof(key_type) + offs.size_compressed();
    }
    uint64_t size_cache() const
    {
      return offs.size_cache();
    }
    
    pos_type insert(const key_type* buf, size_type size, hash_value_type hash)
    {
      const size_type hash_mask = bins.size() - 1;
      const size_type key = hash & hash_mask;
      pos_type i = bins[key];
      for (/**/; i && ! equal_to(i - 1, buf, size); i = nexts[i - 1]);
      if (i)
	return i - 1;
      
      i = nexts.size() + 1;
      
      keys.insert(keys.end(), buf, buf + size);
      offs.push_back(keys.size());
      nexts.push_back(bins[key]);
      
      bins[key] = i;
      return i - 1;
    }
    
    pos_type find(const key_type* buf, size_type size, hash_value_type hash) const
    {
      const size_type hash_mask = bins.size() - 1;
      const size_type key = hash & hash_mask;
      pos_type i = bins[key];
      for (/**/; i && ! equal_to(i - 1, buf, size); i = nexts[i - 1]);
      return i - 1;
    }
  
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      rep["type"] = "succinct-hash";
      dump_file(rep.path("bins"), bins, true);
      dump_file(rep.path("nexts"), nexts, true);
      dump_file(rep.path("keys"), keys, false);
      offs.write(rep.path("offs"));
    }
    
  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data, const bool packed=false) const
    {
      typedef typename _Data::value_type value_type;
      typedef typename Alloc::template rebind<value_type>::other value_alloc_type;
      
      std::auto_ptr<boost::iostreams::filtering_ostream> os(new boost::iostreams::filtering_ostream());
#if BOOST_FILESYSTEM_VERSION == 2
      if (packed)
	os->push(utils::packed_sink<value_type, value_alloc_type>(file));
      else
	os->push(boost::iostreams::file_sink(file.file_string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
#else
      if (packed)
	os->push(utils::packed_sink<value_type, value_alloc_type>(file));
      else
	os->push(boost::iostreams::file_sink(file.string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
#endif
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024) {
	if (! os->write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("error: succinct hash: write");
      }
    }
    
    bool equal_to(pos_type pos, const key_type* buf, size_type size) const
    {
      typename key_set_type::const_iterator first = keys.begin() + offs[pos];
      typename key_set_type::const_iterator last  = keys.begin() + offs[pos + 1];
      
      return (last - first == difference_type(size)) && std::equal(first, last, buf);
    }
    
  private:
    bin_set_type  bins;
    next_set_type nexts;
    key_set_type  keys;
    off_set_type  offs;
  };
  
};

#endif
