// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __SUCCINCT_DB__SUCCINCT_HASH_DB__HPP__
#define __SUCCINCT_DB__SUCCINCT_HASH_DB__HPP__ 1

#include <stdint.h>

#include <vector>
#include <stdexcept>

#include <succinct_db/succinct_hash.hpp>

#include <boost/shared_ptr.hpp>

#include <utils/repository.hpp>
#include <utils/map_file.hpp>
#include <utils/hashmurmur.hpp>

namespace succinctdb
{
  
  template <typename Key, typename Data, typename Hasher=utils::hashmurmur<uint64_t>, typename Alloc=std::allocator<std::pair<Key, Data> > >
  class succinct_hash_db : __succinct_hash_base
  {
  public:
    typedef Key    key_type;
    typedef Data   data_type;
    typedef Data   mapped_type;
    typedef Hasher hasher_type;
    
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
    typedef typename Alloc::template rebind<key_type>::other  key_alloc_type;
    typedef typename Alloc::template rebind<data_type>::other data_alloc_type;
    
    typedef succinct_hash<key_type, key_alloc_type >       succinct_hash_type;
    typedef succinct_hash_mapped<key_type, key_alloc_type> succinct_hash_mapped_type;
    typedef succinct_hash_stream<key_type, key_alloc_type> succinct_hash_stream_type;
    
    typedef std::vector<data_type, data_alloc_type>            data_set_type;
    typedef utils::map_file<data_type, data_alloc_type> data_set_mapped_type;
    
  public:
    succinct_hash_db(size_type bin_size=1024 * 1024 * 4) { open(bin_size); }
    succinct_hash_db(const path_type& path, size_type bin_size=0) { open(path, bin_size); }

  public:
    path_type path() const
    {
      if (__succinct_hash_mapped)
	return __succinct_hash_mapped->path().parent_path();
      else if (__succinct_hash_stream)
	return __succinct_hash_stream->path().parent_path();
      else
	return path_type();
    }

    size_type size() const 
    {
      if (__succinct_hash)
	return __succinct_hash->size();
      else if (__succinct_hash_mapped)
	return __succinct_hash_mapped->size();
      else if (__succinct_hash_stream)
	return __succinct_hash_stream->size();
      else
	return 0;
    }
    
    bool empty() const
    {
      if (__succinct_hash)
	return __succinct_hash->empty();
      else if (__succinct_hash_mapped)
	return __succinct_hash_mapped->empty();
      else if (__succinct_hash_stream)
	return __succinct_hash_stream->emptry();
      else
	return true;
    }

    uint64_t size_bytes() const
    {
      return ((__succinct_hash_mapped ? __succinct_hash_mapped->size_bytes() : uint64_t(0))
	      + (__data_mapped ? __data_mapped->size_bytes() : uint64_t(0)));
    }
    uint64_t size_compressed() const
    {
      return ((__succinct_hash_mapped ? __succinct_hash_mapped->size_compressed() : uint64_t(0))
	      + (__data_mapped ? __data_mapped->size_compressed() : uint64_t(0)));
    }
    uint64_t size_cache() const
    {
      return ((__succinct_hash_mapped ? __succinct_hash_mapped->size_cache() : uint64_t(0))
	      + (__data_mapped ? __data_mapped->size_cache() : uint64_t(0)));
    }
    
  public:
    void close() { clear(); }
    void clear() 
    {
      __succinct_hash.reset();
      __succinct_hash_mapped.reset();
      __succinct_hash_stream.reset();
      
      __data.reset();
      __data_mapped.reset();
      __data_stream.reset();
    }
    
    void open(size_type bin_size)
    {
      clear();
      
      __succinct_hash->reset(new succinct_hash_type(bin_size));
      __data->reset(new data_set_type());
    }

    static bool exists(const path_type& path)
    {
      if (! utils::repository::exists(path)) return false;
      if (! succinct_hash_mapped_type::exists(path / "index")) return false;
      if (! data_set_mapped_type::exists(path / "data")) return false;
      return true;
    }
    
    void open(const path_type& path, size_type bin_size=0)
    {
      typedef utils::repository repository_type;
      
      clear();
      
      if (bin_size > 0) {
	// open stream...
	repository_type rep(path, repository_type::write);
	rep["type"] = "succinct-hash-db";
	
	__succinct_hash_stream.reset(new succinct_hash_stream_type(rep.path("index"), bin_size));
	__data_stream.reset(new boost::iostreams::filtering_ostream());
	__data_stream->push(boost::iostreams::file_sink(rep.path("data").string()), 1024 * 1024);
	__data_stream->exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
      } else {
	// read-only open
	repository_type rep(path, repository_type::read);
	
	__succinct_hash_mapped.reset(new succinct_hash_mapped_type(rep.path("index")));
	__data_mapped.reset(new data_set_mapped_type(rep.path("data")));
      }
    }
    
    void write(const path_type& file) const
    {
      typedef utils::repository repository_type;
      
      if (__succinct_hash) {
	repository_type rep(file, repository_type::write);
	rep["type"] = "succinct-hash-db";
	
	__succinct_hash->write(rep.path("index"));
	dump_file(rep.path("data"), *__data);
      } else if (__succinct_hash_mapped) {
	if (path() != file) {
	  repository_type rep(file, repository_type::write);
	  rep["type"] = "succinct-hash-db";
	  
	  __succinct_hash_mapped->write(rep.path("index"));
	  __data_mapped->write(rep.path("data"));
	}
      }
    }
    
    pos_type insert(const key_type* buf, size_type buf_size, const data_type& data)
    {
      return insert(buf, buf_size, &data);
    }

    pos_type insert(const key_type* buf, size_type buf_size, const data_type* data)
    {
      if (__succinct_hash) {
	const pos_type pos = __succinct_hash->insert(buf, buf_size, __hasher(buf, buf + buf_size, 0));
	if (pos >= __data->size())
	  __data->resize(pos + 1);
	__data->operator[](pos) = *data;
	return pos;
      } else if (__succinct_hash_stream) {
	__data_stream->write((char*) data, sizeof(data_type));
	return __succinct_hash_stream->insert(buf, buf_size, __hasher(buf, buf + buf_size, 0));
      } else
	throw std::runtime_error("we cannot insert...");
    }
    
    pos_type find(const key_type* buf, size_type buf_size) const
    {
      if (__succinct_hash)
	return __succinct_hash->find(buf, buf_size, __hasher(buf, buf + buf_size, 0));
      else if (__succinct_hash_mapped)
	return __succinct_hash_mapped->find(buf, buf_size, __hasher(buf, buf + buf_size, 0));
      else
	throw std::runtime_error("we cannot find...");
    }
    
    data_type operator[](pos_type pos) const
    {
      if (__data)
	return __data->operator[](pos);
      else if (__data_mapped)
	return __data_mapped->operator[](pos);
      else
	throw std::runtime_error("we cannot access...");
    }

  private:
    template <typename _Path, typename _Data>
    inline
    void dump_file(const _Path& file, const _Data& data) const
    {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::file_sink(file.string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
      
      const int64_t file_size = sizeof(typename _Data::value_type) * data.size();
      for (int64_t offset = 0; offset < file_size; offset += 1024 * 1024)
	if (! os.write(((char*) &(*data.begin())) + offset, std::min(int64_t(1024 * 1024), file_size - offset)))
	  throw std::runtime_error("succinct_hash_db::write()");
    }
    
  private:
    hasher_type __hasher;
    
    boost::shared_ptr<succinct_hash_type>        __succinct_hash;
    boost::shared_ptr<succinct_hash_mapped_type> __succinct_hash_mapped;
    boost::shared_ptr<succinct_hash_stream_type> __succinct_hash_stream;
    
    boost::shared_ptr<data_set_type>                       __data;
    boost::shared_ptr<data_set_mapped_type>                __data_mapped;
    boost::shared_ptr<boost::iostreams::filtering_ostream> __data_stream;
  };
  
};

#endif
