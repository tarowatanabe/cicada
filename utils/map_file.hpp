// -*- mode: c++ -*-

#ifndef __UTILS__MAP_FILE__HPP__
#define __UTILS__MAP_FILE__HPP__ 1

#include <stdint.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>

#include <string>
#include <iostream>
#include <stdexcept>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/filesystem.hpp>

#include <utils/config.hpp>
#include <utils/filesystem.hpp>
#include <utils/compress_stream.hpp>

namespace utils
{
  class __map_file_impl
  {
  public:
    typedef std::streamoff off_type;
    typedef size_t         size_type;
    typedef char           byte_type;
    
    typedef boost::filesystem::path path_type;
    
  public:
    __map_file_impl(const std::string& file, const bool writable=false)
      : mmapped(), filesize(), mmap_size(), filename() { open(file, writable); }
    __map_file_impl(const boost::filesystem::path& file, const bool writable=false)
      : mmapped(), filesize(), mmap_size(), filename() { open(file, writable); }
    ~__map_file_impl() { close(); }
    
  public:
    bool is_open() const { return mmapped; }
    
    const void* begin() const { return static_cast<void*>(mmapped); }
    const void* end() const { return static_cast<void*>(mmapped + filesize); }

    void* begin() { return static_cast<void*>(mmapped); }
    void* end() { return static_cast<void*>(mmapped + filesize); }
    
    off_type size() const { return filesize; }
    const boost::filesystem::path& path() const { return filename; }
    
  private:
    void open(const boost::filesystem::path& file, const bool writable=false)
    {
      filename = file;
      filesize = static_cast<off_type>(boost::filesystem::file_size(file));
      modifiable = writable;
      
#if BOOST_FILESYSTEM_VERSION == 2
      int fd = ::open(file.file_string().c_str(), (writable ? O_RDWR | O_NDELAY : O_RDONLY | O_NDELAY));
      if (fd < 0)
	fd = ::open(file.file_string().c_str(), (writable ? O_RDWR : O_RDONLY));
#else
      int fd = ::open(file.string().c_str(), (writable ? O_RDWR | O_NDELAY : O_RDONLY | O_NDELAY));
      if (fd < 0)
	fd = ::open(file.string().c_str(), (writable ? O_RDWR : O_RDONLY));
#endif

      if (fd < 0)
	throw std::runtime_error("map_file::open() open()");
      
      const size_t page_size = getpagesize();
      mmap_size = static_cast<off_type>(std::max(off_type((filesize + page_size - 1) / page_size), off_type(1)) * page_size);
	
      // First, try map_shared
      byte_type* x = static_cast<byte_type*>(::mmap(0, mmap_size, writable ? PROT_WRITE : PROT_READ, MAP_SHARED, fd, 0));
      
      // Second, try map_private
      if (! (x + 1))
	x = static_cast<byte_type*>(::mmap(0, mmap_size, writable ? PROT_WRITE : PROT_READ, MAP_PRIVATE, fd, 0));
	
      // no need to keep file-descriptor
      ::close(fd);
      
      // If successful, use mmap. Otherwise pread|seek/read access...
      if (x + 1)
	mmapped = x;
      else
	throw std::runtime_error(std::string("map_file::open() mmap()") + strerror(errno));
    }

    void close()
    {
      if (mmapped && filesize > 0 && mmap_size > 0)
	::munmap(mmapped, mmap_size);
      
      mmapped = 0;
      filesize = 0;
      mmap_size = 0;
      filename = std::string();
      
      modifiable = false;
    }
      
  private:
    byte_type* mmapped;
    
    off_type filesize;
    off_type mmap_size;
    
    boost::filesystem::path filename;

    bool modifiable;
  };
  
  template <typename _Tp, typename _Alloc=std::allocator<_Tp> >
  class map_file
  {
  private:
    typedef __map_file_impl impl_type;
    
  public:
    typedef _Tp              value_type;
    typedef const _Tp*       const_iterator;
    typedef const_iterator   iterator;
    typedef const _Tp&       const_reference;
    typedef const_reference  reference;

    typedef impl_type::off_type  off_type;
    typedef impl_type::size_type size_type;
    typedef impl_type::byte_type byte_type;
    typedef impl_type::path_type path_type;
    
  public:
    // we will allow copying, but actually, this will not be copied...
    
    map_file(const path_type& file, const bool writable=false) : pimpl(new impl_type(file, writable)) { }
    map_file() {}
    
  public:
    const_reference front() const { return *(begin()); }
    const_reference back() const { return *(end() - 1); }
    const_reference operator[](size_type pos) const { return *(begin() + pos); }
    
    const_iterator begin() const { return static_cast<const_iterator>(pimpl->begin()); }
    const_iterator end()   const { return static_cast<const_iterator>(pimpl->end()); }

    iterator begin() { return static_cast<iterator>(pimpl->begin()); }
    iterator end()   { return static_cast<iterator>(pimpl->end()); }
    
    bool is_open() const { return pimpl && pimpl->is_open(); }
    bool empty() const { return ! is_open() || size() == 0; }
    size_type size() const { return static_cast<size_type>(pimpl->size() / sizeof(value_type)); }
    off_type file_size() const { return pimpl->size(); }

    uint64_t size_bytes() const { return file_size(); }
    uint64_t size_compressed() const { return file_size(); }
    uint64_t size_cache() const { return 0; }
    
    path_type path() const { return pimpl->path(); }

    static bool exists(const path_type& path)
    {
      return boost::filesystem::exists(path);
    }
    
    void open(const std::string& file, const bool writable=false) { pimpl.reset(new impl_type(file, writable)); }
    void open(const path_type& file, const bool writable=false) { pimpl.reset(new impl_type(file, writable)); }
    
    void write(const path_type& file) const
    {
      if (file == path()) return;
      utils::filesystem::copy_file(path(), file);
    }
    
    void close() { pimpl.reset(); }
    void clear() { close(); }
    
    void swap(map_file& x)
    {
      pimpl.swap(x.pimpl);
    }

  private:
    boost::shared_ptr<impl_type> pimpl;
  };
};

namespace std
{
  template <typename T, typename A>
  inline
  void swap(utils::map_file<T,A>& x, utils::map_file<T,A>& y)
  {
    x.swap(y);
  }
};

#endif
