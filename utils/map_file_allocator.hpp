// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MAP_FILE_ALLOCATOR__HPP__
#define __UTILS__MAP_FILE_ALLOCATOR__HPP__ 1

//
// TODO: very strange behavior...
// use of pthread_once to initialize global storage
// use of boost::shared_ptr for safety
//

#include <cstdlib>
#include <cstring>

#include <new>
#include <stdexcept>
#include <iostream>
#include <string>
#include <algorithm>
#include <set>
#include <map>

#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/device/file.hpp>
#include <boost/functional/hash.hpp>

#include <utils/tempfile.hpp>
#include <utils/config.hpp>
#include <utils/spinlock.hpp>
#include <utils/sgi_hash_map.hpp>

namespace utils
{
  class __map_alloc_file
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef char      byte_type;

    typedef boost::filesystem::path path_type;
    
    typedef byte_type* pointer;
    
  private:
    pointer   mapped;
    size_type file_size;
    path_type file;
    
  public:
    __map_alloc_file(const size_t size) : mapped(0), file_size(0), file() { open(size); }
    ~__map_alloc_file() throw() { close(); }
    
  private:
    __map_alloc_file(const __map_alloc_file& x) {}
    __map_alloc_file& operator=(const __map_alloc_file& x) { return *this; }
    
  public:
    pointer begin() { return mapped; }
    pointer end() { return mapped + file_size; }
    
    void close()
    {
      if (mapped && file_size > 0)
	::munmap(mapped, file_size);
      if (! file.empty() && boost::filesystem::exists(file)) {
	boost::filesystem::remove(file);
	utils::tempfile::erase(file);
      }
      
      file_size = 0;
      mapped = 0;
      file = path_type();
    }
    
    void open(const size_type size)
    {
      // create a temporary file
      const path_type tmp_dir = utils::tempfile::tmp_dir();
      path_type file_template = tmp_dir / "map_file_allocator.XXXXXX";
      file = utils::tempfile::file_name(file_template);
      utils::tempfile::insert(file);
      
      // create enough size...
      const size_type alloc_size = ((size + (4096 - 1)) / 4096) * 4096;
      {
	boost::iostreams::filtering_ostream os;
	os.push(boost::iostreams::file_sink(file.file_string()), 1024 * 1024 * 4);
	
	std::vector<byte_type> buffer(4096, 0);
	for (size_type i = 0; i < alloc_size / buffer.size(); ++ i)
	  os.write(&(*buffer.begin()), buffer.size());
      }
      
      file_size = boost::filesystem::file_size(file);
      
      int fd = ::open(file.file_string().c_str(), O_RDWR | O_NDELAY);
      if (fd < 0) {
	boost::filesystem::remove(file);
	utils::tempfile::erase(file);
	throw std::runtime_error("map_file_allocator:: open()");
      }
      
      byte_type* x = (byte_type*) mmap(0, file_size, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
      ::close(fd); // we do not have to keep this...
      if (x + 1 == 0) {
	close();
	throw std::runtime_error("map_file_allocator:: mmap()");
      }
      mapped = x;
    }
  };
  
  struct __map_alloc_impl
  {
    typedef boost::mutex            mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    typedef __map_alloc_file        map_file_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<void*, map_file_type*, boost::hash<void*>, std::equal_to<void*>,
				    std::allocator<std::pair<const void*, map_file_type*> > > mapped_type;
#else
    typedef sgi::hash_map<void*, map_file_type*, boost::hash<void*>, std::equal_to<void*>,
			  std::allocator<std::pair<const void*, map_file_type*> > > mapped_type;
#endif
  
    typedef std::set<void*, std::less<void*>, std::allocator<void*> > removed_type;
    

  public:
    class global_alloc
    {
    public:
      global_alloc() {}
      ~global_alloc() { clear(); }
      
    public:
      void clear()
      {
	lock_type lock(mutex);
	
	for (mapped_type::iterator iter = mapped.begin(); iter != mapped.end(); ++ iter)
	  delete iter->second;
	mapped.clear();
	removed.clear();
      }	

      
      void insert(map_file_type* file)
      {
	{
	  lock_type lock(mutex);
	  
	  mapped.insert(std::make_pair(static_cast<void*>(file->begin()), file));
	}
	prune();
      }
      
      void erase(void* p)
      {
	{
	  lock_type lock(mutex);
	  
	  removed.insert(p);
	}
	prune();
      }
      
      void prune()
      {
	lock_type lock(mutex);
	
	for (removed_type::iterator riter = removed.begin(); riter != removed.end(); /**/) {
	  mapped_type::iterator iter = mapped.find(*riter);
	  if (iter != mapped.end()) {
	    delete iter->second;
	    mapped.erase(iter);
	    removed.erase(riter ++);
	  } else
	    ++ riter;
	}
      }      
      
    private:
      static mutex_type mutex;

    private:
      mapped_type mapped;
      removed_type removed;
    };
    
    mapped_type     mapped;
    utils::tempfile __tempfile;
    static global_alloc __global_alloc;
    
  public:    
    __map_alloc_impl() : mapped() {}
    ~__map_alloc_impl() throw() { clear(); }
    
  private:
    __map_alloc_impl(const __map_alloc_impl& x) {}
    __map_alloc_impl& operator=(const __map_alloc_impl& x) { return *this; }
    
  public:
    
    void clear()
    {
      for (mapped_type::iterator iter = mapped.begin(); iter != mapped.end(); ++ iter)
	__global_alloc.insert(iter->second);
      __global_alloc.prune();
      
      mapped.clear();      
    }
    
    void* allocate(size_t __n) 
    {
      if (__n == 0) return 0;
      
      map_file_type* file = new map_file_type(__n);
      
      mapped.insert(std::make_pair(static_cast<void*>(file->begin()), file));
      
      return static_cast<void*>(file->begin());
    }
    
    void deallocate(void* __p, size_t __n) 
    {
      if (__p) {
	mapped_type::iterator iter = mapped.find(__p);
	if (iter != mapped.end()) {
	  delete iter->second;
	  mapped.erase(iter);
	} else {
	  __global_alloc.erase(__p);
	  __global_alloc.prune();
	}
      }
    }
  };
  
  struct __map_file_allocator_base
  {
    typedef __map_alloc_impl          map_alloc_type;
    typedef __map_file_allocator_base base_type;
    
#ifdef HAVE_TLS
    static __thread map_alloc_type* local_alloc_thread;
#endif
    static boost::thread_specific_ptr<map_alloc_type> local_alloc;
    
    map_alloc_type& map_alloc()
    {
#ifdef HAVE_TLS
      if (! base_type::local_alloc_thread) {
	base_type::local_alloc.reset(new map_alloc_type());
	base_type::local_alloc_thread = local_alloc.get();
      }
      
      return *base_type::local_alloc_thread;
#else
      if (! base_type::local_alloc.get())
	base_type::local_alloc.reset(new map_alloc_type());
      return *base_type::local_alloc;
#endif
    }
    
  };
  
  template <typename _Tp, typename _Alloc=std::allocator<_Tp> >
  struct map_file_allocator : public __map_file_allocator_base,
			      public _Alloc
  {
    
  public:
    typedef size_t     size_type;
    typedef ptrdiff_t  difference_type;
    typedef _Tp*       pointer;
    typedef const _Tp* const_pointer;
    typedef _Tp&       reference;
    typedef const _Tp& const_reference;
    typedef _Tp        value_type;
    
    template<typename _Tp1>
    struct rebind
    { typedef map_file_allocator<_Tp1> other; };
    
  public:
    map_file_allocator() throw() { }
    
    map_file_allocator(const map_file_allocator&) throw() { }
    
    template<typename _Tp1, typename _Alloc1>
    map_file_allocator(const map_file_allocator<_Tp1, _Alloc1>&) throw() { }
    
    ~map_file_allocator() throw() { }
    
  public:
    pointer
    address(reference __x) const { return &__x; }
    
    const_pointer
    address(const_reference __x) const { return &__x; }
    
    size_type
    max_size() const throw() 
    { return size_t(-1) / sizeof(_Tp); }
    
    void 
    construct(pointer __p, const _Tp& __val) 
    { ::new(__p) _Tp(__val); }
    
    void 
    destroy(pointer __p) { __p->~_Tp(); }
    
        
    pointer
    allocate(size_type __n, const void* tmptmp= 0)
    {
      if (__n == 0) return 0;

      if (sizeof(_Tp) * __n < 1024 * 1024)
	return base_allocator().allocate(__n);
      else
	return static_cast<pointer>(map_alloc().allocate(sizeof(_Tp) * __n));
    }
    
    void
    deallocate(pointer __p, size_type __n) {
      if (! __p) return;
      
      if (sizeof(_Tp) * __n < 1024 * 1024)
	base_allocator().deallocate(__p, __n);
      else
	map_alloc().deallocate(static_cast<void*>(__p), sizeof(_Tp) * __n);
    }
    
  private:
    _Alloc& base_allocator()
    {
      return static_cast<_Alloc&>(*this);
    }
    
  };

  template <typename _Tp, typename _Alloc>
  inline bool operator==(const map_file_allocator<_Tp,_Alloc>& x,
                         const map_file_allocator<_Tp,_Alloc>& y)
  { return true; }
  
  template <typename _Tp, typename _Alloc>
  inline bool operator!=(const map_file_allocator<_Tp,_Alloc>& x,
                         const map_file_allocator<_Tp,_Alloc>& y)
  { return false; }

};

#endif
