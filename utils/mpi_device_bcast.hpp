// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_DEVICE_BCAST__HPP__
#define __UTILS__MPI_DEVICE_BCAST__HPP__ 1

// mpi device, but perform bcasting...
//
// remember, this resulted in synchronous access...
//
// we will re-implement this...

#include <string>
#include <algorithm>
#include <vector>
#include <iosfwd>

#include <boost/iostreams/categories.hpp>
#include <boost/shared_ptr.hpp>

#include <mpi.h>
#include <utils/mpi_allocator.hpp>

namespace utils
{
  
  class mpi_device_bcast_source
  {
  public:
    typedef char char_type;
    struct category : public boost::iostreams::source_tag,
		      public boost::iostreams::closable_tag {};
    
    mpi_device_bcast_source(const int rank, const size_t buffer_size=4096) : pimpl(new impl()) { open(rank, buffer_size); }
    mpi_device_bcast_source() : pimpl(new impl()) { }
    
    std::streamsize read(char_type* s, std::streamsize n)
    {
      std::streamsize __n = 0;
      std::streamsize available = n;
      while (available) {
	std::streamsize __read = pimpl->read(s, available);
	if (__read == -1)
	  return (__n > 0 ? __n : -1);
	__n += __read;
	s += __read;
	available -= __read;
      }
      return __n;
    }
    void open(const int rank, const size_t buffer_size=4096) { pimpl->open(rank, buffer_size); }
    bool is_open() const { return pimpl->is_open(); }
    void close() { pimpl->close(); }
    
  private:
    struct impl
    {
      int root;
      volatile unsigned int recv_size;
      volatile size_t buffer_offset;
      std::vector<char_type, std::allocator<char_type> > buffer;
      
      impl() : root(-1), buffer_offset(0) {}
      ~impl() { close(); }
      
      std::streamsize read(char_type* s, std::streamsize n);
      void open(const int rank, const size_t buffer_size=4096);
      bool is_open() const;
      void close();
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  
  class mpi_device_bcast_sink
  {
    public:
    typedef char char_type;
    struct category : public boost::iostreams::sink_tag,
		      public boost::iostreams::closable_tag {};
    
    mpi_device_bcast_sink(const int rank, const size_t buffer_size=4096) : pimpl(new impl()) { open(rank, buffer_size); }
    mpi_device_bcast_sink() : pimpl(new impl()) { }
    
    std::streamsize write(const char_type* s, std::streamsize n)
    {
      std::streamsize __n = 0;
      std::streamsize available = n;
      while (available) {
	const std::streamsize written = pimpl->write(s, available);
	if (written == -1)
	  return (__n > 0 ? __n : -1);
	available -= written;
	s += written;
	__n += written;
      }
      return __n;
    }
    void open(const int rank, const size_t buffer_size=4096) { pimpl->open(rank, buffer_size); }
    bool is_open() const { return pimpl->is_open(); }
    void close() { pimpl->close(); }
    
  private:
    struct impl
    {
      int root;
      volatile unsigned int send_size;
      volatile size_t buffer_offset;
      std::vector<char_type, std::allocator<char_type> > buffer;
      
      impl() : root(-1), buffer_offset(0) {}
      ~impl() { close(); }
      
      std::streamsize write(const char_type* s, std::streamsize n);
      void open(const int rank, const size_t buffer_size=4096);
      bool is_open() const;
      void close();
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  bool mpi_device_bcast_sink::impl::is_open() const
  {
    return root >= 0;
  }
  
  bool mpi_device_bcast_source::impl::is_open() const
  {
    return root >= 0;
  }
  
  
  std::streamsize mpi_device_bcast_sink::impl::write(const char_type* s, std::streamsize n)
  {
    const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - buffer_offset), n);
    std::copy(s, s + copy_size, buffer.begin() + buffer_offset);
    buffer_offset += copy_size;
    
    if (buffer_offset == buffer.size()) {

      send_size = buffer_offset;
      MPI::COMM_WORLD.Bcast((void*) &send_size, 1, MPI::INT, root);
      MPI::COMM_WORLD.Bcast(&(*buffer.begin()), send_size, MPI::CHAR, root);

      buffer_offset = 0;
    }
    return copy_size;
  }
  
  void mpi_device_bcast_sink::impl::close()
  {
    if (buffer_offset > 0) {
      send_size = buffer_offset;
      MPI::COMM_WORLD.Bcast((void*) &send_size, 1, MPI::INT, root);
      MPI::COMM_WORLD.Bcast(&(*buffer.begin()), send_size, MPI::CHAR, root);
    }
    
    if (root >= 0) {
      send_size = 0;
      MPI::COMM_WORLD.Bcast((void*) &send_size, 1, MPI::INT, root);
    }
    
    root = -1;
    buffer_offset = 0;
  }
  
  void mpi_device_bcast_sink::impl::open(const int rank, const size_t buffer_size)
  {
    close();
    
    buffer.reserve(std::max(buffer_size, size_t(64)));
    buffer.resize(std::max(buffer_size, size_t(64)));
    root = rank;
  }
  
  
  std::streamsize mpi_device_bcast_source::impl::read(char_type* s, std::streamsize n)
  {
    if (root < 0) return -1;

    if (buffer_offset == buffer.size()) {
      
      MPI::COMM_WORLD.Bcast((void*) &recv_size, 1, MPI::INT, root);
      
      if (recv_size == 0) {
	close();
	return -1;
      }
      
      buffer.reserve(recv_size);
      buffer.resize(recv_size);
      MPI::COMM_WORLD.Bcast(&(*buffer.begin()), recv_size, MPI::CHAR, root);
      
      buffer_offset = 0;
    }
    
    const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - buffer_offset), n);
    std::copy(buffer.begin() + buffer_offset, buffer.begin() + buffer_offset + copy_size, s);
    buffer_offset += copy_size;
    
    return copy_size;
  };
  
  void mpi_device_bcast_source::impl::close()
  {
    buffer.clear();
    buffer_offset = 0;
    root = -1;
  }
  
  void mpi_device_bcast_source::impl::open(const int rank, const size_t buffer_size)
  {
    close();
    root = rank;
  }

};

#endif
