// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_DEVICE__HPP__
#define __UTILS__MPI_DEVICE__HPP__ 1

// we will re-implement this...

#include <string>
#include <algorithm>
#include <vector>
#include <iosfwd>

#include <boost/iostreams/categories.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/thread.hpp>

#include <mpi.h>

#include <utils/mpi_allocator.hpp>
#include <utils/atomicop.hpp>

namespace utils
{
  template <typename Stream, typename Device>
  inline
  bool mpi_terminate_device(Stream& stream, Device& device)
  {
    bool found = false;
    if (! stream && device && device->test()) {
      if (! device->test_terminate()) {
	if (device->flush() == 0 && device->test())
	  device->terminate();
      } else {
	device->finalize();
	device.reset();
      }
      found = true;
    }
    return found;
  }
  
  template <typename Streams, typename Devices>
  inline
  bool mpi_terminate_devices(Streams& streams, Devices& devices)
  {
    bool found = false;
    
    typename Streams::iterator siter_end = streams.end();
    typename Devices::iterator diter = devices.begin();
    for (typename Streams::iterator siter = streams.begin(); siter != siter_end; ++ siter, ++ diter)
      if (! *siter && *diter && (*diter)->test()) {
	if (! (*diter)->test_terminate()) {
	  if ((*diter)->flush() == 0 && (*diter)->test())
	    (*diter)->terminate();
	} else {
	  (*diter)->finalize();
	  diter->reset();
	}
	found = true;
      }
    return found;
  }
  
  template <typename Streams, typename Devices>
  inline
  bool mpi_flush_device(Streams& stream, Devices& device)
  {
    // do not send too much...
    bool busy = false;
    if (stream && device) {
      device->flush(true);
      busy |=  (! device->test());
    }
    return busy;
  }
  
  template <typename Streams, typename Devices>
  inline
  bool mpi_flush_devices(Streams& streams, Devices& devices)
  {
    // do not send too much...
    // check whether all busy...
    
    bool non_busy = false;
    typename Streams::iterator siter_end = streams.end();
    typename Devices::iterator diter = devices.begin();
    for (typename Streams::iterator siter = streams.begin(); siter != siter_end; ++ siter, ++ diter)
      if (*siter && *diter) {
	(*diter)->flush(true);
	non_busy |= (*diter)->test();
      }
    return ! non_busy;
  }

  

  class mpi_device_source
  {
  public:
    typedef char char_type;
    struct category : public boost::iostreams::source_tag,
		      public boost::iostreams::closable_tag {};
    
    mpi_device_source(MPI::Comm& comm,
		      int rank,
		      int tag,
		      size_t buffer_size=4096)
      : pimpl(new impl()) { open(comm, rank, tag, buffer_size); }
    mpi_device_source(int rank,
		      int tag,
		      size_t buffer_size=4096)
      : pimpl(new impl()) { open(rank, tag, buffer_size); }
    mpi_device_source()
      : pimpl(new impl()) { }
    
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
    void open(int rank, int tag, size_t buffer_size=4096) { pimpl->open(rank, tag, buffer_size); }
    void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096) { pimpl->open(comm, rank, tag, buffer_size); }
    bool is_open() const { return pimpl->is_open(); }
    void close() { pimpl->close(); }
    bool test() const { return pimpl->test(); }
    void wait() const { pimpl->wait(); }
    
  private:
    struct impl
    {
      MPI::Prequest request_size;
      MPI::Prequest request_buffer;
      volatile unsigned int recv_size;
      volatile size_t buffer_offset;
      std::vector<char_type, std::allocator<char_type> > buffer;
      
      impl() {}
      ~impl() { close(); }
      
      std::streamsize read(char_type* s, std::streamsize n);
      void open(int rank, int tag, size_t buffer_size=4096);
      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096);
      bool is_open() const;
      void close();
      bool test() const;
      void wait() const;
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  

  class mpi_device_sink
  {
    public:
    typedef char char_type;
    struct category : public boost::iostreams::sink_tag,
		      public boost::iostreams::closable_tag {};
    
    mpi_device_sink(MPI::Comm& comm,
		    int rank, 
		    int tag, 
		    size_t buffer_size=4096, 
		    bool terminate_on_close=true,
		    bool overcommit=false) 
      : pimpl(new impl()) { open(comm, rank, tag, buffer_size, terminate_on_close, overcommit); }
    mpi_device_sink(int rank,
		    int tag, 
		    size_t buffer_size=4096, 
		    bool terminate_on_close=true,
		    bool overcommit=false) 
      : pimpl(new impl()) { open(rank, tag, buffer_size, terminate_on_close, overcommit); }
    mpi_device_sink() 
      : pimpl(new impl()) { }
    
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
    void open(int rank, int tag, size_t buffer_size=4096, bool terminate_on_close=true, bool overcommit=false)
    {
      pimpl->open(rank, tag, buffer_size, terminate_on_close, overcommit);
    }
    void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096, bool terminate_on_close=true, bool overcommit=false)
    {
      pimpl->open(comm, rank, tag, buffer_size, terminate_on_close, overcommit);
    }
    bool is_open() const { return pimpl->is_open(); }
    void close() { pimpl->close(); }
    bool test() const { return pimpl->test(); }
    void wait() const { pimpl->wait(); }
    std::streamsize flush(const bool if_filled=false) { return pimpl->flush(if_filled); }
    void terminate() { pimpl->terminate(); }
    void finalize() { pimpl->finalize(); }
    bool test_terminate() { return pimpl->test_terminate(); }
    
  private:
    struct impl
    {
      typedef std::vector<char, std::allocator<char_type> > buffer_type;
      
      MPI::Prequest request_size;
      MPI::Prequest request_buffer;
      volatile unsigned int send_size;
      volatile size_t buffer_offset;
      buffer_type buffer;
      buffer_type buffer_overcommit;
      bool terminate_on_close;
      bool overcommit;
      
      impl() {}
      ~impl() { close(); }
      
      std::streamsize write(const char_type* s, std::streamsize n);
      void open(int rank, int tag, size_t buffer_size=4096, bool _terminate_on_close=true, bool _overcommit=false);
      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096, bool _terminate_on_close=true, bool _overcommit=false);
      bool is_open() const;
      void close();
      bool test() const;
      void wait() const;
      std::streamsize flush(const bool if_filled=false);
      void terminate();
      void finalize();
      bool test_terminate();
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  
  void mpi_device_sink::impl::wait() const
  {
    utils::atomicop::memory_barrier();
    
    if (! is_open())
      return;

    if (! const_cast<mpi_device_sink::impl&>(*this).request_size.Test())
      const_cast<mpi_device_sink::impl&>(*this).request_size.Wait();
    
    if (! const_cast<mpi_device_sink::impl&>(*this).request_buffer.Test())
      const_cast<mpi_device_sink::impl&>(*this).request_buffer.Wait();    
  }

  void mpi_device_source::impl::wait() const
  {
    utils::atomicop::memory_barrier();
    
    if (! is_open())
      return;

    if (! const_cast<mpi_device_source::impl&>(*this).request_size.Test())
      const_cast<mpi_device_source::impl&>(*this).request_size.Wait();
    
    if (! const_cast<mpi_device_source::impl&>(*this).request_buffer.Test())
      const_cast<mpi_device_source::impl&>(*this).request_buffer.Wait();    
  }
  
  bool mpi_device_sink::impl::test() const
  {
    utils::atomicop::memory_barrier();

    if (! is_open())
      return true;

    return const_cast<mpi_device_sink::impl&>(*this).request_size.Test() && const_cast<mpi_device_sink::impl&>(*this).request_buffer.Test();
  }
  
  bool mpi_device_source::impl::test() const
  {
    utils::atomicop::memory_barrier();

    if (! is_open()) 
      return true;
    
    return const_cast<mpi_device_source::impl&>(*this).request_size.Test() && const_cast<mpi_device_source::impl&>(*this).request_buffer.Test();
  }
  
  bool mpi_device_sink::impl::is_open() const
  {
    return ! buffer.empty();
  }
  
  bool mpi_device_source::impl::is_open() const
  {
    return ! buffer.empty();
  }
  
  
  std::streamsize mpi_device_sink::impl::write(const char_type* s, std::streamsize n)
  {
    if (! is_open()) return -1;
    
    if (overcommit){
      flush(true);
      
      buffer_overcommit.insert(buffer_overcommit.end(), s, s + n);
      buffer_offset = buffer_overcommit.size();
      
      flush(true);
      
      return n;
    } else {
      
      wait();
      
      const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - buffer_offset), n);
      std::copy(s, s + copy_size, buffer.begin() + buffer_offset);
      buffer_offset += copy_size;
      
      flush(true);
      
      return copy_size;
    }
  }
  
  std::streamsize mpi_device_sink::impl::flush(const bool if_filled)
  {
    if (buffer.empty()) return -1;
    if (! test()) return 0;

    if (overcommit) {
      if (if_filled ? buffer_offset >= buffer.size() : buffer_offset > 0) {
	send_size = std::min(buffer.size(), buffer_overcommit.size());
	std::copy(buffer_overcommit.begin(), buffer_overcommit.begin() + send_size, buffer.begin());
	buffer_overcommit.erase(buffer_overcommit.begin(), buffer_overcommit.begin() + send_size);
	buffer_offset = buffer_overcommit.size();

	if (buffer_offset == 0)
	  buffer_type(buffer_overcommit).swap(buffer_overcommit);
	
	request_size.Start();
	request_buffer.Start();
	return send_size;
      } else
	return 0;
    } else {
      if (if_filled ? buffer_offset == buffer.size() : buffer_offset > 0) {
	send_size = buffer_offset;
	buffer_offset = 0;
	
	request_size.Start();
	request_buffer.Start();
	return send_size;
      } else
	return 0;
    }
  }
  
  bool mpi_device_sink::impl::test_terminate()
  {
    utils::atomicop::memory_barrier();

    return send_size == 0;
  }
  
  void mpi_device_sink::impl::terminate()
  {
    if (buffer.empty()) return;
    
    while (buffer_offset > 0) {
      wait();
      flush();
    }
    
    wait();
    
    send_size = 0;
    
    request_size.Start();
    request_buffer.Start();
  };
  
  void mpi_device_sink::impl::finalize()
  {
    utils::atomicop::memory_barrier();

    if (buffer.empty()) return;
    
    if (send_size != 0)
      terminate();
    
    request_size.Wait();
    request_buffer.Wait();
    
    //request_size.Free();
    //request_buffer.Free();
    
    buffer.clear();
    buffer_overcommit.clear();
    buffer_type(buffer_overcommit).swap(buffer_overcommit);
    
    send_size = 0;
    buffer_offset = 0;
    terminate_on_close = true;
    overcommit = false;
  }
  
  void mpi_device_sink::impl::close()
  {
    if (buffer.empty()) return;
    
    flush();
    
    if (terminate_on_close) {
      terminate();
      
      finalize();
    }
  }
  
  void mpi_device_sink::impl::open(int rank, int tag, size_t buffer_size, bool _terminate_on_close, bool _overcommit)
  {
    open(MPI::COMM_WORLD, rank, tag, buffer_size, _terminate_on_close, _overcommit);
  }

  void mpi_device_sink::impl::open(MPI::Comm& comm, int rank, int tag, size_t buffer_size, bool _terminate_on_close, bool _overcommit)
  {
    typedef unsigned int send_size_type;

    close();
    
    // minimum of 4-bytes...
    buffer.reserve(std::max(buffer_size, size_t(4)));
    buffer.resize(std::max(buffer_size, size_t(4)));
    
    buffer_overcommit.clear();
    buffer_type(buffer_overcommit).swap(buffer_overcommit);
    
    send_size = send_size_type(-1);
    buffer_offset = 0;
    
    request_size = comm.Send_init((void*) &send_size, 1, MPI::INT, rank, (tag << 8) | 0);
    request_buffer = comm.Send_init(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, (tag << 8) | 1);
    
    terminate_on_close = _terminate_on_close;
    overcommit = _overcommit;
  }
  
  std::streamsize mpi_device_source::impl::read(char_type* s, std::streamsize n)
  {
    // try read upto n-bytes...
    if (! is_open()) return -1;
    
    wait();
    
    if (recv_size == 0) {
      close();
      return -1;
    }
    
    const std::streamsize copy_size = std::min(std::streamsize(recv_size - buffer_offset), n);
    std::copy(buffer.begin() + buffer_offset, buffer.begin() + buffer_offset + copy_size, s);
    buffer_offset += copy_size;
    
    if (buffer_offset == recv_size) {
      request_size.Start();
      request_buffer.Start();
      buffer_offset = 0;
    }
    
    return copy_size;
  };
  
  void mpi_device_source::impl::close()
  {
    if (buffer.empty()) return;
    
    //request_size.Free();
    //request_buffer.Free();
    
    recv_size = 0;
    buffer_offset = 0;
  }

  void mpi_device_source::impl::open(int rank, int tag, size_t buffer_size)
  {
    open(MPI::COMM_WORLD, rank, tag, buffer_size);
  }
  
  void mpi_device_source::impl::open(MPI::Comm& comm, int rank, int tag, size_t buffer_size)
  {
    close();
    
    buffer.reserve(std::max(buffer_size, size_t(4)));
    buffer.resize(std::max(buffer_size, size_t(4)));
    
    recv_size = 0;
    buffer_offset = 0;
    
    request_size = comm.Recv_init((void*) &recv_size, 1, MPI::INT, rank, (tag << 8) | 0);
    request_buffer = comm.Recv_init(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, (tag << 8) | 1);
    
    request_size.Start();
    request_buffer.Start();
  }

};

#endif
