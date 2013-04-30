// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
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
      typedef size_t size_type;
      typedef unsigned int stream_size_type;
      typedef std::vector<char_type, utils::mpi_allocator<char_type> > buffer_type;
      
      MPI::Prequest             request;
      volatile stream_size_type recv_size;
      volatile size_type        buffer_offset;
      buffer_type buffer;

      bool ready;
      
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
    std::streamsize committed() const { return pimpl->committed(); }
    void terminate() { pimpl->terminate(); }
    void finalize() { pimpl->finalize(); }
    bool test_terminate() { return pimpl->test_terminate(); }
    
  private:
    struct impl
    {
      typedef size_t size_type;
      typedef unsigned int stream_size_type;
      typedef std::vector<char_type, utils::mpi_allocator<char_type> > buffer_type;
      
      MPI::Prequest             request;
      volatile stream_size_type send_size;
      volatile size_type        buffer_offset;
      buffer_type buffer;
      buffer_type buffer_overcommit;
      bool terminate_on_close;
      bool overcommit;

      bool ready;
      
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
      std::streamsize committed() const;
      void terminate();
      void finalize();
      bool test_terminate();
    };
    
    boost::shared_ptr<impl> pimpl;
  };
  
  
  void mpi_device_sink::impl::wait() const
  {
    const size_t mask = (1 << 6) - 1;
    for (size_t iter = 0; ! test(); ++ iter) {
      if ((iter & mask) == mask) {
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      } else
	boost::thread::yield();
    }    
  }

  void mpi_device_source::impl::wait() const
  {
    const size_t mask = (1 << 6) - 1;
    for (size_t iter = 0; ! test(); ++ iter) {
      if ((iter & mask) == mask) {
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      } else
	boost::thread::yield();
    }

  }
  
  bool mpi_device_sink::impl::test() const
  {
    if (! is_open() || ready)
      return true;

    if (! const_cast<mpi_device_sink::impl&>(*this).request.Test())
      return false;

    const_cast<bool&>(ready) = true;
    
    return true;
  }
  
  bool mpi_device_source::impl::test() const
  {
    if (! is_open() || ready) 
      return true;
    
    if (! const_cast<mpi_device_source::impl&>(*this).request.Test())
      return false;
    
    std::copy(buffer.end() - sizeof(stream_size_type), buffer.end(), (char_type*) &recv_size);
    const_cast<bool&>(ready) = true;
    
    return true;
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
      
      const std::streamsize copy_size = std::min(std::streamsize(buffer.size() - sizeof(stream_size_type) - buffer_offset), n);
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
      if (if_filled ? buffer_offset >= buffer.size() - sizeof(stream_size_type) : buffer_offset > 0) {
	
	send_size = std::min(buffer.size() - sizeof(stream_size_type), buffer_overcommit.size());
	
	std::copy(buffer_overcommit.begin(), buffer_overcommit.begin() + send_size, buffer.begin());
	
	// if the capacity is twice as large as the buffer size, shrink...
	
	buffer_overcommit.erase(buffer_overcommit.begin(), buffer_overcommit.begin() + send_size);
	if (buffer_overcommit.capacity() > (buffer_overcommit.size() << 1))
	  buffer_type(buffer_overcommit).swap(buffer_overcommit);
	
	buffer_offset = buffer_overcommit.size();
	
	std::copy((char_type*) &send_size, ((char_type*) &send_size) + sizeof(stream_size_type), buffer.end() - sizeof(stream_size_type));
	
	request.Start();
	ready = false;
	
	return send_size;
      } else
	return 0;
    } else {
      if (if_filled ? buffer_offset == buffer.size() - sizeof(stream_size_type) : buffer_offset > 0) {
	send_size = buffer_offset;
	buffer_offset = 0;
	
	std::copy((char_type*) &send_size, ((char_type*) &send_size) + sizeof(stream_size_type), buffer.end() - sizeof(stream_size_type));
	
	request.Start();
	ready = false;
	
	return send_size;
      } else
	return 0;
    }
  }

  std::streamsize mpi_device_sink::impl::committed() const
  {
    return buffer_overcommit.size();
  }
  
  bool mpi_device_sink::impl::test_terminate()
  {
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
    
    std::copy((char_type*) &send_size, ((char_type*) &send_size) + sizeof(stream_size_type), buffer.end() - sizeof(stream_size_type));
    
    request.Start();
    ready = false;
  };
  
  void mpi_device_sink::impl::finalize()
  {
    if (buffer.empty()) return;
    
    if (send_size != 0)
      terminate();
    
    wait();
    
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
    close();
    
    // minimum of 4-bytes...
    buffer.reserve(std::max(buffer_size, size_t(128) + sizeof(stream_size_type)));
    buffer.resize(std::max(buffer_size, size_t(128) + sizeof(stream_size_type)));
    
    buffer_overcommit.clear();
    buffer_type(buffer_overcommit).swap(buffer_overcommit);
    
    send_size = stream_size_type(-1);
    buffer_offset = 0;
    
    request = comm.Send_init(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, tag);
    ready = true;
    
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
      request.Start();
      ready = false;
      
      buffer_offset = 0;
    }
    
    return copy_size;
  };
  
  void mpi_device_source::impl::close()
  {
    if (buffer.empty()) return;
    
    buffer.clear();
    
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
    
    buffer.reserve(std::max(buffer_size, size_t(128) + sizeof(stream_size_type)));
    buffer.resize(std::max(buffer_size, size_t(128) + sizeof(stream_size_type)));
    
    recv_size = 0;
    buffer_offset = 0;
    
    request = comm.Recv_init(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, tag);
    
    request.Start();
    ready = false;
  }

};

#endif
