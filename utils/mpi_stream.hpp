// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_STREAM__HPP__
#define __UTILS__MPI_STREAM__HPP__ 1

//
// basically, we maintain two streams, one to keep track of data size
// the other, data stream, by streaming data in a stream, meaning that larger data is split into chunk.
//

#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>

#include <boost/thread.hpp>

#include <mpi.h>

#include <utils/atomicop.hpp>

namespace utils
{
  struct __basic_mpi_stream_base
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    static const int tag_shift = 8;
    
    enum {
      tag_ack    = 0,
      tag_size   = 1,
      tag_buffer = 2,
      tag_ready  = 3,
    };
  };

  template <typename Alloc=std::allocator<char> >
  class basic_mpi_ostream
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  public:
    basic_mpi_ostream(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096) : pimpl(new impl()) { open(comm, rank, tag, buffer_size); }
    basic_mpi_ostream(int rank, int tag, size_t buffer_size=4096) : pimpl(new impl()) { open(rank, tag, buffer_size); }

  public:
    void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096) { pimpl->open(comm, rank, tag, buffer_size); }
    void open(int rank, int tag, size_t buffer_size=4096) { pimpl->open(MPI::COMM_WORLD, rank, tag, buffer_size); }
    
    basic_mpi_ostream& write(const std::string& data) { pimpl->write(data); return *this; }
    void close() { pimpl->close(); }
    
    bool test() { return pimpl->test(); }
    void wait() { pimpl->wait(); };
   
    void terminate() { pimpl->terminate(); }
    bool terminated() { return pimpl->terminated(); }

    operator bool() const { return pimpl->is_open(); }

  private:
    struct impl : public __basic_mpi_stream_base
    {
      typedef std::vector<char, Alloc> buffer_type;

      impl() : comm(0) {}
      ~impl() { close(); }

      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096);
      void close();

      void write(const std::string& data);

      bool test();
      void wait();
      
      void terminate();
      bool terminated();

      bool is_open() const;

      MPI::Comm* comm;
      int        rank;
      int        tag;
      
      buffer_type  buffer;
      volatile int buffer_size;
      volatile int ack;
      volatile int ack_expected;
      volatile int state;
      
      MPI::Prequest request_ack;
      MPI::Prequest request_size;
      MPI::Request  request_buffer;
    };

    boost::shared_ptr<impl> pimpl;
  };

  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::write(const std::string& data)
  {
    wait();
    
    buffer_size = data.size();
    buffer.clear();
    buffer.insert(buffer.end(), data.begin(), data.end());
    
    request_size.Start();
    state = tag_size;
    
    if (! buffer.empty())
      request_buffer = comm->Isend(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, (tag << tag_shift) | tag_buffer);
    
    ++ ack_expected;
    // request_ack.Start();
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::terminate()
  {
    wait();
    
    if (! terminated()) {
      buffer_size = -1;
      request_size.Start();
      state = tag_size;
    }
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::open(MPI::Comm& __comm, int __rank, int __tag, size_t __buffer_size)
  {
    close();
    
    comm = &__comm;
    rank = __rank;
    tag  = __tag;
    
    buffer.clear();
    buffer_size = 0;
    ack = 0;
    ack_expected = 0;
    
    request_ack  = comm->Recv_init(const_cast<int*>(&ack), 1, MPI::INT, rank, (tag << tag_shift) | tag_ack);
    request_size = comm->Send_init(const_cast<int*>(&buffer_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_size);
    
    request_ack.Start();
    state = tag_ack;
  }

  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::close()
  {
    if (comm) {
      while (! terminated()) {
	terminate();
	boost::thread::yield();
      }
      wait();
      
      //request_ack.Free();
      //request_size.Free();
      //request_buffer.Free();
    }
    
    comm = 0;
    rank = -1;
    tag = -1;
    state = tag_ready;
    
    buffer.clear();
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream<Alloc>::impl::terminated()
  {
    return test() && buffer_size < 0;
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream<Alloc>::impl::is_open() const
  {
    return comm;
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream<Alloc>::impl::test()
  {
    utils::atomicop::memory_barrier();
    
    switch (state) {
    case tag_size:
      if (! request_size.Test()) return false;
      
      if (buffer_size < 0) {	
	state = tag_ready;
	return true;
      } else if (buffer_size == 0) {
	request_ack.Start();
	state = tag_ack;
	if (! request_ack.Test()) return false;
	state = tag_ready;
	return true;
      } else
	state = tag_buffer;
    case tag_buffer:
      if (! request_buffer.Test()) return false;
      request_ack.Start();
      state = tag_ack;
    case tag_ack:
      if (! request_ack.Test()) return false;
      state = tag_ready;
    default:
      return true;
    }
  }
  
  template <typename Alloc>
  void basic_mpi_ostream<Alloc>::impl::wait()
  {
    const size_type mask = (1 << 6) - 1;
    for (size_type iter = 0; ! test(); ++ iter) {
      if ((iter & mask) == mask) {
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      } else
	boost::thread::yield();
    }
  }



  template <typename Alloc=std::allocator<char> >
  class basic_mpi_istream
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  public:
    basic_mpi_istream(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096, bool no_ready=false)
      : pimpl(new impl()) { open(comm, rank, tag, buffer_size, no_ready); }
    basic_mpi_istream(int rank, int tag, size_t buffer_size=4096, bool no_ready=false)
      : pimpl(new impl()) { open(rank, tag, buffer_size, no_ready); }
    
  public:
    void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096, bool no_ready=false) { pimpl->open(comm, rank, tag, buffer_size, no_ready); }
    void open(int rank, int tag, size_t buffer_size=4096, bool no_ready=false) { pimpl->open(MPI::COMM_WORLD, rank, tag, buffer_size, no_ready); }

    basic_mpi_istream& read(std::string& data) { pimpl->read(data); return *this; }
    void close() { pimpl->close(); }

    bool test() { return pimpl->test(); }
    void wait() { pimpl->wait(); };
    void ready() { pimpl->ready(); }

    operator bool() const { return pimpl->is_open(); }

  private: 
    struct impl : public __basic_mpi_stream_base
    {
      typedef std::vector<char, Alloc> buffer_type;

      impl() : comm(0) {}
      ~impl() { close(); }

      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096, bool __no_ready=false);
      void close();

      void read(std::string& data);

      bool test();
      void wait();
      void ready();

      bool is_open() const;

      MPI::Comm* comm;
      int        rank;
      int        tag;

      buffer_type  buffer;
      volatile int buffer_size;
      volatile int ack;
      volatile int state;
      
      MPI::Prequest request_ack;
      MPI::Prequest request_size;
      MPI::Request  request_buffer;
      
      bool no_ready;
    };

    boost::shared_ptr<impl> pimpl;
  };


  template <typename Alloc>
  void basic_mpi_istream<Alloc>::impl::open(MPI::Comm& __comm, int __rank, int __tag, size_t __buffer_size, bool __no_ready)
  {
    close();
    
    comm = &__comm;
    rank = __rank;
    tag  = __tag;
    
    no_ready = __no_ready;
    
    buffer.clear();
    buffer_size = -1;
    ack = 0;
    
    request_ack  = comm->Send_init(const_cast<int*>(&ack), 1, MPI::INT, rank, (tag << tag_shift) | tag_ack);
    request_size = comm->Recv_init(const_cast<int*>(&buffer_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_size);

    request_ack.Start();
    state = tag_ack;
    
    test();
  }

  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::close()
  {
    if (comm) {
      wait();
      
      //request_ack.Free();
      //request_size.Free();
      //request_buffer.Free();
    }
    
    comm = 0;
    rank = -1;
    tag = -1;
    state = tag_ready;
    
    buffer.clear();
    buffer_size = -1;
    
    no_ready = false;
  }

  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::read(std::string& data)
  {
    data.clear();
    
    wait();
    
    if (buffer_size < 0)
      close();
    else {
      data.insert(data.end(), buffer.begin(), buffer.end());
      buffer.clear();
      
      if (! no_ready) {
	++ ack;
	request_ack.Start();
	state = tag_ack;
	
	test();
      }
    }
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::ready()
  {
    if (! no_ready) return;
    
    ++ ack;
    request_ack.Start();
    state = tag_ack;
    
    test();
  }


  template <typename Alloc>
  inline
  bool basic_mpi_istream<Alloc>::impl::is_open() const
  {
    return comm;
  }
  
  template <typename Alloc>
  inline
  bool basic_mpi_istream<Alloc>::impl::test()
  {
    utils::atomicop::memory_barrier();
    
    switch (state) {
    case tag_ack:
      if (! request_ack.Test()) return false;
      request_size.Start();
      state = tag_size;
    case tag_size:
      if (! request_size.Test()) return false;
      if (buffer_size <= 0) {
	state = tag_ready;
	return true;
      } else {
	buffer.resize(buffer_size);
	request_buffer = comm->Irecv(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, (tag << tag_shift) | tag_buffer);
	state = tag_buffer;
      }
    case tag_buffer:
      if (! request_buffer.Test()) return false;
      state = tag_ready;
    default:
      return true;
    }
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::wait()
  {
    const size_type mask = (1 << 6) - 1;
    for (size_type iter = 0; ! test(); ++ iter) {
      if ((iter & mask) == mask) {
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      } else
	boost::thread::yield();
    }
  }
  
  // basic_mpi_{i,o}stream
  typedef basic_mpi_ostream<> mpi_ostream;
  typedef basic_mpi_istream<> mpi_istream;
};


#endif
