// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_STREAM_SIMPLE__HPP__
#define __UTILS__MPI_STREAM_SIMPLE__HPP__ 1

//
// no ack-variant of mpi-stream!
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
  struct __basic_mpi_stream_simple_base
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
  class basic_mpi_ostream_simple
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  public:
    basic_mpi_ostream_simple(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096) : pimpl(new impl()) { open(comm, rank, tag, buffer_size); }
    basic_mpi_ostream_simple(int rank, int tag, size_t buffer_size=4096) : pimpl(new impl()) { open(rank, tag, buffer_size); }

  public:
    void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096) { pimpl->open(comm, rank, tag, buffer_size); }
    void open(int rank, int tag, size_t buffer_size=4096) { pimpl->open(MPI::COMM_WORLD, rank, tag, buffer_size); }
    
    basic_mpi_ostream_simple& write(const std::string& data) { pimpl->write(data); return *this; }
    void close() { pimpl->close(); }
    
    bool test() { return pimpl->test(); }
    void wait() { pimpl->wait(); };
   
    void terminate() { pimpl->terminate(); }
    bool terminated() { return pimpl->terminated(); }

    operator bool() const { return pimpl->is_open(); }

  private:
    struct impl : public __basic_mpi_stream_simple_base
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
      volatile int state;
      
      MPI::Prequest request_size;
      MPI::Request  request_buffer;
    };

    boost::shared_ptr<impl> pimpl;
  };

  template <typename Alloc>
  inline
  void basic_mpi_ostream_simple<Alloc>::impl::write(const std::string& data)
  {
    wait();
    
    buffer_size = data.size();
    buffer.clear();
    buffer.insert(buffer.end(), data.begin(), data.end());
    
    request_size.Start();
    state = tag_size;
    
    if (! buffer.empty())
      request_buffer = comm->Isend(&(*buffer.begin()), buffer.size(), MPI::CHAR, rank, (tag << tag_shift) | tag_buffer);
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_ostream_simple<Alloc>::impl::terminate()
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
  void basic_mpi_ostream_simple<Alloc>::impl::open(MPI::Comm& __comm, int __rank, int __tag, size_t __buffer_size)
  {
    close();
    
    comm = &__comm;
    rank = __rank;
    tag  = __tag;
    
    buffer.clear();
    buffer_size = 0;
    
    request_size = comm->Send_init(const_cast<int*>(&buffer_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_size);
    state = tag_ready;
  }

  template <typename Alloc>
  inline
  void basic_mpi_ostream_simple<Alloc>::impl::close()
  {
    if (comm) {
      while (! terminated()) {
	terminate();
	boost::thread::yield();
      }
      wait();
      
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
  bool basic_mpi_ostream_simple<Alloc>::impl::terminated()
  {
    return test() && buffer_size < 0;
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream_simple<Alloc>::impl::is_open() const
  {
    return comm;
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream_simple<Alloc>::impl::test()
  {
    utils::atomicop::memory_barrier();
    
    switch (state) {
    case tag_size:
      if (! request_size.Test()) return false;
      if (buffer_size <= 0) {
	state = tag_ready;
	return true;
      } else
	state = tag_buffer;
    case tag_buffer:
      if (! request_buffer.Test()) return false;
      state = tag_ready;
    default:
      return true;
    }
  }
  
  template <typename Alloc>
  void basic_mpi_ostream_simple<Alloc>::impl::wait()
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
  class basic_mpi_istream_simple
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
  public:
    basic_mpi_istream_simple(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096)
      : pimpl(new impl()) { open(comm, rank, tag, buffer_size); }
    basic_mpi_istream_simple(int rank, int tag, size_t buffer_size=4096)
      : pimpl(new impl()) { open(rank, tag, buffer_size); }
    
  public:
    void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096) { pimpl->open(comm, rank, tag, buffer_size); }
    void open(int rank, int tag, size_t buffer_size=4096) { pimpl->open(MPI::COMM_WORLD, rank, tag, buffer_size); }

    basic_mpi_istream_simple& read(std::string& data) { pimpl->read(data); return *this; }
    void close() { pimpl->close(); }

    bool test() { return pimpl->test(); }
    void wait() { pimpl->wait(); };

    operator bool() const { return pimpl->is_open(); }

  private: 
    struct impl : public __basic_mpi_stream_simple_base
    {
      typedef std::vector<char, Alloc> buffer_type;

      impl() : comm(0) {}
      ~impl() { close(); }

      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096);
      void close();

      void read(std::string& data);

      bool test();
      void wait();

      bool is_open() const;

      MPI::Comm* comm;
      int        rank;
      int        tag;

      buffer_type  buffer;
      volatile int buffer_size;
      volatile int state;
      
      MPI::Prequest request_size;
      MPI::Request  request_buffer;
    };

    boost::shared_ptr<impl> pimpl;
  };


  template <typename Alloc>
  void basic_mpi_istream_simple<Alloc>::impl::open(MPI::Comm& __comm, int __rank, int __tag, size_t __buffer_size)
  {
    close();
    
    comm = &__comm;
    rank = __rank;
    tag  = __tag;
    
    buffer.clear();
    buffer_size = -1;
    
    request_size = comm->Recv_init(const_cast<int*>(&buffer_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_size);
    request_size.Start();
    state = tag_size;
  }

  template <typename Alloc>
  inline
  void basic_mpi_istream_simple<Alloc>::impl::close()
  {
    if (comm) {
      wait();
      
      //request_size.Free();
      //request_buffer.Free();
    }
    
    comm = 0;
    rank = -1;
    tag = -1;
    state = tag_ready;
    
    buffer.clear();
    buffer_size = -1;
  }

  template <typename Alloc>
  inline
  void basic_mpi_istream_simple<Alloc>::impl::read(std::string& data)
  {
    data.clear();
    
    wait();
    
    if (buffer_size < 0)
      close();
    else {
      data.insert(data.end(), buffer.begin(), buffer.end());
      buffer.clear();
      
      request_size.Start();
      state = tag_size;
    }
  }


  template <typename Alloc>
  inline
  bool basic_mpi_istream_simple<Alloc>::impl::is_open() const
  {
    return comm;
  }
  
  template <typename Alloc>
  inline
  bool basic_mpi_istream_simple<Alloc>::impl::test()
  {
    utils::atomicop::memory_barrier();
    
    switch (state) {
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
  void basic_mpi_istream_simple<Alloc>::impl::wait()
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
  
  // basic_mpi_{i,o}stream_simple
  typedef basic_mpi_ostream_simple<> mpi_ostream_simple;
  typedef basic_mpi_istream_simple<> mpi_istream_simple;
};


#endif
