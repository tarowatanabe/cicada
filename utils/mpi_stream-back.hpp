// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_STREAM__HPP__
#define __UTILS__MPI_STREAM__HPP__ 1

//
// basically, we maintain two streams, one to keep track of data size
// the other, data stream, by streaming data in a stream, meaning that larger data is split into chunk.
//

#include <string>
#include <vector>
#include <algorithm>

#include <boost/thread.hpp>

#include <mpi.h>

namespace utils
{
  struct __basic_mpi_stream_base
  {
    static const int tag_shift = 8;
    
    static const int tag_ack         = 0;
    static const int tag_size        = 1;
    static const int tag_buffer      = 2;
    static const int tag_buffer_size = 3;
  };

  template <typename Alloc=std::allocator<char> >
  class basic_mpi_ostream
  {
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
    size_t flush() { return pimpl->flush(); }
   
    void terminate() { pimpl->terminate(); }
    bool terminated() { return pimpl->terminated(); }

    operator bool() const { return pimpl->is_open(); }

  private:
    struct impl : public __basic_mpi_stream_base
    {
      typedef std::vector<char, Alloc> buffer_type;

      impl() {}
      ~impl() { close(); }

      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096);
      void close();

      void write(const std::string& data);

      bool test();
      void wait();
      size_t flush();
      
      void terminate();
      bool terminated();

      bool is_open() const;
      
      buffer_type  buffer;
      volatile int buffer_size;

      buffer_type  buffer_send;
      volatile int buffer_send_size;

      MPI::Prequest request_ack;
      MPI::Prequest request_size;

      MPI::Prequest request_buffer;
      MPI::Prequest request_buffer_size;
    };

    boost::shared_ptr<impl> pimpl;
  };

  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::write(const std::string& data)
  {
    wait();
    
    buffer_size = data.size();
    buffer.insert(buffer.end(), data.begin(), data.end());
    
    request_size.Start();
    request_ack.Start();
    
    flush();
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::terminate()
  {
    wait();
    
    if (buffer_size >= 0 && flush() == 0 && request_buffer_size.Test() && request_buffer.Test()) {
      buffer_size = -1;
      request_size.Start();
    }
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::open(MPI::Comm& comm, int rank, int tag, size_t __buffer_size)
  {
    close();

    buffer.clear();
    buffer_size = 0;

    buffer_send.reserve(__buffer_size);
    buffer_send.resize(__buffer_size);
    buffer_send_size = 0;

    request_ack  = comm.Recv_init(0, 0, MPI::INT, rank, (tag << tag_shift) | tag_ack);
    request_size = comm.Send_init(const_cast<int*>(&buffer_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_size);
    
    request_buffer      = comm.Send_init(&(*buffer_send.begin()), buffer_send.size(), MPI::CHAR, rank, (tag << tag_shift) | tag_buffer);
    request_buffer_size = comm.Send_init(const_cast<int*>(&buffer_send_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_buffer_size);

    request_ack.Start();
  }

  template <typename Alloc>
  inline
  void basic_mpi_ostream<Alloc>::impl::close()
  {
    request_size.Cancel();
    request_ack.Cancel();
    request_buffer.Cancel();
    request_buffer_size.Cancel();
    
    buffer.clear();
    buffer_send.clear();
  }

  template <typename Alloc>
  inline
  size_t basic_mpi_ostream<Alloc>::impl::flush()
  {
    // returns remaining buffer sizes...

    const bool test_size   = request_buffer_size.Test();
    const bool test_buffer = request_buffer.Test();
    
    if (! test_size || ! test_buffer)
      return buffer.size() + buffer_send_size;
    
    if (buffer.empty()) return 0;
    
    buffer_send_size = std::min(buffer_send.size(), buffer.size());
    std::copy(buffer.begin(), buffer.begin() + buffer_send_size, buffer_send.begin());
    buffer.erase(buffer.begin(), buffer.begin() + buffer_send_size);
    
    request_buffer_size.Start();
    request_buffer.Start();
    
    return buffer.size() + buffer_send_size;
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream<Alloc>::impl::terminated()
  {
    return test() && flush() == 0 && buffer_size < 0 && request_buffer_size.Test() && request_buffer.Test();
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream<Alloc>::impl::is_open() const
  {
    return ! terminated();
  }
  
  template <typename Alloc>
  bool basic_mpi_ostream<Alloc>::impl::test()
  {
    const bool test_size = request_size.Test();
    const bool test_ack  = request_ack.Test();

    flush();
    
    return test_size && test_ack;
  }
  
  template <typename Alloc>
  void basic_mpi_ostream<Alloc>::impl::wait()
  {
    for (;;) {
      const bool test_size = request_size.Test();
      const bool test_ack  = request_ack.Test();
      
      flush();

      if (test_size && test_ack)
	return;
      else
	boost::thread::yield();
    }
  }



  template <typename Alloc=std::allocator<char> >
  class basic_mpi_istream
  {
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
    size_t fill() { return pimpl->fill(); }

    operator bool() const { return pimpl->is_open(); }

  private: 
    struct impl : public __basic_mpi_stream_base
    {
      typedef std::vector<char, Alloc> buffer_type;

      impl() {}
      ~impl() { close(); }

      void open(MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096, bool __no_ready=false);
      void close();

      void read(std::string& data);

      bool test();
      void wait();
      size_t fill();
      void ready();

      bool is_open() const;

      buffer_type  buffer;
      volatile int buffer_size;

      buffer_type  buffer_recv;
      volatile int buffer_recv_size;

      MPI::Prequest request_ack;
      MPI::Prequest request_size;

      MPI::Prequest request_buffer;
      MPI::Prequest request_buffer_size;

      bool no_ready;
    };

    boost::shared_ptr<impl> pimpl;
  };


  template <typename Alloc>
  void basic_mpi_istream<Alloc>::impl::open(MPI::Comm& comm, int rank, int tag, size_t __buffer_size, bool __no_ready)
  {
    close();

    no_ready = __no_ready;

    buffer.clear();
    buffer_size = -1;

    buffer_recv.reserve(__buffer_size);
    buffer_recv.resize(__buffer_size);
    buffer_recv_size = 0;

    request_ack  = comm.Send_init(0, 0, MPI::INT, rank, (tag << tag_shift) | tag_ack);
    request_size = comm.Recv_init(const_cast<int*>(&buffer_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_size);

    request_buffer      = comm.Recv_init(&(*buffer_recv.begin()), buffer_recv.size(), MPI::CHAR, rank, (tag << tag_shift) | tag_buffer);
    request_buffer_size = comm.Recv_init(const_cast<int*>(&buffer_recv_size), 1, MPI::INT, rank, (tag << tag_shift) | tag_buffer_size);
    
    request_ack.Start();
    request_size.Start();
    request_buffer.Start();
    request_buffer_size.Start();
  }

  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::close()
  {
    request_size.Cancel();
    request_ack.Cancel();
    request_buffer.Cancel();
    request_buffer_size.Cancel();

    buffer.clear();
    buffer_recv.clear();
    
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
      data.insert(data.end(), buffer.begin(), buffer.begin() + buffer_size);
      buffer.erase(buffer.begin(), buffer.begin() + buffer_size);

      if (! no_ready) {
	request_size.Start();
	request_ack.Start();
	
	fill();
      }
    }
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::ready()
  {
    if (! no_ready) return;

    for (;;) {
      const bool test_size = request_size.Test();
      const bool test_ack  = request_ack.Test();

      fill();
      
      if (test_size && test_ack)
	break;
      else
	boost::thread::yield();
    }
    
    request_size.Start();
    request_ack.Start();
    
    fill();
  }

  template <typename Alloc>
  inline
  size_t basic_mpi_istream<Alloc>::impl::fill()
  {
    const bool test_size   = request_buffer_size.Test();
    const bool test_buffer = request_buffer.Test();

    if (! test_size || ! test_buffer)
      return buffer.size();
    
    const int copy_size = buffer_recv_size;
    buffer.insert(buffer.end(), buffer_recv.begin(), buffer_recv.begin() + copy_size);
    
    request_buffer_size.Start();
    request_buffer.Start();
    
    return buffer.size() + copy_size;
  }

  template <typename Alloc>
  inline
  bool basic_mpi_istream<Alloc>::impl::is_open() const
  {
    return ! buffer_recv.empty();
  }

  template <typename Alloc>
  inline
  bool basic_mpi_istream<Alloc>::impl::test()
  {
    const bool test_size = request_size.Test();
    const bool test_ack  = request_ack.Test();
    
    fill();
    
    return (test_size && test_ack && (buffer_size <= 0 || buffer.size() >= buffer_size));
  }
  
  template <typename Alloc>
  inline
  void basic_mpi_istream<Alloc>::impl::wait()
  {
    for (;;) {
      const bool test_size = request_size.Test();
      const bool test_ack  = request_ack.Test();
      
      fill();

      if (test_size && test_ack && (buffer_size <= 0 || buffer.size() >= buffer_size))
	return;
      else
	boost::thread::yield();
    }
  }


  // basic_mpi_{i,o}stream
  typedef basic_mpi_ostream<> mpi_ostream;
  typedef basic_mpi_istream<> mpi_istream;
};


#endif
