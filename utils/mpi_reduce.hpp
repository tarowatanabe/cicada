// -*- mode: c++ -*-
//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_REDUCE__HPP__
#define __UTILS__MPI_REDUCE__HPP__ 1

#include <mpi.h>

#include <vector>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/shared_ptr.hpp>

#include <utils/mpi_device.hpp>

namespace utils
{
  
  namespace impl
  {
    template <typename Reducer, typename RankIterator>
    void mpi_reduce_recv(Reducer& reducer, MPI::Comm& comm, RankIterator first, RankIterator last, int tag, size_t buffer_size);
    
    template <typename Reducer>
    void mpi_reduce_recv(Reducer& reducer, MPI::Comm& comm, int rank, int tag, size_t buffer_size);
    
    template <typename Reducer>
    void mpi_reduce_send(const Reducer& reducer, MPI::Comm& comm, int rank, int tag, size_t buffer_size);
  };

  template <typename Reducer>
  inline
  void mpi_reduce(const Reducer& reducer, int rank, int tag, size_t buffer_size=4096)
  {
    mpi_reduce(reducer, MPI::COMM_WORLD, rank, tag, buffer_size);
  }
  
  template <typename Reducer>
  inline
  void mpi_reduce(const Reducer& reducer, MPI::Comm& comm, int rank, int tag, size_t buffer_size=4096)
  {
    typedef std::vector<int, std::allocator<int> > rank_set_type;
    
    const int mpi_rank = comm.Get_rank();
    const int mpi_size = comm.Get_size();
    
    // we will "swap" rank == 0 and rank!
    const int mpi_rank_logical = (mpi_rank == rank ? 0 : (mpi_rank == 0 ? rank : mpi_rank));
    
    rank_set_type ranks;
    int merge_size = mpi_size;
    
    while (merge_size > 1 && mpi_rank_logical < merge_size) {
      const int reduce_size = (merge_size / 2 == 0 ? 1 : merge_size / 2);
      
      if (mpi_rank_logical < reduce_size) {
	ranks.clear();
	
	for (int i = reduce_size; i < merge_size; ++ i)
	  if (i % raduce_size == mpi_rank_logical)
	    ranks.push_back(i == rank ? 0 : (i == 0 ? rank : i));
	
	if (ranks.empty()) continue;
	
	// receiving
	if (ranks.size() == 1)
	  impl::mpi_reduce_recv(const_cast<Reducer&>(reducer), comm, ranks.front(), tag, buffer_size);
	else
	  impl::mpi_reduce_recv(const_cast<Reducer&>(reducer), comm, ranks.begin(), ranks.end(), tag, buffer_size);
      } else {
	const int mpi_rank_reduce_logical = mpi_rank_logical % reduce_size;
	const int mpi_rank_reduce = (mpi_rank_reduce_logical == rank ? 0 : (mpi_rank_reduce_logical == 0 ? rank : mpi_rank_reduce_logical));
	
	// sending...
	impl::mpi_reduce_send(reducer, comm, mpi_rank_reduce, tag, buffer_size);
      }
      
      merge_size = reduce_size;
    }
  }
  
  namespace impl
  {
    template <typename Reducer, typename RankIterator>
    inline
    void mpi_reduce_recv(Reducer& reducer, MPI::Comm& comm, RankIterator first, RankIterator last, int tag, size_t buffer_size)
    {
      typedef utils::mpi_device_source            device_type;
      typedef boost::iostreams::filtering_istream stream_type;
      
      typedef boost::shared_ptr<device_type> device_ptr_type;
      typedef boost::shared_ptr<stream_type> stream_ptr_type;
      
      typedef std::vector<device_ptr_type, std::allocator<device_ptr_type> > device_ptr_set_type;
      typedef std::vector<stream_ptr_type, std::allocator<stream_ptr_type> > stream_ptr_set_type;
      
      device_ptr_set_type device;
      stream_ptr_set_type stream;
      
      for (/**/; first != last; ++ first) {
	device.push_back(device_ptr_type(new device_type(comm, *first, tag, buffer_size)));
	stream.push_back(stream_ptr_type(new stream_type()));
	
	stream.back()->push(boost::iostreams::zlib_decompressor());
	stream.back()->push(*device.back());
      }
      
      int non_found_iter = 0;
      
      while (1) {
	bool found = false;
	
	for (size_t i = 0; i != device.size(); ++ i)
	  while (stream[i] && device[i] && device[i]->test()) {
	    if (! reducer(*stream[i])) {
	      stream[i].reset();
	      device[i].reset();
	    }
	    
	    found = true;
	  }
	
	if (std::count(device.begin(), device.end(), device_ptr_type()) == static_cast<int>(device.size())) break;
	
	if (! found) {
	  boost::thread::yield();
	  ++ non_found_iter;
	} else
	  non_found_iter = 0;
	
	if (non_found_iter >= 64) {
	  struct timespec tm;
	  tm.tv_sec = 0;
	  tm.tv_nsec = 2000001; // above 2ms
	  nanosleep(&tm, NULL);
	  
	  non_found_iter = 0;
	}
      }
    }
    
    
    template <typename Reducer>
    inline
    void mpi_reduce_recv(Reducer& reducer, MPI::Comm& comm, int rank, int tag, size_t buffer_size)
    {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::zlib_decompressor());
      is.push(utils::mpi_device_source(comm, rank, tag, buffer_size));
      
      while (is)
	reducer(is);
    }
    

    template <typename Reducer>
    inline
    void mpi_reduce_send(const Reducer& reducer, MPI::Comm& comm, int rank, int tag, size_t buffer_size)
    {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::zlib_compressor());
      os.push(utils::mpi_device_sink(comm, rank, tag, buffer_size));
      
      reducer(os);
    }
    
  };
  
};


#endif
