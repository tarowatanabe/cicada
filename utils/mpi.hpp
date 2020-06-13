// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI__HPP__
#define __UTILS__MPI__HPP__ 1

#include <memory>

#include <mpi.h>

namespace utils
{
  struct mpi_world
  {
    mpi_world(int argc, char** argv) {  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &__provided); }
    ~mpi_world() {  MPI_Finalize(); }
    
    const int& provided() const { return __provided; }

    int __provided;
  };
  
  struct mpi_intercomm
  {
    mpi_intercomm(const MPI_Comm& _comm) : comm(_comm) {}
    ~mpi_intercomm() { MPI_Comm_free(&comm); }
    
    operator MPI_Comm&() { return comm; }
    
    MPI_Comm comm;
  };

  struct mpi_comm {
    mpi_comm(const MPI_Comm& _comm) : comm(_comm) {}
    mpi_comm() : comm(MPI_COMM_WORLD) {}

    operator MPI_Comm&() { return comm; }

    int rank() const {
      int rank = 0;
      MPI_Comm_rank(comm, &rank);
      return rank;
    }

    int size() const {
      int size = 0;
      MPI_Comm_size(comm, &size);
      return size;
    }

    MPI_Comm comm;
  };

  struct mpi_request
  {
    void start() {
      MPI_Start(&request);
    }

    void free() {
      MPI_Request_free(&request);
    }

    bool test() const {
      int flag = false;
      MPI_Status status;
      MPI_Test(&const_cast<MPI_Request&>(request), &flag, &status);
      return flag;
    }

    void wait() const {
      MPI_Status status;
      MPI_Wait(&const_cast<MPI_Request&>(request), &status);
    }

    MPI_Request request;
  };

  inline void mpi_abort(const int errorcode) {
    MPI_Abort(MPI_COMM_WORLD, errorcode);
  }

};

#endif
