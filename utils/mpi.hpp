// -*- mode: c++ -*-

#ifndef __UTILS__MPI__HPP__
#define __UTILS__MPI__HPP__ 1

#include <mpi.h>

namespace utils
{
  struct mpi_world
  {
    mpi_world(int argc, char** argv) {  __provided = MPI::Init_thread(argc, argv, MPI::THREAD_SERIALIZED); }
    ~mpi_world() {  MPI::Finalize(); }
    
    const int& provided() const { return __provided; }

    int __provided;
  };
  
  struct mpi_intercomm
  {
    mpi_intercomm(const MPI::Intercomm& _comm) : comm(_comm) {}
    ~mpi_intercomm() { comm.Free(); }
    
    operator MPI::Intercomm&() { return comm; }
    
    MPI::Intercomm comm;
  };
};

#endif
