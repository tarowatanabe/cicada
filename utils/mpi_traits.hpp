// -*- mode: c++ -*-

#ifndef __UTILS__MPI_TRAITS__HPP__
#define __UTILS__MPI_TRAITS__HPP__ 1

#include <stdint.h>

#include <mpi.h>

namespace utils
{
  template <typename Tp>
  struct mpi_traits {};
  
  template <>
  struct mpi_traits<double>
  {
    typedef double value_type;
    static inline MPI::Datatype data_type() { return MPI::DOUBLE; }
  };
  
  template <>
  struct mpi_traits<float>
  {
    typedef float value_type;
    static inline MPI::Datatype data_type() { return MPI::FLOAT; }
  };

  template <>
  struct mpi_traits<int8_t>
  {
    typedef int8_t value_type;
    static inline MPI::Datatype data_type() { return MPI::CHAR; }
  };
  
  template <>
  struct mpi_traits<uint8_t>
  {
    typedef uint8_t value_type;
    static inline MPI::Datatype data_type() { return MPI::UNSIGNED_CHAR; }
  };
  
  template <>
  struct mpi_traits<int16_t>
  {
    typedef int16_t value_type;
    static inline MPI::Datatype data_type() { return MPI::SHORT; }
  };
  
  template <>
  struct mpi_traits<uint16_t>
  {
    typedef uint16_t value_type;
    static inline MPI::Datatype data_type() { return MPI::UNSIGNED_SHORT; }
  };
  
  template <>
  struct mpi_traits<int32_t>
  {
    typedef int32_t value_type;
    static inline MPI::Datatype data_type() { return MPI::INT; }
  };
  
  template <>
  struct mpi_traits<uint32_t>
  {
    typedef uint32_t value_type;
    static inline MPI::Datatype data_type() { return MPI::UNSIGNED; }
  };
  
};

#endif
