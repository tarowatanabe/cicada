// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__MPI_TRAITS__HPP__
#define __UTILS__MPI_TRAITS__HPP__ 1

#include <stdint.h>

#include <mpi.h>

namespace utils
{
  template <typename Tp>
  struct mpi_traits {};

  template <>
  struct mpi_traits<long double>
  {
    typedef long double value_type;
    static inline MPI_Datatype data_type() { return MPI_LONG_DOUBLE; }
  };
  
  template <>
  struct mpi_traits<double>
  {
    typedef double value_type;
    static inline MPI_Datatype data_type() { return MPI_DOUBLE; }
  };
  
  template <>
  struct mpi_traits<float>
  {
    typedef float value_type;
    static inline MPI_Datatype data_type() { return MPI_FLOAT; }
  };

  template <>
  struct mpi_traits<bool>
  {
    typedef bool value_type;
    static inline MPI_Datatype data_type() { return MPI_BOOL; }
  };

  template <>
  struct mpi_traits<char>
  {
    typedef char value_type;
    static inline MPI_Datatype data_type() { return MPI_CHAR; }
  };
  
  template <>
  struct mpi_traits<unsigned char>
  {
    typedef unsigned char value_type;
    static inline MPI_Datatype data_type() { return MPI_UNSIGNED_CHAR; }
  };
  
  template <>
  struct mpi_traits<short>
  {
    typedef short value_type;
    static inline MPI_Datatype data_type() { return MPI_SHORT; }
  };
  
  template <>
  struct mpi_traits<unsigned short>
  {
    typedef unsigned short value_type;
    static inline MPI_Datatype data_type() { return MPI_UNSIGNED_SHORT; }
  };
  
  template <>
  struct mpi_traits<int>
  {
    typedef int value_type;
    static inline MPI_Datatype data_type() { return MPI_INT; }
  };
  
  template <>
  struct mpi_traits<unsigned int>
  {
    typedef unsigned int value_type;
    static inline MPI_Datatype data_type() { return MPI_UNSIGNED; }
  };

  template <>
  struct mpi_traits<long>
  {
    typedef long value_type;
    static inline MPI_Datatype data_type() { return MPI_LONG; }
  };
  
  template <>
  struct mpi_traits<unsigned long>
  {
    typedef unsigned long value_type;
    static inline MPI_Datatype data_type() { return MPI_UNSIGNED_LONG; }
  };
  
};

#endif
