// -*- mode: c++ -*-
//
//  Copyright(C) 2009-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __UTILS__ATOMICOP__HPP__
#define __UTILS__ATOMICOP__HPP__ 1

#include <cstddef>

#include <utils/config.hpp>

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __APPLE__
#include <libkern/OSAtomic.h>
#endif

namespace utils
{
  namespace atomicop
  {
    inline
    void memory_barrier()
    {
#if defined(_WIN32)
      ::MemoryBarrier();
#elif defined(__APPLE__)
      OSMemoryBarrier();
#elif defined(__GNUC__) && defined(__x86_64__)
      __asm__ __volatile__("mfence" ::: "memory");
#elif defined(__GNUC__)
      __sync_synchronize();
#else
#   warning "no memory barrier implemented for this platform"
#endif
    }
    
    template <size_t _Size>
    struct __struct_add_and_fetch{};
    

    template <>
    struct __struct_add_and_fetch<4> 
    {
      template <typename Tp>
      static inline
      Tp result(volatile Tp* ptr, Tp addend)
      {
	return (Tp) __result((volatile int32_t*) ptr, (int32_t) addend);
      }

      static inline
      int32_t __result(volatile int32_t* ptr, int32_t addend)
      {
#if defined(_WIN32)
	return ::InterlockedExchangeAdd((void*)ptr, addend) + addend;
#elif defined(__APPLE__)
	return OSAtomicAdd32Barrier(addend, (int32_t*) ptr);
#elif defined(__GNUC__)
	return __sync_add_and_fetch(ptr, addend);
#else	//fallback, slow
#pragma message("slow add_and_fetch_32")
	int32_t res;
	{
	  res = *ptr;
	  *(ptr) += addend;
	}
	return res + addend;
#endif
      }
    };

    template <>
    struct __struct_add_and_fetch<8>
    {
      template <typename Tp>
      static inline
      Tp result(volatile Tp* ptr, Tp addend)
      {
	return (Tp) __result((volatile int64_t*) ptr, (int64_t) addend);
      }

      static inline
      int64_t __result(volatile int64_t* ptr, int64_t addend)
      {
#if defined(_WIN32)
	return ::InterlockedExchangeAdd64((void*)ptr, addend) + addend;
#elif defined(__APPLE__)
	return OSAtomicAdd64Barrier(addend, (int64_t*) ptr);
#elif defined(__GNUC__)
	return __sync_add_and_fetch(ptr, addend);
#else	//fallback, slow
#pragma message("slow fetch_and_add_64")
	int64_t res;
	{
	  res = *ptr;
	  *(ptr) += addend;
	}
	return res + addend;
#endif
      }
    };

    template <typename Tp>
    inline
    Tp add_and_fetch(volatile Tp* ptr, Tp addend)
    {
      return __struct_add_and_fetch<sizeof(Tp)>::result(ptr, addend);
    }
    
    template <typename Tp>
    inline
    Tp add_and_fetch(volatile Tp& ref, Tp addend)
    {
      return __struct_add_and_fetch<sizeof(Tp)>::result(&ref, addend);
    }
    
    
    template <size_t _Size>
    struct __struct_fetch_and_add{};
    
    template <>
    struct __struct_fetch_and_add<4> 
    {
      template <typename Tp>
      static inline
      Tp result(volatile Tp* ptr, Tp addend)
      {
	return (Tp) __result((volatile int32_t*) ptr, (int32_t) addend);
      }

      static inline
      int32_t __result(volatile int32_t* ptr, int32_t addend)
      {
#if defined(_WIN32)
	return ::InterlockedExchangeAdd((void*)ptr, addend);
#elif defined(__APPLE__)
	return OSAtomicAdd32Barrier(addend, (int32_t*) ptr) - addend;
#elif defined(__GNUC__)
	return __sync_fetch_and_add(ptr, addend);
#else	//fallback, slow
#pragma message("slow fetch_and_add_32")
	int32_t res;
	{
	  res = *ptr;
	  *(ptr) += addend;
	}
	return res;
#endif
      }
    };
    
    template <>
    struct __struct_fetch_and_add<8>
    {
      template <typename Tp>
      static inline
      Tp result(volatile Tp* ptr, Tp addend)
      {
	return (Tp) __result((volatile int64_t*) ptr, (int64_t) addend);
      }

      static inline
      int64_t __result(volatile int64_t* ptr, int64_t addend)
      {
#if defined(_WIN32)
	return ::InterlockedExchangeAdd64((void*)ptr, addend);
#elif defined(__APPLE__)
	return OSAtomicAdd64Barrier(addend, (int64_t*) ptr) - addend;
#elif defined(__GNUC__)
	return __sync_fetch_and_add(ptr, addend);
#else	//fallback, slow
#pragma message("slow fetch_and_add_64")
	int64_t res;
	{
	  res = *ptr;
	  *(ptr) += addend;
	}
	return res;
#endif
      }
    };

    template <typename Tp>
    inline
    Tp fetch_and_add(volatile Tp* ptr, Tp addend)
    {
      return __struct_fetch_and_add<sizeof(Tp)>::result(ptr, addend);
    }
    
    template <typename Tp>
    inline
    Tp fetch_and_add(volatile Tp& ref, Tp addend)
    {
      return __struct_fetch_and_add<sizeof(Tp)>::result(&ref, addend);
    }

    template <size_t ByteSize>
    struct __struct_compare_and_swap {};
    
    template <>
    struct __struct_compare_and_swap<4>
    {
      template <typename Tp>
      static inline
      bool result(volatile Tp* ptr, Tp comparand, Tp replacement)
      {
	return __result((volatile int32_t*) ptr, (int32_t) comparand, (int32_t) replacement);
      }
      
      
      static inline
      bool __result(volatile int32_t* ptr, int32_t comparand, int32_t replacement)
      {
#if defined(_WIN32)
	return ::InterlockedCompareExchange((void*)ptr, replacement, comparand) == comparand;
#elif defined(__APPLE__)
	return OSAtomicCompareAndSwap32Barrier(comparand, replacement, (int32_t*) ptr);
#elif defined(__GNUC__)
	return __sync_bool_compare_and_swap(ptr, comparand, replacement);
#else
#pragma message("slow compare_and_swap_32")
	bool res = false;
	{
	  if(*ptr == comparand) {
	    *ptr = replacement;
	    res = true;
	  }
	}
	return res;
#endif
      }
    };
    
    template <>
    struct __struct_compare_and_swap<8>
    {
      template <typename Tp>
      static inline
      bool result(volatile Tp* ptr, Tp comparand, Tp replacement)
      {
	return __result((volatile int64_t*) ptr, (int64_t) comparand, (int64_t) replacement);
      }
      
      static inline
      bool __result(volatile int64_t* ptr, int64_t comparand, int64_t replacement)
      {
#if defined(_WIN32) || defined(_WIN64)
	return ::InterlockedCompareExchange64((void*)ptr, replacement, comparand) == comparand;
#elif defined(__APPLE__)
	return OSAtomicCompareAndSwap64Barrier(comparand, replacement, (int64_t*) ptr);
#elif defined(__GNUC__)
	return __sync_bool_compare_and_swap(ptr, comparand, replacement);
#else
#pragma message("slow compare_and_swap_64")
	bool res = false;
	{
	  if(*ptr == comparand) {
	    *ptr = replacement;
	    res = true;
	  }
	}
	return res;
#endif
      }
    };

    template <typename Tp>
    inline
    bool compare_and_swap(volatile Tp* ptr, Tp comparand, Tp replacement)
    {
      return __struct_compare_and_swap<sizeof(Tp)>::result(ptr, comparand, replacement);
    }
    
    template <typename Tp>
    inline
    bool compare_and_swap(volatile Tp& ref, Tp comparand, Tp replacement)
    {
      return __struct_compare_and_swap<sizeof(Tp)>::result(&ref, comparand, replacement);
    }

    template <size_t Size>
    struct __struct_ptr_cast {};
    
    template <>
    struct __struct_ptr_cast<4> {
      template <typename Tp>
      static inline
      volatile int32_t* result(volatile Tp** ptr) { return (volatile int32_t*) ptr; }
      
      template <typename Tp>
      static inline
      int32_t result(volatile Tp* ptr) { return (int32_t) ptr; }
    };
    
    template <>
    struct __struct_ptr_cast<8> {
      template <typename Tp>
      static inline
      volatile int64_t* result(volatile Tp** ptr) { return (volatile int64_t*) ptr; }
      
      template <typename Tp>
      static inline
      int64_t result(volatile Tp* ptr) { return (int64_t) ptr; }
    };

    template <typename Tp>
    inline
    bool compare_and_swap(volatile Tp** ptr, Tp* comparand, Tp* replacement)
    {
      return __struct_compare_and_swap<sizeof(Tp*)>::result(__struct_ptr_cast<sizeof(Tp*)>::result(ptr),
							    __struct_ptr_cast<sizeof(Tp*)>::result(comparand),
							    __struct_ptr_cast<sizeof(Tp*)>::result(replacement));
    }
    
  };
};

#endif
