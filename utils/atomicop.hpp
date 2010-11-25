// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
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

#if defined(__SUNPRO_CC) && defined(__sparc)
#include <sys/atomic.h>
#endif


#if defined(_MSC_VER)
#include <Windows.h>
#include <intrin.h>
#endif


namespace utils
{
  namespace atomicop
  {
    inline
    void memory_barrier()
    {
#if defined(__GNUC__) && ( (__GNUC__ > 4) || ((__GNUC__ >= 4) && (__GNUC_MINOR__ >= 1)) )
      __sync_synchronize();
#elif defined(_MSC_VER) && (_MSC_VER >= 1300)
      _ReadWriteBarrier();
#elif defined(__APPLE__)
      OSMemoryBarrier();
#elif defined(AO_HAVE_nop_full)
      AO_nop_full(); // for libatomic_ops from GC...
#else
#   warning "no memory barrier implemented for this platform"
#endif
    }
    
    
    
#if defined(__ICC) 
    template<typename must_be_int = int>
    inline
    int32_t __faa32(int32_t* x, int32_t inc)
    {
      asm volatile ( 
		    "lock xadd %0,%1" 
		    : "=r" (inc), "=m" (*x) 
		    : "0" (inc) 
		    : "memory");
      return inc;
    }
#if defined(__x86_64)
    template<typename must_be_int = int>
    inline
    int64_t __faa64(int64_t* x, int64_t inc)
    {
      asm volatile ( 
		    "lock xadd %0,%1" 
		    : "=r" (inc), "=m" (*x) 
		    : "0" (inc) 
		    : "memory");
      return inc;
    }
#endif
#endif /* __ICC */
    
    
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
#if defined(__ICC)
	return _InterlockedExchangeAdd((void*)ptr, addend);
#elif defined(__ECC)
	return _InterlockedExchangeAdd((void*)ptr, addend);
#elif defined(__ICL) || defined(_MSC_VER)
	return _InterlockedExchangeAdd(reinterpret_cast<volatile long*>(ptr), addend);
#elif defined(__APPLE__)
	return OSAtomicAdd32Barrier(addend, (int32_t*) ptr);
#elif defined(__GNUC__)
	return __sync_fetch_and_add(ptr, addend);
#elif defined(__SUNPRO_CC) && defined(__sparc)
	volatile int32 before, after;
	do
	  {
	    before = *ptr;
	    after = before + addend;
	  } while(atomic_cas_32((volatile unsigned int*)ptr, before, after) != before);
	return before;
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
#if defined(__ICC) && defined(__x86_64)
	return __faa64<int>((int64*)ptr, addend);
#elif defined(__ECC)
	return _InterlockedExchangeAdd64((void*)ptr, addend);
#elif defined(__ICL) || defined(_MSC_VER)
#ifndef _WIN64
	assert(false);	//not available in this case
	return 0;
#else
	return _InterlockedExchangeAdd64(ptr, addend);
#endif
	//#elif defined(_WIN32)
	//  return _InterlockedExchangeAdd64(reinterpret_cast<volatile long long*>(ptr), addend);
#elif defined(__APPLE__)
	return OSAtomicAdd64Barrier(addend, (int64_t*) ptr);
#elif defined(__GNUC__) && defined(__x86_64)
	return __sync_fetch_and_add(ptr, addend);
#elif defined(__GNUC__) && defined(__i386) && (defined(__i686) || defined(__pentium4) || defined(__athlon))
	return __sync_fetch_and_add(ptr, addend);
#elif defined(__SUNPRO_CC) && defined(__sparc)
	volatile int64 before, after;
	do
	  {
	    before = *ptr;
	    after = before + addend;
	  } while(atomic_cas_64((volatile unsigned long long*)ptr, before, after) != before);
	return before;
#else	//fallback, slow
#if defined(__GNUC__) && defined(__i386)
#warning "please compile with -march=i686 or better"
#endif
#pragma message("slow fetch_and_add_64")
	int64 res;
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
    
#if defined(__ICC)

    template<typename must_be_int = int>
    inline 
    int32_t __cas32(volatile int32_t* ptr, int32_t old, int32_t nw)
    {
      int32_t before;
      __asm__ __volatile__("lock; cmpxchgl %1,%2"
			   : "=a"(before)
			   : "q"(nw), "m"(*(volatile long long*)(ptr)), "0"(old)
			   : "memory");
      return before;
    }
#if defined(__x86_64)
    template<typename must_be_int = int>
    inline 
    int64_t __cas64(volatile int64_t *ptr, int64_t old, int64_t nw)
    {
      int64_t before;
      __asm__ __volatile__("lock; cmpxchgq %1,%2"
			   : "=a"(before)
			   : "q"(nw), "m"(*(volatile long long*)(ptr)), "0"(old)
			   : "memory");
      return before;
    }
#endif
#endif

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
#if defined(__ICC)	//x86 version
	return _InterlockedCompareExchange((void*)ptr, replacement, comparand) == comparand;
#elif defined(__ECC)	//IA-64 version
	return _InterlockedCompareExchange((void*)ptr, replacement, comparand) == comparand;
#elif defined(__ICL) || defined(_MSC_VER)
	return _InterlockedCompareExchange(reinterpret_cast<volatile long*>(ptr), replacement, comparand) == comparand;
#elif defined(__APPLE__)
	return OSAtomicCompareAndSwap32Barrier(comparand, replacement, (int32_t*) ptr);
#elif defined(__GNUC__)
	return __sync_bool_compare_and_swap(ptr, comparand, replacement);
#elif defined(__SUNPRO_CC) && defined(__sparc)
	return atomic_cas_32((volatile unsigned int*)ptr, comparand, replacement) == comparand;
#elif defined(AO_HAVE_compare_and_swap_full)
	return AO_compare_and_swap_full(reinterpret_cast<volatile AO_t*>(ptr),
					reinterpret_cast<AO_t>(comparand),
					reinterpret_cast<AO_t>(replacement));
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
#if defined(__ICC) && defined(__x86_64)	//x86 version
	return __cas64<int>(ptr, comparand, replacement) == comparand;
#elif defined(__ECC)	//IA-64 version
	return _InterlockedCompareExchange64((void*)ptr, replacement, comparand) == comparand;
#elif defined(__ICL) || defined(_MSC_VER)
#ifndef _WIN64
	assert(false);	//not available in this case
	return 0;
#else
	return _InterlockedCompareExchange64(ptr, replacement, comparand) == comparand;
#endif
	//#elif defined(_WIN32)
	//	return _InterlockedCompareExchange64(reinterpret_cast<volatile __int64*>(ptr), replacement, comparand) == comparand;
#elif defined(__APPLE__)
	return OSAtomicCompareAndSwap64Barrier(comparand, replacement, (int64_t*) ptr);
#elif defined(__GNUC__) && defined(__x86_64)
	return __sync_bool_compare_and_swap(ptr, comparand, replacement);
#elif defined(__GNUC__) && defined(__i386) && (defined(__i686) || defined(__pentium4) || defined(__athlon))
	return __sync_bool_compare_and_swap(ptr, comparand, replacement);
#elif defined(__SUNPRO_CC) && defined(__sparc)
	return atomic_cas_64((volatile unsigned long long*)ptr, comparand, replacement) == comparand;
#else
#if defined(__GNUC__) && defined(__i386)
#warning "please compile with -march=i686 or better"
#endif
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
