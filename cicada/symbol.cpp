//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "symbol.hpp"

#include "vocab.hpp"


namespace cicada
{
    
  Symbol::mutex_type    Symbol::__mutex;

  class __symbol_set_instance
  {
  public:
    __symbol_set_instance() : instance(0) {}
    ~__symbol_set_instance() { if (instance) delete instance; }
    
    Symbol::symbol_set_type* instance;
  };
  
  static boost::once_flag      __symbols_once = BOOST_ONCE_INIT;
  static __symbol_set_instance __symbols_instance;
  
  static void __symbols_init()
  {
    __symbols_instance.instance = new Symbol::symbol_set_type();
  }
  
  Symbol::symbol_set_type& Symbol::__symbols()
  {
    boost::call_once(__symbols_once, __symbols_init);

    return *__symbols_instance.instance;
  }
  
  Symbol::symbol_map_type& Symbol::__symbol_maps()
  {
#ifdef HAVE_TLS
    static __thread symbol_map_type* __maps_tls = 0;
    static boost::thread_specific_ptr<symbol_map_type> __maps;
      
    if (! __maps_tls) {
      __maps.reset(new symbol_map_type());
      __maps->reserve(allocated());
	
      __maps_tls = __maps.get();
    }
      
    return *__maps_tls;
#else
    static boost::thread_specific_ptr<symbol_map_type> __maps;
      
    if (! __maps.get()) {
      __maps.reset(new symbol_map_type());
      __maps->reserve(allocated());
    }
      
    return *__maps;
#endif
  }

  Symbol::index_map_type& Symbol::__index_maps()
  {
#ifdef HAVE_TLS
    static __thread index_map_type* __maps_tls = 0;
    static boost::thread_specific_ptr<index_map_type> __maps;
      
    if (! __maps_tls) {
      __maps.reset(new index_map_type());
      __maps->reserve(allocated());
      
      __maps_tls = __maps.get();
    }
      
    return *__maps_tls;
#else
    static boost::thread_specific_ptr<index_map_type> __maps;
      
    if (! __maps.get()) {
      __maps.reset(new index_map_type());
      __maps->reserve(allocated());
    }
      
    return *__maps;
#endif
  }

  Symbol::non_terminal_map_type& Symbol::__non_terminal_maps()
  {
#ifdef HAVE_TLS
    static __thread non_terminal_map_type* __maps_tls = 0;
    static boost::thread_specific_ptr<non_terminal_map_type> __maps;
      
    if (! __maps_tls) {
      __maps.reset(new non_terminal_map_type());
      __maps->reserve(allocated());
	
      __maps_tls = __maps.get();
    }
      
    return *__maps_tls;
#else
    static boost::thread_specific_ptr<non_terminal_map_type> __maps;
      
    if (! __maps.get()) {
      __maps.reset(new non_terminal_map_type());
      __maps->reserve(allocated());
    }
      
    return *__maps;
#endif
  }

  Symbol::non_terminal_symbol_map_type& Symbol::__non_terminal_symbol_maps()
  {
#ifdef HAVE_TLS
    static __thread non_terminal_symbol_map_type* __maps_tls = 0;
    static boost::thread_specific_ptr<non_terminal_symbol_map_type> __maps;
      
    if (! __maps_tls) {
      __maps.reset(new non_terminal_symbol_map_type());
      __maps->reserve(allocated());
	
      __maps_tls = __maps.get();
    }
      
    return *__maps_tls;
#else
    static boost::thread_specific_ptr<non_terminal_symbol_map_type> __maps;
      
    if (! __maps.get()) {
      __maps.reset(new non_terminal_symbol_map_type());
      __maps->reserve(allocated());
    }
      
    return *__maps;
#endif
  }


  void Symbol::write(const path_type& path)
  {
    lock_type lock(__mutex);
    
    Vocab vocab(path, std::max(size_type(__symbols().size() / 2), size_type(1024)));
    symbol_set_type::const_iterator siter_end = __symbols().end();
    for (symbol_set_type::const_iterator siter = __symbols().begin(); siter != siter_end; ++ siter)
      vocab.insert(*siter);
  }
  
};
