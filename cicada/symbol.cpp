//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utils/config.hpp>

#include "symbol.hpp"

#include "vocab.hpp"


namespace cicada
{
  struct SymbolImpl
  {
    typedef Symbol::symbol_set_type              symbol_set_type;
    typedef Symbol::symbol_map_type              symbol_map_type;
    typedef Symbol::index_map_type               index_map_type;
    typedef Symbol::non_terminal_map_type        non_terminal_map_type;
    typedef Symbol::non_terminal_symbol_map_type non_terminal_symbol_map_type;
    
    static boost::once_flag once;

    static Symbol::symbol_set_type* symbols; 
   
#ifdef HAVE_TLS
    static __thread symbol_map_type*              symbol_maps_tls;
    static __thread index_map_type*               index_maps_tls;
    static __thread non_terminal_map_type*        non_terminal_maps_tls;
    static __thread non_terminal_symbol_map_type* non_terminal_symbol_maps_tls;
#endif
    
    static boost::thread_specific_ptr<symbol_map_type>              symbol_maps;
    static boost::thread_specific_ptr<index_map_type>               index_maps;
    static boost::thread_specific_ptr<non_terminal_map_type>        non_terminal_maps;
    static boost::thread_specific_ptr<non_terminal_symbol_map_type> non_terminal_symbol_maps;

    static void initialize()
    {
      symbols = new symbol_set_type();
    }
  };
  
  
  boost::once_flag SymbolImpl::once = BOOST_ONCE_INIT;
  
  SymbolImpl::symbol_set_type* SymbolImpl::symbols = 0;
  
#ifdef HAVE_TLS
  __thread SymbolImpl::symbol_map_type*              SymbolImpl::symbol_maps_tls = 0;
  __thread SymbolImpl::index_map_type*               SymbolImpl::index_maps_tls = 0;
  __thread SymbolImpl::non_terminal_map_type*        SymbolImpl::non_terminal_maps_tls = 0;
  __thread SymbolImpl::non_terminal_symbol_map_type* SymbolImpl::non_terminal_symbol_maps_tls = 0;
#endif
  
  boost::thread_specific_ptr<SymbolImpl::symbol_map_type>              SymbolImpl::symbol_maps;
  boost::thread_specific_ptr<SymbolImpl::index_map_type>               SymbolImpl::index_maps;
  boost::thread_specific_ptr<SymbolImpl::non_terminal_map_type>        SymbolImpl::non_terminal_maps;
  boost::thread_specific_ptr<SymbolImpl::non_terminal_symbol_map_type> SymbolImpl::non_terminal_symbol_maps;
  
  
  Symbol::mutex_type    Symbol::__mutex;
  
  Symbol::symbol_set_type& Symbol::__symbols()
  {
    boost::call_once(SymbolImpl::once, SymbolImpl::initialize);
    
    return *SymbolImpl::symbols;
  }
  

  Symbol::symbol_map_type& Symbol::__symbol_maps()
  {
#ifdef HAVE_TLS
    if (! SymbolImpl::symbol_maps_tls) {
      SymbolImpl::symbol_maps.reset(new symbol_map_type());
      SymbolImpl::symbol_maps->reserve(allocated());
      
      SymbolImpl::symbol_maps_tls = SymbolImpl::symbol_maps.get();
    }
    
    return *SymbolImpl::symbol_maps_tls;
#else
    if (! SymbolImpl::symbol_maps.get()) {
      SymbolImpl::symbol_maps.reset(new symbol_map_type());
      SymbolImpl::symbol_maps->reserve(allocated());
    }
    
    return *SymbolImpl::symbol_maps;
#endif
  }

  
  Symbol::index_map_type& Symbol::__index_maps()
  {
#ifdef HAVE_TLS
    if (! SymbolImpl::index_maps_tls) {
      SymbolImpl::index_maps.reset(new index_map_type());
      SymbolImpl::index_maps->reserve(allocated());
      
      SymbolImpl::index_maps_tls = SymbolImpl::index_maps.get();
    }
    
    return *SymbolImpl::index_maps_tls;
#else
    if (! SymbolImpl::index_maps.get()) {
      SymbolImpl::index_maps.reset(new index_map_type());
      SymbolImpl::index_maps->reserve(allocated());
    }
    
    return *SymbolImpl::index_maps;
#endif
  }
  
  Symbol::non_terminal_map_type& Symbol::__non_terminal_maps()
  {
#ifdef HAVE_TLS
    if (! SymbolImpl::non_terminal_maps_tls) {
      SymbolImpl::non_terminal_maps.reset(new non_terminal_map_type());
      SymbolImpl::non_terminal_maps->reserve(allocated());
      
      SymbolImpl::non_terminal_maps_tls = SymbolImpl::non_terminal_maps.get();
    }
    
    return *SymbolImpl::non_terminal_maps_tls;
#else
    if (! SymbolImpl::non_terminal_maps.get()) {
      SymbolImpl::non_terminal_maps.reset(new non_terminal_map_type());
      SymbolImpl::non_terminal_maps->reserve(allocated());
    }
    
    return *SymbolImpl::non_terminal_maps;
#endif
  }
  
  Symbol::non_terminal_symbol_map_type& Symbol::__non_terminal_symbol_maps()
  {
#ifdef HAVE_TLS
    if (! SymbolImpl::non_terminal_symbol_maps_tls) {
      SymbolImpl::non_terminal_symbol_maps.reset(new non_terminal_symbol_map_type());
      SymbolImpl::non_terminal_symbol_maps->reserve(allocated());
      
      SymbolImpl::non_terminal_symbol_maps_tls = SymbolImpl::non_terminal_symbol_maps.get();
    }
    
    return *SymbolImpl::non_terminal_symbol_maps_tls;
#else
    if (! SymbolImpl::non_terminal_symbol_maps.get()) {
      SymbolImpl::non_terminal_symbol_maps.reset(new non_terminal_symbol_map_type());
      SymbolImpl::non_terminal_symbol_maps->reserve(allocated());
    }
    
    return *SymbolImpl::non_terminal_symbol_maps;
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
