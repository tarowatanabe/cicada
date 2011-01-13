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
    typedef Symbol::symbol_map_type               symbol_map_type;
    typedef Symbol::index_map_type                index_map_type;
    typedef Symbol::non_terminal_map_type         non_terminal_map_type;
    typedef Symbol::non_terminal_symbol_map_type  non_terminal_symbol_map_type;
  };
  
  Symbol::mutex_type    Symbol::__mutex;
  
  Symbol::symbol_set_type& Symbol::__symbols()
  {
    static symbol_set_type symbols;
    
    return symbols;
  }
  

#ifdef HAVE_TLS
  static __thread SymbolImpl::symbol_map_type*              symbol_maps_tls = 0;
  static __thread SymbolImpl::index_map_type*               index_maps_tls = 0;
  static __thread SymbolImpl::non_terminal_map_type*        non_terminal_maps_tls = 0;
  static __thread SymbolImpl::non_terminal_symbol_map_type* non_terminal_symbol_maps_tls = 0;
#endif
  static boost::thread_specific_ptr<SymbolImpl::symbol_map_type>              symbol_maps;
  static boost::thread_specific_ptr<SymbolImpl::index_map_type>               index_maps;
  static boost::thread_specific_ptr<SymbolImpl::non_terminal_map_type>        non_terminal_maps;
  static boost::thread_specific_ptr<SymbolImpl::non_terminal_symbol_map_type> non_terminal_symbol_maps;

  Symbol::symbol_map_type& Symbol::__symbol_maps()
  {
#ifdef HAVE_TLS
    if (! symbol_maps_tls) {
      symbol_maps.reset(new symbol_map_type());
      symbol_maps->reserve(allocated());
      symbol_maps_tls = symbol_maps.get();
    }
    
    return *symbol_maps_tls;
#else
    if (! symbol_maps.get()) {
      symbol_maps.reset(new symbol_map_type());
      symbol_maps->reserve(allocated());
    }
    
    return *symbol_maps;
#endif
  }

  
  Symbol::index_map_type& Symbol::__index_maps()
  {
#ifdef HAVE_TLS
    if (! index_maps_tls) {
      index_maps.reset(new index_map_type());
      index_maps_tls = index_maps.get();
    }
    
    return *index_maps_tls;
#else
    if (! index_maps.get())
      index_maps.reset(new index_map_type());
    
    return *index_maps;
#endif
  }
  
  Symbol::non_terminal_map_type& Symbol::__non_terminal_maps()
  {
#ifdef HAVE_TLS
    if (! non_terminal_maps_tls) {
      non_terminal_maps.reset(new non_terminal_map_type());
      non_terminal_maps_tls = non_terminal_maps.get();
    }
    
    return *non_terminal_maps_tls;
#else
    if (! non_terminal_maps.get())
      non_terminal_maps.reset(new non_terminal_map_type());
    
    return *non_terminal_maps;
#endif
  }
  
  Symbol::non_terminal_symbol_map_type& Symbol::__non_terminal_symbol_maps()
  {
#ifdef HAVE_TLS
    if (! non_terminal_symbol_maps_tls) {
      non_terminal_symbol_maps.reset(new non_terminal_symbol_map_type());
      non_terminal_symbol_maps_tls = non_terminal_symbol_maps.get();
    }
    
    return *non_terminal_symbol_maps_tls;
#else
    if (! non_terminal_symbol_maps.get())
      non_terminal_symbol_maps.reset(new non_terminal_symbol_map_type());
    
    return *non_terminal_symbol_maps;
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
