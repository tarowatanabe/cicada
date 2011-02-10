//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iterator>

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <utils/config.hpp>
#include <utils/thread_specific_ptr.hpp>
#include <utils/lexical_cast.hpp>

#include "symbol.hpp"

#include "vocab.hpp"


namespace cicada
{
  struct SymbolImpl
  {
    typedef Symbol::symbol_map_type               symbol_map_type;
    
    typedef std::vector<int, std::allocator<int> >   index_map_type;
    typedef std::vector<bool, std::allocator<bool> > non_terminal_map_type;
    typedef std::vector<Symbol::id_type, std::allocator<Symbol::id_type> > non_terminal_symbol_map_type;
  };
  
  Symbol::mutex_type    Symbol::__mutex_index;
  Symbol::mutex_type    Symbol::__mutex_data;
  

#ifdef HAVE_TLS
  static __thread SymbolImpl::symbol_map_type*              symbol_maps_tls = 0;
  static __thread SymbolImpl::index_map_type*               index_maps_tls = 0;
  static __thread SymbolImpl::non_terminal_map_type*        non_terminal_maps_tls = 0;
  static __thread SymbolImpl::non_terminal_symbol_map_type* non_terminal_symbol_maps_tls = 0;

  static boost::thread_specific_ptr<SymbolImpl::symbol_map_type>              symbol_maps;
  static boost::thread_specific_ptr<SymbolImpl::index_map_type>               index_maps;
  static boost::thread_specific_ptr<SymbolImpl::non_terminal_map_type>        non_terminal_maps;
  static boost::thread_specific_ptr<SymbolImpl::non_terminal_symbol_map_type> non_terminal_symbol_maps;
#else
  static utils::thread_specific_ptr<SymbolImpl::symbol_map_type>              symbol_maps;
  static utils::thread_specific_ptr<SymbolImpl::index_map_type>               index_maps;
  static utils::thread_specific_ptr<SymbolImpl::non_terminal_map_type>        non_terminal_maps;
  static utils::thread_specific_ptr<SymbolImpl::non_terminal_symbol_map_type> non_terminal_symbol_maps;
#endif


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

  Symbol::piece_type Symbol::non_terminal_strip() const
  {
    if (! is_non_terminal())
      return *this;
    
    piece_type stripped(static_cast<const std::string&>(non_terminal()));
    
    return stripped.substr(1, stripped.size() - 2);
  }

  
  int Symbol::non_terminal_index() const
  {
    if (! is_non_terminal())
      return 0;
    
#ifdef HAVE_TLS
    if (! index_maps_tls) {
      index_maps.reset(new SymbolImpl::index_map_type());
      index_maps_tls = index_maps.get();
    }
    
    SymbolImpl::index_map_type& maps = *index_maps_tls;
#else
    if (! index_maps.get())
      index_maps.reset(new SymbolImpl::index_map_type());
    
    SymbolImpl::index_map_type& maps = *index_maps;
#endif
    
    if (__id >= maps.size())
      maps.resize(__id + 1, -1);

    if (maps[__id] < 0) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      const symbol_type& word = symbol();
	
      // at leas we have [<char>,<digit>]
      maps[__id] = 0;

      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();
      
      qi::parse(iter, iter_end,
		qi::lit('[')
		>> +(standard::char_ - ',' - ']') >> ',' >> qi::int_ [ phoenix::ref(maps[__id]) = qi::_1 ]
		>> qi::lit(']'));
    }
      
    return maps[__id];
  }
  
  bool Symbol::is_non_terminal() const
  {
#ifdef HAVE_TLS
    if (! non_terminal_maps_tls) {
      non_terminal_maps.reset(new SymbolImp::non_terminal_map_type());
      non_terminal_maps_tls = non_terminal_maps.get();
    }
    
    SymbolImpl::non_terminal_map_type& maps =  *non_terminal_maps_tls;
#else
    if (! non_terminal_maps.get())
      non_terminal_maps.reset(new SymbolImpl::non_terminal_map_type());
    
    SymbolImpl::non_terminal_map_type& maps =  *non_terminal_maps;
#endif
    
    const size_type scan_pos = (__id << 1);
    const size_type flag_pos = (__id << 1) + 1;
    
    if (flag_pos >= maps.size())
      maps.resize(flag_pos + 1, false);
    
    if (! maps[scan_pos]) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      const symbol_type& word = symbol();
      
      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();
      
      maps[scan_pos] = true;
      maps[flag_pos] = qi::parse(iter, iter_end, '[' >> +(standard::char_ - ',' - ']') >> -(',' >> qi::int_) >> ']') && iter == iter_end;
    }
    
    return maps[flag_pos];
  }

  Symbol Symbol::non_terminal(const int index) const
  {
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;

    if (! is_non_terminal())
      return *this;
    if (index <= 0)
      return non_terminal();
    
    std::string generated;
    std::back_insert_iterator<std::string> iter(generated);
    
    karma::generate(iter, '[' << +(standard::char_) << ',' << karma::int_ << ']', non_terminal_strip(), index);
    
    return Symbol(generated);
  }
  
  Symbol Symbol::non_terminal() const
  {
    if (! is_non_terminal())
      return *this;

#ifdef HAVE_TLS
    if (! non_terminal_symbol_maps_tls) {
      non_terminal_symbol_maps.reset(new SymbolImpl::non_terminal_symbol_map_type());
      non_terminal_symbol_maps_tls = non_terminal_symbol_maps.get();
    }
    
    SymbolImpl::non_terminal_symbol_map_type& maps = *non_terminal_symbol_maps_tls;
#else
    if (! non_terminal_symbol_maps.get())
      non_terminal_symbol_maps.reset(new SymbolImpl::non_terminal_symbol_map_type());
    
    SymbolImpl::non_terminal_symbol_map_type& maps = *non_terminal_symbol_maps;
#endif
    
    if (__id >= maps.size())
      maps.resize(__id + 1, id_type(-1));
      
    if (maps[__id] == id_type(-1)) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;

      const symbol_type& word = symbol();
      
      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();

      maps[__id] = __id;
      
      std::string label;
      if (qi::parse(iter, iter_end,
		    qi::lit('[')
		    >> +(standard::char_ - ',' - ']') [ phoenix::ref(label) = qi::_1 ] >> -(',' >> qi::int_)
		    >> qi::lit(']')))
	maps[__id] = Symbol('[' + label + ']').id();
    }
      
    return maps[__id];
  }

  void Symbol::write(const path_type& path)
  {
    lock_type lock(__mutex_data);
    
    Vocab vocab(path, std::max(size_type(__symbols().size() / 2), size_type(1024)));
    symbol_set_type::const_iterator siter_end = __symbols().end();
    for (symbol_set_type::const_iterator siter = __symbols().begin(); siter != siter_end; ++ siter)
      vocab.insert(*siter);
  }
  
};
