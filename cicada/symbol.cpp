//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iterator>

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/xpressive/xpressive.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/config.hpp>
#include <utils/thread_specific_ptr.hpp>

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
    typedef std::vector<Symbol::id_type, std::allocator<Symbol::id_type> > pos_symbol_map_type;
    typedef std::vector<Symbol::id_type, std::allocator<Symbol::id_type> > terminal_symbol_map_type;
    typedef std::vector<Symbol::id_type, std::allocator<Symbol::id_type> > coarse_symbol_map_type;
  };
  
  Symbol::mutex_type    Symbol::__mutex_index;
  Symbol::mutex_type    Symbol::__mutex_data;
  

#ifdef HAVE_TLS
  static __thread SymbolImpl::symbol_map_type*              symbol_maps_tls = 0;
  static __thread SymbolImpl::index_map_type*               index_maps_tls = 0;
  static __thread SymbolImpl::non_terminal_map_type*        non_terminal_maps_tls = 0;
  static __thread SymbolImpl::non_terminal_symbol_map_type* non_terminal_symbol_maps_tls = 0;
  static __thread SymbolImpl::pos_symbol_map_type*          pos_symbol_maps_tls = 0;
  static __thread SymbolImpl::terminal_symbol_map_type*     terminal_symbol_maps_tls = 0;
  static __thread SymbolImpl::coarse_symbol_map_type*       coarse_symbol_maps_tls = 0;

  static boost::thread_specific_ptr<SymbolImpl::symbol_map_type>              symbol_maps;
  static boost::thread_specific_ptr<SymbolImpl::index_map_type>               index_maps;
  static boost::thread_specific_ptr<SymbolImpl::non_terminal_map_type>        non_terminal_maps;
  static boost::thread_specific_ptr<SymbolImpl::non_terminal_symbol_map_type> non_terminal_symbol_maps;
  static boost::thread_specific_ptr<SymbolImpl::pos_symbol_map_type>          pos_symbol_maps;
  static boost::thread_specific_ptr<SymbolImpl::terminal_symbol_map_type>     terminal_symbol_maps;
  static boost::thread_specific_ptr<SymbolImpl::coarse_symbol_map_type>       coarse_symbol_maps;
#else
  static utils::thread_specific_ptr<SymbolImpl::symbol_map_type>              symbol_maps;
  static utils::thread_specific_ptr<SymbolImpl::index_map_type>               index_maps;
  static utils::thread_specific_ptr<SymbolImpl::non_terminal_map_type>        non_terminal_maps;
  static utils::thread_specific_ptr<SymbolImpl::non_terminal_symbol_map_type> non_terminal_symbol_maps;
  static utils::thread_specific_ptr<SymbolImpl::pos_symbol_map_type>          pos_symbol_maps;
  static utils::thread_specific_ptr<SymbolImpl::terminal_symbol_map_type>     terminal_symbol_maps;
  static utils::thread_specific_ptr<SymbolImpl::coarse_symbol_map_type>       coarse_symbol_maps;
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
      
      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();
      
      if (! qi::parse(iter, iter_end,
		      '[' >> +(standard::char_ - ',' - ']') >> ',' >> qi::int_ [ phoenix::ref(maps[__id]) = qi::_1 ] >> ']'))
	maps[__id] = 0;
    }
      
    return maps[__id];
  }
  
  bool Symbol::is_non_terminal() const
  {
#ifdef HAVE_TLS
    if (! non_terminal_maps_tls) {
      non_terminal_maps.reset(new SymbolImpl::non_terminal_map_type());
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

    if (! is_non_terminal())
      return *this;
    if (index <= 0)
      return non_terminal();
    
    std::string generated;
    std::back_insert_iterator<std::string> iter(generated);
    
    karma::generate(iter, '[' << standard::string << ',' << karma::int_ << ']', non_terminal_strip(), index);
    
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

      const symbol_type& word = symbol();

      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();

      maps[__id] = __id;
      
      std::string label;
      int index;
      if (qi::parse(iter, iter_end,
		    qi::lit('[')
		    >> (+(standard::char_ - ',' - ']')) >> -(',' >> qi::int_)
		    >> qi::lit(']'), label, index))
	maps[__id] = Symbol('[' + label + ']').id();
    }
      
    return maps[__id];
  }

  Symbol Symbol::pos() const
  {
    if (! is_terminal()) return Symbol();
    
#ifdef HAVE_TLS
    if (! pos_symbol_maps_tls) {
      pos_symbol_maps.reset(new SymbolImpl::pos_symbol_map_type());
      pos_symbol_maps_tls = pos_symbol_maps.get();
    }
    
    SymbolImpl::pos_symbol_map_type& maps = *pos_symbol_maps_tls;
#else
    if (! pos_symbol_maps.get())
      pos_symbol_maps.reset(new SymbolImpl::pos_symbol_map_type());
    
    SymbolImpl::pos_symbol_map_type& maps = *pos_symbol_maps;
#endif

    if (__id >= maps.size())
      maps.resize(__id + 1, id_type(-1));
    
    if (maps[__id] == id_type(-1)) {
      namespace xpressive = boost::xpressive;
      
      typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
      typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
      
      static pregex re = (+(~xpressive::_s)) >> (xpressive::set= '|', '/') >> (xpressive::s1= '[' >> -+(~(xpressive::set= ']')) >> ']');
      
      pmatch what;
      if (xpressive::regex_match(utils::piece(symbol()), what, re))
	maps[__id] = Symbol(what[1]).id();
      else
	maps[__id] = Symbol().id();
    }
    return maps[__id];
  }
  
  Symbol Symbol::terminal() const
  {
    if (! is_terminal()) return *this;
    
#ifdef HAVE_TLS
    if (! terminal_symbol_maps_tls) {
      terminal_symbol_maps.reset(new SymbolImpl::terminal_symbol_map_type());
      terminal_symbol_maps_tls = terminal_symbol_maps.get();
    }
    
    SymbolImpl::terminal_symbol_map_type& maps = *terminal_symbol_maps_tls;
#else
    if (! terminal_symbol_maps.get())
      terminal_symbol_maps.reset(new SymbolImpl::terminal_symbol_map_type());
    
    SymbolImpl::terminal_symbol_map_type& maps = *terminal_symbol_maps;
#endif
    
    if (__id >= maps.size())
      maps.resize(__id + 1, id_type(-1));

    if (maps[__id] == id_type(-1)) {
      namespace xpressive = boost::xpressive;
      
      typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
      typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
      
      static pregex re = (xpressive::s1= +(~xpressive::_s)) >> (xpressive::set= '|', '/') >> ('[' >> -+(~(xpressive::set= ']')) >> ']');
      
      pmatch what;
      if (xpressive::regex_match(utils::piece(symbol()), what, re))
	maps[__id] = Symbol(what[1]).id();
      else
	maps[__id] = __id;
    }
    return maps[__id];
  }

  bool Symbol::binarized() const
  {
    if (! is_non_terminal()) return false;
    
    return non_terminal_strip().find('^') != piece_type::npos();
  }

  Symbol Symbol::annotation(const int pos, const bool bit) const
  {
    if (! is_non_terminal()) return *this;
    
    namespace xpressive = boost::xpressive;
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
    typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
    
    static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (xpressive::s2= -+xpressive::_d);
    
    const int __non_terminal_index = non_terminal_index();
    const piece_type piece = non_terminal_strip();
    const int mask = 1 << pos;
    
    std::string generated;
    std::back_insert_iterator<std::string> iter(generated);
    
    pmatch what;
    if (xpressive::regex_match(piece, what, re)) {
      const int value = (utils::lexical_cast<int>(what[2]) & (~mask)) | (-bit & mask);
      karma::generate(iter, '[' << standard::string << '@' << karma::int_ << ']', utils::piece(what[1]), value);
    } else
      karma::generate(iter, '[' << standard::string << '@' << karma::int_ << ']', piece, -bit & mask);
    
    return Symbol(generated).non_terminal(__non_terminal_index);
  }
  
  Symbol Symbol::coarse(const int pos) const
  {
    if (! is_non_terminal()) return *this;
    
    // even coarser!
    if (pos < 0) return coarse();
    
    namespace xpressive = boost::xpressive;
    namespace qi = boost::spirit::qi;
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
      
    typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
    typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
    
    static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (xpressive::s2= -+xpressive::_d);
    
    const int __non_terminal_index = non_terminal_index();
    const piece_type piece = non_terminal_strip();
    const int mask = (1 << pos) - 1;

    std::string generated;
    std::back_insert_iterator<std::string> iter(generated);
    
    pmatch what;
    if (xpressive::regex_match(piece, what, re)) {
      const int value = (utils::lexical_cast<int>(what[2]) & mask);
      
      karma::generate(iter, '[' << standard::string << '@' << karma::int_ << ']', utils::piece(what[1]), value);
    } else
      karma::generate(iter, '[' << standard::string << '@' << karma::int_ << ']', piece, 0);
    
    return Symbol(generated).non_terminal(__non_terminal_index);
  }

  Symbol Symbol::coarse() const
  {
    if (! is_non_terminal()) return *this;

#ifdef HAVE_TLS
    if (! coarse_symbol_maps_tls) {
      coarse_symbol_maps.reset(new SymbolImpl::coarse_symbol_map_type());
      coarse_symbol_maps_tls = coarse_symbol_maps.get();
    }
    
    SymbolImpl::coarse_symbol_map_type& maps = *coarse_symbol_maps_tls;
#else
    if (! coarse_symbol_maps.get())
      coarse_symbol_maps.reset(new SymbolImpl::coarse_symbol_map_type());
    
    SymbolImpl::coarse_symbol_map_type& maps = *coarse_symbol_maps;
#endif
    
    if (__id >= maps.size())
      maps.resize(__id + 1, id_type(-1));
    
    if (maps[__id] == id_type(-1)) {
      namespace xpressive = boost::xpressive;
      
      typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
      typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
      
      static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (-+xpressive::_d);
      
      const int __non_terminal_index = non_terminal_index();
      const piece_type piece = non_terminal_strip();
      
      pmatch what;
      if (xpressive::regex_match(piece, what, re))
	maps[__id] = Symbol('[' + what[1] + ']').non_terminal(__non_terminal_index).id();
      else
	maps[__id] = __id;
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
