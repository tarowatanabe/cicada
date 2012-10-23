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
#include <utils/simple_vector.hpp>

#include "symbol.hpp"

#include "vocab.hpp"


namespace cicada
{
  struct SymbolImpl
  {
    typedef Symbol::symbol_map_type               symbol_map_type;
    typedef Symbol::id_type id_type;
    typedef Symbol::mutex_type mutex_type;
    
    typedef utils::indexed_set<id_type, utils::hashmurmur<size_t>, std::equal_to<id_type>, std::allocator<id_type> > non_terminal_set_type;
    
    typedef std::vector<int, std::allocator<int> >   index_map_type;
    typedef std::vector<bool, std::allocator<bool> > non_terminal_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > non_terminal_id_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > non_terminal_symbol_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > pos_symbol_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > terminal_symbol_map_type;
    typedef std::vector<id_type, std::allocator<id_type> > coarse_symbol_map_type;
    
    typedef utils::simple_vector<id_type, std::allocator<id_type> > id_set_type;
    typedef std::vector<id_set_type, std::allocator<id_set_type> >  coarser_symbol_map_type;
    
    symbol_map_type              symbol_maps;
    index_map_type               index_maps;
    non_terminal_map_type        non_terminal_maps;
    non_terminal_id_map_type     non_terminal_id_maps;
    non_terminal_symbol_map_type non_terminal_symbol_maps;
    pos_symbol_map_type          pos_symbol_maps;
    terminal_symbol_map_type     terminal_symbol_maps;
    coarse_symbol_map_type       coarse_symbol_maps;
    coarser_symbol_map_type      coarser_symbol_maps;
  };
  
  Symbol::ticket_type    Symbol::__mutex;

  static SymbolImpl::mutex_type            __non_terminal_mutex;
  static SymbolImpl::non_terminal_set_type __non_terminal_map;

  namespace symbol_impl
  {
#ifdef HAVE_TLS
    static __thread SymbolImpl*                   impl_tls = 0;
    static utils::thread_specific_ptr<SymbolImpl> impl;
#else
    static utils::thread_specific_ptr<SymbolImpl> impl;
#endif
    
    static SymbolImpl& instance()
    {
#ifdef HAVE_TLS
      if (! impl_tls) {
	impl.reset(new SymbolImpl());
	impl_tls = impl.get();
      }
      
      return *impl_tls;
#else
      if (! impl.get())
	impl.reset(new SymbolImpl());
      return *impl;
#endif
    }
  };



  Symbol::symbol_map_type& Symbol::__symbol_maps()
  {
    return symbol_impl::instance().symbol_maps;
  }

  Symbol::piece_type Symbol::sgml_tag() const
  {
    if (! is_sgml_tag())
      return Symbol::piece_type();

    const symbol_type& __str = static_cast<const symbol_type&>(*this);
    
    if (is_empty_tag())
      return piece_type(__str.begin() + 1, __str.end() - 2);
    else if (is_end_tag())
      return piece_type(__str.begin() + 2, __str.end() - 1);
    else
      return piece_type(__str.begin() + 1, __str.end() - 1);
  }

  bool Symbol::is_sgml_tag() const
  {
    const symbol_type& __str = static_cast<const symbol_type&>(*this);
    const size_type str_size = __str.size();
    
    return str_size >= 3 && __str[0] == '<' && __str[str_size - 1] == '>';
  }
  
  bool Symbol::is_start_tag() const
  {
    const symbol_type& __str = static_cast<const symbol_type&>(*this);
    const size_type str_size = __str.size();
    
    return (str_size >= 3
	    && __str[0] == '<'
	    && __str[1] != '/'
	    && __str[str_size - 2] != '/'
	    && __str[str_size - 1] == '>');
  }

  bool Symbol::is_end_tag() const
  {
    const symbol_type& __str = static_cast<const symbol_type&>(*this);
    const size_type str_size = __str.size();
    
    return (str_size >= 4
	    && __str[0] == '<'
	    && __str[1] == '/'
	    && __str[str_size - 2] != '/'
	    && __str[str_size - 1] == '>');
  }

  bool Symbol::is_empty_tag() const
  {
    const symbol_type& __str = static_cast<const symbol_type&>(*this);
    const size_type str_size = __str.size();
    
    return (str_size >= 4
	    && __str[0] == '<'
	    && __str[1] != '/'
	    && __str[str_size - 2] == '/'
	    && __str[str_size - 1] == '>');
  }

  bool Symbol::is_non_terminal() const
  {
    SymbolImpl::non_terminal_map_type& maps =  symbol_impl::instance().non_terminal_maps;
    
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
  
  Symbol::id_type Symbol::non_terminal_id() const
  {
    if (! is_non_terminal()) return id_type(-1);
    
    SymbolImpl::non_terminal_id_map_type& maps = symbol_impl::instance().non_terminal_id_maps;

    if (__id >= maps.size()) {
      const size_type size = __id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2, id_type(-1));
    }
    
    if (maps[__id] == id_type(-1)) {
      mutex_type::scoped_lock lock(__non_terminal_mutex);
      
      SymbolImpl::non_terminal_set_type::iterator iter = __non_terminal_map.insert(__id).first;
      
      maps[__id] = iter - __non_terminal_map.begin();
    }
    return maps[__id];
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
    
    SymbolImpl::index_map_type& maps = symbol_impl::instance().index_maps;

    const id_type __non_terminal_id = non_terminal_id();
    
    if (__non_terminal_id >= maps.size()) {
      const size_type size = __non_terminal_id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2, -1);
    }

    if (maps[__non_terminal_id] < 0) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      const symbol_type& word = symbol();
      
      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();
      
      if (! qi::parse(iter, iter_end,
		      '[' >> +(standard::char_ - ',' - ']') >> ',' >> qi::int_ [ phoenix::ref(maps[__non_terminal_id]) = qi::_1 ] >> ']'))
	maps[__non_terminal_id] = 0;
    }
      
    return maps[__non_terminal_id];
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
    
    SymbolImpl::non_terminal_symbol_map_type& maps = symbol_impl::instance().non_terminal_symbol_maps;

    const id_type __non_terminal_id = non_terminal_id();
    
    if (__non_terminal_id >= maps.size()) {
      const size_type size = __non_terminal_id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2, id_type(-1));
    }
      
    if (maps[__non_terminal_id] == id_type(-1)) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      const symbol_type& word = symbol();

      symbol_type::const_iterator iter = word.begin();
      symbol_type::const_iterator iter_end = word.end();

      maps[__non_terminal_id] = __id;
      
      std::string label;
      int index;
      if (qi::parse(iter, iter_end,
		    qi::lit('[')
		    >> (+(standard::char_ - ',' - ']')) >> -(',' >> qi::int_)
		    >> qi::lit(']'), label, index))
	maps[__non_terminal_id] = Symbol('[' + label + ']').id();
    }
      
    return maps[__non_terminal_id];
  }

  Symbol Symbol::pos() const
  {
    if (! is_terminal()) return Symbol();
    
    SymbolImpl::pos_symbol_map_type& maps = symbol_impl::instance().pos_symbol_maps;
    
    if (__id >= maps.size()) {
      const size_type size = __id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2, id_type(-1));
    }
    
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
    
    SymbolImpl::terminal_symbol_map_type& maps = symbol_impl::instance().terminal_symbol_maps;
    
    if (__id >= maps.size()) {
      const size_type size = __id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2, id_type(-1));
    }

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


  bool Symbol::annotated() const
  {
    if (! is_non_terminal()) return false;
    
    return non_terminal_strip().find('@') != piece_type::npos();
  }

  Symbol Symbol::annotate(const int pos, const bool bit) const
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
    
    SymbolImpl::coarser_symbol_map_type& maps = symbol_impl::instance().coarser_symbol_maps;

    const id_type __non_terminal_id = non_terminal_id();

    if (__non_terminal_id >= maps.size()) {
      const size_type size = __non_terminal_id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2);
    }

    if (pos >= static_cast<int>(maps[__non_terminal_id].size()))
      maps[__non_terminal_id].resize(pos + 1, id_type(-1));

    if (maps[__non_terminal_id][pos] == id_type(-1)) {
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
      
      maps[__non_terminal_id][pos] = Symbol(generated).non_terminal(__non_terminal_index).id();
    }
    return maps[__non_terminal_id][pos];
  }

  Symbol Symbol::coarse() const
  {
    if (! is_non_terminal()) return *this;
    
    SymbolImpl::coarse_symbol_map_type& maps = symbol_impl::instance().coarse_symbol_maps;

    const id_type __non_terminal_id = non_terminal_id();
    
    if (__non_terminal_id >= maps.size()) {
      const size_type size = __non_terminal_id + 1;
      const size_type power2 = utils::bithack::branch(utils::bithack::is_power2(size),
						      size,
						      size_type(utils::bithack::next_largest_power2(size)));
      maps.reserve(power2);
      maps.resize(power2, id_type(-1));
    }
    
    if (maps[__non_terminal_id] == id_type(-1)) {
      namespace xpressive = boost::xpressive;
      
      typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
      typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
      
      static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (-+xpressive::_d);
      
      const int __non_terminal_index = non_terminal_index();
      const piece_type piece = non_terminal_strip();
      
      pmatch what;
      if (xpressive::regex_match(piece, what, re))
	maps[__non_terminal_id] = Symbol('[' + what[1] + ']').non_terminal(__non_terminal_index).id();
      else
	maps[__non_terminal_id] = __id;
    }
    return maps[__non_terminal_id];
  }
  
  void Symbol::write(const path_type& path)
  {
    ticket_type::scoped_reader_lock lock(__mutex);
    
    Vocab vocab(path, std::max(size_type(__symbols().size() / 2), size_type(1024)));
    symbol_set_type::const_iterator siter_end = __symbols().end();
    for (symbol_set_type::const_iterator siter = __symbols().begin(); siter != siter_end; ++ siter)
      vocab.insert(*siter);
  }
  
};
