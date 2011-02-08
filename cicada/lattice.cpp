//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE 1
#define PHOENIX_THREADSAFE 1

#include "lattice.hpp"
#include "vocab.hpp"

#include "unite.hpp"

#include <boost/numeric/conversion/bounds.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/thread_specific_ptr.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Lattice::arc_type,
			  (cicada::Lattice::symbol_type, label)
			  (cicada::Lattice::feature_set_type, features)
			  (int, distance)
			  )

namespace cicada
{

  void Lattice::unite(const Lattice& x)
  {
    cicada::unite(*this, x);
  }
  
  void Lattice::initialize_distance()
  {
    const int pos_infinity = boost::numeric::bounds<int>::highest();
    const int neg_infinity = boost::numeric::bounds<int>::lowest();
    
    const size_type dist_size = size() + 1;
    
    // edge-cost is inifinity if no-path
    dist_short.clear();
    dist_short.reserve(dist_size, dist_size);
    dist_short.resize(dist_size, dist_size, pos_infinity);

    dist_long.clear();
    dist_long.reserve(dist_size, dist_size);
    dist_long.resize(dist_size, dist_size, neg_infinity);
    
    // edge-cost dist(i, i) = 0
    for (size_t i = 0; i != dist_size; ++ i) {
      dist_short(i, i) = 0;
      dist_long(i, i)  = 0;
    }

    // we will assume lattice structure to compute shortest and longest distance...
    typedef std::set<int, std::less<int>, std::allocator<int> > closure_type;
    typedef std::vector<closure_type, std::allocator<closure_type> > closure_set_type;
    
    closure_set_type closure(lattice.size() + 1);
    
    for (int node = lattice.size() - 1; node >= 0; -- node) {
      arc_set_type::const_iterator aiter_end = lattice[node].end();
      for (arc_set_type::const_iterator aiter = lattice[node].begin(); aiter != aiter_end; ++ aiter) {
	const int last = node + aiter->distance;
	const int dist = aiter->label != Vocab::EPSILON;
	
	dist_short(node, last) = utils::bithack::min(dist_short(node, last), dist);
	dist_long(node, last)  = utils::bithack::max(dist_long(node, last),  dist);
      }
      
      for (arc_set_type::const_iterator aiter = lattice[node].begin(); aiter != aiter_end; ++ aiter) {
	const int last = node + aiter->distance;
	
	closure[node].insert(last);
	closure_type::const_iterator citer_end = closure[last].end();
	for (closure_type::const_iterator citer = closure[last].begin(); citer != citer_end; ++ citer) {
	  
	  dist_short(node, *citer) = utils::bithack::min(dist_short(node, *citer), dist_short(node, last) + dist_short(last, *citer));
	  dist_long(node, *citer)  = utils::bithack::max(dist_long(node, *citer),  dist_long(node, last)  + dist_long(last, *citer));
	  
	  closure[node].insert(*citer);
	}
      }
    }
  }
  
  template <typename Iterator>
  struct lattice_grammar_parser : boost::spirit::qi::grammar<Iterator, Lattice::lattice_type(), boost::spirit::standard::space_type>
  {
    lattice_grammar_parser() : lattice_grammar_parser::base_type(lattice_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
    
      jlf_escape_char.add
	("\\\"", '\"')
	("\\\\", '\\')
	("\\/", '/')
	("\\b", '\b')
	("\\f", '\f')
	("\\n", '\n')
	("\\r", '\r')
	("\\t", '\t');

      plf_escape_char.add
	("\\\"", '\"')
	("\\\'", '\'')
	("\\\\", '\\')
	("\\a", '\a')
	("\\b", '\b')
	("\\f", '\f')
	("\\n", '\n')
	("\\r", '\r')
	("\\t", '\t')
	("\\v", '\v');
      
      jlf_label_double_quote %= '"' >> qi::lexeme[*(jlf_escape_char | ~standard::char_('"'))] >> '"';
      plf_label_double_quote %= '"' >> qi::lexeme[*(plf_escape_char | ~standard::char_('"'))] >> '"';
      plf_label_single_quote %= '\'' >> qi::lexeme[*(plf_escape_char | ~standard::char_('\''))] >> '\'';
      
      jlf_lattice_score %= jlf_label_double_quote >> ':' >> qi::double_;
      plf_lattice_score %= qi::attr("lattice-cost") >> qi::double_;
      
      jlf_lattice_arc %= '[' >> jlf_label_double_quote >> ',' >> '{' >> -(jlf_lattice_score % ',') >> '}' >> ',' >> qi::int_ >> ']';
      plf_lattice_arc %= '(' >> (plf_label_double_quote | plf_label_single_quote) >> ',' >> qi::repeat(1)[plf_lattice_score] >> ',' >> qi::int_ >> ')';
      
      jlf_lattice_set %= '[' >> -(jlf_lattice_arc % ',') >> ']';
      plf_lattice_set %= '(' >> +(plf_lattice_arc >> ',') >> ')';

      lattice_grammar %= qi::hold['(' >> *(plf_lattice_set >> ',') >> ')'] | ('[' >> -(jlf_lattice_set % ',') >> ']');
    }
    
    typedef boost::spirit::standard::space_type space_type;

    boost::spirit::qi::symbols<char, char> jlf_escape_char;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type> jlf_label_double_quote;
    
    boost::spirit::qi::rule<Iterator, std::pair<std::string, double >(), space_type> jlf_lattice_score;
    boost::spirit::qi::rule<Iterator, Lattice::arc_type(), space_type>               jlf_lattice_arc;
    boost::spirit::qi::rule<Iterator, Lattice::arc_set_type(), space_type>           jlf_lattice_set;
    
    boost::spirit::qi::symbols<char, char> plf_escape_char;
    
    boost::spirit::qi::rule<Iterator, std::string(), space_type> plf_label_double_quote;
    boost::spirit::qi::rule<Iterator, std::string(), space_type> plf_label_single_quote;
    
    boost::spirit::qi::rule<Iterator, std::pair<std::string, double >(), space_type> plf_lattice_score;
    boost::spirit::qi::rule<Iterator, Lattice::arc_type(), space_type>               plf_lattice_arc;
    boost::spirit::qi::rule<Iterator, Lattice::arc_set_type(), space_type>           plf_lattice_set;
    
    boost::spirit::qi::rule<Iterator, Lattice::lattice_type(), space_type>           lattice_grammar;
  };

  void Lattice::assign(const utils::piece& x)
  {
    std::string::const_iterator iter(x.begin());
    std::string::const_iterator end(x.end());

    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("LATTICE format parsing failed...");
  }

  namespace lattice_grammar_parser_impl
  {
    typedef lattice_grammar_parser<std::string::const_iterator > grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif

    static grammar_type& instance()
    {
#ifdef HAVE_TLS
      if (! __grammar_tls) {
	__grammar.reset(new grammar_type());
	__grammar_tls = __grammar.get();
      }
      
      return *__grammar_tls;
#else
      if (! __grammar.get())
	__grammar.reset(new grammar_type());
      
      return *__grammar;
#endif
    }
  };
  
  bool Lattice::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    // empty lattice...
    if (iter == end) return true;
    
    std::string::const_iterator iter_back = iter;
    
    if (qi::phrase_parse(iter, end, lattice_grammar_parser_impl::instance(), standard::space, lattice)) {
      initialize_distance();
      return true;
    } else {
      clear();
      
      // fallback to sentence...
      iter = iter_back;
      Sentence sentence;
      
      if (sentence.assign(iter, end)) {
	Sentence::const_iterator iter_end = sentence.end();
	for (Sentence::const_iterator iter = sentence.begin(); iter != iter_end; ++ iter)
	  lattice.push_back(arc_set_type(1, arc_type(*iter)));
	
	return true;
      } else {
	clear();
	return false;
      }
    }
  }
  
  std::istream& operator>>(std::istream& is, Lattice& x)
  {
    std::string line;
    
    x.clear();
    if (std::getline(is, line) && ! line.empty())
      x.assign(line);
    
    return is;
  }

  // generator...
  template <typename Iterator>
  struct lattice_grammar_generator : boost::spirit::karma::grammar<Iterator, Lattice::lattice_type()>
  {
    lattice_grammar_generator() : lattice_grammar_generator::base_type(lattice_grammar)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      // json grammar...
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('/', "\\/")
	('\b', "\\b")
	('\f', "\\f")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t");
      
      label_double_quote %= '\"' << *(escape_char | ~standard::char_('\"')) << '\"';

      lattice_score %= label_double_quote << ": " << karma::double_;
      
      lattice_arc %= '[' << label_double_quote << ", " << '{' << -(lattice_score % ',') << '}' << ", " << karma::int_ << ']';
      lattice_set %= '[' << -(lattice_arc % ", ") << ']';
      lattice_grammar %= '[' << -(lattice_set % ", ") << ']';
    }
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    
    boost::spirit::karma::rule<Iterator, std::string()> label_double_quote;
    boost::spirit::karma::rule<Iterator, std::pair<std::string, double >()> lattice_score;
    boost::spirit::karma::rule<Iterator, Lattice::arc_type()> lattice_arc;
    boost::spirit::karma::rule<Iterator, Lattice::arc_set_type()> lattice_set;
    boost::spirit::karma::rule<Iterator, Lattice::lattice_type()> lattice_grammar;
  };

  namespace lattice_grammar_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef lattice_grammar_generator<iterator_type> grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
#else
    static utils::thread_specific_ptr<grammar_type > __grammar;
#endif
    
    static grammar_type& instance()
    {
#ifdef HAVE_TLS
      if (! __grammar_tls) {
	__grammar.reset(new grammar_type());
	__grammar_tls = __grammar.get();
      }
      
      return *__grammar_tls;
#else
      if (! __grammar.get())
	__grammar.reset(new grammar_type());
      
      return *__grammar;
#endif
    }
  };


  std::ostream& operator<<(std::ostream& os, const Lattice& x)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    lattice_grammar_generator_impl::iterator_type iter(os);
    
    if (! karma::generate(iter, lattice_grammar_generator_impl::instance(), x.lattice))
      throw std::runtime_error("failed lattice generation!");

    return os;
  }
  
};
