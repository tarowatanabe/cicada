#define BOOST_SPIRIT_THREADSAFE 1
#define PHOENIX_THREADSAFE 1

#include "lattice.hpp"
#include "vocab.hpp"

#include <boost/numeric/conversion/bounds.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/thread.hpp>

#include "utils/config.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Lattice::arc_type,
			  (cicada::Lattice::symbol_type, label)
			  (cicada::Lattice::feature_set_type, features)
			  (int, distance)
			  )

namespace cicada
{
  
  void Lattice::initialize_distance()
  {
    const size_type dist_size = size() + 1;

    const int pos_infinity = boost::numeric::bounds<int>::highest();
    const int neg_infinity = boost::numeric::bounds<int>::lowest();
    
    // edge-cost is inifinity if no-path
    dist_short.clear();
    dist_short.reserve(dist_size, dist_size);
    dist_short.resize(dist_size, dist_size, pos_infinity);

    dist_long.clear();
    dist_long.reserve(dist_size, dist_size);
    dist_long.resize(dist_size, dist_size, neg_infinity);

#if 0
    // we will assume lattice structure to compute shortest and longest distance...
    typedef std::set<int, std::less<int>, std::allocator<int> > closure_type;
    typedef std::vector<closure_type, std::allocator<closure_type> > closure_set_type;
    
    closure_set_type closure(lattice.size() + 1);
    
    for (int node = lattice.size() - 1; node >= 0; -- node) {
      arc_set_type::const_iterator aiter_end = lattice[node].end();
      for (arc_set_type::const_iterator aiter = lattice[node].begin(); aiter != aiter_end; ++ aiter) {
	const int last = node + aiter->distance;
	
	if (aiter->label != vocab_type::EPSILON) {
	  closure[node].insert(last);
	  
	  closure_type::const_iterator citer_end = closure[last].end();
	  for (closure_type::const_iterator citer = closure[last].begin(); citer != citer_end; ++ citer)
	    closure[node].insert(*citer);
	} else {
	  
	  
	  
	}
      }
    }
#endif
    
#if 1
    // edge-cost for dist(i, j)
    for (size_t i = 0; i != lattice.size(); ++ i)
      for (size_t j = 0; j != lattice[i].size(); ++ j) {
	dist_short(i, i + lattice[i][j].distance) = (lattice[i][j].label != Vocab::EPSILON);
	dist_long(i, i + lattice[i][j].distance)  = (lattice[i][j].label != Vocab::EPSILON);
      }
    
    // edge-cost dist(i, i) = 0
    for (size_t i = 0; i != dist_size; ++ i) {
      dist_short(i, i) = 0;
      dist_long(i, i)  = 0;
    }
    
    // Floyd-Warshall algorithm to compute shortest path
    for (size_t k = 0; k != dist_size; ++ k) 
      for (size_t i = 0; i != dist_size; ++ i)
	for (size_t j = 0; j != dist_size; ++ j) {
	  if (dist_short(i, k) != pos_infinity && dist_short(k, j) != pos_infinity)
	    dist_short(i, j) = utils::bithack::min(dist_short(i, j), dist_short(i, k) + dist_short(k, j));
	    
	  if (dist_long(i, k) != neg_infinity && dist_long(k, j) != neg_infinity)
	    dist_long(i, j) = utils::bithack::max(dist_long(i, j), dist_long(i, k) + dist_long(k, j));
	}
#endif
  }
  
  template <typename Iterator>
  struct lattice_grammar_parser : boost::spirit::qi::grammar<Iterator, Lattice::lattice_type(), boost::spirit::standard::space_type>
  {
    lattice_grammar_parser() : lattice_grammar_parser::base_type(lattice_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
    
      using qi::phrase_parse;
      using qi::lexeme;
      using qi::repeat;
      using qi::hold;
      using qi::lit;
      using qi::attr;
      using standard::char_;
      using qi::double_;
      using qi::int_;
      
      using namespace qi::labels;
      
      using standard::space;

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
      
      jlf_label_double_quote %= '"' >> lexeme[*(jlf_escape_char | ~char_('"'))] >> '"';
      plf_label_double_quote %= '"' >> lexeme[*(plf_escape_char | ~char_('"'))] >> '"';
      plf_label_single_quote %= '\'' >> lexeme[*(plf_escape_char | ~char_('\''))] >> '\'';
      
      jlf_lattice_score %= jlf_label_double_quote >> ':' >> double_;
      plf_lattice_score %= attr("lattice-cost") >> double_;
      
      jlf_lattice_arc %= '[' >> jlf_label_double_quote >> ',' >> '{' >> -(jlf_lattice_score % ',') >> '}' >> ',' >> int_ >> ']';
      plf_lattice_arc %= '(' >> (plf_label_double_quote | plf_label_single_quote) >> ',' >> repeat(1)[plf_lattice_score] >> ',' >> int_ >> ')';
      
      jlf_lattice_set %= '[' >> -(jlf_lattice_arc % ',') >> ']';
      plf_lattice_set %= '(' >> +(plf_lattice_arc >> ',') >> ')';

      lattice_grammar %= hold[lit('(') >> *(plf_lattice_set >> ',') >> lit(')')] | (lit('[') >> -(jlf_lattice_set % ',') >> lit(']'));
    }
    

    boost::spirit::qi::symbols<char, char> jlf_escape_char;
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> jlf_label_double_quote;
    
    boost::spirit::qi::rule<Iterator, std::pair<std::string, double >(), boost::spirit::standard::space_type> jlf_lattice_score;
    boost::spirit::qi::rule<Iterator, Lattice::arc_type(), boost::spirit::standard::space_type>               jlf_lattice_arc;
    boost::spirit::qi::rule<Iterator, Lattice::arc_set_type(), boost::spirit::standard::space_type>           jlf_lattice_set;
    
    boost::spirit::qi::symbols<char, char> plf_escape_char;
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> plf_label_double_quote;
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> plf_label_single_quote;
    
    boost::spirit::qi::rule<Iterator, std::pair<std::string, double >(), boost::spirit::standard::space_type> plf_lattice_score;
    boost::spirit::qi::rule<Iterator, Lattice::arc_type(), boost::spirit::standard::space_type>               plf_lattice_arc;
    boost::spirit::qi::rule<Iterator, Lattice::arc_set_type(), boost::spirit::standard::space_type>           plf_lattice_set;
    
    boost::spirit::qi::rule<Iterator, Lattice::lattice_type(), boost::spirit::standard::space_type>           lattice_grammar;
  };

  void Lattice::assign(const std::string& x)
  {
    std::string::const_iterator iter = x.begin();
    std::string::const_iterator end = x.end();

    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("LATTICE format parsing failed...");
  }
  
  bool Lattice::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    typedef lattice_grammar_parser<std::string::const_iterator > grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
    
    if (! __grammar_tls) {
      __grammar.reset(new grammar_type());
      __grammar_tls = __grammar.get();
    }
    
    grammar_type& grammar = *__grammar_tls;
#else
    static boost::thread_specific_ptr<grammar_type > __grammar;
    if (! __grammar.get())
      __grammar.reset(new grammar_type());
    
    grammar_type& grammar = *__grammar;
#endif
    
    clear();
    
    // empty lattice...
    if (iter == end) return true;
    
    std::string::const_iterator iter_back = iter;
    
    if (boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, lattice)) {
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
      namespace phoenix = boost::phoenix;
    
      using karma::lit;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      
      using namespace karma::labels;
      
      using standard::space;

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
      
      label_double_quote %= '\"' << *(escape_char | ~char_('\"')) << '\"';

      lattice_score %= label_double_quote << ": " << double_;
      
      lattice_arc %= '[' << label_double_quote << ", " << '{' << -(lattice_score % ',') << '}' << ", " << int_ << ']';
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


  std::ostream& operator<<(std::ostream& os, const Lattice& x)
  {
    typedef std::ostream_iterator<char> iterator_type;

    typedef lattice_grammar_generator<iterator_type> grammar_type;
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static boost::thread_specific_ptr<grammar_type > __grammar;
    
    if (! __grammar_tls) {
      __grammar.reset(new grammar_type());
      __grammar_tls = __grammar.get();
    }
    
    grammar_type& grammar = *__grammar_tls;
#else
    static boost::thread_specific_ptr<grammar_type > __grammar;
    if (! __grammar.get())
      __grammar.reset(new grammar_type());
    
    grammar_type& grammar = *__grammar;
#endif
    
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, x.lattice))
      throw std::runtime_error("failed lattice generation!");

    return os;
  }
  
};
