
#include "lattice.hpp"

#define BOOST_SPIRIT_THREADSAFE

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

    const difference_type infinity = boost::numeric::bounds<difference_type>::highest();
    
    // edge-cost is inifinity if no-path
    dist.clear();
    dist.reserve(dist_size, dist_size);
    dist.resize(dist_size, dist_size, infinity);
    
    // edge-cost for dist(i, j)
    for (int i = 0; i < lattice.size(); ++ i)
      for (int j = 0; j < lattice[i].size(); ++ j)
	dist(i, i + lattice[i][j].distance) = 1;
    
    // edge-cost dist(i, i) = 0
    for (int i = 0; i < dist_size; ++ i)
      dist(i, i) = 0;
    
    // Floyd-Warshall algorithm to compute shortest path
    for (int k = 0; k < dist_size; ++ k) 
      for (int i = 0; i < dist_size; ++ i)
	for (int j = 0; j < dist_size; ++ j)
	  if (dist(i, k) != infinity && dist(k, j) != infinity)
	    dist(i, j) = utils::bithack::min(dist(i, j), dist(i, k) + dist(k, j));
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
      using qi::omit;
      using qi::repeat;
      using qi::lit;
      using qi::inf;
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
      
      jlf_lattice_set %= '[' >> (plf_lattice_arc % ',') >> ']';
      plf_lattice_set %= '(' >> +(plf_lattice_arc >> ',') >> ')';
      
      lattice_grammar %= (lit('(') >> *(plf_lattice_set >> ',') >> lit(')') | '[' >> -(jlf_lattice_set % ',') >> ']'); 
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
    
    hypergraph_parser<grammar_type>& grammar = *__grammar;
#endif

    clear();
    if (x.empty()) return;
    
    std::string::const_iterator iter = x.begin();
    std::string::const_iterator end = x.end();
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, lattice);
    if (! result || iter != end)
      throw std::runtime_error("LATTICE format parsing failed...");
    
    initialize_distance();
  }
  
  std::istream& operator>>(std::istream& is, Lattice& x)
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
    
    hypergraph_parser<grammar_type>& grammar = *__grammar;
#endif

    std::string line;
    
    x.clear();
    if (std::getline(is, line) && ! line.empty()) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, x.lattice);
      if (result && iter == end)
	x.initialize_distance();
      else
	x.clear();
    }
    
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
    
      using karma::omit;
      using karma::repeat;
      using karma::lit;
      using karma::inf;
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
      lattice_set %= '[' << (lattice_arc % ", ") << ']';
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
    
    hypergraph_parser<grammar_type>& grammar = *__grammar;
#endif
    
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, x.lattice))
      throw std::runtime_error("failed generation!");

    return os;
  }
  
};
