#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/numeric/conversion/bounds.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/array.hpp>

#include <string>
#include <iostream>
#include <vector>

typedef std::string phrase_type;
typedef std::vector<double, std::allocator<double> > scores_type;
typedef boost::fusion::tuple<phrase_type, phrase_type, scores_type> phrase_pair_type;

template <typename Iterator>
struct phrase_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_type(), boost::spirit::standard::space_type>
{
  phrase_pair_parser() : phrase_pair_parser::base_type(phrase_pair)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    using qi::phrase_parse;
    using qi::lexeme;
    using qi::hold;
    using qi::repeat;
    using standard::char_;
    using qi::double_;
    using standard::space;
    
    phrase %= lexeme[+(char_ - (space >> "|||" >> space))];
    scores %= repeat(5)[double_];
    phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> scores;
  }
  
  boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
  boost::spirit::qi::rule<Iterator, scores_type(), boost::spirit::standard::space_type> scores;
  boost::spirit::qi::rule<Iterator, phrase_pair_type(), boost::spirit::standard::space_type> phrase_pair;
};

int main(int argc, char** argv)
{
  // we will split inputs into "lhs ||| rhs ||| scores"
  
  phrase_pair_parser<std::string::const_iterator> parser;
  
  std::string line;
  while (std::getline(std::cin, line)) {
    phrase_pair_type phrase_pair;
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();

    const bool result = boost::spirit::qi::phrase_parse(iter, end, parser, boost::spirit::standard::space, phrase_pair);
    if (! result || iter != end)
      std::cerr << "invalid input: " <<  line << " pos: " << (iter - line.begin()) << std::endl;

    std::cout << boost::fusion::get<0>(phrase_pair)
	      << " --- "
	      << boost::fusion::get<1>(phrase_pair)
	      << std::endl;
  }
  
  
  
}
