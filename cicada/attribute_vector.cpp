//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/thread.hpp>

#include "attribute_vector.hpp"

namespace cicada
{
  
  // format????
  // attribute value
  //
  // {"key": value, "key": value, ...}

  typedef std::pair<std::string, AttributeVector::data_type> attribute_parsed_type;
  typedef std::vector<attribute_parsed_type, std::allocator<attribute_parsed_type> > attribute_set_parsed_type;

  template <typename Iterator>
  struct attribute_vector_parser : boost::spirit::qi::grammar<Iterator, attribute_set_parsed_type(), boost::spirit::standard::space_type>
  {
    attribute_vector_parser() : attribute_vector_parser::base_type(attributes)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::hold;
      using qi::lexeme;
      using standard::char_;
      using standard::space;
      using qi::double_;

      escape_char.add
	("\\\"", '\"')
	("\\\\", '\\')
	("\\/", '/')
	("\\b", '\b')
	("\\f", '\f')
	("\\n", '\n')
	("\\r", '\r')
	("\\t", '\t')
	("\\u0020", ' ');
      
      key %= ('\"' >> lexeme[*(escape_char | (char_ - '\"' - space))] >> '\"');
      data_value %= ('\"' >> lexeme[*(escape_char | (char_ - '\"'))] >> '\"');
      data %= data_value | double_dot | int64_;
      
      attribute %= key >> ':' >> data;
      attributes %= '{' >> -(attribute % ',') >> '}';
    }
    
    typedef boost::spirit::standard::space_type space_type;
    
    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1> int64_;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > double_dot;
    
    boost::spirit::qi::symbols<char, char> escape_char;
    boost::spirit::qi::rule<Iterator, std::string(), space_type>                key;
    boost::spirit::qi::rule<Iterator, std::string(), space_type>                data_value;
    boost::spirit::qi::rule<Iterator, AttributeVector::data_type(), space_type> data;
    boost::spirit::qi::rule<Iterator, attribute_parsed_type(), space_type>      attribute;
    boost::spirit::qi::rule<Iterator, attribute_set_parsed_type(), space_type>  attributes;
  };

  template <typename Iterator>
  struct attribute_vector_generator : boost::spirit::karma::grammar<Iterator, AttributeVector::attribute_vector_type()>
  {
    attribute_vector_generator() : attribute_vector_generator::base_type(attributes)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::buffer;
      using standard::char_;
      using karma::double_;
      
      escape_char.add
	('\\', "\\\\")
	('\"', "\\\"")
	('/', "\\/")
	('\b', "\\b")
	('\f', "\\f")
	('\n', "\\n")
	('\r', "\\r")
	('\t', "\\t")
	(' ', "\\u0020");
      
      key %= ('\"' << +(escape_char | ~char_('\"')) << '\"');
      data %= int64_ | double10 | key;
      
      attribute %= key << ':' << data;
      attributes %= '{' << -(attribute % ',') << '}';
    }
    
    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 10;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double10;
    
    boost::spirit::karma::int_generator<AttributeVector::int_type, 10, false> int64_;
    
    boost::spirit::karma::symbols<char, const char*> escape_char;
    boost::spirit::karma::rule<Iterator, std::string()> key;
    boost::spirit::karma::rule<Iterator, AttributeVector::data_type()> data;
    boost::spirit::karma::rule<Iterator, AttributeVector::value_type()> attribute;
    boost::spirit::karma::rule<Iterator, AttributeVector::attribute_vector_type()> attributes;
  };

  bool AttributeVector::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    typedef attribute_vector_parser<std::string::const_iterator > grammar_type;
    
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

    // empty attribute vector...
    if (iter == end) return true;

    attribute_set_parsed_type parsed;
    
    if (boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, parsed)) {
      __values.insert(parsed.begin(), parsed.end());
      return true;
    } else
      return false;
  }

  void AttributeVector::assign(const std::string& x)
  {
    std::string::const_iterator iter = x.begin();
    std::string::const_iterator end  = x.end();

    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("attribute-vector format parsing failed...");
  }

  
  std::ostream& operator<<(std::ostream& os, const AttributeVector& x)
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef attribute_vector_generator<iterator_type> grammar_type;    
 
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
    
    if (! boost::spirit::karma::generate(iter, grammar, x.__values))
      throw std::runtime_error("failed lattice generation!");

    return os;
  }

  std::istream& operator>>(std::istream& is, AttributeVector& x)
  {
    x.clear();
    
    std::string token;
    if (is >> token)
      x.assign(token);
    
    return is;
  }
};
