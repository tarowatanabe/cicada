//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/thread.hpp>

#include "attribute_vector.hpp"

#include "utils/thread_specific_ptr.hpp"
#include "utils/json_string_parser.hpp"
#include "utils/json_string_generator.hpp"

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
      
      data %= data_string | data_double | data_int;
      
      attribute %= key >> ':' >> data;
      attributes %= '{' >> -(attribute % ',') >> '}';
    }
    
    typedef boost::spirit::standard::space_type space_type;
    
    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1>                      data_int;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > data_double;
    utils::json_string_parser<Iterator>                                                      data_string;
    
    utils::json_string_parser<Iterator>                                         key;
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
      
      data %= data_int | data_double | data_string;
      
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
    
    boost::spirit::karma::real_generator<double, real_precision>              data_double;
    boost::spirit::karma::int_generator<AttributeVector::int_type, 10, false> data_int;
    utils::json_string_generator<Iterator, true>                              data_string;

    utils::json_string_generator<Iterator, true>                        key;
    boost::spirit::karma::rule<Iterator, AttributeVector::data_type()>  data;
    boost::spirit::karma::rule<Iterator, AttributeVector::value_type()> attribute;
    boost::spirit::karma::rule<Iterator, AttributeVector::attribute_vector_type()> attributes;
  };
  
  namespace attribute_vector_parser_impl
  {
    typedef attribute_vector_parser<utils::piece::const_iterator > grammar_type;

#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static utils::thread_specific_ptr<grammar_type > __grammar;
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


  bool AttributeVector::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    clear();
    
    // empty attribute vector...
    if (iter == end) return true;
    
    const char* citer_begin = &(*iter);
    const char* citer       = &(*iter);
    const char* citer_end   = &(*end);
    const bool result = assign(citer, citer_end);
    
    iter += citer - citer_begin;
    
    return result;
  }

  bool AttributeVector::assign(utils::piece::const_iterator& iter, utils::piece::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    // empty attribute vector...
    if (iter == end) return true;
    
    attribute_set_parsed_type parsed;
    
    if (qi::phrase_parse(iter, end, attribute_vector_parser_impl::instance(), standard::space, parsed)) {
      __values.insert(parsed.begin(), parsed.end());
      return true;
    } else
      return false;
  }

  void AttributeVector::assign(const utils::piece& x)
  {
    utils::piece::const_iterator iter(x.begin());
    utils::piece::const_iterator end(x.end());

    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("attribute-vector format parsing failed...");
  }

  namespace attribute_vector_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef attribute_vector_generator<iterator_type> grammar_type;    

#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static utils::thread_specific_ptr<grammar_type > __grammar;
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

  
  std::ostream& operator<<(std::ostream& os, const AttributeVector& x)
  {
    namespace karma = boost::spirit::karma;
    
    attribute_vector_generator_impl::iterator_type iter(os);
    
    if (! karma::generate(iter, attribute_vector_generator_impl::instance(), x.__values))
      throw std::runtime_error("failed attribute vector generation!");

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

  typedef AttributeVector::data_type attribute_data_type;

  template <typename Iterator>
  struct attribute_data_parser : boost::spirit::qi::grammar<Iterator, attribute_data_type(), boost::spirit::standard::space_type>
  {
    attribute_data_parser() : attribute_data_parser::base_type(data)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      data %= data_string | data_double | data_int;
    }
    
    typedef boost::spirit::standard::space_type space_type;

    utils::json_string_parser<Iterator>                                                      data_string;
    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1>                      data_int;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > data_double;
    
    boost::spirit::qi::rule<Iterator, AttributeVector::data_type(), space_type> data;
  };

  template <typename Iterator>
  struct attribute_data_generator : boost::spirit::karma::grammar<Iterator, attribute_data_type()>
  {
    attribute_data_generator() : attribute_data_generator::base_type(data)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      data %= data_int | data_double | data_string;
    }
    
    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 10;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision>              data_double;
    boost::spirit::karma::int_generator<AttributeVector::int_type, 10, false> data_int;
    utils::json_string_generator<Iterator, true>                              data_string;
    
    boost::spirit::karma::rule<Iterator, AttributeVector::data_type()> data;
  };

  namespace attribute_data_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef attribute_data_generator<iterator_type> grammar_type;    
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static utils::thread_specific_ptr<grammar_type > __grammar;
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

  std::ostream& operator<<(std::ostream& os, const AttributeVector::data_type& x)
  {
    namespace karma = boost::spirit::karma;
    
    attribute_data_generator_impl::iterator_type iter(os);
    
    if (! karma::generate(iter, attribute_data_generator_impl::instance(), x))
      throw std::runtime_error("failed attribute vector data generation!");
    
    return os;
  }

  namespace attribute_data_parser_impl
  {
    typedef attribute_data_parser<std::string::const_iterator > grammar_type;
    
#ifdef HAVE_TLS
    static __thread grammar_type* __grammar_tls = 0;
    static utils::thread_specific_ptr<grammar_type > __grammar;
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

  std::istream& operator>>(std::istream& is, AttributeVector::data_type& x)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    x = attribute_data_type();
    
    std::string token;
    if (is >> token) {
      std::string::const_iterator iter = token.begin();
      std::string::const_iterator end  = token.end();
      
      if (qi::phrase_parse(iter, end, attribute_data_parser_impl::instance(), standard::space, x))
	if (iter != end)
	  throw std::runtime_error("failed parsing attribute vector data");
    }

    return is;
  }

  
};
