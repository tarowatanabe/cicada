//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include "ngram_count_set.hpp"

#include "utils/config.hpp"
#include "utils/thread_specific_ptr.hpp"
#include "utils/json_string_parser.hpp"
#include "utils/json_string_generator.hpp"
#include "utils/getline.hpp"

namespace cicada
{
  typedef NGramCountSet::ngram_type ngram_type;
  typedef NGramCountSet::count_type count_type;
  
  typedef utils::unordered_map<ngram_type, count_type, boost::hash<ngram_type>, std::equal_to<ngram_type>,
			       std::allocator<std::pair<const ngram_type, count_type> > >::type ngram_set_type;

  template <typename Iterator>
  struct ngram_count_set_parser : boost::spirit::qi::grammar<Iterator, ngram_set_type(), boost::spirit::standard::space_type>
  {
    typedef std::pair<ngram_type, count_type> value_type;
    typedef NGramCountSet::word_type word_type;
    typedef std::vector<word_type, std::allocator<word_type> > sequence_type;
    
    ngram_count_set_parser() : ngram_count_set_parser::base_type(ngrams)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      ngram %= '[' >> -(word % ',') >> ']';
      value = (ngram >> ':' >> qi::double_) [qi::_val = phoenix::construct<value_type>(phoenix::construct<ngram_type>(phoenix::begin(qi::_1),
														      phoenix::end(qi::_1)),
										       qi::_2)];
      ngrams %= '{' >> -(value % ',') >> '}';
    }
    
    typedef boost::spirit::standard::space_type space_type;
    
    utils::json_string_parser<Iterator> word;
    boost::spirit::qi::rule<Iterator, sequence_type(), space_type> ngram;
    boost::spirit::qi::rule<Iterator, value_type(), space_type> value;
    boost::spirit::qi::rule<Iterator, ngram_set_type(), space_type> ngrams;
  };
  
  template <typename Iterator>
  struct ngram_count_set_generator : boost::spirit::karma::grammar<Iterator, ngram_set_type()>
  {
    typedef ngram_set_type::value_type value_type;
    
    ngram_count_set_generator() : ngram_count_set_generator::base_type(ngrams)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      value  %= (karma::lit('[') << -(word % ',') << karma::lit(']')) << karma::lit(':') << count;
      ngrams %= karma::lit('{') << -(value % ',') << karma::lit('}');
    }
    
    struct real_precision : boost::spirit::karma::real_policies<count_type>
    {
      static unsigned int precision(count_type) 
      { 
        return 10;
      }
    };
    
    boost::spirit::karma::real_generator<count_type, real_precision> count;
    utils::json_string_generator<Iterator>                           word;
    boost::spirit::karma::rule<Iterator, value_type()>               value;
    boost::spirit::karma::rule<Iterator, ngram_set_type()>           ngrams;
  };
  
  namespace ngram_count_set_parser_impl
  {
    typedef utils::piece::const_iterator iterator_type;
    typedef ngram_count_set_parser<iterator_type> grammar_type;
    
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
  
  namespace ngram_count_set_generator_impl
  {
    typedef std::ostream_iterator<char> iterator_type;
    typedef ngram_count_set_generator<iterator_type> grammar_type;
    
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

  void NGramCountSet::assign(const utils::piece& x)
  {
    utils::piece::const_iterator iter(x.begin());
    utils::piece::const_iterator end(x.end());
    
    const bool result = assign(iter, end);
    if (! result || iter != end)
      throw std::runtime_error("NGramCountSet format parsing failed...");
  }

  bool NGramCountSet::assign(std::string::const_iterator& iter, std::string::const_iterator end)
  {
    const char* citer_begin = &(*iter);
    const char* citer       = &(*iter);
    const char* citer_end   = &(*end);
    
    const bool result = assign(citer, citer_end);
    
    iter += citer - citer_begin;
    
    return result;
  }

  bool NGramCountSet::assign(utils::piece::const_iterator& iter, utils::piece::const_iterator end)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    clear();
    
    if (iter == end) return true;
    
    return qi::phrase_parse(iter, end, ngram_count_set_parser_impl::instance(), standard::space, ngrams);
  }
  
  std::ostream& operator<<(std::ostream& os, const NGramCountSet& x)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    ngram_count_set_generator_impl::iterator_type iter(os);
    
    if (! karma::generate(iter, ngram_count_set_generator_impl::instance(), x.ngrams))
     throw std::runtime_error("failed ngram count set generation!");
    
    return os;
  }
  
  std::istream& operator>>(std::istream& is, NGramCountSet& x)
  {
    std::string line;
    
    if (utils::getline(is, line))
      x.assign(line);
    
    return is;
  }
  
  
};
