//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE 1
#define PHOENIX_THREADSAFE 1

#include "lattice.hpp"
#include "vocab.hpp"

#include "unite.hpp"

#include <boost/numeric/conversion/bounds.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_statement.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/thread.hpp>

#include "utils/array_power2.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/config.hpp"
#include "utils/thread_specific_ptr.hpp"
#include "utils/python_string_parser.hpp"
#include "utils/json_string_parser.hpp"
#include "utils/json_string_generator.hpp"
#include "utils/getline.hpp"

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Lattice::arc_type,
			  (cicada::Lattice::symbol_type, label)
			  (cicada::Lattice::feature_set_type, features)
			  (cicada::Lattice::attribute_set_type, attributes)
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
    
    num_edge = 0;
    for (int node = lattice.size() - 1; node >= 0; -- node) {
      num_edge += lattice[node].size();
      
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

    typedef Lattice::feature_set_type   feature_set_type;
    typedef Lattice::attribute_set_type attribute_set_type;
    
    lattice_grammar_parser() : lattice_grammar_parser::base_type(lattice_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;

      attribute_data = (data_string [qi::_val = data_string_cache(qi::_1)]
			| data_double [qi::_val = qi::_1]
			| data_int [qi::_val = qi::_1]);
      attribute  %= attribute_key >> ':' >> attribute_data;
      attributes %= qi::hold[qi::lit(',') >> '{' >> -(attribute % ',') >> '}'] | qi::eps;
      
      jlf_lattice_score %= jlf_label >> ':' >> qi::double_;
      plf_lattice_score %= qi::attr("lattice-cost") >> qi::double_;

      jlf_lattice_scores %= '{' >> -(jlf_lattice_score % ',') >> '}';
      
      jlf_lattice_arc = (('[' >> jlf_label
			  >> ',' >> jlf_lattice_scores
			  >> attributes
			  >> ',' >> qi::int_
			   >> ']')
			 [qi::_val = phoenix::construct<Lattice::arc_type>(qi::_1,
									   qi::_2,
									   phoenix::construct<attribute_set_type>(phoenix::begin(qi::_3),
														  phoenix::end(qi::_3)),
									   qi::_4)]);
      plf_lattice_arc = (('('
			  >> plf_label
			  >> ',' >> qi::repeat(1)[plf_lattice_score]
			  >> ',' >> qi::int_ >> ')')
			 [qi::_val = phoenix::construct<Lattice::arc_type>(qi::_1,
									   phoenix::construct<feature_set_type>(phoenix::begin(qi::_2),
														phoenix::end(qi::_2)),
									   qi::_3)]);
      
      jlf_lattice_set %= '[' >> -(jlf_lattice_arc % ',') >> ']';
      plf_lattice_set %= '(' >> +(plf_lattice_arc >> ',') >> ')';

      lattice_grammar %= qi::hold['(' >> *(plf_lattice_set >> ',') >> ')'] | ('[' >> -(jlf_lattice_set % ',') >> ']');
    }
    
    typedef boost::spirit::standard::space_type space_type;
    
    typedef attribute_set_type::data_type attribute_data_type;
    typedef std::pair<std::string, attribute_data_type> attribute_parsed_type;
    typedef std::vector<attribute_parsed_type> attribute_parsed_set_type;

    struct attribute_string_cache_type : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      typedef std::string string_type;
      
      typedef utils::array_power2<string_type, 1024 * 4, std::allocator<string_type> > string_set_type;
      
      template <typename >
      struct result { typedef const string_type& type; };
      
      const string_type& operator()(const string_type& x) const
      {
	const size_t pos = hasher_type::operator()(x.begin(), x.end(), 0) & (caches.size() - 1);
	
	string_type& ret = const_cast<string_type&>(caches[pos]);
	if (ret != x)
	  ret = x;
	
	return ret;
      }
      
      string_set_type caches;
    };
    
    boost::spirit::qi::int_parser<attribute_set_type::int_type, 10, 1, -1>                   data_int;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > data_double;
    utils::json_string_parser<Iterator>                                                      data_string;
    boost::phoenix::function<attribute_string_cache_type>                                    data_string_cache;
    
    utils::json_string_parser<Iterator>                                            attribute_key;
    boost::spirit::qi::rule<Iterator, attribute_set_type::data_type(), space_type> attribute_data;
    
    boost::spirit::qi::rule<Iterator, attribute_parsed_type(), space_type> attribute;
    boost::spirit::qi::rule<Iterator, attribute_parsed_set_type(), space_type> attributes;

    
    utils::json_string_parser<Iterator>   jlf_label;
    utils::python_string_parser<Iterator> plf_label;
    
    boost::spirit::qi::rule<Iterator, std::pair<std::string, double >(), space_type> jlf_lattice_score;
    boost::spirit::qi::rule<Iterator, feature_set_type(), space_type>                jlf_lattice_scores;
    boost::spirit::qi::rule<Iterator, Lattice::arc_type(), space_type>               jlf_lattice_arc;
    boost::spirit::qi::rule<Iterator, Lattice::arc_set_type(), space_type>           jlf_lattice_set;
    
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
    if (utils::getline(is, line) && ! line.empty())
      x.assign(line);
    
    return is;
  }

  // generator...
  template <typename Iterator>
  struct lattice_grammar_generator : boost::spirit::karma::grammar<Iterator, Lattice::lattice_type()>
  {
    typedef Lattice::feature_set_type   feature_set_type;
    typedef Lattice::attribute_set_type attribute_set_type;
    
    lattice_grammar_generator() : lattice_grammar_generator::base_type(lattice_grammar)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      attribute_data %= data_int | data_double | data_string;
      attribute %= attribute_key << ':' << attribute_data;
      
      // json grammar...
      lattice_score %= label_double_quote << ':' << double10;
      
      lattice_arc %= ('['
		      << label_double_quote
		      << karma::lit(", ") << '{' << -(lattice_score % ", ") << '}'
		      << -karma::buffer[karma::lit(", ") << '{' << (attribute % ", ") << '}']
		      << karma::lit(", ") << karma::int_ << ']');
      lattice_set %= '[' << -(lattice_arc % ", ") << ']';
      lattice_grammar %= '[' << -(lattice_set % ", ") << ']';
    }
    
    utils::json_string_generator<Iterator> label_double_quote;
    
    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 10;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision>                 data_double;
    boost::spirit::karma::int_generator<attribute_set_type::int_type, 10, false> data_int;
    utils::json_string_generator<Iterator>                                       data_string;
    
    utils::json_string_generator<Iterator>                                 attribute_key;
    boost::spirit::karma::rule<Iterator, attribute_set_type::data_type()>  attribute_data;
    boost::spirit::karma::rule<Iterator, attribute_set_type::value_type()> attribute;

    boost::spirit::karma::real_generator<double, real_precision> double10;

    boost::spirit::karma::rule<Iterator, feature_set_type::value_type()> lattice_score;
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
