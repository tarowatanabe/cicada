#ifndef __CICADA__FILTER_EXTRACT_IMPL__HPP__
#define __CICADA__FILTER_EXTRACT_IMPL__HPP__ 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <string>
#include <vector>
#include <iostream>

#include <utils/hashmurmur.hpp>


struct RootCount
{
  typedef std::string label_type;
  typedef double count_type;
  typedef std::vector<count_type, std::allocator<count_type> > counts_type;
  
  label_type  label;
  counts_type counts;
  
  double observed_joint;
  double observed;
  
  RootCount() : label(), counts() {}
  RootCount(const label_type& __label) : label(__label), counts() {}
  
  void clear()
  {
    label.clear();
    counts.clear();
    observed_joint = 0;
    observed = 0;
  }
  
  friend
  size_t  hash_value(RootCount const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    return hasher_type()(x.label.begin(), x.label.end(), 0);
  }
  
  friend
  bool operator==(const RootCount& x, const RootCount& y) 
  {
    return x.label == y.label;
  }
  
  friend
  bool operator!=(const RootCount& x, const RootCount& y) 
  {
    return x.label != y.label;
  }
  
  friend
  bool operator<(const RootCount& x, const RootCount& y)
  {
    return x.label < y.label;
  }

  friend
  bool operator>(const RootCount& x, const RootCount& y)
  {
    return y < x;
  }
};

struct PhrasePair
{
  typedef std::string phrase_type;
  typedef double count_type;
  typedef std::vector<count_type, std::allocator<count_type> > counts_type;
  
  phrase_type source;
  phrase_type target;
  counts_type counts;
  counts_type counts_source;
  counts_type counts_target;
  double observed_source;
  double observed_target;
  double lexicon_target_source;
  double lexicon_source_target;
  
  PhrasePair()
    : source(), target(), counts(), counts_source(), counts_target() {}

  void clear()
  {
    source.clear();
    target.clear();
    counts.clear();
    counts_source.clear();
    counts_target.clear();
    observed_source = 0;
    observed_target = 0;
    lexicon_target_source = 0;
    lexicon_source_target = 0;
  }

  friend
  size_t  hash_value(PhrasePair const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
   
    return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
  }


  friend
  bool operator==(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source == y.source && x.target == y.target;
  }
  
  friend
  bool operator!=(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source != y.source || x.target != y.target;
  }
  
  friend
  bool operator<(const PhrasePair& x, const PhrasePair& y)
  {
    return (x.source < y.source || (!(y.source < x.source) && x.target < y.target));
  }

  friend
  bool operator>(const PhrasePair& x, const PhrasePair& y)
  {
    return y < x;
  }

};

BOOST_FUSION_ADAPT_STRUCT(
			  RootCount,
			  (RootCount::label_type,  label)
			  (RootCount::counts_type, counts)
			  (RootCount::count_type,  observed_joint)
			  (RootCount::count_type, observed)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePair,
			  (PhrasePair::phrase_type, source)
			  (PhrasePair::phrase_type, target)
			  (PhrasePair::counts_type, counts)
			  (PhrasePair::counts_type, counts_source)
			  (PhrasePair::counts_type, counts_target)
			  (PhrasePair::count_type,  observed_source)
			  (PhrasePair::count_type,  observed_target)
			  (PhrasePair::count_type,  lexicon_target_source)
			  (PhrasePair::count_type,  lexicon_source_target)
			  )

struct RootCountParser
{
  typedef RootCount root_count_type;
  
  typedef root_count_type::label_type  label_type;
  typedef root_count_type::counts_type counts_type;

  RootCountParser() : grammar() {}
  RootCountParser(const RootCountParser& x) : grammar() {}
  
  template <typename Iterator>
  struct root_count_parser : boost::spirit::qi::grammar<Iterator, root_count_type(), boost::spirit::standard::space_type>
  {
    root_count_parser() : root_count_parser::base_type(root_count)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::phrase_parse;
      using qi::lexeme;
      using qi::hold;
      using standard::char_;
      using qi::double_;
      using standard::space;
      
      label %= lexeme[+(char_ - (space >> "|||" >> space))];
      counts %= +double_;
      root_count %= label >> "|||" >> counts >> "|||" >> double_ >> double_;
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> label;
    boost::spirit::qi::rule<Iterator, counts_type(), boost::spirit::standard::space_type> counts;
    boost::spirit::qi::rule<Iterator, root_count_type(), boost::spirit::standard::space_type> root_count;
  };

  bool operator()(std::istream& is, root_count_type& root_count)
  {
    root_count.clear();
    
    std::string line;
    if (! getline(is, line)) return false;
    
    return operator()(line, root_count);
  }
  
  bool operator()(const std::string& line, root_count_type& root_count)
  {
    root_count.clear();
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    const bool result =  boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, root_count);
    return result && iter == end;
  }
  
  root_count_parser<std::string::const_iterator> grammar;
};

struct PhrasePairParser
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  PhrasePairParser() : grammar() {}
  PhrasePairParser(const PhrasePairParser& x) : grammar() {}

  template <typename Iterator>
  struct phrase_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_type(), boost::spirit::standard::space_type>
  {
    phrase_pair_parser() : phrase_pair_parser::base_type(phrase_pair)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::lexeme;
      using qi::lit;
      using qi::hold;
      using qi::repeat;
      using qi::double_;
      using qi::int_;
      
      using standard::char_;
      using standard::space;
      
      phrase %= lexeme[+(char_ - (space >> "|||" >> space))];
      counts %= +double_;
      
      phrase_pair %= (phrase
		      >> "|||" >> phrase
		      >> "|||" >> counts
		      >> "|||" >> counts
		      >> "|||" >> counts
		      >> "|||" >> double_ >> double_
		      >> "|||" >> double_ >> double_);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, counts_type(), boost::spirit::standard::space_type> counts;
    boost::spirit::qi::rule<Iterator, phrase_pair_type(), boost::spirit::standard::space_type> phrase_pair;
  };

  bool operator()(std::istream& is, phrase_pair_type& phrase_pair)
  {
    phrase_pair.clear();
    
    std::string line;
    if (! getline(is, line)) return false;
    
    return operator()(line, phrase_pair);
  }
  
  bool operator()(const std::string& line, phrase_pair_type& phrase_pair)
  {
    phrase_pair.clear();
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, phrase_pair);
    
    if (result && iter == end)
      return true;
    else {
      std::cerr << "WARNING: parsing failed: " << line << std::endl;
      return false;
    }
  }
  
  phrase_pair_parser<std::string::const_iterator> grammar;
};

struct PhrasePairGenerator
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  PhrasePairGenerator() : grammar() {}
  PhrasePairGenerator(const PhrasePairGenerator& x) : grammar() {}

  
  template <typename Iterator>
  struct phrase_pair_generator : boost::spirit::karma::grammar<Iterator, phrase_pair_type()>
  {
    phrase_pair_generator() : phrase_pair_generator::base_type(phrase_pair)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::repeat;
      using standard::char_;
      using karma::int_;
      using standard::space;
      
      phrase %= +char_;
      counts %= double20 % ' ';
      phrase_pair %= (phrase
		      << " ||| " << phrase
		      << " ||| " << counts
		      << " ||| " << counts
		      << " ||| " << counts
		      << " ||| " << double20 << ' ' << double20
		      << " ||| " << double20 << ' ' << double20);
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, std::string()> phrase;
    boost::spirit::karma::rule<Iterator, counts_type()> counts;
    boost::spirit::karma::rule<Iterator, phrase_pair_type()> phrase_pair;
  };

  typedef std::ostream_iterator<char> iterator_type;
  
  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair)
  {
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, phrase_pair))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  phrase_pair_generator<iterator_type> grammar;
};




struct ExtractRootPhrase
{
  const std::string& operator()(const std::string& phrase) const
  {
    static const std::string __label("[x]");
    return __label;
  }
};

struct ExtractRootSCFG
{
  // extract the first word...
  std::string operator()(const std::string& phrase) const
  {
    using boost::spirit::standard::char_;
    using boost::spirit::standard::space;
    
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();

    std::string label;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end,
							boost::spirit::qi::lexeme[+(char_ - space)],
							space,
							label);
    if (! result)
      throw std::runtime_error("no label?");
    
    return label;
  }
};

struct ExtractPhraseSCFG
{
  // extract the non-first word...
  std::pair<std::string, std::string> operator()(const std::string& phrase) const
  {
    using boost::spirit::standard::char_;
    using boost::spirit::standard::space;
    using boost::spirit::qi::lexeme;
    
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();
    
    std::pair<std::string, std::string> label_pair;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end,
							lexeme[+(char_ - space)] >> lexeme[+(char_ - space) >> *char_],
							space,
							label_pair);
    if (! result || iter != end)
      throw std::runtime_error("no label?");
    
    return label_pair;
  }
};

struct ExtractRootGHKM
{
  // extract the first word...
  std::string operator()(const std::string& phrase) const
  {
    using boost::spirit::standard::char_;
    using boost::spirit::standard::space;
    
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();
    
    std::string label;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end,
							boost::spirit::qi::lexeme[+((char_ - space - '(') | "\\(")],
							space,
							label);
    if (! result)
      throw std::runtime_error("no label?");
    
    return label;
  }
};

#endif

