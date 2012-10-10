//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXTRACT_SCORE_IMPL__HPP__
#define __CICADA__EXTRACT_SCORE_IMPL__HPP__ 1

#include <unistd.h>
#include <cstring>

#include <memory>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>
#include <boost/tokenizer.hpp>
#include <boost/range.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <utility>
#include <queue>
#include <deque>
#include <sstream>
#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/alignment.hpp>
#include <cicada/tree_rule.hpp>

#include <utils/piece.hpp>
#include <utils/space_separator.hpp>
#include <utils/tempfile.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/compress_stream.hpp>
#include <utils/array_power2.hpp>
#include <utils/repository.hpp>
#include <utils/map_file.hpp>
#include <utils/alloc_vector.hpp>
#include <utils/unordered_set.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/base64.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/simple_vector.hpp>
#include <utils/double_base64_parser.hpp>
#include <utils/double_base64_generator.hpp>
#include <utils/dense_hash_map.hpp>
#include <utils/map_file_allocator.hpp>
#include <utils/group_aligned_code.hpp>

class RootCount
{
public:
  typedef std::string label_type;
  typedef utils::simple_vector<double, std::allocator<double> > counts_type;

  label_type  label;
  counts_type counts;

  double observed;

  RootCount(const label_type& __label) : label(__label), counts(), observed(0) {}
  RootCount() : label(), counts(), observed(0) {}

  void clear()
  {
    label.clear();
    counts.clear();
    observed = 0;
  }

  void swap(RootCount& x)
  {
    label.swap(x.label);
    counts.swap(x.counts);

    std::swap(observed, x.observed);
  }


  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    const size_t size_max = utils::bithack::max(counts.size(), size_t(std::distance(first, last)));

    counts.resize(size_max, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  std::ostream& operator<<(std::ostream& os, const RootCount& x)
  {
    os << x.label << " ||| ";
    std::copy(x.counts.begin(), x.counts.end(), std::ostream_iterator<double>(os, " "));
    os << "||| " << x.observed;

    return os;
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

class PhrasePair
{
public:
  typedef std::string phrase_type;
  typedef utils::simple_vector<double, std::allocator<double> > counts_type;
  typedef cicada::Alignment alignment_type;
  typedef alignment_type::point_type point_type;

  phrase_type    source;
  phrase_type    target;
  alignment_type alignment;
  counts_type    counts;

  PhrasePair() : source(), target(), alignment(), counts() {}

  void clear()
  {
    source.clear();
    target.clear();
    alignment.clear();
    counts.clear();
  }

  void swap(PhrasePair& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    alignment.swap(x.alignment);
    counts.swap(x.counts);
  }

  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    const size_t size_max = utils::bithack::max(counts.size(), size_t(std::distance(first, last)));

    counts.resize(size_max, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  size_t  hash_value(PhrasePair const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    return hasher_type()(x.source.begin(), x.source.end(),
			 hasher_type()(x.target.begin(), x.target.end(),
				       hasher_type()(x.alignment.begin(), x.alignment.end(), 0)));
  }

  friend
  bool operator==(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
  }

  friend
  bool operator!=(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
  }

  friend
  bool operator<(const PhrasePair& x, const PhrasePair& y)
  {
    return (x.source < y.source
	    || (!(y.source < x.source)
		&& (x.target < y.target
		    || (!(y.target < x.target)
			&& x.alignment < y.alignment))));
  }

  friend
  bool operator>(const PhrasePair& x, const PhrasePair& y)
  {
    return y < x;
  }

  friend
  std::ostream& operator<<(std::ostream& os, const PhrasePair& x)
  {
    os << x.source
       << " ||| " << x.target
       << " ||| " << x.alignment
       << " |||";

    PhrasePair::counts_type::const_iterator citer_end = x.counts.end();
    for (PhrasePair::counts_type::const_iterator citer = x.counts.begin(); citer != citer_end; ++ citer)
      os << ' ' << *citer;

    return os;
  }

};

class PhrasePairSimple
{
public:
  typedef PhrasePair phrase_pair_type;
  typedef phrase_pair_type::phrase_type phrase_type;
  typedef phrase_pair_type::counts_type counts_type;

  phrase_type    source;
  phrase_type    target;
  counts_type    counts;

  PhrasePairSimple() : source(), target(), counts() {}
  PhrasePairSimple(const phrase_type& __source, const phrase_type& __target, const counts_type& __counts)
    : source(__source), target(__target), counts(__counts) {}

  void clear()
  {
    source.clear();
    target.clear();
    counts.clear();
  }

  void swap(PhrasePairSimple& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    counts.swap(x.counts);
  }

  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    const size_t size_max = utils::bithack::max(counts.size(), size_t(std::distance(first, last)));

    counts.resize(size_max, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  size_t  hash_value(PhrasePairSimple const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
  }


  friend
  bool operator==(const PhrasePairSimple& x, const PhrasePairSimple& y) 
  {
    return x.source == y.source && x.target == y.target;
  }

  friend
  bool operator!=(const PhrasePairSimple& x, const PhrasePairSimple& y) 
  {
    return x.source != y.source || x.target != y.target;
  }

  friend
  bool operator<(const PhrasePairSimple& x, const PhrasePairSimple& y)
  {
    return (x.source < y.source || (!(y.source < x.source) && x.target < y.target));
  }

  friend
  bool operator>(const PhrasePairSimple& x, const PhrasePairSimple& y)
  {
    return y < x;
  }
};

class PhraseCount
{
public:
  typedef PhrasePair phrase_pair_type;
  typedef phrase_pair_type::phrase_type phrase_type;
  typedef phrase_pair_type::counts_type counts_type;

  phrase_type    phrase;
  counts_type    counts;

  PhraseCount() : phrase(), counts() {}
  PhraseCount(const phrase_type& __phrase, const counts_type& __counts)
    : phrase(__phrase), counts(__counts) {}

  void clear()
  {
    phrase.clear();
    counts.clear();
  }

  void swap(PhraseCount& x)
  {
    phrase.swap(x.phrase);
    counts.swap(x.counts);
  }

  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    const size_t size_max = utils::bithack::max(counts.size(), size_t(std::distance(first, last)));

    counts.resize(size_max, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  size_t  hash_value(PhraseCount const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    return hasher_type()(x.phrase.begin(), x.phrase.end(), 0);
  }


  friend
  bool operator==(const PhraseCount& x, const PhraseCount& y) 
  {
    return x.phrase == y.phrase;
  }

  friend
  bool operator!=(const PhraseCount& x, const PhraseCount& y) 
  {
    return x.phrase != y.phrase;
  }

  friend
  bool operator<(const PhraseCount& x, const PhraseCount& y)
  {
    return x.phrase < y.phrase;
  }

  friend
  bool operator>(const PhraseCount& x, const PhraseCount& y)
  {
    return y < x;
  }
};


namespace std
{
  inline
  void swap(RootCount& x, RootCount& y)
  {
    x.swap(y);
  }

  inline
  void swap(PhrasePair& x, PhrasePair& y)
  {
    x.swap(y);
  }

  inline
  void swap(PhrasePairSimple& x, PhrasePairSimple& y)
  {
    x.swap(y);
  }

  inline
  void swap(PhraseCount& x, PhraseCount& y)
  {
    x.swap(y);
  }

};

BOOST_FUSION_ADAPT_STRUCT(
			  RootCount,
			  (RootCount::label_type, label)
			  (RootCount::counts_type, counts)
			  (double, observed)
			  )

BOOST_FUSION_ADAPT_STRUCT(PhrasePair::point_type,
			  (PhrasePair::alignment_type::index_type, source)
			  (PhrasePair::alignment_type::index_type, target)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePair,
			  (PhrasePair::phrase_type, source)
			  (PhrasePair::phrase_type, target)
			  (PhrasePair::alignment_type, alignment)
			  (PhrasePair::counts_type, counts)
			  )
BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePairSimple,
			  (PhrasePairSimple::phrase_type, source)
			  (PhrasePairSimple::phrase_type, target)
			  (PhrasePairSimple::counts_type, counts)
			  )
BOOST_FUSION_ADAPT_STRUCT(
			  PhraseCount,
			  (PhrasePairSimple::phrase_type, phrase)
			  (PhrasePairSimple::counts_type, counts)
			  )

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
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();

    std::string label;

    const bool result = qi::phrase_parse(iter, end,
					 qi::lexeme[+(standard::char_ - standard::space)],
					 standard::space,
					 label);
    if (! result)
      throw std::runtime_error("no label?");

    return label;
  }
};

struct ExtractRootGHKM
{
  // extract the first word...
  std::string operator()(const std::string& phrase) const
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();

    std::string label;

    const bool result = qi::phrase_parse(iter, end,
					 qi::lexeme[+((standard::char_ - standard::space - '(') | "\\(")],
					 standard::space,
					 label);
    if (! result)
      throw std::runtime_error("no label?");

    return label;
  }
};


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

      label %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];

      count %= 'B' >> count_base64 | qi::double_;

      counts %= +count;
      root_count %= label >> "|||" >> counts >> "|||" >> count;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> label;
    utils::double_base64_parser<Iterator> count_base64;
    boost::spirit::qi::rule<Iterator, double(), boost::spirit::standard::space_type> count;
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


struct RootCountGenerator
{
  typedef RootCount root_count_type;

  typedef root_count_type::label_type  label_type;
  typedef root_count_type::counts_type counts_type;

  std::ostream& operator()(std::ostream& os, const root_count_type& root_count) const
  {
    os << root_count.label << " |||";
    counts_type::const_iterator citer_end = root_count.counts.end();
    for (counts_type::const_iterator citer = root_count.counts.begin(); citer != citer_end; ++ citer) {
      os << " B";
      utils::encode_base64(*citer, std::ostream_iterator<char>(os));
    }

    os << " |||";
    os << " B";
    utils::encode_base64(root_count.observed, std::ostream_iterator<char>(os));

    return os;
  }
};

struct PhrasePairParser
{
  typedef PhrasePair phrase_pair_type;

  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
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

      phrase %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];
      point %= qi::int_ >> '-' >> qi::int_;
      alignment %= *point;
      
      counts %= +('B' >> count_base64 | qi::double_);
      phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> alignment >> "|||" >> counts;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    
    boost::spirit::qi::rule<Iterator, alignment_type::point_type(), boost::spirit::standard::space_type> point;
    boost::spirit::qi::rule<Iterator, alignment_type(), boost::spirit::standard::space_type> alignment;
    
    utils::double_base64_parser<Iterator> count_base64;
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
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::counts_type    counts_type;

  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair) const
  {
    os << phrase_pair.source
       << " ||| " << phrase_pair.target
       << " ||| " << phrase_pair.alignment
       << " |||";

    counts_type::const_iterator citer_end = phrase_pair.counts.end();
    for (counts_type::const_iterator citer = phrase_pair.counts.begin(); citer != citer_end; ++ citer) {
      os << " B";
      utils::encode_base64(*citer, std::ostream_iterator<char>(os));
    }

    return os;
  }
};


struct PhrasePairSimpleParser
{
  typedef PhrasePairSimple phrase_pair_type;

  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;

  PhrasePairSimpleParser() : grammar() {}
  PhrasePairSimpleParser(const PhrasePairSimpleParser& x) : grammar() {}

  template <typename Iterator>
  struct phrase_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_type(), boost::spirit::standard::space_type>
  {
    phrase_pair_parser() : phrase_pair_parser::base_type(phrase_pair)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      phrase %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];

      counts %= +('B' >> count_base64 | qi::double_);
      phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> counts;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    utils::double_base64_parser<Iterator> count_base64;
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

struct PhrasePairSimpleGenerator
{
  typedef PhrasePairSimple phrase_pair_type;

  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;

  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair) const
  {
    os << phrase_pair.source << " ||| " << phrase_pair.target << " |||";

    counts_type::const_iterator citer_end = phrase_pair.counts.end();
    for (counts_type::const_iterator citer = phrase_pair.counts.begin(); citer != citer_end; ++ citer) {
      os << " B";
      utils::encode_base64(*citer, std::ostream_iterator<char>(os));
    }

    return os;
  }
};


struct PhraseCountParser
{
  typedef PhraseCount phrase_count_type;
  
  typedef phrase_count_type::phrase_type    phrase_type;
  typedef phrase_count_type::counts_type    counts_type;
  
  PhraseCountParser() : grammar() {}
  PhraseCountParser(const PhraseCountParser& x) : grammar() {}
  
  template <typename Iterator>
  struct phrase_count_parser : boost::spirit::qi::grammar<Iterator, phrase_count_type(), boost::spirit::standard::space_type>
  {
    phrase_count_parser() : phrase_count_parser::base_type(phrase_count)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      phrase %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];
      
      counts %= +('B' >> count_base64 | qi::double_);
      phrase_count %= phrase >> "|||" >> counts;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    utils::double_base64_parser<Iterator> count_base64;
    boost::spirit::qi::rule<Iterator, counts_type(), boost::spirit::standard::space_type> counts;
    boost::spirit::qi::rule<Iterator, phrase_count_type(), boost::spirit::standard::space_type> phrase_count;
  };

  bool operator()(std::istream& is, phrase_count_type& phrase_count)
  {
    phrase_count.clear();
    
    std::string line;
    if (! getline(is, line)) return false;
    
    return operator()(line, phrase_count);
  }
  
  bool operator()(const std::string& line, phrase_count_type& phrase_count)
  {
    phrase_count.clear();
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, phrase_count);
    
    if (result && iter == end)
      return true;
    else {
      std::cerr << "WARNING: parsing failed: " << line << std::endl;
      return false;
    }
  }
  
  phrase_count_parser<std::string::const_iterator> grammar;
};

struct PhraseCountGenerator
{
  typedef PhraseCount phrase_count_type;
  
  typedef phrase_count_type::phrase_type    phrase_type;
  typedef phrase_count_type::counts_type    counts_type;

  std::ostream& operator()(std::ostream& os, const phrase_count_type& phrase_count) const
  {
    os << phrase_count.phrase << " |||";

    counts_type::const_iterator citer_end = phrase_count.counts.end();
    for (counts_type::const_iterator citer = phrase_count.counts.begin(); citer != citer_end; ++ citer) {
      os << " B";
      utils::encode_base64(*citer, std::ostream_iterator<char>(os));
    }

    return os;
  }
};

// compute observation counts for source-side
struct PhrasePairSource
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;

  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef RootCount        root_count_type;
  typedef PhrasePair       phrase_pair_type;
  typedef PhraseCount      phrase_count_type;
  typedef PhrasePairSimple simple_type;

  typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
  
  typedef utils::lockfree_list_queue<simple_type, std::allocator<simple_type> > queue_type;
  
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;
};

struct PhrasePairSourceMapper
{
  typedef PhrasePairSource map_reduce_type;

  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::simple_type      simple_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef PhrasePairParser phrase_pair_parser_type;
  
  hasher_type             hasher;
  phrase_pair_parser_type phrase_pair_parser;

  path_set_type paths;
  queue_ptr_set_type& queues;
  double max_malloc;
  int    debug;
  
  PhrasePairSourceMapper(const path_set_type& __paths,
			 queue_ptr_set_type& __queues,
			 const double __max_malloc,
			 const int __debug)
    : paths(__paths),
      queues(__queues),
      max_malloc(__max_malloc),
      debug(__debug) {}
  
    template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }

    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    phrase_pair_type phrase_pair;

    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_pair_parser(line, phrase_pair)) continue;
      
      if (counts.empty() || counts.back().source != phrase_pair.source)
	counts.push_back(simple_type(phrase_pair.source, phrase_pair.target, phrase_pair.counts));
      else if (counts.back().target != phrase_pair.target) {
	phrase_pair.source = counts.back().source;
	counts.push_back(simple_type(phrase_pair.source, phrase_pair.target, phrase_pair.counts));
      } else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }
  
  inline
  int loop_sleep(bool found, int non_found_iter)
  {
    if (! found) {
      boost::thread::yield();
      ++ non_found_iter;
    } else
      non_found_iter = 0;

    if (non_found_iter >= 64) {
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);

      non_found_iter = 0;
    }
    return non_found_iter;
  }
  
  void operator()()
  {
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;

    typedef std::deque<simple_type, std::allocator<simple_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;

    pqueue_type pqueue;

    istream_ptr_set_type   istreams(paths.size());
    buffer_stream_set_type buffer_streams(paths.size());

    size_t pos = 0;
    path_set_type::const_iterator piter_end = paths.end();
    for (path_set_type::const_iterator piter = paths.begin(); piter != piter_end; ++ piter, ++ pos) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->string());

      istreams[pos].reset(new istream_type(*piter, 1024 * 1024));
      
      buffer_stream_type* buffer_stream = &buffer_streams[pos];
      buffer_stream->second = &(*istreams[pos]);
      
      read_phrase_pair(*istreams[pos], buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    simple_type counts;

    int iter = 0;
    const int iteration_mask = (1 << 4) - 1;
    const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
    bool malloc_full = false;
    int non_found_iter = 0;
    
    while (! pqueue.empty()) {
      buffer_stream_type* buffer_stream(pqueue.top());
      pqueue.pop();
      
      simple_type& curr = buffer_stream->first.front();
      
      if (counts != curr) {
	if (! counts.counts.empty()) {
	  const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
	  queues[shard]->push_swap(counts);
	}
	
	if ((iter & iteration_mask) == iteration_mask)
	  malloc_full = (utils::malloc_stats::used() > malloc_threshold);
	
	++ iter;
	
	non_found_iter = loop_sleep(! malloc_full, non_found_iter);

	counts.swap(curr);
      } else
	counts.increment(curr.counts.begin(), curr.counts.end());
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! counts.counts.empty()) {
      const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
      queues[shard]->push_swap(counts);
    }
    
    // termination...
    std::vector<bool, std::allocator<bool> > terminated(queues.size(), false);
    counts.clear();
    
    while (1) {
      bool found = false;
      
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (! terminated[shard] && queues[shard]->push_swap(counts, true)) {
	  counts.clear();
	  
	  terminated[shard] = true;
	  found = true;
	}
      
      if (std::count(terminated.begin(), terminated.end(), true) == static_cast<int>(terminated.size())) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
};

template <typename ExtractRoot>
struct PhrasePairSourceReducer
{
  typedef PhrasePairSource map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::phrase_count_type phrase_count_type;
  typedef map_reduce_type::simple_type       simple_type;
  
  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;

  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef PhraseCountParser    phrase_count_parser_type;
  typedef PhraseCountGenerator phrase_count_generator_type;
  
  typedef ExtractRoot extract_root_type;
  
  phrase_count_parser_type    parser;
  phrase_count_generator_type generator;

  queue_ptr_set_type& queues;
  
  path_type      prefix;
  path_type&     path;
  
  root_count_set_type& joint_counts;
  root_count_set_type& root_counts;
  
  extract_root_type extract_root;
  
  double              max_malloc;
  int                 debug;
  
  PhrasePairSourceReducer(queue_ptr_set_type&  __queues,
			  const path_type&     __prefix,
			  path_type&           __path,
			  root_count_set_type& __joint_counts,
			  root_count_set_type& __root_counts,
			  const double __max_malloc,
			  const int    __debug)
    : queues(__queues),
      prefix(__prefix),
      path(__path),
      joint_counts(__joint_counts),
      root_counts(__root_counts),
      max_malloc(__max_malloc),
      debug(__debug) {}
  
  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first > y->first;
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first > y->first;
    }
  };
  
  void operator()()
  {
    typedef std::pair<simple_type, queue_type*> buffer_queue_type;
    typedef std::vector<buffer_queue_type*, std::allocator<buffer_queue_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_queue_type*, pqueue_base_type, greater_buffer<buffer_queue_type> > pqueue_type;
    
    pqueue_type pqueue;
    std::vector<buffer_queue_type, std::allocator<buffer_queue_type> > buffer_queues(queues.size());
    
    {
      size_t pos = 0;
      for (queue_ptr_set_type::iterator qiter = queues.begin(); qiter != queues.end(); ++ qiter, ++ pos) {
	queue_ptr_type& queue = *qiter;
	
	buffer_queue_type* buffer_queue = &buffer_queues[pos];
	
	queue->pop_swap(buffer_queue->first);
	buffer_queue->second = &(*queue);
	
	if (! buffer_queue->first.source.empty())
	  pqueue.push(buffer_queue);
      }
    }
    
    {
      const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.source.XXXXXX");
      utils::tempfile::insert(counts_file_tmp);
      const path_type counts_file = counts_file_tmp.string() + ".gz";
      utils::tempfile::insert(counts_file);
      
      path = counts_file;
    }
    
    utils::compress_ostream os(path, 1024 * 1024);
    simple_type counts;
    size_type observed = 0;
    
    root_count_set_type::iterator jiter;
    root_count_set_type::iterator riter;
    
    while (! pqueue.empty()) {
      buffer_queue_type* buffer_queue(pqueue.top());
      pqueue.pop();
      
      simple_type& curr = buffer_queue->first;
      
      if (curr.source != counts.source) {
	if (observed) {
	  counts.counts.push_back(observed);
	  
	  generator(os, phrase_count_type(counts.source, counts.counts)) << '\n';
	}
	
	// increment root_observed(lhs+rhs)
	jiter = joint_counts.insert(extract_root(curr.source)+extract_root(curr.target)).first;
	const_cast<root_count_type&>(*jiter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*jiter).observed += 1;
	
	// increment root_observed(lhs)
	riter = root_counts.insert(extract_root(curr.source)).first;
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*riter).observed += 1;
	
	counts.swap(curr);
	observed = 1;
      } else if (curr.target != counts.target) {
	// increment root_observed(lhs+rhs)
	jiter = joint_counts.insert(extract_root(curr.source)+extract_root(curr.target)).first;
	const_cast<root_count_type&>(*jiter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*jiter).observed += 1;
	
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	
	counts.target.swap(curr.target);
	counts.increment(curr.counts.begin(), curr.counts.end());
	
	observed += 1;
      } else {
	const_cast<root_count_type&>(*jiter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	
	counts.increment(curr.counts.begin(), curr.counts.end());
      }
      
      buffer_queue->second->pop_swap(buffer_queue->first);
      if (! buffer_queue->first.source.empty())
	pqueue.push(buffer_queue);
    }
    
    if (observed) {
      counts.counts.push_back(observed);
      
      generator(os, phrase_count_type(counts.source, counts.counts)) << '\n';
    }
    
  }
};

// modify counts... simply map from source to target, and collect counts
struct PhrasePairReverse
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;

  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef RootCount          root_count_type;
  typedef PhrasePair         phrase_pair_type;
  typedef PhrasePairSimple simple_type;

  typedef utils::chunk_vector<simple_type, 4096 / sizeof(simple_type), std::allocator<simple_type> >  simple_set_type;

  typedef utils::lockfree_list_queue<simple_type, std::allocator<simple_type> > queue_type;

  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;


};

struct PhrasePairReverseMapper
{
  typedef PhrasePairReverse map_reduce_type;

  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;

  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::phrase_pair_type phrase_pair_type;

  typedef map_reduce_type::simple_type     simple_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef PhrasePairParser phrase_pair_parser_type;

  hasher_type             hasher;
  phrase_pair_parser_type phrase_pair_parser;

  path_set_type paths;
  queue_ptr_set_type& queues;
  double max_malloc;
  int    debug;

  PhrasePairReverseMapper(const path_set_type& __paths,
			  queue_ptr_set_type& __queues,
			  const double __max_malloc,
			  const int __debug)
    : paths(__paths),
      queues(__queues),
      max_malloc(__max_malloc),
      debug(__debug) {}

  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }

    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    phrase_pair_type phrase_pair;

    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_pair_parser(line, phrase_pair)) continue;

      if (counts.empty() || counts.back().source != phrase_pair.source)
	counts.push_back(simple_type(phrase_pair.source, phrase_pair.target, phrase_pair.counts));
      else if (counts.back().target != phrase_pair.target) {
	phrase_pair.source = counts.back().source;
	counts.push_back(simple_type(phrase_pair.source, phrase_pair.target, phrase_pair.counts));
      } else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }

  inline
  int loop_sleep(bool found, int non_found_iter)
  {
    if (! found) {
      boost::thread::yield();
      ++ non_found_iter;
    } else
      non_found_iter = 0;

    if (non_found_iter >= 64) {
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);

      non_found_iter = 0;
    }
    return non_found_iter;
  }

  void operator()()
  {
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;

    typedef std::deque<simple_type, std::allocator<simple_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;

    pqueue_type pqueue;

    istream_ptr_set_type   istreams(paths.size());
    buffer_stream_set_type buffer_streams(paths.size());

    size_t pos = 0;
    path_set_type::const_iterator piter_end = paths.end();
    for (path_set_type::const_iterator piter = paths.begin(); piter != piter_end; ++ piter, ++ pos) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->string());

      istreams[pos].reset(new istream_type(*piter, 1024 * 1024));
      
      buffer_stream_type* buffer_stream = &buffer_streams[pos];
      buffer_stream->second = &(*istreams[pos]);
      
      read_phrase_pair(*istreams[pos], buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    simple_type counts;

    int iter = 0;
    const int iteration_mask = (1 << 4) - 1;
    const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
    bool malloc_full = false;
    int non_found_iter = 0;
    
    while (! pqueue.empty()) {
      buffer_stream_type* buffer_stream(pqueue.top());
      pqueue.pop();
      
      simple_type& curr = buffer_stream->first.front();
      
      if (counts != curr) {
	if (! counts.counts.empty()) {
	  // swap source and target!
	  counts.source.swap(counts.target);
	  
	  const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
	  queues[shard]->push_swap(counts);
	}
	
	if ((iter & iteration_mask) == iteration_mask)
	  malloc_full = (utils::malloc_stats::used() > malloc_threshold);
	
	++ iter;
	
	non_found_iter = loop_sleep(! malloc_full, non_found_iter);

	counts.swap(curr);
      } else
	counts.increment(curr.counts.begin(), curr.counts.end());
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! counts.counts.empty()) {
      // swap source and target!
      counts.source.swap(counts.target);
      
      const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
      queues[shard]->push_swap(counts);
    }
    
     // termination...
    std::vector<bool, std::allocator<bool> > terminated(queues.size(), false);
    counts.clear();
    
    while (1) {
      bool found = false;
      
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (! terminated[shard] && queues[shard]->push_swap(counts, true)) {
	  counts.clear();
	  
	  terminated[shard] = true;
	  found = true;
	}
      
      if (std::count(terminated.begin(), terminated.end(), true) == static_cast<int>(terminated.size())) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
};

struct PhrasePairReverseReducer
{
  typedef PhrasePairReverse map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::simple_type     simple_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef PhrasePairSimpleParser    simple_parser_type;
  typedef PhrasePairSimpleGenerator simple_generator_type;

  typedef utils::unordered_set<simple_type, boost::hash<simple_type>, std::equal_to<simple_type>,
			       std::allocator<simple_type> >::type simple_unique_type;
  
  simple_parser_type    parser;
  simple_generator_type generator;

  queue_type&    queue;
  path_type      prefix;
  path_set_type& paths;
  
  int            shard_size;
  double         max_malloc;
  int            debug;
  
  PhrasePairReverseReducer(queue_type&    __queue,
			   const path_type& __prefix,
			   path_set_type& __paths,
			   const int      __shard_size,
			   const double   __max_malloc,
			   const int      __debug)
    : queue(__queue),
      prefix(__prefix),
      paths(__paths),
      shard_size(__shard_size),
      max_malloc(__max_malloc),
      debug(__debug) {}

  struct less_file_size
  {
    bool operator()(const path_type& x, const path_type& y) const
    {
      return boost::filesystem::file_size(x) < boost::filesystem::file_size(y);
    }
  };


  // merge counts from two streams into os..
  void merge_counts(std::istream& is1, std::istream& is2, std::ostream& os)
  {
    simple_type simple1;
    simple_type simple2;
    
    bool parsed1 = parser(is1, simple1);
    bool parsed2 = parser(is2, simple2);
    
    while (parsed1 && parsed2) {
      if (simple1 < simple2) {
	generator(os, simple1) << '\n';
	parsed1 = parser(is1, simple1);
      } else if (simple2 < simple1) {
	generator(os, simple2) << '\n';
	parsed2 = parser(is2, simple2);
      } else {
	simple1.increment(simple2.counts.begin(), simple2.counts.end());
	generator(os, simple1) << '\n';
	
	parsed1 = parser(is1, simple1);
	parsed2 = parser(is2, simple2);
      }
    }
    
    // dump remaining...
    while (parsed1) {
      generator(os, simple1) << '\n';
      parsed1 = parser(is1, simple1);
    }
    
    while (parsed2) {
      generator(os, simple2) << '\n';
      parsed2 = parser(is2, simple2);
    }
  }
  
  
  // merge from smallest files...
  void merge_counts(path_set_type& paths)
  {
    if (paths.size() <= 128) return;

    while (paths.size() > 128) {
      
      // sort according to the file-size...
      std::sort(paths.begin(), paths.end(), less_file_size());
      
      const path_type file1 = paths.front();
      paths.erase(paths.begin());
      
      const path_type file2 = paths.front();
      paths.erase(paths.begin());
      
      const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.reversed.XXXXXX");
      utils::tempfile::insert(counts_file_tmp);
      const path_type counts_file = counts_file_tmp.string() + ".gz";
      utils::tempfile::insert(counts_file);
      
      paths.push_back(counts_file);
      {
	utils::compress_istream is1(file1, 1024 * 1024);
	utils::compress_istream is2(file2, 1024 * 1024);
	
	utils::compress_ostream os(counts_file, 1024 * 1024);
	os.exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
	
	merge_counts(is1, is2, os);
      }
      
      boost::filesystem::remove(file1);
      boost::filesystem::remove(file2);
      
      utils::tempfile::erase(file1);
      utils::tempfile::erase(file2);
    }
  }
  

  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };
  
  void dump_counts(path_set_type& paths, const simple_unique_type& counts)
  {
    // sort...!
    typedef std::vector<const simple_type*, std::allocator<const simple_type*> > sorted_type;
    
    // sorting...
    sorted_type sorted(counts.size());
    {
      sorted_type::iterator siter = sorted.begin();
      simple_unique_type::const_iterator citer_end = counts.end();
      for (simple_unique_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    std::sort(sorted.begin(), sorted.end(), less_ptr<simple_type>());
    
    // tempfile...
    const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.reversed.XXXXXX");
    utils::tempfile::insert(counts_file_tmp);
    const path_type counts_file = counts_file_tmp.string() + ".gz";
    utils::tempfile::insert(counts_file);
    
    paths.push_back(counts_file);

    // final dump!
    utils::compress_ostream os(counts_file, 1024 * 1024);
    os.exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      generator(os, *(*siter)) << '\n';
  }
  
  void operator()()
  {
    simple_type        reversed;
    simple_unique_type counts;

    int num_termination = 0;
    
    size_type min_counts_size = 0;
    const size_type iteration_mask = (1 << 5) - 1;
    const size_type malloc_threshold = size_type(max_malloc * 1024 * 1024 * 1024);
    
    for (size_type iteration = 0; /**/; ++ iteration) {
      reversed.clear();
      queue.pop_swap(reversed);
      
      if (reversed.source.empty()) {
	++ num_termination;
	
	if (num_termination == shard_size)
	  break;
	else
	  continue;
      }
      
      std::pair<simple_unique_type::iterator, bool> result = counts.insert(reversed);
      if (! result.second)
	const_cast<simple_type&>(*result.first).increment(reversed.counts.begin(), reversed.counts.end());
      
      if (((iteration & iteration_mask) == iteration_mask)
	  && ! counts.empty()
	  && (! min_counts_size || counts.size() > min_counts_size)
	  && (utils::malloc_stats::used() > malloc_threshold)) {
	if (! min_counts_size)
	  min_counts_size = counts.size() >> 1;
	    
	dump_counts(paths, counts);
	counts.clear();
	simple_unique_type(counts).swap(counts);
      }
    }
    
    if (! counts.empty()) {
      dump_counts(paths, counts);
      counts.clear();
      simple_unique_type(counts).swap(counts);
    }
    
    merge_counts(paths);
  }
};

struct PhrasePairTarget
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef double count_type;
  
  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;
  
  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
  
  typedef PhrasePair         phrase_pair_type;
  typedef phrase_pair_type::phrase_type phrase_type;
  typedef RootCount          root_count_type;
  typedef PhrasePairSimple simple_type;
  
  typedef utils::chunk_vector<simple_type, 4096 / sizeof(simple_type), std::allocator<simple_type> > simple_set_type;

  typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
  
  typedef utils::lockfree_list_queue<simple_type, std::allocator<simple_type> > queue_type;
  
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;
};

template <typename ExtractRoot>
struct PhrasePairTargetMapper
{
  typedef PhrasePairTarget map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;
  
  typedef map_reduce_type::count_type  count_type;
  
  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::phrase_type      phrase_type;
  typedef map_reduce_type::simple_type      simple_type;

  typedef map_reduce_type::simple_set_type simple_set_type;
    
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;

#if 1
  struct PhraseSet
  {
    typedef uint32_t length_type;
    typedef char     char_type;

    typedef std::vector<char_type, std::allocator<char_type> > buffer_type;
    typedef std::vector<length_type, std::allocator<length_type> > lengths_type;
    
    struct const_iterator
    {
      const_iterator() {}
      const_iterator(typename buffer_type::const_iterator  __biter,
		     typename buffer_type::const_iterator  __biter_end,
		     typename lengths_type::const_iterator __liter,
		     typename lengths_type::const_iterator __liter_end)
	: biter(__biter), biter_end(__biter_end),
	  liter(__liter), liter_end(__liter_end)
      {
	operator++();
      }

      const std::string& operator*() const
      {
	return curr;
      }
      
      const_iterator& operator++()
      {
	if (biter == biter_end || liter == liter_end) {
	  curr.clear();
	  return *this;
	}

	// replace!
	curr.replace(curr.begin() + *liter, curr.end(), biter, biter + *(liter + 1));
	biter += *(liter + 1);
	liter += 2;
	
	return *this;
      }
      
      friend
      bool operator==(const const_iterator& x, const const_iterator& y)
      {
	return x.curr == y.curr;
      }
      
      friend
      bool operator!=(const const_iterator& x, const const_iterator& y)
      {
	return x.curr != y.curr;
      }
      
      typename buffer_type::const_iterator  biter;
      typename buffer_type::const_iterator  biter_end;
      typename lengths_type::const_iterator liter;
      typename lengths_type::const_iterator liter_end;
      
      std::string curr;
    };
    
    void clear()
    {
      buffer.clear();
      lengths.clear();
      last.clear();
    }
    
    bool empty() const { return lengths.empty(); }
    size_t size() const { return lengths.size() >> 1; }

    void swap(PhraseSet& x)
    {
      buffer.swap(x.buffer);
      lengths.swap(x.lengths);
      last.swap(x.last);
    }
    
    void push_back(const std::string& x)
    {
      size_type pos = 0;
      const size_type pos_last = utils::bithack::min(last.size(), x.size());
      for (/**/; pos != pos_last && x[pos] == last[pos]; ++ pos) {}
      
      buffer.insert(buffer.end(), x.begin() + pos, x.end());
      lengths.push_back(pos);
      lengths.push_back(x.size() - pos);
      last = x;
    }

    const_iterator begin() const { return const_iterator(buffer.begin(), buffer.end(),
							 lengths.begin(), lengths.end()); }
    const_iterator end() const { return const_iterator(); }
    
    buffer_type  buffer;
    lengths_type lengths;
    std::string  last;
  };
#endif
  
#if 0
  struct PhraseSet
  {
    typedef uint32_t length_type;
    typedef char     char_type;

    typedef std::vector<char_type, std::allocator<char_type> > buffer_type;
    typedef std::vector<length_type, std::allocator<length_type> > lengths_type;
    
    struct const_iterator
    {
      const_iterator(typename buffer_type::const_iterator  __biter,
		     typename lengths_type::const_iterator __liter)
	: biter(__biter), liter(__liter) {}

      std::string operator*() const
      {
	return std::string(biter, biter + *liter);
      }
      
      const_iterator& operator++()
      {
	biter += *liter;
	++ liter;
	return *this;
      }
      
      friend
      bool operator==(const const_iterator& x, const const_iterator& y)
      {
	return x.biter == y.biter && x.liter == y.liter;
      }
      
      friend
      bool operator!=(const const_iterator& x, const const_iterator& y)
      {
	return x.biter != y.biter || x.liter != y.liter;
      }
      
      typename buffer_type::const_iterator  biter;
      typename lengths_type::const_iterator liter;
    };
    
    void clear()
    {
      buffer.clear();
      lengths.clear();
    }
    
    bool empty() const { return lengths.empty(); }
    size_t size() const { return lengths.size(); }

    void swap(PhraseSet& x)
    {
      buffer.swap(x.buffer);
      lengths.swap(x.lengths);
    }
    
    void push_back(const std::string& x)
    {
      buffer.insert(buffer.end(), x.begin(), x.end());
      lengths.push_back(x.size());
    }

    const_iterator begin() const { return const_iterator(buffer.begin(), lengths.begin()); }
    const_iterator end() const { return const_iterator(buffer.end(), lengths.end()); }
    
    buffer_type  buffer;
    lengths_type lengths;
  };
#endif

  typedef PhraseSet phrase_set_type;

  typedef PhrasePairSimpleParser    simple_parser_type;
  typedef PhrasePairSimpleGenerator simple_generator_type;
  
  typedef ExtractRoot extract_root_type;
  
  hasher_type hasher;
  
  extract_root_type extract_root;
  
  path_set_type paths;
  queue_ptr_set_type& queues;
  root_count_set_type& root_counts;
  
  double max_malloc;
  int    debug;
  
  PhrasePairTargetMapper(const path_set_type& __paths,
			 queue_ptr_set_type& __queues,
			 root_count_set_type& __root_counts,
			 const double __max_malloc,
			 const int __debug)
    : paths(__paths),
      queues(__queues),
      root_counts(__root_counts),
      max_malloc(__max_malloc),
      debug(__debug) {}

  
  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };

  simple_parser_type    parser;
  simple_generator_type generator;

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    simple_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! parser(line, phrase_pair)) continue;
      
      if (counts.empty() || counts.back().source != phrase_pair.source)
	counts.push_back(phrase_pair);
      else if (counts.back().target != phrase_pair.target) {
	phrase_pair.source = counts.back().source;
	counts.push_back(phrase_pair);
      } else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }
  
  void operator()()
  {
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<simple_type, std::allocator<simple_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;

    pqueue_type            pqueue;
    istream_ptr_set_type   istreams(paths.size());
    buffer_stream_set_type buffer_streams(paths.size());
    
    size_t pos = 0;
    for (path_set_type::const_iterator piter = paths.begin(); piter != paths.end(); ++ piter, ++ pos) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->string());
      
      istreams[pos].reset(new istream_type(*piter, 1024 * 1024));
      
      buffer_stream_type* buffer_stream = &buffer_streams[pos];
      buffer_stream->second = &(*istreams[pos]);
      
      read_phrase_pair(*istreams[pos], buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    phrase_set_type phrases;
    simple_type     counts;
    count_type      observed(0);
    
    root_count_set_type::iterator riter;
    
    int iter = 0;
    const int iteration_mask = (1 << 4) - 1;
    const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
    bool malloc_full = false;
    int non_found_iter = 0;

    while (! pqueue.empty()) {
      buffer_stream_type* buffer_stream(pqueue.top());
      pqueue.pop();
      
      simple_type& curr = buffer_stream->first.front();
      
      if (curr.source != counts.source) {
	
	if (! phrases.empty()) {
	  if (observed != phrases.size())
	    throw std::runtime_error("invalid observation count");

	  counts.counts.push_back(observed);
	  
	  typename phrase_set_type::const_iterator piter_end = phrases.end();
	  for (typename phrase_set_type::const_iterator piter = phrases.begin(); piter != piter_end; ++ piter) {
	    const std::string& phrase = *piter;
	    
	    const int shard = hasher(phrase.begin(), phrase.end(), 0) % queues.size();
	    
	    queues[shard]->push(simple_type(phrase, counts.source, counts.counts));
	  }

	  phrases.clear();
	}
	
	if ((iter & iteration_mask) == iteration_mask) {
	  malloc_full = (utils::malloc_stats::used() > malloc_threshold);
	  
	  if (malloc_full) {
	    phrase_set_type(phrases).swap(phrases);
	    //phrases.shrink();
	    
	    malloc_full = (utils::malloc_stats::used() > malloc_threshold);
	  }
	}
	
	++ iter;
	
	non_found_iter = loop_sleep(! malloc_full, non_found_iter);
	
	counts.swap(curr);
	
	phrases.push_back(counts.target);
	observed = 1;
	
	// increment root_observed(lhs)
	riter = root_counts.insert(extract_root(counts.source)).first;
	const_cast<root_count_type&>(*riter).increment(counts.counts.begin(), counts.counts.end());
	const_cast<root_count_type&>(*riter).observed += 1;
	
      } else if (curr.target != counts.target) {
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	
	counts.target.swap(curr.target);
	counts.increment(curr.counts.begin(), curr.counts.end());
	
	phrases.push_back(counts.target);
	observed += 1;
      } else {
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	counts.increment(curr.counts.begin(), curr.counts.end());
      }
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! phrases.empty()) {
      if (observed != phrases.size())
	throw std::runtime_error("invalid observation count");
      
      counts.counts.push_back(observed);
      
      typename phrase_set_type::const_iterator piter_end = phrases.end();
      for (typename phrase_set_type::const_iterator piter = phrases.begin(); piter != piter_end; ++ piter) {
	const std::string& phrase = *piter;
	
	const int shard = hasher(phrase.begin(), phrase.end(), 0) % queues.size();
	
	queues[shard]->push(simple_type(phrase, counts.source, counts.counts));
      }
      
      phrases.clear();
      phrase_set_type(phrases).swap(phrases);
      //phrases.shrink();
    }
    
    //
    // map-merged files by inversing
    //
    
    std::vector<bool, std::allocator<bool> > terminated(queues.size(), false);
    counts.clear();

    for (;;) {
      bool found = false;
      
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (! terminated[shard] && queues[shard]->push_swap(counts, true)) {
	  counts.clear();
	  
	  terminated[shard] = true;
	  found = true;
	}
      
      if (std::count(terminated.begin(), terminated.end(), true) == static_cast<int>(queues.size())) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
  
  inline
  int loop_sleep(bool found, int non_found_iter)
  {
    if (! found) {
      boost::thread::yield();
      ++ non_found_iter;
    } else
      non_found_iter = 0;
    
    if (non_found_iter >= 64) {
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    
      non_found_iter = 0;
    }
    return non_found_iter;
  }
};


struct PhrasePairTargetReducer
{
  typedef PhrasePairTarget map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;
  
  typedef map_reduce_type::count_type  count_type;
  
  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::simple_type    simple_type;

  typedef map_reduce_type::simple_set_type simple_set_type;
    
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;

  typedef PhrasePairSimpleParser    simple_parser_type;
  typedef PhrasePairSimpleGenerator simple_generator_type;

  simple_parser_type    parser;
  simple_generator_type generator;
  
  queue_type&    queue;
  path_type      prefix;
  path_set_type& paths;
  
  int            shard_size;
  double         max_malloc;
  int            debug;
  
  PhrasePairTargetReducer(queue_type&    __queue,
			  const path_type& __prefix,
			  path_set_type& __paths,
			  const int      __shard_size,
			  const double   __max_malloc,
			  const int      __debug)
    : queue(__queue),
      prefix(__prefix),
      paths(__paths),
      shard_size(__shard_size),
      max_malloc(__max_malloc),
      debug(__debug) {}
  

  struct less_file_size
  {
    bool operator()(const path_type& x, const path_type& y) const
    {
      return boost::filesystem::file_size(x) < boost::filesystem::file_size(y);
    }
  };


  // merge counts from two streams into os..
  void merge_counts(std::istream& is1, std::istream& is2, std::ostream& os)
  {
    simple_type simple1;
    simple_type simple2;
    
    bool parsed1 = parser(is1, simple1);
    bool parsed2 = parser(is2, simple2);
    
    while (parsed1 && parsed2) {
      if (simple1 < simple2) {
	generator(os, simple1) << '\n';
	parsed1 = parser(is1, simple1);
      } else if (simple2 < simple1) {
	generator(os, simple2) << '\n';
	parsed2 = parser(is2, simple2);
      } else {
	simple1.increment(simple2.counts.begin(), simple2.counts.end());
	generator(os, simple1) << '\n';
	
	parsed1 = parser(is1, simple1);
	parsed2 = parser(is2, simple2);
      }
    }
    
    // dump remaining...
    while (parsed1) {
      generator(os, simple1) << '\n';
      parsed1 = parser(is1, simple1);
    }
    
    while (parsed2) {
      generator(os, simple2) << '\n';
      parsed2 = parser(is2, simple2);
    }
  }
  
  // merge from smallest files...
  void merge_counts(path_set_type& paths)
  {
    if (paths.size() <= 128) return;

    while (paths.size() > 128) {
      
      // sort according to the file-size...
      std::sort(paths.begin(), paths.end(), less_file_size());
      
      const path_type file1 = paths.front();
      paths.erase(paths.begin());
      
      const path_type file2 = paths.front();
      paths.erase(paths.begin());
      
      const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.target.XXXXXX");
      utils::tempfile::insert(counts_file_tmp);
      const path_type counts_file = counts_file_tmp.string() + ".gz";
      utils::tempfile::insert(counts_file);
      
      paths.push_back(counts_file);
      {
	utils::compress_istream is1(file1, 1024 * 1024);
	utils::compress_istream is2(file2, 1024 * 1024);
	
	utils::compress_ostream os(counts_file, 1024 * 1024);
	os.exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
	
	merge_counts(is1, is2, os);
      }
      
      boost::filesystem::remove(file1);
      boost::filesystem::remove(file2);
      
      utils::tempfile::erase(file1);
      utils::tempfile::erase(file2);
    }
  }
  

  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };

  void dump_counts(path_set_type& paths, const simple_set_type& counts)
  {
    // sort...!
    typedef std::vector<const simple_type*, std::allocator<const simple_type*> > sorted_type;
    
    // sorting...
    sorted_type sorted(counts.size());
    {
      sorted_type::iterator siter = sorted.begin();
      simple_set_type::const_iterator citer_end = counts.end();
      for (simple_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    std::sort(sorted.begin(), sorted.end(), less_ptr<simple_type>());
    
    // tempfile...
    const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.target.XXXXXX");
    utils::tempfile::insert(counts_file_tmp);
    const path_type counts_file = counts_file_tmp.string() + ".gz";
    utils::tempfile::insert(counts_file);
    
    paths.push_back(counts_file);

    // final dump!
    utils::compress_ostream os(counts_file, 1024 * 1024);
    os.exceptions(std::ostream::eofbit | std::ostream::failbit | std::ostream::badbit);
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      generator(os, *(*siter)) << '\n';
  }

  void operator()()
  {
    // we already know that the entries are uniqued through the previous map-reduce!
    
    simple_type     target;
    simple_set_type counts;
    
    int num_termination = 0;
    
    size_type min_counts_size = 0;
    const size_type iteration_mask = (1 << 5) - 1;
    const size_type malloc_threshold = size_type(max_malloc * 1024 * 1024 * 1024);
    
    for (size_type iteration = 0; /**/; ++ iteration) {
      target.clear();
      queue.pop_swap(target);
      
      if (target.source.empty()) {
	++ num_termination;
	
	if (num_termination == shard_size)
	  break;
	else
	  continue;
      }
      
      counts.push_back(target);
      
      if (((iteration & iteration_mask) == iteration_mask)
	  && ! counts.empty()
	  && (! min_counts_size || counts.size() > min_counts_size)
	  && (utils::malloc_stats::used() > malloc_threshold)) {
	if (! min_counts_size)
	  min_counts_size = counts.size() >> 1;
	
	dump_counts(paths, counts);
	counts.clear();
      }
    }
    
    if (! counts.empty()) {
      dump_counts(paths, counts);
      counts.clear();
    }
    
    merge_counts(paths);
  }
};

// final scoring...
struct PhrasePairScore
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef double count_type;
  
  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;

  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef PhrasePair       phrase_pair_type;
  typedef PhrasePairSimple simple_type;
  typedef RootCount        root_count_type;
  typedef PhraseCount      phrase_count_type;
  
  typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
  
  typedef utils::lockfree_list_queue<phrase_pair_type, std::allocator<phrase_pair_type> > queue_type;
  
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;
  
};

struct PhrasePairScoreMapper
{
  typedef PhrasePairScore map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::count_type  count_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::simple_type    simple_type;
    
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef PhrasePairParser phrase_pair_parser_type;
  
  hasher_type             hasher;
  phrase_pair_parser_type phrase_pair_parser;
  
  path_set_type paths;
  queue_ptr_set_type& queues;
  double max_malloc;
  int debug;
  
  PhrasePairScoreMapper(const path_set_type& __paths,
			queue_ptr_set_type& __queues,
			const double __max_malloc,
			const int __debug)
    : paths(__paths), queues(__queues), max_malloc(__max_malloc), debug(__debug) {}
  
  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    phrase_pair_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_pair_parser(line, phrase_pair)) continue;

      if (counts.empty() || counts.back().source != phrase_pair.source)
	counts.push_back(phrase_pair);
      else if (counts.back().target != phrase_pair.target) {
	phrase_pair.source = counts.back().source;
	counts.push_back(phrase_pair);
      } else if (counts.back().alignment != phrase_pair.alignment) {
	phrase_pair.source = counts.back().source;
	phrase_pair.target = counts.back().target;
	counts.push_back(phrase_pair);
      } else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }

  inline
  int loop_sleep(bool found, int non_found_iter)
  {
    if (! found) {
      boost::thread::yield();
      ++ non_found_iter;
    } else
      non_found_iter = 0;
  
    if (non_found_iter >= 64) {
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    
      non_found_iter = 0;
    }
    return non_found_iter;
  }
  
  void operator()()
  {
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<phrase_pair_type, std::allocator<phrase_pair_type> > buffer_type;
    
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;
    
    pqueue_type            pqueue;
    istream_ptr_set_type   istreams(paths.size());
    buffer_stream_set_type buffer_streams(paths.size());
    
    size_t pos = 0;
    path_set_type::const_iterator piter_end = paths.end();
    for (path_set_type::const_iterator piter = paths.begin(); piter != piter_end; ++ piter, ++ pos) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->string());
      
      istreams[pos].reset(new istream_type(*piter, 1024 * 1024));
      
      buffer_stream_type* buffer_stream = &buffer_streams[pos];
      buffer_stream->second = &(*istreams[pos]);
      
      read_phrase_pair(*istreams[pos], buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    phrase_pair_type counts;
    
    int iter = 0;
    const int iteration_mask = (1 << 4) - 1;
    const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
    bool malloc_full = false;
    int non_found_iter = 0;

    while (! pqueue.empty()) {
      buffer_stream_type* buffer_stream(pqueue.top());
      pqueue.pop();

      phrase_pair_type& curr = buffer_stream->first.front();
      
      if (counts != curr) {
	if (! counts.counts.empty()) {
	  if (debug >= 5)
	    std::cerr << "score count mapper send" << std::endl;
	  
	  const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
	  queues[shard]->push_swap(counts);
	}
	
	
	if ((iter & iteration_mask) == iteration_mask)
	  malloc_full = (utils::malloc_stats::used() > malloc_threshold);
	
	++ iter;
	
	non_found_iter = loop_sleep(! malloc_full, non_found_iter);
	
	counts.swap(curr);
      } else
	counts.increment(curr.counts.begin(), curr.counts.end());
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! counts.counts.empty()) {
      if (debug >= 5)
	std::cerr << "score count mapper send" << std::endl;

      const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
      queues[shard]->push_swap(counts);
    }
    
    // termination...
    // we will terminate asynchronously...
    std::vector<bool, std::allocator<bool> > terminated(queues.size(), false);
    counts.clear();
    
    while (1) {
      bool found = false;
      
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (! terminated[shard] && queues[shard]->push_swap(counts, true)) {
	  counts.clear();
	  
	  terminated[shard] = true;
	  found = true;
	}
      
      if (std::count(terminated.begin(), terminated.end(), true) == static_cast<int>(terminated.size())) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
};

struct PhrasePairScoreReducer
{
  typedef PhrasePairScore map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::count_type  count_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type  phrase_pair_type;
  typedef map_reduce_type::phrase_count_type phrase_count_type;
  typedef map_reduce_type::simple_type       simple_type;

  typedef phrase_pair_type::alignment_type alignment_type;
  
  typedef std::vector<phrase_pair_type, std::allocator<phrase_pair_type> > phrase_pair_set_type;
  typedef std::vector<simple_type, std::allocator<simple_type> > simple_set_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;  

  typedef PhrasePairSimpleParser    simple_parser_type;
  typedef PhraseCountParser         phrase_parser_type;
  
  const path_type&     path_source;
  const path_set_type& path_targets;
  queue_ptr_set_type& queues;
  std::ostream& os;
  int debug;
  
  PhrasePairScoreReducer(const path_type& __path_source,
			 const path_set_type& __path_targets,
			 queue_ptr_set_type& __queues,
			 std::ostream& __os,
			 int __debug)
    : path_source(__path_source),
      path_targets(__path_targets),
      queues(__queues),
      os(__os),
      debug(__debug) {}
  
  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first > y->first;
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first > y->first;
    }
  };

  template <typename Tp>
  struct greater_stream
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };

  struct real_precision : boost::spirit::karma::real_policies<count_type>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };
  
  boost::spirit::karma::real_generator<count_type, real_precision> double20;
  
  void dump_phrase_pair(const phrase_pair_set_type& counts, const phrase_count_type& source, const simple_type& target)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    // counts are grouped by source/target
    
    const phrase_count_type::counts_type& counts_source = source.counts;
    const simple_type::counts_type&       counts_target = target.counts;

    os << target.source << " ||| " << target.target << " |||";
    
    std::ostream_iterator<char> oiter(os);
    
    phrase_pair_type counts_pair;
    
    phrase_pair_set_type::const_iterator citer_end = counts.end();
    for (phrase_pair_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      
      if (citer->source != source.phrase || citer->source != target.source)
	throw std::runtime_error("different source?");
      if (citer->target != target.target)
	throw std::runtime_error("different target?");
      
      if (! karma::generate(oiter, (karma::lit(' ') << karma::lit('(')
				    << -((karma::int_ << '-' << karma::int_) % ' ')
				    << karma::lit(')')),
			    citer->alignment))
	throw std::runtime_error("alignment generation failed");
      
      counts_pair.increment(citer->counts.begin(), citer->counts.end());
    }
    
    // cont(LHS RHS)
    if (! karma::generate(oiter, " ||| " << -(double20 % ' '), counts_pair.counts))
      throw std::runtime_error("generation failed");
    
    // count(LHS)
    if (! karma::generate(oiter, " ||| " << -(double20 % ' '),
			  boost::make_iterator_range(counts_source.begin(), counts_source.end() - 1)))
      throw std::runtime_error("generation failed");
    
    // count(RHS)
    if (! karma::generate(oiter, " ||| " << -(double20 % ' '),
			  boost::make_iterator_range(counts_target.begin(), counts_target.end() - 1)))
      throw std::runtime_error("generation failed");
    
    // observed(LHS) observed(RHS)
    if (! karma::generate(oiter, " ||| " << double20 << ' ' << double20 << '\n', counts_source.back(), counts_target.back()))
      throw std::runtime_error("generation failed");
  }
  
  simple_parser_type simple_parser;
  phrase_parser_type phrase_parser;

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    simple_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! simple_parser(line, phrase_pair)) continue;
      
      if (counts.empty() || counts.back().source != phrase_pair.source)
	counts.push_back(phrase_pair);
      else if (counts.back().target != phrase_pair.target) {
	phrase_pair.source = counts.back().source;
	counts.push_back(phrase_pair);
      } else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }
  
  template <typename Counts>
  void read_phrase(std::istream& is, Counts& counts)
  {
    std::string line;
    phrase_count_type phrase;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_parser(line, phrase)) continue;
      
      if (counts.empty() || counts.back().phrase != phrase.phrase)
	counts.push_back(phrase);
      else
	counts.back().increment(phrase.counts.begin(), phrase.counts.end());
    }
  }
  
  void operator()()
  {
    typedef std::pair<phrase_pair_type, queue_type*> buffer_queue_type;
    typedef std::vector<buffer_queue_type*, std::allocator<buffer_queue_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_queue_type*, pqueue_base_type, greater_buffer<buffer_queue_type> > pqueue_type;
    
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<simple_type, std::allocator<simple_type> > simple_buffer_type;
    typedef std::pair<simple_buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > simple_pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, simple_pqueue_base_type, greater_stream<buffer_stream_type> > simple_pqueue_type;

    typedef std::deque<phrase_count_type, std::allocator<phrase_count_type> > phrase_buffer_type;
    
    os.precision(20);
    
    pqueue_type pqueue;
    std::vector<buffer_queue_type, std::allocator<buffer_queue_type> > buffer_queues(queues.size());
    
    simple_pqueue_type     queue_target;
    istream_ptr_set_type   istreams(path_targets.size());
    buffer_stream_set_type buffer_streams(path_targets.size());
    
    {
      size_t pos = 0;
      for (queue_ptr_set_type::iterator qiter = queues.begin(); qiter != queues.end(); ++ qiter, ++ pos) {
	queue_ptr_type& queue = *qiter;
	
	buffer_queue_type* buffer_queue = &buffer_queues[pos];
	
	queue->pop_swap(buffer_queue->first);
	buffer_queue->second = &(*queue);
	
	if (! buffer_queue->first.source.empty())
	  pqueue.push(buffer_queue);
      }
    }

    {
      size_t pos = 0;
      for (path_set_type::const_iterator piter = path_targets.begin(); piter != path_targets.end(); ++ piter, ++ pos) {
	if (! boost::filesystem::exists(*piter))
	  throw std::runtime_error("no file? " + piter->string());
	
	istreams[pos].reset(new istream_type(*piter, 1024 * 1024));
	
	buffer_stream_type* buffer_stream = &buffer_streams[pos];
	buffer_stream->second = &(*istreams[pos]);
	
	read_phrase_pair(*istreams[pos], buffer_stream->first);
	
	if (! buffer_stream->first.empty())
	  queue_target.push(buffer_stream);
      }
    }

    utils::compress_istream is_source(path_source, 1024 * 1024);
    phrase_buffer_type buffer_source;
    
    read_phrase(is_source, buffer_source);
    
    phrase_pair_set_type counts;
    phrase_count_type    source;
    simple_type          target;
    
    while (! pqueue.empty()) {
      buffer_queue_type* buffer_queue(pqueue.top());
      pqueue.pop();

      phrase_pair_type& curr = buffer_queue->first;
      
      if (counts.empty() || counts.back().source != curr.source) {
	
	if (! counts.empty()) {
	  if (debug >= 4)
	    std::cerr << "score count reducer: " << counts.size() << std::endl;
	  
	  dump_phrase_pair(counts, source, target);
	  
	  counts.clear();
	}
	
	counts.push_back(curr);
	
	// next source..
	if (buffer_source.empty())
	  throw std::runtime_error("source counts and phrase pair do not match");
	
	source.swap(buffer_source.front());
	buffer_source.pop_front();
	
	if (buffer_source.empty())
	  read_phrase(is_source, buffer_source);
	
	// next target...
	if (queue_target.empty())
	  throw std::runtime_error("target counts and phrase pair do not match");
	
	buffer_stream_type* buffer_stream(queue_target.top());
	queue_target.pop();
	
	target.swap(buffer_stream->first.front());
	
	buffer_stream->first.pop_front();
	
	if (buffer_stream->first.empty())
	  read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
	
	if (! buffer_stream->first.empty())
	  queue_target.push(buffer_stream);
	
	if (counts.back().source != target.source)
	  throw std::runtime_error("source mismatch? " + counts.back().source + " target: " + target.source);
	
	if (counts.back().target != target.target)
	  throw std::runtime_error("target mismatch? " + counts.back().target + " target: " + target.target);
	
      } else if (counts.back().target != curr.target) {
	if (! counts.empty()) {
	  if (debug >= 4)
	    std::cerr << "score count reducer: " << counts.size() << std::endl;
	  
	  dump_phrase_pair(counts, source, target);
	  
	  counts.clear();
	}
	
	counts.push_back(curr);
	
	// next target...
	if (queue_target.empty())
	  throw std::runtime_error("simple counts and phrase pair do not match: sharing the same source");
	
	buffer_stream_type* buffer_stream(queue_target.top());
	queue_target.pop();
	
	target.swap(buffer_stream->first.front());
	
	buffer_stream->first.pop_front();
	
	if (buffer_stream->first.empty())
	  read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
	
	if (! buffer_stream->first.empty())
	  queue_target.push(buffer_stream);
	
	if (counts.back().source != target.source)
	  throw std::runtime_error("source mismatch? " + counts.back().source + " simple: " + target.source);
	
	if (counts.back().target != target.target)
	  throw std::runtime_error("target mismatch? " + counts.back().target + " simple: " + target.target);
	
      } else if (counts.back().alignment != curr.alignment)
	counts.push_back(curr);
      else
	counts.back().increment(curr.counts.begin(), curr.counts.end());
      
      buffer_queue->second->pop_swap(buffer_queue->first);
      if (! buffer_queue->first.source.empty())
	pqueue.push(buffer_queue);
    }
    
    if (! queue_target.empty())
      throw std::runtime_error("queue is empty but queue-target is not empty!");
    
    if (! counts.empty()) {
      if (debug >= 4)
	std::cerr << "score count reducer: " << counts.size() << std::endl;
      
      dump_phrase_pair(counts, source, target);
    }
  }
  
};

#endif
