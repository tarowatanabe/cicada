//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include <utils/sgi_hash_set.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/base64.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/chunk_vector.hpp>

#include <google/dense_hash_map>

class RootCount
{
public:
  typedef std::string label_type;
  typedef std::vector<double, std::allocator<double> > counts_type;

  label_type  label;
  counts_type counts;

  double observed_joint;
  double observed;

  RootCount(const label_type& __label) : label(__label), counts(), observed_joint(0), observed(0) {}
  RootCount() : label(), counts(), observed_joint(0), observed(0) {}

  void clear()
  {
    label.clear();
    counts.clear();
    observed_joint = 0;
    observed = 0;
  }

  void swap(RootCount& x)
  {
    label.swap(x.label);
    counts.swap(x.counts);

    std::swap(observed_joint, x.observed_joint);
    std::swap(observed, x.observed);
  }


  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    const size_t size_max = utils::bithack::max(counts.size(), size_t(std::distance(first, last)));

    counts.reserve(size_max);
    counts.resize(size_max, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  std::ostream& operator<<(std::ostream& os, const RootCount& x)
  {
    os << x.label << " ||| ";
    std::copy(x.counts.begin(), x.counts.end(), std::ostream_iterator<double>(os, " "));
    os << "||| " << x.observed_joint << ' ' << x.observed;

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
  typedef std::vector<double, std::allocator<double> > counts_type;
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

    counts.reserve(size_max);
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

class PhrasePairModified
{
public:
  typedef PhrasePair phrase_pair_type;
  typedef phrase_pair_type::phrase_type phrase_type;
  typedef phrase_pair_type::counts_type counts_type;

  phrase_type    source;
  phrase_type    target;
  counts_type    counts;

  PhrasePairModified() : source(), target(), counts() {}
  PhrasePairModified(const phrase_type& __source, const phrase_type& __target, const counts_type& __counts)
    : source(__source), target(__target), counts(__counts) {}

  void clear()
  {
    source.clear();
    target.clear();
    counts.clear();
  }

  void swap(PhrasePairModified& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    counts.swap(x.counts);
  }

  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    const size_t size_max = utils::bithack::max(counts.size(), size_t(std::distance(first, last)));

    counts.reserve(size_max);
    counts.resize(size_max, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  size_t  hash_value(PhrasePairModified const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;

    return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
  }


  friend
  bool operator==(const PhrasePairModified& x, const PhrasePairModified& y) 
  {
    return x.source == y.source && x.target == y.target;
  }

  friend
  bool operator!=(const PhrasePairModified& x, const PhrasePairModified& y) 
  {
    return x.source != y.source || x.target != y.target;
  }

  friend
  bool operator<(const PhrasePairModified& x, const PhrasePairModified& y)
  {
    return (x.source < y.source || (!(y.source < x.source) && x.target < y.target));
  }

  friend
  bool operator>(const PhrasePairModified& x, const PhrasePairModified& y)
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
  void swap(PhrasePairModified& x, PhrasePairModified& y)
  {
    x.swap(y);
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  RootCount,
			  (RootCount::label_type, label)
			  (RootCount::counts_type, counts)
			  (double, observed_joint)
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
			  PhrasePairModified,
			  (PhrasePairModified::phrase_type, source)
			  (PhrasePairModified::phrase_type, target)
			  (PhrasePairModified::counts_type, counts)
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

struct LexiconModel
{
  typedef cicada::Symbol word_type;
  typedef cicada::Vocab vocab_type;

  typedef boost::filesystem::path path_type;



  struct table_type
  {
    typedef google::dense_hash_map<word_type, double, boost::hash<word_type> , std::equal_to<word_type> > __table_type;

    table_type() : table() { table.set_empty_key(word_type()); }

    __table_type table;
  };

  typedef utils::alloc_vector<table_type, std::allocator<table_type> > table_set_type;

  LexiconModel() : smooth(1e-7), tables() {}
  LexiconModel(const path_type& path) : smooth(), tables() { open(path); }

  void open(const path_type& path)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    typedef boost::fusion::tuple<std::string, std::string, double> parsed_type;

    qi::rule<std::string::const_iterator, std::string(), standard::space_type> word;
    qi::rule<std::string::const_iterator, parsed_type(), standard::space_type> lexicon;

    word %= qi::lexeme[+(standard::char_ - standard::space)];
    lexicon %= word >> word >> qi::double_;

    smooth = std::numeric_limits<double>::infinity();
    tables.clear();

    utils::compress_istream is(path, 1024 * 1024);
    std::string line;
    parsed_type parsed;

    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();

      boost::fusion::get<0>(parsed).clear();
      boost::fusion::get<1>(parsed).clear();

      const bool result = qi::phrase_parse(iter, end, lexicon, standard::space, parsed);
      if (! result || iter != end) continue;


      tables[word_type(boost::fusion::get<1>(parsed)).id()].table[boost::fusion::get<0>(parsed)] = boost::fusion::get<2>(parsed);

      smooth = std::min(smooth, boost::fusion::get<2>(parsed));
    }

  }

  double operator()(const word_type& source, const word_type& target) const
  {
    if (tables.exists(source.id())) {
      const table_type& table = tables[source.id()];

      table_type::__table_type::const_iterator iter = table.table.find(target);
      return (iter != table.table.end() ? iter->second : smooth);
    } else 
      return smooth;
  }

  double smooth;
  table_set_type tables;
};

struct LexiconBase
{
  typedef PhrasePair::point_type     point_type;
  typedef PhrasePair::alignment_type alignment_type;

  typedef cicada::Sentence sentence_type;

  typedef std::vector<int, std::allocator<int> > align_set_type;
  typedef std::vector<align_set_type, std::allocator<align_set_type> > align_map_type;

  typedef LexiconModel lexicon_model_type;

  LexiconBase(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : lexicon_source_target(__lexicon_source_target),
      lexicon_target_source(__lexicon_target_source) {}

  const lexicon_model_type& lexicon_source_target;
  const lexicon_model_type& lexicon_target_source;

  align_map_type aligns_source_impl;
  align_map_type aligns_target_impl;

  sentence_type source;
  sentence_type target;

  std::pair<double, double> operator()(const alignment_type& alignment) const
  {
    const size_t source_size = source.size();
    const size_t target_size = target.size();

    double lex_source_target = 1.0;
    double lex_target_source = 1.0;

    align_map_type& aligns_source = const_cast<align_map_type&>(aligns_source_impl);
    align_map_type& aligns_target = const_cast<align_map_type&>(aligns_target_impl);

    aligns_source.clear();
    aligns_target.clear();

    aligns_source.resize(source_size);
    aligns_target.resize(target_size);

    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      aligns_source[aiter->source].push_back(aiter->target);
      aligns_target[aiter->target].push_back(aiter->source);
    }

    for (size_t trg = 0; trg != target_size; ++ trg) 
      if (target[trg].is_terminal()) {
	if (aligns_target[trg].empty())
	  lex_source_target *= lexicon_source_target(lexicon_model_type::vocab_type::EPSILON, target[trg]);
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_target[trg].end();
	  for (align_set_type::const_iterator aiter = aligns_target[trg].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_source_target(source[*aiter], target[trg]);

	  lex_source_target *= score / aligns_target[trg].size();
	}
      }

    for (size_t src = 0; src != source_size; ++ src) 
      if (source[src].is_terminal()) {
	if (aligns_source[src].empty())
	  lex_target_source *= lexicon_target_source(lexicon_model_type::vocab_type::EPSILON, source[src]);
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_source[src].end();
	  for (align_set_type::const_iterator aiter = aligns_source[src].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_target_source(target[*aiter], source[src]);

	  lex_target_source *= score / aligns_source[src].size();
	}
      }

    return std::make_pair(lex_source_target, lex_target_source);
  }
};

struct LexiconPhrase : public LexiconBase
{

  LexiconPhrase(const lexicon_model_type& __lexicon_source_target,
		const lexicon_model_type& __lexicon_target_source)
    : LexiconBase(__lexicon_source_target, __lexicon_target_source) {}

  void assign_source(const std::string& phrase)
  {
    source.assign(phrase);
  }

  void assign_target(const std::string& phrase)
  {
    target.assign(phrase);
  }

};

struct LexiconSCFG : public LexiconBase
{

  LexiconSCFG(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : LexiconBase(__lexicon_source_target, __lexicon_target_source) {}

  void assign_source(const std::string& phrase)
  {
    source.assign(phrase);
    source.erase(source.begin());
  }

  void assign_target(const std::string& phrase)
  {
    target.assign(phrase);
    target.erase(target.begin());
  }
};

struct LexiconGHKM : public LexiconBase
{

  LexiconGHKM(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : LexiconBase(__lexicon_source_target, __lexicon_target_source) {}

  void assign_source(const utils::piece& phrase)
  {
    source.clear();
    cicada::TreeRule(phrase).frontier(std::back_inserter(source));
  }

  void assign_target(const utils::piece& phrase)
  {
    target.clear();
    cicada::TreeRule(phrase).frontier(std::back_inserter(target));
  }
};


struct RootCountParser
{
  typedef RootCount root_count_type;

  typedef root_count_type::label_type  label_type;
  typedef root_count_type::counts_type counts_type;

  RootCountParser() : grammar() {}
  RootCountParser(const RootCountParser& x) : grammar() {}

  class double_base64_type : public std::string
  {
  public:
    operator double() const { return utils::decode_base64<double>(static_cast<const std::string&>(*this)); }
  };

  template <typename Iterator>
  struct root_count_parser : boost::spirit::qi::grammar<Iterator, root_count_type(), boost::spirit::standard::space_type>
  {
    root_count_parser() : root_count_parser::base_type(root_count)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      label %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];

      token %= qi::lexeme[+(standard::char_ - standard::space)];
      count_base64 %= token;
      count %= 'B' >> count_base64 | qi::double_;

      counts %= +count;
      root_count %= label >> "|||" >> counts >> "|||" >> count >> count;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> label;
    boost::spirit::qi::rule<Iterator, double_base64_type(), boost::spirit::standard::space_type> token;
    boost::spirit::qi::rule<Iterator, double(), boost::spirit::standard::space_type> count_base64;
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
    utils::encode_base64(root_count.observed_joint, std::ostream_iterator<char>(os));
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

  class double_base64_type : public std::string
  {
  public:
    operator double() const { return utils::decode_base64<double>(static_cast<const std::string&>(*this)); }
  };

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
      
      token %= qi::lexeme[+(standard::char_ - standard::space)];
      count_base64 %= token;

      counts %= +('B' >> count_base64 | qi::double_);
      phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> alignment >> "|||" >> counts;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;

    boost::spirit::qi::rule<Iterator, alignment_type::point_type(), boost::spirit::standard::space_type> point;
    boost::spirit::qi::rule<Iterator, alignment_type(), boost::spirit::standard::space_type> alignment;

    boost::spirit::qi::rule<Iterator, double_base64_type(), boost::spirit::standard::space_type> token;
    boost::spirit::qi::rule<Iterator, double(), boost::spirit::standard::space_type> count_base64;
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


struct PhrasePairModifiedParser
{
  typedef PhrasePairModified phrase_pair_type;

  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;

  PhrasePairModifiedParser() : grammar() {}
  PhrasePairModifiedParser(const PhrasePairModifiedParser& x) : grammar() {}

  class double_base64_type : public std::string
  {
  public:
    operator double() const { return utils::decode_base64<double>(static_cast<const std::string&>(*this)); }
  };

  template <typename Iterator>
  struct phrase_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_type(), boost::spirit::standard::space_type>
  {
    phrase_pair_parser() : phrase_pair_parser::base_type(phrase_pair)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      phrase %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];

      token %= qi::lexeme[+(standard::char_ - standard::space)];
      count_base64 %= token;

      counts %= +('B' >> count_base64 | qi::double_);
      phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> counts;
    }

    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, double_base64_type(), boost::spirit::standard::space_type> token;
    boost::spirit::qi::rule<Iterator, double(), boost::spirit::standard::space_type> count_base64;
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

struct PhrasePairModifiedGenerator
{
  typedef PhrasePairModified phrase_pair_type;

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


// modify counts... simply map from source to target, and collect counts
struct PhrasePairModify
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;

  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef RootCount          root_count_type;
  typedef PhrasePair         phrase_pair_type;
  typedef PhrasePairModified modified_type;

  typedef utils::chunk_vector<modified_type, 4096 / sizeof(modified_type), std::allocator<modified_type> >  modified_set_type;

  typedef utils::lockfree_list_queue<modified_type, std::allocator<modified_type> > queue_type;

  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;


};

struct PhrasePairModifyMapper
{
  typedef PhrasePairModify map_reduce_type;

  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;

  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::phrase_pair_type phrase_pair_type;

  typedef map_reduce_type::modified_type     modified_type;
  typedef map_reduce_type::modified_set_type modified_set_type;

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

  PhrasePairModifyMapper(const path_set_type& __paths,
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
	counts.push_back(modified_type(phrase_pair.source, phrase_pair.target, phrase_pair.counts));
      else if (counts.back().target != phrase_pair.target) {
	phrase_pair.source = counts.back().source;
	counts.push_back(modified_type(phrase_pair.source, phrase_pair.target, phrase_pair.counts));
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

    if (non_found_iter >= 16) {
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

    typedef std::deque<modified_type, std::allocator<modified_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;

    typedef std::vector<modified_set_type, std::allocator<modified_set_type> > modified_map_type;

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
    
    modified_type counts;
    
    while (! pqueue.empty()) {
      buffer_stream_type* buffer_stream(pqueue.top());
      pqueue.pop();
      
      modified_type& curr = buffer_stream->first.front();
      
      if (counts != curr) {
	if (! counts.counts.empty()) {
	  // swap source and target!
	  counts.source.swap(counts.target);
	  
	  const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
	  queues[shard]->push_swap(counts);
	}
	
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

struct PhrasePairModifyReducer
{
  typedef PhrasePairModify map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::modified_type     modified_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef PhrasePairModifiedParser    modified_parser_type;
  typedef PhrasePairModifiedGenerator modified_generator_type;


#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<modified_type, boost::hash<modified_type>, std::equal_to<modified_type>,
				  std::allocator<modified_type> > modified_unique_type;
#else
typedef sgi::hash_set<modified_type, boost::hash<modified_type>, std::equal_to<modified_type>,
		      std::allocator<modified_type> > modified_unique_type;
#endif
  
  modified_parser_type    parser;
  modified_generator_type generator;

  queue_type&    queue;
  path_type      prefix;
  path_set_type& paths;
  
  int            shard_size;
  double         max_malloc;
  int            debug;
  
  PhrasePairModifyReducer(queue_type&    __queue,
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
    modified_type modified1;
    modified_type modified2;
    
    bool parsed1 = parser(is1, modified1);
    bool parsed2 = parser(is2, modified2);
    
    while (parsed1 && parsed2) {
      if (modified1 < modified2) {
	generator(os, modified1) << '\n';
	parsed1 = parser(is1, modified1);
      } else if (modified2 < modified1) {
	generator(os, modified2) << '\n';
	parsed2 = parser(is2, modified2);
      } else {
	modified1.increment(modified2.counts.begin(), modified2.counts.end());
	generator(os, modified1) << '\n';
	
	parsed1 = parser(is1, modified1);
	parsed2 = parser(is2, modified2);
      }
    }
    
    // dump remaining...
    while (parsed1) {
      generator(os, modified1) << '\n';
      parsed1 = parser(is1, modified1);
    }
    
    while (parsed2) {
      generator(os, modified2) << '\n';
      parsed2 = parser(is2, modified2);
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
      
      const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.modified.XXXXXX");
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
  
  void dump_counts(path_set_type& paths, const modified_unique_type& counts)
  {
    // sort...!
    typedef std::vector<const modified_type*, std::allocator<const modified_type*> > sorted_type;
    
    // sorting...
    sorted_type sorted(counts.size());
    {
      sorted_type::iterator siter = sorted.begin();
      modified_unique_type::const_iterator citer_end = counts.end();
      for (modified_unique_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    std::sort(sorted.begin(), sorted.end(), less_ptr<modified_type>());
    
    // tempfile...
    const path_type counts_file_tmp = utils::tempfile::file_name(prefix / "cicada.extract.modified.XXXXXX");
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
    modified_type    modified;
    modified_unique_type counts;

    int num_termination = 0;
    
    size_type min_counts_size = 0;
    const size_type iteration_mask = (1 << 5) - 1;
    const size_type malloc_threshold = size_type(max_malloc * 1024 * 1024 * 1024);
    
    for (size_type iteration = 0; /**/; ++ iteration) {
      modified.clear();
      queue.pop_swap(modified);
      
      if (modified.source.empty()) {
	++ num_termination;
	
	if (num_termination == shard_size)
	  break;
	else
	  continue;
      }
      
      std::pair<modified_unique_type::iterator, bool> result = counts.insert(modified);
      if (! result.second)
	const_cast<modified_type&>(*result.first).increment(modified.counts.begin(), modified.counts.end());
      
      if (((iteration & iteration_mask) == iteration_mask)
	  && ! counts.empty()
	  && (! min_counts_size || counts.size() > min_counts_size)
	  && (utils::malloc_stats::used() > malloc_threshold)) {
	if (! min_counts_size)
	  min_counts_size = counts.size() >> 1;
	    
	dump_counts(paths, counts);
	counts.clear();
	modified_unique_type(counts).swap(counts);
      }
    }
    
    if (! counts.empty()) {
      dump_counts(paths, counts);
      counts.clear();
      modified_unique_type(counts).swap(counts);
    }
    
    merge_counts(paths);
  }
};

struct PhrasePairReverse
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef double count_type;
  
  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;
  
  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
  
  typedef PhrasePair         phrase_pair_type;
  typedef RootCount          root_count_type;
  typedef PhrasePairModified modified_type;
  
  typedef utils::chunk_vector<modified_type, 4096 / sizeof(modified_type), std::allocator<modified_type> > modified_set_type;

  typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
  
  typedef utils::lockfree_list_queue<modified_type, std::allocator<modified_type> > queue_type;
  
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;
};

template <typename ExtractRoot>
struct PhrasePairReverseMapper
{
  typedef PhrasePairReverse map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;
  
  typedef map_reduce_type::count_type  count_type;
  
  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::modified_type    modified_type;

  typedef map_reduce_type::modified_set_type modified_set_type;
    
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;

  typedef PhrasePairModifiedParser    modified_parser_type;
  typedef PhrasePairModifiedGenerator modified_generator_type;
  
  typedef ExtractRoot extract_root_type;
  
  hasher_type hasher;
  
  extract_root_type extract_root;
  
  path_set_type paths;
  queue_ptr_set_type& queues;
  root_count_set_type& root_counts;
  
  double max_malloc;
  int    debug;
  
  PhrasePairReverseMapper(const path_set_type& __paths,
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

  modified_parser_type    parser;
  modified_generator_type generator;

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    modified_type phrase_pair;
    
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
    
    typedef std::deque<modified_type, std::allocator<modified_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;

    typedef std::vector<modified_set_type, std::allocator<modified_set_type> > modified_map_type;
    
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
    
    modified_set_type counts;
    modified_type     modified;
    count_type        observed(0);
    
    root_count_set_type::iterator riter;
    
    while (! pqueue.empty()) {
      buffer_stream_type* buffer_stream(pqueue.top());
      pqueue.pop();
      
      modified_type& curr = buffer_stream->first.front();
      
      if (curr.source != modified.source) {
	
	if (! counts.empty()) {
	  // dump counts... but we use the counts from modified and additional observed...
	  // this observed is the same as counts.size()!
	  
	  modified.counts.push_back(observed);
	  
	  modified_set_type::iterator citer_end = counts.end();
	  for (modified_set_type::iterator citer = counts.begin(); citer != citer_end; ++ citer) {
	    citer->source.swap(citer->target);
	    citer->counts = modified.counts;
	    
	    const int shard = hasher(citer->source.begin(), citer->source.end(), 0) % queues.size();
	    queues[shard]->push_swap(*citer);
	  }
	  
	  counts.clear();
	}
	
	counts.push_back(curr);
	modified.swap(curr);
	observed = 1;
	
	// increment root_observed(lhs)
	riter = root_counts.insert(extract_root(modified.source)).first;
	const_cast<root_count_type&>(*riter).increment(modified.counts.begin(), modified.counts.end());
	const_cast<root_count_type&>(*riter).observed_joint += 1;
	const_cast<root_count_type&>(*riter).observed += 1;
	
      } else if (curr.target != modified.target) {
	// increment root_observed(lhs, rhs)
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*riter).observed_joint += 1;
	
	counts.push_back(curr);
	modified.target.swap(curr.target);
	modified.increment(curr.counts.begin(), curr.counts.end());
	observed += 1;
      } else {
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	modified.increment(curr.counts.begin(), curr.counts.end());
      }
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! counts.empty()) {
      modified.counts.push_back(observed);
      
      modified_set_type::iterator citer_end = counts.end();
      for (modified_set_type::iterator citer = counts.begin(); citer != citer_end; ++ citer) {
	citer->source.swap(citer->target);
	citer->counts = modified.counts;
	
	const int shard = hasher(citer->source.begin(), citer->source.end(), 0) % queues.size();
	queues[shard]->push_swap(*citer);
      }
      
      counts.clear();
    }

    std::vector<bool, std::allocator<bool> > terminated(queues.size(), false);
    modified.clear();
    
    for (;;) {
      bool found = false;
      
      if (! terminated[shard] && queues[shard]->push_swap(modified, true)) {
	modified.clear();
	
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
    
    if (non_found_iter >= 16) {
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    
      non_found_iter = 0;
    }
    return non_found_iter;
  }
};


struct PhrasePairReverseReducer
{
  typedef PhrasePairReverse map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;
  
  typedef map_reduce_type::count_type  count_type;
  
  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::modified_type    modified_type;

  typedef map_reduce_type::modified_set_type modified_set_type;
    
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;

  typedef PhrasePairModifiedParser    modified_parser_type;
  typedef PhrasePairModifiedGenerator modified_generator_type;

  modified_parser_type    parser;
  modified_generator_type generator;
  
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
    modified_type modified1;
    modified_type modified2;
    
    bool parsed1 = parser(is1, modified1);
    bool parsed2 = parser(is2, modified2);
    
    while (parsed1 && parsed2) {
      if (modified1 < modified2) {
	generator(os, modified1) << '\n';
	parsed1 = parser(is1, modified1);
      } else if (modified2 < modified1) {
	generator(os, modified2) << '\n';
	parsed2 = parser(is2, modified2);
      } else {
	modified1.increment(modified2.counts.begin(), modified2.counts.end());
	generator(os, modified1) << '\n';
	
	parsed1 = parser(is1, modified1);
	parsed2 = parser(is2, modified2);
      }
    }
    
    // dump remaining...
    while (parsed1) {
      generator(os, modified1) << '\n';
      parsed1 = parser(is1, modified1);
    }
    
    while (parsed2) {
      generator(os, modified2) << '\n';
      parsed2 = parser(is2, modified2);
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

  void dump_counts(path_set_type& paths, const modified_set_type& counts)
  {
    // sort...!
    typedef std::vector<const modified_type*, std::allocator<const modified_type*> > sorted_type;
    
    // sorting...
    sorted_type sorted(counts.size());
    {
      sorted_type::iterator siter = sorted.begin();
      modified_set_type::const_iterator citer_end = counts.end();
      for (modified_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    std::sort(sorted.begin(), sorted.end(), less_ptr<modified_type>());
    
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
    // we already know that the entries are uniqued through the previous map-reduce!
    
    modified_type modified;
    modified_set_type counts;
    
    int num_termination = 0;
    
    size_type min_counts_size = 0;
    const size_type iteration_mask = (1 << 5) - 1;
    const size_type malloc_threshold = size_type(max_malloc * 1024 * 1024 * 1024);
    
    for (size_type iteration = 0; /**/; ++ iteration) {
      modified.clear();
      queue.pop_swap(modified);
      
      if (modified.source.empty()) {
	++ num_termination;
	
	if (num_termination == shard_size)
	  break;
	else
	  continue;
      }
      
      counts.push_back(modified);
      
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

  typedef PhrasePair         phrase_pair_type;
  typedef PhrasePairModified modified_type;
  typedef RootCount          root_count_type;
  
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
  typedef map_reduce_type::modified_type    modified_type;
    
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
  
    if (non_found_iter >= 16) {
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
    
    while (1) {
      bool found = false;
      
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (! terminated[shard]) {
	  counts.clear();
	  
	  if (queues[shard]->push_swap(counts, true)) {
	    counts.clear();
	    terminated[shard] = true;
	    found = true;
	  }
	}
      
      if (std::count(terminated.begin(), terminated.end(), true) == static_cast<int>(terminated.size())) break;
      
      non_found_iter = loop_sleep(found, non_found_iter);
    }
  }
};

template <typename ExtractRoot, typename Lexicon>
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
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::modified_type    modified_type;
  
  typedef std::vector<phrase_pair_type, std::allocator<phrase_pair_type> > phrase_pair_set_type;
  typedef std::vector<modified_type, std::allocator<modified_type> > modified_set_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;  

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;

  typedef PhrasePairModifiedParser    modified_parser_type;
  
  typedef ExtractRoot extract_root_type;
  typedef Lexicon     lexicon_type;
  
  root_count_set_type& root_counts;
  
  extract_root_type extract_root;
  lexicon_type      lexicon;
  
  const path_set_type& paths;
  queue_ptr_set_type& queues;
  std::ostream& os;
  int debug;
  
  PhrasePairScoreReducer(root_count_set_type& __root_counts,
			 const extract_root_type& __extract_root,
			 const lexicon_type& __lexicon,
			 const path_set_type& __paths,
			 queue_ptr_set_type& __queues,
			 std::ostream& __os,
			 int __debug)
    : root_counts(__root_counts),
      extract_root(__extract_root),
      lexicon(__lexicon),
      paths(__paths),
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

  struct real_precision : boost::spirit::karma::real_policies<count_type>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };
  
  boost::spirit::karma::real_generator<count_type, real_precision> double20;
  
  
  void dump_phrase_pair(const phrase_pair_set_type& counts, const modified_set_type& modified)
  {
    typedef std::vector<phrase_pair_set_type::const_iterator, std::allocator<phrase_pair_set_type::const_iterator> > iterators_type;

    // counts are grouped by source
    
    //
    // compute source counts and merged counts (ignoring difference of word alignment)
    //
    phrase_pair_type counts_source;
    count_type       observed_source(0);
    
    iterators_type iterators;
    
    phrase_pair_set_type::const_iterator citer_end = counts.end();
    for (phrase_pair_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      counts_source.increment(citer->counts.begin(), citer->counts.end());
      
      if (counts_source.target != citer->target) {
	iterators.push_back(citer);
	
	observed_source += 1;
	counts_source.target = citer->target;
      }
    }
    iterators.push_back(citer_end);

    if (modified.size() != iterators.size() - 1)
      throw std::runtime_error("source/target size mismatch?");
    
    phrase_pair_type counts_pair;

    lexicon.assign_source(counts.front().source);
    
    modified_set_type::const_iterator miter = modified.begin();
    
    iterators_type::const_iterator iiter_end = iterators.end() - 1;
    for (iterators_type::const_iterator iiter = iterators.begin(); iiter != iiter_end; ++ iiter, ++ miter) {
      phrase_pair_set_type::const_iterator first = *iiter;
      phrase_pair_set_type::const_iterator last  = *(iiter + 1);
      
      // compute lexical weights from [first, last)... we will take "max"
      
      if (first->source != miter->source)
	throw std::runtime_error("different source?") ;
      if (first->target != miter->target)
	throw std::runtime_error("different source?") ;
      
      counts_pair.clear();
      
      lexicon.assign_target(first->target);
      
      double lex_source_target = 0.0;
      double lex_target_source = 0.0;
      
      for (phrase_pair_set_type::const_iterator iter = first; iter != last; ++ iter) {
	const std::pair<double, double> lex = lexicon(iter->alignment);
	
	lex_source_target = std::max(lex_source_target, lex.first);
	lex_target_source = std::max(lex_target_source, lex.second);
	
	counts_pair.increment(iter->counts.begin(), iter->counts.end());
      }
      
      const modified_type::counts_type& counts_target = miter->counts;
      
      os << first->source << " ||| " << first->target;
      
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      typedef std::ostream_iterator<char> iterator_type;

      iterator_type iter(os);
      
      // cont(LHS RHS)
      if (! karma::generate(iter, " ||| " << -(double20 % ' '), counts_pair.counts))
	throw std::runtime_error("generation failed");
      
      // count(LHS)
      if (! karma::generate(iter, " ||| " << -(double20 % ' '), counts_source.counts))
	throw std::runtime_error("generation failed");
      
      // count(RHS)
      if (! karma::generate(iter, " ||| " << -(double20 % ' '), boost::make_iterator_range(counts_target.begin(), counts_target.end() - 1)))
	throw std::runtime_error("generation failed");
      
      // observed(LHS) observed(RHS)
      if (! karma::generate(iter, " ||| " << double20 << ' ' << double20, observed_source, counts_target.back()))
	throw std::runtime_error("generation failed");
      
      // lex(rhs | lhs) lex(rhs | lhs)
      if (! karma::generate(iter, " ||| " << double20 << ' ' << double20 << '\n', lex_target_source, lex_source_target))
	throw std::runtime_error("generation failed");
    }
  }

  modified_parser_type    modified_parser;

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    modified_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! modified_parser(line, phrase_pair)) continue;
      
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
    typedef std::pair<phrase_pair_type, queue_type*> buffer_queue_type;
    typedef std::vector<buffer_queue_type*, std::allocator<buffer_queue_type*> > pqueue_base_type;
    typedef std::priority_queue<buffer_queue_type*, pqueue_base_type, greater_buffer<buffer_queue_type> > pqueue_type;
    
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<modified_type, std::allocator<modified_type> > modified_buffer_type;
    typedef std::pair<modified_buffer_type, istream_type*> buffer_stream_type;
    typedef std::vector<buffer_stream_type, std::allocator<buffer_stream_type> > buffer_stream_set_type;
    
    typedef std::vector<buffer_stream_type*, std::allocator<buffer_stream_type*> > modified_pqueue_base_type;
    typedef std::priority_queue<buffer_stream_type*, modified_pqueue_base_type, greater_buffer<buffer_stream_type> > modified_pqueue_type;
    
    os.precision(20);
    
    pqueue_type pqueue;
    std::vector<buffer_queue_type, std::allocator<buffer_queue_type> > buffer_queues(queues.size());

    modified_pqueue_type   mqueue;
    istream_ptr_set_type   istreams(paths.size());
    buffer_stream_set_type buffer_streams(paths.size());
    
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
      for (path_set_type::const_iterator piter = paths.begin(); piter != paths.end(); ++ piter, ++ pos) {
	if (! boost::filesystem::exists(*piter))
	  throw std::runtime_error("no file? " + piter->string());
	
	istreams[pos].reset(new istream_type(*piter, 1024 * 1024));
	
	buffer_stream_type* buffer_stream = &buffer_streams[pos];
	buffer_stream->second = &(*istreams[pos]);
	
	read_phrase_pair(*istreams[pos], buffer_stream->first);
	
	if (! buffer_stream->first.empty())
	  mqueue.push(buffer_stream);
      }
    }
    
    phrase_pair_set_type counts;
    modified_set_type    modified;

    root_count_set_type::iterator riter;
    
    while (! pqueue.empty()) {
      buffer_queue_type* buffer_queue(pqueue.top());
      pqueue.pop();

      phrase_pair_type& curr = buffer_queue->first;
      
      if (counts.empty() || counts.back().source != curr.source) {
	if (! counts.empty()) {
	  if (debug >= 4)
	    std::cerr << "score count reducer: " << counts.size() << std::endl;
	  
	  dump_phrase_pair(counts, modified);
	  
	  counts.clear();
	  modified.clear();
	  
	  phrase_pair_set_type(counts).swap(counts);
	  modified_set_type(modified).swap(modified);
	}
	
	// increment root_observed(lhs)
	riter = root_counts.insert(extract_root(curr.source)).first;
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*riter).observed_joint += 1;
	const_cast<root_count_type&>(*riter).observed += 1;
	
	counts.push_back(curr);
	
	if (mqueue.empty())
	  throw std::runtime_error("modified counts and phrase pair do not match");
	
	buffer_stream_type* buffer_stream(mqueue.top());
	mqueue.pop();
	
	modified.push_back(buffer_stream->first.front());
	
	buffer_stream->first.pop_front();
	
	if (buffer_stream->first.empty())
	  read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
	
	if (! buffer_stream->first.empty())
	  mqueue.push(buffer_stream);

	if (counts.back().source != modified.back().source)
	  throw std::runtime_error("source mismatch? " + counts.back().source + " modified: " + modified.back().source);
	
	if (counts.back().target != modified.back().target)
	  throw std::runtime_error("target mismatch? " + counts.back().target + " modified: " + modified.back().target);
	
      } else if (counts.back().target != curr.target) {
	// increment root_observed(lhs, rhs)
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	const_cast<root_count_type&>(*riter).observed_joint += 1;
	
	counts.push_back(curr);
	
	if (mqueue.empty())
	  throw std::runtime_error("modified counts and phrase pair do not match: sharing the same source");
	
	buffer_stream_type* buffer_stream(mqueue.top());
	mqueue.pop();
	
	modified.push_back(buffer_stream->first.front());
	
	buffer_stream->first.pop_front();
	
	if (buffer_stream->first.empty())
	  read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
	
	if (! buffer_stream->first.empty())
	  mqueue.push(buffer_stream);
	
	if (counts.back().source != modified.back().source)
	  throw std::runtime_error("source mismatch? " + counts.back().source + " modified: " + modified.back().source);
	
	if (counts.back().target != modified.back().target)
	  throw std::runtime_error("target mismatch? " + counts.back().target + " modified: " + modified.back().target);
	
      } else if (counts.back().alignment != curr.alignment) {
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	counts.push_back(curr);
      } else {
	const_cast<root_count_type&>(*riter).increment(curr.counts.begin(), curr.counts.end());
	counts.back().increment(curr.counts.begin(), curr.counts.end());
      }
      
      buffer_queue->second->pop_swap(buffer_queue->first);
      if (! buffer_queue->first.source.empty())
	pqueue.push(buffer_queue);
    }
    
    if (! mqueue.empty())
      throw std::runtime_error("queue is empty but modified-queue is not empty!");
    
    if (! counts.empty()) {
      if (debug >= 4)
	std::cerr << "score count reducer: " << counts.size() << std::endl;
      
      dump_phrase_pair(counts, modified);
    }
  }
  
};

#endif
