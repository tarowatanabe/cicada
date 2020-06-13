//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXTRACT_SCFG_IMPL__HPP__
#define __CICADA__EXTRACT_SCFG_IMPL__HPP__ 1

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>

#include <cstring>
#include <string>
#include <vector>

#include <boost/array.hpp>

#include "cicada/sentence.hpp"
#include "cicada/alignment.hpp"
#include "cicada/span_vector.hpp"
#include "cicada/vocab.hpp"

#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/compact_set.hpp"
#include "utils/chart.hpp"
#include "utils/indexed_set.hpp"

#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/chunk_vector.hpp>
#include <utils/hashmurmur3.hpp>
#include <utils/getline.hpp>

struct Bitext
{
  typedef cicada::Sentence   sentence_type;
  typedef cicada::Alignment  alignment_type;
  typedef cicada::SpanVector span_set_type;
  
  sentence_type  source;
  sentence_type  target;
  alignment_type alignment;
  span_set_type  spans_source;
  span_set_type  spans_target;
  
  Bitext() : source(), target(), alignment(), spans_source(), spans_target() {}
  Bitext(const sentence_type& __source,
	 const sentence_type& __target,
	 const alignment_type& __alignment)
    : source(__source), target(__target), alignment(__alignment), spans_source(), spans_target() {}
  
  void swap(Bitext& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    alignment.swap(x.alignment);
    spans_source.swap(x.spans_source);
    spans_target.swap(x.spans_target);
  }

  void clear()
  {
    source.clear();
    target.clear();
    alignment.clear();
    spans_source.clear();
    spans_target.clear();
  }
  
  friend
  std::istream& operator>>(std::istream& is, Bitext& bitext)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    bitext.clear();
    std::string line;
    if (utils::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end  = line.end();
      
      if((!bitext.source.assign(iter, end))
	 || (!qi::phrase_parse(iter, end, qi::lit("|||"), standard::space))
	 || (!bitext.target.assign(iter, end))
	 || (!qi::phrase_parse(iter, end, qi::lit("|||"), standard::space))
	 || (!bitext.alignment.assign(iter, end))
	 || (!qi::phrase_parse(iter, end, qi::lit("|||"), standard::space))
	 || (!bitext.spans_source.assign(iter, end))
	 || (!qi::phrase_parse(iter, end, qi::lit("|||"), standard::space))
	 || (!bitext.spans_target.assign(iter, end))
	 || iter != end)
	bitext.clear();
    }
    return is;
  }
  
  friend
  std::ostream& operator<<(std::ostream& os, const Bitext& bitext)
  {
    os << bitext.source
       << " ||| " << bitext.target
       << " ||| " << bitext.alignment
       << " ||| " << bitext.spans_source
       << " ||| " << bitext.spans_target;
    return os;
  }
  
};

struct RulePair
{
  typedef std::string phrase_type;
  typedef cicada::Alignment alignment_type;
  typedef double count_type;
  typedef int64_t frequency_type;
  typedef boost::array<frequency_type, 3> freqs_type;
  
  const phrase_type*    source;
  const phrase_type*    target;
  const alignment_type* alignment;
  count_type     count;
  freqs_type     freqs;

  RulePair() : source(0), target(0), alignment(0), count(0), freqs() {}
  RulePair(const phrase_type& __source,
	   const phrase_type& __target,
	   const alignment_type& __alignment,
	   const count_type& __count)
    : source(&__source), target(&__target), alignment(&__alignment), count(__count), freqs() {}

  friend
  size_t hash_value(RulePair const& x)
  {
    typedef utils::hashmurmur3<size_t> hasher_type;
    
    return hasher_type()(x.source, hasher_type()(x.target, (uintptr_t) x.alignment));
  }

  friend
  bool operator==(const RulePair& x, const RulePair& y) 
  {
    return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
  }
  
  friend
  bool operator!=(const RulePair& x, const RulePair& y) 
  {
    return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
  }
  
  friend
  bool operator<(const RulePair& x, const RulePair& y)
  {
    return (RulePair::compare(x.source, y.source)
	    || (!RulePair::compare(y.source, x.source)
		&& (RulePair::compare(x.target, y.target)
		    || (!RulePair::compare(y.target, x.target)
			&& RulePair::compare(x.alignment, y.alignment)))));
  }
  
  friend
  bool operator>(const RulePair& x, const RulePair& y)
  {
    return y < x;
  }
  
  template <typename Tp>
  static inline
  bool compare(const Tp* x, const Tp* y)
  {
    return (x && y && *x < *y) || (!x && y);
  }
};

BOOST_FUSION_ADAPT_STRUCT(RulePair::alignment_type::point_type,
			  (RulePair::alignment_type::index_type, source)
			  (RulePair::alignment_type::index_type, target)
			  )

struct RulePairGenerator
{
  typedef RulePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::frequency_type frequency_type;
  typedef phrase_pair_type::count_type     count_type;
  
  RulePairGenerator() {}
  
  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair) const
  {
    typedef std::ostream_iterator<char> iterator_type;
  
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    iterator_type iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string
			  << " ||| " << -((karma::int_ << '-' << karma::int_) % ' ')
			  << " ||| " << double20 << ' ' << (frequency % ' '),
			  *phrase_pair.source,
			  *phrase_pair.target,
			  *phrase_pair.alignment,
			  phrase_pair.count,
			  phrase_pair.freqs))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return 20;
    }
  };
  
  boost::spirit::karma::real_generator<double, real_precision> double20;
  boost::spirit::karma::int_generator<frequency_type> frequency;
};

namespace std
{
  inline
  void swap(Bitext& x, Bitext& y)
  {
    x.swap(y);
  }
};

struct ExtractSCFG
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::Sentence   sentence_type;
  typedef cicada::Alignment  alignment_type;
  typedef cicada::SpanVector span_set_type;
  typedef cicada::Symbol     symbol_type;
  typedef cicada::Vocab      vocab_type;

  typedef std::pair<int, int> span_type;
  
  struct span_pair_type
  {
    span_type source;
    span_type target;
    
    span_pair_type() : source(), target() {}
    span_pair_type(const span_type& __source, const span_type& __target) : source(__source), target(__target) {}

    friend
    size_t hash_value(span_pair_type const& x)
    {
      return utils::hashmurmur3<size_t>()(x, 0);
    }
    
    friend
    bool operator==(const span_pair_type& x, const span_pair_type& y)
    {
      return x.source == y.source && x.target == y.target;
    }
    
    friend
    bool operator!=(const span_pair_type& x, const span_pair_type& y)
    {
      return x.source != y.source || x.target != y.target;
    }

    friend
    bool operator<(const span_pair_type& x, const span_pair_type& y)
    {
      return x.source < y.source || (!(y.source < x.source) && x.target < y.target);
    }
    
    friend
    bool operator>(const span_pair_type& x, const span_pair_type& y)
    {
      return y < x;
    }
    
  };

  typedef RulePair rule_pair_type;
  typedef rule_pair_type::phrase_type phrase_type;
  typedef rule_pair_type::count_type  count_type;

  struct rule_compact_set_type
  {
    struct string_hash : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      size_t operator()(const std::string& x) const
      {
	return hasher_type::operator()(x.begin(), x.end(), 0);
      }
    };
    
    typedef utils::indexed_set<phrase_type, string_hash, std::equal_to<phrase_type>, 
			       std::allocator<phrase_type> > phrase_set_type;
    
    typedef phrase_set_type::index_type      index_type;
    typedef phrase_set_type::size_type       size_type;
    typedef phrase_set_type::difference_type difference_type;

    void clear()
    {
      phrases.clear();
    }

    void swap(rule_compact_set_type& x)
    {
      phrases.swap(x.phrases);
    }
    
    index_type insert(const phrase_type& x)
    {
      phrase_set_type::iterator iter = phrases.insert(x).first;
      return iter - phrases.begin();
    }
    
    const phrase_type& operator[](index_type x) const { return phrases[x]; }
    
    size_type size() const { return phrases.size(); }
    bool empty() const { return phrases.empty(); }
    
    phrase_set_type phrases;
  };

  struct alignment_set_type
  {
    typedef rule_pair_type::alignment_type alignment_type;

    typedef utils::indexed_set<alignment_type,
			       boost::hash<alignment_type>,
			       std::equal_to<alignment_type>,
			       std::allocator<alignment_type> > align_set_type;
    
    typedef align_set_type::index_type      index_type;
    typedef align_set_type::size_type       size_type;
    typedef align_set_type::difference_type difference_type;

    void clear()
    {
      aligns.clear();
    }
    
    void swap(alignment_set_type& x)
    {
      aligns.swap(x.aligns);
    }
    
    index_type insert(const alignment_type& x)
    {
      align_set_type::iterator iter = aligns.insert(x).first;
      return iter - aligns.begin();
    }
    
    const alignment_type& operator[](index_type x) const { return aligns[x]; }
    
    size_type size() const { return aligns.size(); }
    bool empty() const { return aligns.empty(); }
    
    align_set_type aligns;
  };

  struct RulePairCompact
  {
    rule_compact_set_type::index_type source;
    rule_compact_set_type::index_type target;
    alignment_set_type::index_type    alignment;
    count_type                        count;
    
    RulePairCompact()
      : source(rule_compact_set_type::index_type(-1)),
	target(rule_compact_set_type::index_type(-1)),
	alignment(alignment_set_type::index_type(-1)),
	count(0) {}

    friend
    size_t hash_value(RulePairCompact const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.target, hasher_type()(x.alignment, x.source));
    }
    
    friend
    bool operator==(const RulePairCompact& x, const RulePairCompact& y) 
    {
      return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
    }
    
    friend
    bool operator!=(const RulePairCompact& x, const RulePairCompact& y)
    {
      return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
    }
  };
  
  typedef RulePairCompact rule_pair_compact_type;

  struct rule_pair_compact_unassigned
  {
    rule_pair_compact_type operator()() const
    {
      return rule_pair_compact_type();
    }
  };
  
  typedef utils::compact_set<rule_pair_compact_type,
			     rule_pair_compact_unassigned, rule_pair_compact_unassigned,
			     boost::hash<rule_pair_compact_type>, std::equal_to<rule_pair_compact_type>,
			     std::allocator<rule_pair_compact_type> > rule_pair_compact_set_type;


  struct rule_pair_set_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef rule_pair_type::alignment_type alignment_type;
    
    struct string_hash : public utils::hashmurmur3<size_t>
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      size_t operator()(const std::string& x) const
      {
	return hasher_type::operator()(x.begin(), x.end(), 0);
      }
    };
    
    typedef utils::indexed_set<phrase_type, string_hash, std::equal_to<phrase_type>, 
			       std::allocator<phrase_type> > phrase_set_type;

    typedef utils::indexed_set<alignment_type,
			       boost::hash<alignment_type>,
			       std::equal_to<alignment_type>,
			       std::allocator<alignment_type> > alignment_set_type;

    struct rule_pair_unassigned
    {
      rule_pair_type operator()() const
      {
	return rule_pair_type();
      }
    };
    
    typedef utils::compact_set<rule_pair_type,
			       rule_pair_unassigned, rule_pair_unassigned,
			       boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
			       std::allocator<rule_pair_type> > pair_set_type;
    
    typedef pair_set_type::const_iterator const_iterator;
    typedef pair_set_type::iterator       iterator;
    
    bool empty() const { return rules.empty(); }
    size_type size() const { return rules.size(); }

    const_iterator begin() const { return rules.begin(); }
    const_iterator end() const { return rules.end(); }

    std::pair<iterator, bool> insert(const rule_pair_type& x)
    {
      return rules.insert(rule_pair_type(*sources.insert(*x.source).first,
					 *targets.insert(*x.target).first,
					 *alignments.insert(*x.alignment).first,
					 x.count));
    }
    
    void clear()
    {
      rules.clear();
      sources.clear();
      targets.clear();
      alignments.clear();
    }

    void swap(rule_pair_set_type& x)
    {
      rules.swap(x.rules);
      sources.swap(x.sources);
      targets.swap(x.targets);
      alignments.swap(x.alignments);
    }
    
    pair_set_type        rules;
    phrase_set_type      sources;
    phrase_set_type      targets;
    alignment_set_type   alignments;
  };  

  typedef std::vector<bool, std::allocator<bool> > unique_set_type;
  
  typedef std::pair<rule_compact_set_type::index_type, rule_compact_set_type::index_type> unique_pair_type;
  
  struct unique_pair_unassigned
  {
    unique_pair_type operator()() const
    {
      return unique_pair_type(rule_compact_set_type::index_type(-1),
			      rule_compact_set_type::index_type(-1));
    }
  };
  
  typedef utils::compact_set<unique_pair_type,
			     unique_pair_unassigned, unique_pair_unassigned,
			     utils::hashmurmur3<size_t>, std::equal_to<unique_pair_type>,
			     std::allocator<unique_pair_type> > unique_pair_set_type;
  
  typedef utils::chart<span_type, std::allocator<span_type> >          span_chart_type;
  typedef std::vector<int, std::allocator<int> >                       alignment_count_set_type;
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  
  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > alignment_multiple_type;
  
  typedef utils::chart<phrase_type, std::allocator<phrase_type> >      phrase_chart_type;
  
  typedef std::pair<span_type, span_type> span_parent_type;
  
  ExtractSCFG(const int __max_length,
	      const int __max_fertility,
	      const int __max_span_source,
	      const int __max_span_target,
	      const int __min_hole_source,
	      const int __min_hole_target,
	      const int __max_rank,
	      const int __max_scope,
	      const bool __exhaustive,
	      const bool __constrained,
	      const bool __exclude,
	      const bool __sentential,
	      const bool __inverse)
    : max_length(__max_length),
      max_fertility(__max_fertility),
      max_span_source(__max_span_source),
      max_span_target(__max_span_target),
      min_hole_source(__min_hole_source),
      min_hole_target(__min_hole_target),
      max_rank(__max_rank),
      max_scope(__max_scope),
      exhaustive(__exhaustive),
      constrained(__constrained),
      exclude(__exclude),
      sentential(__sentential),
      inverse(__inverse) {}
  
  int max_length;
  int max_fertility;
  int max_span_source;
  int max_span_target;
  int min_hole_source;
  int min_hole_target;
  int max_rank;
  int max_scope;
  bool exhaustive;
  bool constrained;
  bool exclude;
  bool sentential;
  bool inverse;
  
  alignment_multiple_type alignment_source_target;
  alignment_multiple_type alignment_target_source;
  
  span_chart_type span_source_chart;
  span_chart_type span_target_chart;
  
  alignment_count_set_type alignment_count_source;
  alignment_count_set_type alignment_count_target;
  
  span_pair_set_type        spans;
  span_pair_set_type        spans_unique;

  struct ExtractCategory
  {
    symbol_type operator()(const span_pair_type& spans) const
    {
      return vocab_type::X;
    }

    symbol_type operator()(const std::string& ) const
    {
      return vocab_type::X;
    }
  };
  
  struct ExtractCategoryBase
  {
    typedef utils::chart<symbol_type, std::allocator<symbol_type> > label_chart_type;
    
    typedef utils::unordered_map<span_type, symbol_type, utils::hashmurmur3<size_t>, std::equal_to<span_type>,
				 std::allocator<std::pair<const span_type, symbol_type> > >::type label_map_type;
    
    label_map_type   label_map;
    label_chart_type label_chart;

    ExtractCategoryBase(const sentence_type& sentence, const span_set_type& spans)
      : label_map(), label_chart(sentence.size() + 1, vocab_type::EMPTY)
    {
      span_set_type::const_iterator siter_end = spans.end();
      for (span_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)
	if (! siter->label.empty()) {
	  label_map[std::make_pair(siter->first, siter->last)] = siter->label;
	  label_chart(siter->first, siter->last) = siter->label;
	}
    }

    std::string operator()(const std::string& x) const
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
    
      std::string::const_iterator iter = x.begin();
      std::string::const_iterator end  = x.end();

      std::string label;
      
      const bool result = qi::phrase_parse(iter, end,
					   qi::lexeme[+(standard::char_ - standard::space)],
					   standard::space,
					   label);
      if (! result)
	throw std::runtime_error("no label?");
      
      return label;
    }
    
    symbol_type extract(const span_type& span) const
    {
      symbol_type& label = const_cast<symbol_type&>(label_chart(span.first, span.second));
      
      if (label != vocab_type::EMPTY)
	return label;
      
      // note that we have already inserted exact match into label-chart...
      
      // try binary combination...
      for (int middle = span.first + 1; middle < span.second; ++ middle) {
	label_map_type::const_iterator piter = label_map.find(std::make_pair(span.first, middle));
	label_map_type::const_iterator niter = label_map.find(std::make_pair(middle, span.second));
	
	if (piter != label_map.end() && niter != label_map.end()) {
	  label = '[' + piter->second.non_terminal_strip() + '+' + niter->second.non_terminal_strip() + ']';
	  return label;
	}
      }
      
      // try right-subsitution...
      for (int last_super = span.second + 1; last_super < static_cast<int>(label_chart.size()); ++ last_super) {
	label_map_type::const_iterator siter = label_map.find(std::make_pair(span.first,  last_super));
	label_map_type::const_iterator riter = label_map.find(std::make_pair(span.second, last_super));
	
	if (siter != label_map.end() && riter != label_map.end()) {
	  label = '[' + siter->second.non_terminal_strip() + '/' + riter->second.non_terminal_strip() + ']';
	  return label;
	}
      }
      
      // try left-substitution...
      for (int first_super = span.first - 1; first_super >= 0; -- first_super) {
	label_map_type::const_iterator siter = label_map.find(std::make_pair(first_super, span.second));
	label_map_type::const_iterator liter = label_map.find(std::make_pair(first_super, span.first));
	
	if (siter != label_map.end() && liter != label_map.end()) { 
	  label = '[' + siter->second.non_terminal_strip() + '\\' + liter->second.non_terminal_strip() + ']';
	  return label;
	}
      }
      
      // try tripple combination...
      for (int middle1 = span.first + 1; middle1 < span.second; ++ middle1)
	for (int middle2 = middle1 + 1; middle2 < span.second; ++ middle2) {
	  label_map_type::const_iterator iter1 = label_map.find(std::make_pair(span.first, middle1));
	  label_map_type::const_iterator iter2 = label_map.find(std::make_pair(middle1, middle2));
	  label_map_type::const_iterator iter3 = label_map.find(std::make_pair(middle2, span.second));
	  
	  if (iter1 != label_map.end() && iter2 != label_map.end() && iter3 != label_map.end()) {
	    label = '[' + iter1->second.non_terminal_strip() + '+' + iter2->second.non_terminal_strip() + '+' + iter3->second.non_terminal_strip() + ']';
	    return label;
	  }
	}
      
      // try longest left and longest right
      {
	label_map_type::const_iterator liter = label_map.end();
	for (int last_left = span.second - 1; span.first < last_left && liter == label_map.end(); -- last_left)
	  liter = label_map.find(std::make_pair(span.first, last_left));
	
	label_map_type::const_iterator riter = label_map.end();
	for (int first_right = span.first + 1; first_right < span.second && riter == label_map.end(); ++ first_right)
	  riter = label_map.find(std::make_pair(first_right, span.second));
	
	if (liter != label_map.end() && riter != label_map.end()) {
	  label = '[' + liter->second.non_terminal_strip() + ".." + riter->second.non_terminal_strip() + ']';
	  return label;
	}
      }

      // fallback...
      return vocab_type::X;
    }
  };

  struct ExtractCategorySource : public ExtractCategoryBase
  {
    ExtractCategorySource(const sentence_type& sentence, const span_set_type& spans)
      : ExtractCategoryBase(sentence, spans) {}

    symbol_type operator()(const span_pair_type& spans) const
    {
      return extract(spans.source);
    }

    std::string operator()(const std::string& x) const
    {
      return ExtractCategoryBase::operator()(x);
    }
  };
  
  struct ExtractCategoryTarget : public ExtractCategoryBase
  {
    ExtractCategoryTarget(const sentence_type& sentence, const span_set_type& spans)
      : ExtractCategoryBase(sentence, spans) {}
    
    symbol_type operator()(const span_pair_type& spans) const
    {
      return extract(spans.target);
    }
    
    std::string operator()(const std::string& x) const
    {
      return ExtractCategoryBase::operator()(x);
    }
  };

  rule_pair_compact_set_type rule_pairs_local;
  rule_compact_set_type      rules_source;
  rule_compact_set_type      rules_target;
  alignment_set_type         rules_alignment;
  
  unique_set_type            uniques_source;
  unique_set_type            uniques_target;
  unique_pair_set_type       uniques_pair;
  
  template <typename Dumper>
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  const alignment_type& alignment,
		  const span_set_type& spans_source,
		  const span_set_type& spans_target,
		  rule_pair_set_type& rule_pairs,
		  const Dumper& dumper)
  {
    compute_spans(source, target, alignment);

    if (sentential) {
      if (! spans_source.empty())
	extract_sentential_rules(source, target, alignment, rule_pairs, ExtractCategorySource(source, spans_source), dumper);
      else if (! spans_target.empty())
	extract_sentential_rules(source, target, alignment, rule_pairs, ExtractCategoryTarget(target, spans_target), dumper);
      else
	extract_sentential_rules(source, target, alignment, rule_pairs, ExtractCategory(), dumper);
    } else {
      if (! spans_source.empty())
	extract_rules(source, target, alignment, rule_pairs, ExtractCategorySource(source, spans_source), dumper);
      else if (! spans_target.empty())
	extract_rules(source, target, alignment, rule_pairs, ExtractCategoryTarget(target, spans_target), dumper);
      else
	extract_rules(source, target, alignment, rule_pairs, ExtractCategory(), dumper);
    }
  }
  
  template <typename Category, typename Dumper>
  void extract_sentential_rules(const sentence_type& source,
				const sentence_type& target,
				const alignment_type& alignment,
				rule_pair_set_type& rule_pairs,
				const Category& category,
				const Dumper& dumper)
  {
    typedef utils::chunk_vector<rule_pair_compact_type, 4096 / sizeof(rule_pair_compact_type), std::allocator<rule_pair_compact_type> > rule_pair_compact_list_type;
    
    const int source_length = source.size();
    const int target_length = target.size();

    const span_pair_type span(span_type(0, source_length), span_type(0, target_length));
    
    rule_pair_compact_list_type rule_pair_list;
    rule_pair_compact_type      rule_pair;
    
    rule_pairs_local.clear();
    rules_source.clear();
    rules_target.clear();
    rules_alignment.clear();
    
    if (! exclude) // phrasal rule
      if (max_length <= 0 || (source_length <= max_length && target_length <= max_length))
	if (max_fertility <= 0 || fertility(source_length, target_length) < max_fertility) {
	  extract_rule(source, target, span, category, rule_pair, true);
	  
	  rule_pair_type& rule_pair_new = const_cast<rule_pair_type&>(*(rule_pairs.insert(rule_pair_type(rules_source[rule_pair.source],
													 rules_target[rule_pair.target],
													 rules_alignment[rule_pair.alignment],
													 rule_pair.count)).first));
	  
	  rule_pair_new.count += 1;
	  rule_pair_new.freqs[0] += 1;
	  rule_pair_new.freqs[1] += 1;
	  rule_pair_new.freqs[2] += 1;
	}
    
    const int source_count = alignment_count_source[span.source.second] - alignment_count_source[span.source.first];
    //const int target_count = alignment_count_target[span.target.second] - alignment_count_target[span.target.first];
    
    span_pair_set_type::const_iterator niter_begin = (exhaustive ? spans.begin() : spans_unique.begin());
    span_pair_set_type::const_iterator niter_end   = (exhaustive ? spans.end()   : spans_unique.end());

    // first non-terminal...
    if (max_rank >= 1)
      for (span_pair_set_type::const_iterator niter1 = niter_begin; niter1 != niter_end; ++ niter1) 
	if (span != *niter1
	    && (min_hole_source <= 1 || (niter1->source.second - niter1->source.first) >= min_hole_source)
	    && (min_hole_target <= 1 || (niter1->target.second - niter1->target.first) >= min_hole_target)
	    && is_parent(span.source, niter1->source)
	    && is_parent(span.target, niter1->target)) {
	
	  const int source_count1 = alignment_count_source[niter1->source.second] - alignment_count_source[niter1->source.first];
	  //const int target_count1 = alignment_count_target[niter1->target.second] - alignment_count_target[niter1->target.first];
	
	  const int source_length1 = source_length - (niter1->source.second - niter1->source.first);
	  const int target_length1 = target_length - (niter1->target.second - niter1->target.first);
	
	  if (source_count1 != source_count)
	    if (max_length <= 0 || (source_length1 <= max_length && target_length1 <= max_length))
	      if (max_fertility <= 0 || fertility(source_length1, target_length1) < max_fertility) {
		extract_rule(source, target, span, *niter1, category, rule_pair, true);
		rule_pair_list.push_back(rule_pair);
	      }
	
	  // second non-terminal...
	  if (max_rank >= 2)
	    for (span_pair_set_type::const_iterator niter2 = niter1 + 1; niter2 != niter_end; ++ niter2) 
	      if (span != *niter2
		  && (min_hole_source <= 1 || (niter2->source.second - niter2->source.first) >= min_hole_source)
		  && (min_hole_target <= 1 || (niter2->target.second - niter2->target.first) >= min_hole_target)
		  && is_parent(span.source, niter2->source)
		  && is_parent(span.target, niter2->target)
		  && is_disjoint(niter1->source, niter2->source)
		  && is_disjoint(niter1->target, niter2->target)
		  && ! is_adjacent(niter1->source, niter2->source)) {
	    
		const int source_count2 = alignment_count_source[niter2->source.second] - alignment_count_source[niter2->source.first];
		//const int target_count2 = alignment_count_target[niter2->target.second] - alignment_count_target[niter2->target.first];
	    
		const int source_length2 = source_length1 - (niter2->source.second - niter2->source.first);
		const int target_length2 = target_length1 - (niter2->target.second - niter2->target.first);
	    
		if (source_count != source_count1 + source_count2)
		  if (max_length <= 0 || (source_length2 <= max_length && target_length2 <= max_length))
		    if (max_fertility <= 0 || fertility(source_length2, target_length2) < max_fertility) {
		      extract_rule(source, target, span, *niter1, *niter2, category, rule_pair, true);
		      rule_pair_list.push_back(rule_pair);
		    }
		
		// third non-temrinal...
		if (max_rank >= 3)
		  for (span_pair_set_type::const_iterator niter3 = niter2 + 1; niter3 != niter_end; ++ niter3) 
		    if (span != *niter3
			&& (min_hole_source <= 1 || (niter3->source.second - niter3->source.first) >= min_hole_source)
			&& (min_hole_target <= 1 || (niter3->target.second - niter3->target.first) >= min_hole_target)
			&& is_parent(span.source, niter3->source)
			&& is_parent(span.target, niter3->target)
			&& is_disjoint(niter1->source, niter3->source)
			&& is_disjoint(niter1->target, niter3->target)
			&& is_disjoint(niter2->source, niter3->source)
			&& is_disjoint(niter2->target, niter3->target)
			&& ! is_adjacent(niter1->source, niter3->source)
			&& ! is_adjacent(niter2->source, niter3->source)) {
		  
		      const int source_count3 = alignment_count_source[niter3->source.second] - alignment_count_source[niter3->source.first];
		      //const int target_count3 = alignment_count_target[niter3->target.second] - alignment_count_target[niter3->target.first];
		  
		      const int source_length3 = source_length2 - (niter3->source.second - niter3->source.first);
		      const int target_length3 = target_length2 - (niter3->target.second - niter3->target.first);
		  
		      if (source_count != source_count1 + source_count2 + source_count3)
			if (max_length <= 0 || (source_length3 <= max_length && target_length3 <= max_length))
			  if (max_fertility <= 0 || fertility(source_length3, target_length3) < max_fertility) {
			    extract_rule(source, target, span, *niter1, *niter2, *niter3, category, rule_pair, true);
			    rule_pair_list.push_back(rule_pair);
			  }

		      // forth non-temrinal...
		      if (max_rank >= 4)
			  for (span_pair_set_type::const_iterator niter4 = niter3 + 1; niter4 != niter_end; ++ niter4) 
			    if (span != *niter4
				&& (min_hole_source <= 1 || (niter4->source.second - niter4->source.first) >= min_hole_source)
				&& (min_hole_target <= 1 || (niter4->target.second - niter4->target.first) >= min_hole_target)
				&& is_parent(span.source, niter4->source)
				&& is_parent(span.target, niter4->target)
				&& is_disjoint(niter1->source, niter4->source)
				&& is_disjoint(niter1->target, niter4->target)
				&& is_disjoint(niter2->source, niter4->source)
				&& is_disjoint(niter2->target, niter4->target)
				&& is_disjoint(niter3->source, niter4->source)
				&& is_disjoint(niter3->target, niter4->target)
				&& ! is_adjacent(niter1->source, niter4->source)
				&& ! is_adjacent(niter2->source, niter4->source)
				&& ! is_adjacent(niter3->source, niter4->source)) {
			      
			      
			      const int source_count4 = alignment_count_source[niter4->source.second] - alignment_count_source[niter4->source.first];
			      //const int target_count4 = alignment_count_target[niter4->target.second] - alignment_count_target[niter4->target.first];
			      
			      const int source_length4 = source_length3 - (niter4->source.second - niter4->source.first);
			      const int target_length4 = target_length3 - (niter4->target.second - niter4->target.first);
			      
			      if (source_count != source_count1 + source_count2 + source_count3 + source_count4)
				if (max_length <= 0 || (source_length4 <= max_length && target_length4 <= max_length))
				  if (max_fertility <= 0 || fertility(source_length4, target_length4) < max_fertility) {
				    extract_rule(source, target, span, *niter1, *niter2, *niter3, *niter4, category, rule_pair, true);
				    rule_pair_list.push_back(rule_pair);
				  }
			    }
		    }
	      }
	}
    
    if (! rule_pair_list.empty()) {
      const double count = 1.0 / rule_pair_list.size();
      
      rule_pair_compact_list_type::const_iterator riter_end = rule_pair_list.end();
      for (rule_pair_compact_list_type::const_iterator riter = rule_pair_list.begin(); riter != riter_end; ++ riter)
	const_cast<rule_pair_compact_type&>(*(rule_pairs_local.insert(*riter).first)).count += count;
    }
    
    uniques_source.clear();
    uniques_target.clear();
    uniques_source.resize(rules_source.size(), true);
    uniques_target.resize(rules_target.size(), true);

    uniques_pair.clear();
    
    size_t iter = 0;
    size_t iter_mask = (1024 * 1024 * 4) - 1;

    rule_pair_compact_set_type::const_iterator riter_end = rule_pairs_local.end();
    for (rule_pair_compact_set_type::const_iterator riter = rule_pairs_local.begin(); riter != riter_end; ++ riter, ++ iter) {
      // uncover phrasal representation!
      const bool unique_source = uniques_source[riter->source];
      const bool unique_target = uniques_target[riter->target];

      uniques_source[riter->source] = false;
      uniques_target[riter->target] = false;
      
      std::pair<rule_pair_set_type::iterator, bool> result = rule_pairs.insert(rule_pair_type(rules_source[riter->source],
											      rules_target[riter->target],
											      rules_alignment[riter->alignment],
											      riter->count));
      
      rule_pair_type& rule_pair = const_cast<rule_pair_type&>(*result.first);
      
      if (! result.second)
	rule_pair.count += riter->count;
      
      rule_pair.freqs[0] += uniques_pair.insert(std::make_pair(riter->source, riter->target)).second;;
      rule_pair.freqs[1] += unique_source;
      rule_pair.freqs[2] += unique_target;

      // dumping...
      if ((iter & iter_mask) == iter_mask)
	dumper(rule_pairs);
      
      //rule_pairs_local.erase(riter ++);
    }
    
    rule_pairs_local.clear();
    rules_source.clear();
    rules_target.clear();
    rules_alignment.clear();

    uniques_source.clear();
    uniques_target.clear();
    uniques_pair.clear();
    
    dumper(rule_pairs);
    
    if (rule_pairs.empty()) {
      rule_pair_compact_set_type(rule_pairs_local).swap(rule_pairs_local);
      rule_compact_set_type(rules_source).swap(rules_source);
      rule_compact_set_type(rules_target).swap(rules_target);
      alignment_set_type(rules_alignment).swap(rules_alignment);
      
      unique_set_type(uniques_source).swap(uniques_source);
      unique_set_type(uniques_target).swap(uniques_target);
      unique_pair_set_type(uniques_pair).swap(uniques_pair);
    }
  }

  template <typename Category, typename Dumper>
  void extract_rules(const sentence_type& source,
		     const sentence_type& target,
		     const alignment_type& alignment,
		     rule_pair_set_type& rule_pairs,
		     const Category& category,
		     const Dumper& dumper)
  {
    typedef utils::chunk_vector<rule_pair_compact_type, 4096 / sizeof(rule_pair_compact_type), std::allocator<rule_pair_compact_type> > rule_pair_compact_list_type;

    rule_pair_compact_list_type rule_pair_list;
    rule_pair_compact_type      rule_pair;

    rule_pairs_local.clear();
    rules_source.clear();
    rules_target.clear();
    rules_alignment.clear();
    
    if (! exclude) { // phrase extraction
      span_pair_set_type::const_iterator iter_end = spans.end();
      for (span_pair_set_type::const_iterator iter = spans.begin(); iter != iter_end; ++ iter) {
	const int source_length = iter->source.second - iter->source.first;
	const int target_length = iter->target.second - iter->target.first;
	
	if (max_length <= 0 || (source_length <= max_length && target_length <= max_length))
	  if (max_fertility <= 0 || fertility(source_length, target_length) < max_fertility) {
	    extract_rule(source, target, *iter, category, rule_pair);
	    const_cast<rule_pair_compact_type&>(*(rule_pairs_local.insert(rule_pair).first)).count += 1;
	  }
      }
    }

    span_pair_set_type::const_iterator iter_begin = (constrained ? spans_unique.begin() : spans.begin());
    span_pair_set_type::const_iterator iter_end   = (constrained ? spans_unique.end()   : spans.end());
    
    for (span_pair_set_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter) {
      const int source_length = iter->source.second - iter->source.first;
      const int target_length = iter->target.second - iter->target.first;
      
      if (max_span_source > 0 && source_length > max_span_source) continue;
      if (max_span_target > 0 && target_length > max_span_target) continue;
      
      const int source_count = alignment_count_source[iter->source.second] - alignment_count_source[iter->source.first];
      //const int target_count = alignment_count_target[iter->target.second] - alignment_count_target[iter->target.first];
      
      rule_pair_list.clear();
      
      span_pair_set_type::const_iterator niter_begin = (exhaustive ? spans.begin() : spans_unique.begin());
      span_pair_set_type::const_iterator niter_end   = (exhaustive ? spans.end()   : spans_unique.end());
      
      // first non-terminal...
      if (max_rank >= 1)
	for (span_pair_set_type::const_iterator niter1 = niter_begin; niter1 != niter_end; ++ niter1) 
	  if (*iter != *niter1
	      && (min_hole_source <= 1 || (niter1->source.second - niter1->source.first) >= min_hole_source)
	      && (min_hole_target <= 1 || (niter1->target.second - niter1->target.first) >= min_hole_target)
	      && is_parent(iter->source, niter1->source)
	      && is_parent(iter->target, niter1->target)) {
	  
	    const int source_count1 = alignment_count_source[niter1->source.second] - alignment_count_source[niter1->source.first];
	    //const int target_count1 = alignment_count_target[niter1->target.second] - alignment_count_target[niter1->target.first];
	  
	    const int source_length1 = source_length - (niter1->source.second - niter1->source.first);
	    const int target_length1 = target_length - (niter1->target.second - niter1->target.first);
	  
	    if (source_count1 != source_count)
	      if (max_length <= 0 || (source_length1 <= max_length && target_length1 <= max_length))
		if (max_fertility <= 0 || fertility(source_length1, target_length1) < max_fertility) {
		  extract_rule(source, target, *iter, *niter1, category, rule_pair);
		  rule_pair_list.push_back(rule_pair);
		}
	  
	    // second non-terminal...
	    if (max_rank >= 2)
	      for (span_pair_set_type::const_iterator niter2 = niter1 + 1; niter2 != niter_end; ++ niter2) 
		if (*iter != *niter2
		    && (min_hole_source <= 1 || (niter2->source.second - niter2->source.first) >= min_hole_source)
		    && (min_hole_target <= 1 || (niter2->target.second - niter2->target.first) >= min_hole_target)
		    && is_parent(iter->source, niter2->source)
		    && is_parent(iter->target, niter2->target)
		    && is_disjoint(niter1->source, niter2->source)
		    && is_disjoint(niter1->target, niter2->target)
		    && ! is_adjacent(niter1->source, niter2->source)) {
	      
		  const int source_count2 = alignment_count_source[niter2->source.second] - alignment_count_source[niter2->source.first];
		  //const int target_count2 = alignment_count_target[niter2->target.second] - alignment_count_target[niter2->target.first];
	      
		  const int source_length2 = source_length1 - (niter2->source.second - niter2->source.first);
		  const int target_length2 = target_length1 - (niter2->target.second - niter2->target.first);
	      
		  if (source_count != source_count1 + source_count2)
		    if (max_length <= 0 || (source_length2 <= max_length && target_length2 <= max_length))
		      if (max_fertility <= 0 || fertility(source_length2, target_length2) < max_fertility) {
			extract_rule(source, target, *iter, *niter1, *niter2, category, rule_pair);
			rule_pair_list.push_back(rule_pair);
		      }
	      
		  // third non-temrinal...
		  if (max_rank >= 3)
		    for (span_pair_set_type::const_iterator niter3 = niter2 + 1; niter3 != niter_end; ++ niter3) 
		      if (*iter != *niter3
			  && (min_hole_source <= 1 || (niter3->source.second - niter3->source.first) >= min_hole_source)
			  && (min_hole_target <= 1 || (niter3->target.second - niter3->target.first) >= min_hole_target)
			  && is_parent(iter->source, niter3->source)
			  && is_parent(iter->target, niter3->target)
			  && is_disjoint(niter1->source, niter3->source)
			  && is_disjoint(niter1->target, niter3->target)
			  && is_disjoint(niter2->source, niter3->source)
			  && is_disjoint(niter2->target, niter3->target)
			  && ! is_adjacent(niter1->source, niter3->source)
			  && ! is_adjacent(niter2->source, niter3->source)) {
		    
			const int source_count3 = alignment_count_source[niter3->source.second] - alignment_count_source[niter3->source.first];
			//const int target_count3 = alignment_count_target[niter3->target.second] - alignment_count_target[niter3->target.first];
		    
			const int source_length3 = source_length2 - (niter3->source.second - niter3->source.first);
			const int target_length3 = target_length2 - (niter3->target.second - niter3->target.first);
		    
			if (source_count != source_count1 + source_count2 + source_count3)
			  if (max_length <= 0 || (source_length3 <= max_length && target_length3 <= max_length))
			    if (max_fertility <= 0 || fertility(source_length3, target_length3) < max_fertility) {
			      extract_rule(source, target, *iter, *niter1, *niter2, *niter3, category, rule_pair);
			      rule_pair_list.push_back(rule_pair);
			    }
			
			// forth non-temrinal...
			if (max_rank >= 4)
			  for (span_pair_set_type::const_iterator niter4 = niter3 + 1; niter4 != niter_end; ++ niter4) 
			    if (*iter != *niter4
				&& (min_hole_source <= 1 || (niter4->source.second - niter4->source.first) >= min_hole_source)
				&& (min_hole_target <= 1 || (niter4->target.second - niter4->target.first) >= min_hole_target)
				&& is_parent(iter->source, niter4->source)
				&& is_parent(iter->target, niter4->target)
				&& is_disjoint(niter1->source, niter4->source)
				&& is_disjoint(niter1->target, niter4->target)
				&& is_disjoint(niter2->source, niter4->source)
				&& is_disjoint(niter2->target, niter4->target)
				&& is_disjoint(niter3->source, niter4->source)
				&& is_disjoint(niter3->target, niter4->target)
				&& ! is_adjacent(niter1->source, niter4->source)
				&& ! is_adjacent(niter2->source, niter4->source)
				&& ! is_adjacent(niter3->source, niter4->source)) {
			      
			      
			      const int source_count4 = alignment_count_source[niter4->source.second] - alignment_count_source[niter4->source.first];
			      //const int target_count4 = alignment_count_target[niter4->target.second] - alignment_count_target[niter4->target.first];
			      
			      const int source_length4 = source_length3 - (niter4->source.second - niter4->source.first);
			      const int target_length4 = target_length3 - (niter4->target.second - niter4->target.first);
			      
			      if (source_count != source_count1 + source_count2 + source_count3 + source_count4)
				if (max_length <= 0 || (source_length4 <= max_length && target_length4 <= max_length))
				  if (max_fertility <= 0 || fertility(source_length4, target_length4) < max_fertility) {
				    extract_rule(source, target, *iter, *niter1, *niter2, *niter3, *niter4, category, rule_pair);
				    rule_pair_list.push_back(rule_pair);
				  }
			    }
		      }
		}
	  }
      
      // add into rule-pairs...
      if (! rule_pair_list.empty()) {
	const double count = 1.0 / rule_pair_list.size();
	
	rule_pair_compact_list_type::const_iterator riter_end = rule_pair_list.end();
	for (rule_pair_compact_list_type::const_iterator riter = rule_pair_list.begin(); riter != riter_end; ++ riter)
	  const_cast<rule_pair_compact_type&>(*(rule_pairs_local.insert(*riter).first)).count += count;
	
	rule_pair_list.clear();
      }
    }
    
    uniques_source.clear();
    uniques_target.clear();
    uniques_source.resize(rules_source.size(), true);
    uniques_target.resize(rules_target.size(), true);

    uniques_pair.clear();

    size_t iter = 0;
    size_t iter_mask = (1024 * 1024 * 4) - 1;

    rule_pair_compact_set_type::const_iterator riter_end = rule_pairs_local.end();
    for (rule_pair_compact_set_type::const_iterator riter = rule_pairs_local.begin(); riter != riter_end; ++ riter, ++ iter) {
      // uncover phrasal representation!
      const bool unique_source = uniques_source[riter->source];
      const bool unique_target = uniques_target[riter->target];

      uniques_source[riter->source] = false;
      uniques_target[riter->target] = false;
      
      std::pair<rule_pair_set_type::iterator, bool> result = rule_pairs.insert(rule_pair_type(rules_source[riter->source],
											      rules_target[riter->target],
											      rules_alignment[riter->alignment],
											      riter->count));
      
      rule_pair_type& rule_pair = const_cast<rule_pair_type&>(*result.first);
      
      if (! result.second)
	rule_pair.count += riter->count;
      
      rule_pair.freqs[0] += uniques_pair.insert(std::make_pair(riter->source, riter->target)).second;;
      rule_pair.freqs[1] += unique_source;
      rule_pair.freqs[2] += unique_target;
      
      // dumping...
      if ((iter & iter_mask) == iter_mask)
	dumper(rule_pairs);

      //rule_pairs_local.erase(riter ++);
    }
    
    rule_pairs_local.clear();
    rules_source.clear();
    rules_target.clear();
    rules_alignment.clear();
    
    uniques_source.clear();
    uniques_target.clear();
    uniques_pair.clear();
    
    dumper(rule_pairs);
    
    if (rule_pairs.empty()) {
      rule_pair_compact_set_type(rule_pairs_local).swap(rule_pairs_local);
      rule_compact_set_type(rules_source).swap(rules_source);
      rule_compact_set_type(rules_target).swap(rules_target);
      alignment_set_type(rules_alignment).swap(rules_alignment);
      
      unique_set_type(uniques_source).swap(uniques_source);
      unique_set_type(uniques_target).swap(uniques_target);
      unique_pair_set_type(uniques_pair).swap(uniques_pair);
    }
  }

  struct Builder
  {
    typedef std::vector<char, std::allocator<char> > buffer_type;
    
    typedef buffer_type::iterator       iterator;
    typedef buffer_type::const_iterator const_iterator;
    
    friend
    Builder& operator<<(Builder& builder, const Builder& x)
    {
      builder.buffer.insert(builder.buffer.end(), x.buffer.begin(), x.buffer.end());
      return builder;
    }

    friend
    Builder& operator<<(Builder& builder, const symbol_type& x)
    {
      builder.buffer.insert(builder.buffer.end(), x.begin(), x.end());
      return builder;
    }

    friend
    Builder& operator<<(Builder& builder, const std::string& x)
    {
      builder.buffer.insert(builder.buffer.end(), x.begin(), x.end());
      return builder;
    }

    friend
    Builder& operator<<(Builder& builder, const char* x)
    {
      builder.buffer.insert(builder.buffer.end(), x, x + std::strlen(x));
      return builder;
    }
    
    template <size_t N>
    friend
    Builder& operator<<(Builder& builder, const char (&x)[N])
    {
      builder.buffer.insert(builder.buffer.end(), x, x + N);
      return builder;
    }

    friend
    Builder& operator<<(Builder& builder, const char x)
    {
      builder.buffer.push_back(x);
      return builder;
    }
    
    void clear() { buffer.clear(); }

    operator std::string() const { return std::string(buffer.begin(), buffer.end()); }
    
  private:
    buffer_type buffer;
  };

  Builder builder;
  
  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const Category& category,
		    rule_pair_compact_type& rule_pair,
		    const bool sentential=false)
  {
    const symbol_type& lhs = (sentential ? vocab_type::S : category(spans));

    builder.clear();
    builder << lhs;
    for (int src = spans.source.first; src != spans.source.second; ++ src)
      builder << ' ' << source[src];
    rule_pair.source = rules_source.insert(builder);

    builder.clear();
    builder << lhs;
    for (int trg = spans.target.first; trg != spans.target.second; ++ trg)
      builder << ' ' << target[trg];
    rule_pair.target = rules_target.insert(builder);
    
    rule_pair_type::alignment_type alignment;
    for (int src = spans.source.first; src != spans.source.second; ++ src) {
      point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
      point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
      
      for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	alignment.push_back(std::make_pair(src - spans.source.first, *aiter - spans.target.first));
    }

    rule_pair.alignment = rules_alignment.insert(alignment);
  }

  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const span_pair_type& spans_nt1,
		    const Category& category,
		    rule_pair_compact_type& rule_pair,
		    const bool sentential=false)
  {
    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    const symbol_type nt1 = symbol_type(category(spans_nt1)).non_terminal(1);

    builder.clear();
    builder << lhs;
    for (int src = spans.source.first; src != spans_nt1.source.first; ++ src)
      builder << ' ' << source[src];
    builder << ' ' << nt1;
    for (int src = spans_nt1.source.second; src != spans.source.second; ++ src)
      builder << ' ' << source[src];
    rule_pair.source = rules_source.insert(builder);
    
    builder.clear();
    builder << lhs;
    for (int trg = spans.target.first; trg != spans_nt1.target.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << nt1;
    for (int trg = spans_nt1.target.second; trg != spans.target.second; ++ trg)
      builder << ' ' << target[trg];
    rule_pair.target = rules_target.insert(builder);

    const int nt1_source_size = spans_nt1.source.second - spans_nt1.source.first;
    const int nt1_target_size = spans_nt1.target.second - spans_nt1.target.first;
    
    rule_pair_type::alignment_type alignment;
    for (int src = spans.source.first; src != spans.source.second; ++ src) 
      if (is_out_of_span(spans_nt1.source, src)) {
	point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	
	for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (is_out_of_span(spans_nt1.target, *aiter)) {
	    const int mask_source = - (src >= spans_nt1.source.second);
	    const int mask_target = - (*aiter >= spans_nt1.target.second);
	    
	    // we shift -1, since we have to take into account the <x1> token...
	    const int shift_source = (mask_source & (nt1_source_size - 1)) + spans.source.first;
	    const int shift_target = (mask_target & (nt1_target_size - 1)) + spans.target.first;
	    
	    //const int shift_source = (src >= spans_nt1.source.second ? spans_nt1.source.second - spans_nt1.source.first - 1 : 0) + ...;
	    //const int shift_target = (trg >= spans_nt1.target.second ? spans_nt1.target.second - spans_nt1.target.first - 1 : 0) + ...;
	    
	    alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }

    rule_pair.alignment = rules_alignment.insert(alignment);
  }
  
  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const span_pair_type& __spans_nt1,
		    const span_pair_type& __spans_nt2,
		    const Category& category,
		    rule_pair_compact_type& rule_pair,
		    const bool sentential=false)
  {
    // we will reorder wrt source side spans...
    const bool __is_ordered = is_ordered(__spans_nt1.source, __spans_nt2.source);
    const span_pair_type& spans_nt1 = (__is_ordered ? __spans_nt1 : __spans_nt2);
    const span_pair_type& spans_nt2 = (__is_ordered ? __spans_nt2 : __spans_nt1);
    
    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    const symbol_type nt1 = symbol_type(category(spans_nt1)).non_terminal(1);
    const symbol_type nt2 = symbol_type(category(spans_nt2)).non_terminal(2);
    
    builder.clear();
    builder << lhs;
    for (int src = spans.source.first; src != spans_nt1.source.first; ++ src)
      builder << ' ' << source[src];
    builder << ' ' << nt1;
    for (int src = spans_nt1.source.second; src != spans_nt2.source.first; ++ src)
      builder << ' ' << source[src];
    builder << ' ' << nt2;
    for (int src = spans_nt2.source.second; src != spans.source.second; ++ src)
      builder << ' ' << source[src];
    rule_pair.source = rules_source.insert(builder);
    
    builder.clear();
    builder << lhs;
    if (is_ordered(spans_nt1.target, spans_nt2.target)) {
      for (int trg = spans.target.first; trg != spans_nt1.target.first; ++ trg)
	builder << ' ' << target[trg];
      builder << ' ' << nt1;
      for (int trg = spans_nt1.target.second; trg != spans_nt2.target.first; ++ trg)
	builder << ' ' << target[trg];
      builder << ' ' << nt2;
      for (int trg = spans_nt2.target.second; trg != spans.target.second; ++ trg)
	builder << ' ' << target[trg];
    } else {
      for (int trg = spans.target.first; trg != spans_nt2.target.first; ++ trg)
	builder << ' ' << target[trg];
      builder << ' ' << nt2;
      for (int trg = spans_nt2.target.second; trg != spans_nt1.target.first; ++ trg)
	builder << ' ' << target[trg];
      builder << ' ' << nt1;
      for (int trg = spans_nt1.target.second; trg != spans.target.second; ++ trg)
	builder << ' ' << target[trg];
    }
    rule_pair.target = rules_target.insert(builder);
    
    const int nt1_source_size = spans_nt1.source.second - spans_nt1.source.first;
    const int nt2_source_size = spans_nt2.source.second - spans_nt2.source.first;
    const int nt1_target_size = spans_nt1.target.second - spans_nt1.target.first;
    const int nt2_target_size = spans_nt2.target.second - spans_nt2.target.first;

    rule_pair_type::alignment_type alignment;
    for (int src = spans.source.first; src != spans.source.second; ++ src) 
      if (is_out_of_span(spans_nt1.source, src) && is_out_of_span(spans_nt2.source, src)) {
	point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	
	for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (is_out_of_span(spans_nt1.target, *aiter) && is_out_of_span(spans_nt2.target, *aiter)) {
	    const int mask_source1 = - (src >= spans_nt1.source.second);
	    const int mask_source2 = - (src >= spans_nt2.source.second);
	    
	    const int mask_target1 = - (*aiter >= spans_nt1.target.second);
	    const int mask_target2 = - (*aiter >= spans_nt2.target.second);
	    
	    // we shift -1, since we have to take into account the <x1> token... (also for <x2>)
	    const int shift_source = (spans.source.first
				      + (mask_source1 & (nt1_source_size - 1))
				      + (mask_source2 & (nt2_source_size - 1)));
	    const int shift_target = (spans.target.first
				      + (mask_target1 & (nt1_target_size - 1))
				      + (mask_target2 & (nt2_target_size - 1)));
	    
	    alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }

    rule_pair.alignment = rules_alignment.insert(alignment);
  }

  struct less_source
  {
    template <typename Tp>
    bool operator()(const Tp& x, const Tp& y) const
    {
      return x.source < y.source;
    }
  };
  
  struct less_first
  {
    template <typename Tp>
    bool operator()(const Tp& x, const Tp& y) const
    {
      return x.first < y.first;
    }
  };
  
  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const span_pair_type& __spans_nt1,
		    const span_pair_type& __spans_nt2,
		    const span_pair_type& __spans_nt3,
		    const Category& category,
		    rule_pair_compact_type& rule_pair,
		    const bool sentential=false)
  {
    typedef std::pair<span_type, symbol_type> span_category_type;
    
    // sort by source-side span...
    
    boost::array<span_pair_type, 3> spans_nt;
    spans_nt[0] = __spans_nt1;
    spans_nt[1] = __spans_nt2;
    spans_nt[2] = __spans_nt3;
    
    std::sort(spans_nt.begin(), spans_nt.end(), less_source());
    
    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    
    boost::array<symbol_type, 3> cat;
    cat[0] = symbol_type(category(spans_nt[0])).non_terminal(1);
    cat[1] = symbol_type(category(spans_nt[1])).non_terminal(2);
    cat[2] = symbol_type(category(spans_nt[2])).non_terminal(3);
    
    builder.clear();
    builder << lhs;
    for (int src = spans.source.first; src != spans_nt[0].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[0];
    for (int src = spans_nt[0].source.second; src != spans_nt[1].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[1];
    for (int src = spans_nt[1].source.second; src != spans_nt[2].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[2];
    for (int src = spans_nt[2].source.second; src != spans.source.second; ++ src)
      builder << ' '  << source[src];
    rule_pair.source = rules_source.insert(builder);
    
    // sort by target-side span with category,...
    boost::array<span_category_type, 3> spans_cat;
    spans_cat[0] = std::make_pair(spans_nt[0].target, cat[0]);
    spans_cat[1] = std::make_pair(spans_nt[1].target, cat[1]);
    spans_cat[2] = std::make_pair(spans_nt[2].target, cat[2]);
    
    std::sort(spans_cat.begin(), spans_cat.end(), less_first());
    
    builder.clear();
    builder << lhs;
    for (int trg = spans.target.first; trg != spans_cat[0].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[0].second;
    for (int trg = spans_cat[0].first.second; trg != spans_cat[1].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[1].second;
    for (int trg = spans_cat[1].first.second; trg != spans_cat[2].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[2].second;
    for (int trg = spans_cat[2].first.second; trg != spans.target.second; ++ trg)
      builder << ' ' << target[trg];
    rule_pair.target = rules_target.insert(builder);
    
    const int nt1_source_size = spans_nt[0].source.second - spans_nt[0].source.first;
    const int nt2_source_size = spans_nt[1].source.second - spans_nt[1].source.first;
    const int nt3_source_size = spans_nt[2].source.second - spans_nt[2].source.first;
    
    const int nt1_target_size = spans_nt[0].target.second - spans_nt[0].target.first;
    const int nt2_target_size = spans_nt[1].target.second - spans_nt[1].target.first;
    const int nt3_target_size = spans_nt[2].target.second - spans_nt[2].target.first;
    
    rule_pair_type::alignment_type alignment;
    for (int src = spans.source.first; src != spans.source.second; ++ src)
      if (is_out_of_span(spans_nt[0].source, src)
	  && is_out_of_span(spans_nt[1].source, src)
	  && is_out_of_span(spans_nt[2].source, src)) {
	point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	
	for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (is_out_of_span(spans_nt[0].target, *aiter)
	      && is_out_of_span(spans_nt[1].target, *aiter)
	      && is_out_of_span(spans_nt[2].target, *aiter)) {
	    
	    const int mask_source1 = - (src >= spans_nt[0].source.second);
	    const int mask_source2 = - (src >= spans_nt[1].source.second);
	    const int mask_source3 = - (src >= spans_nt[2].source.second);
	    
	    const int mask_target1 = - (*aiter >= spans_nt[0].target.second);
	    const int mask_target2 = - (*aiter >= spans_nt[1].target.second);
	    const int mask_target3 = - (*aiter >= spans_nt[2].target.second);
	    
	    // we shift -1, since we have to take into account the <x1> token... (also for <x2>)
	    const int shift_source = (spans.source.first
				      + (mask_source1 & (nt1_source_size - 1))
				      + (mask_source2 & (nt2_source_size - 1))
				      + (mask_source3 & (nt3_source_size - 1)));
	    const int shift_target = (spans.target.first
				      + (mask_target1 & (nt1_target_size - 1))
				      + (mask_target2 & (nt2_target_size - 1))
				      + (mask_target3 & (nt3_target_size - 1)));
	    
	    alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }
    
    rule_pair.alignment = rules_alignment.insert(alignment);
  }
  
  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const span_pair_type& __spans_nt1,
		    const span_pair_type& __spans_nt2,
		    const span_pair_type& __spans_nt3,
		    const span_pair_type& __spans_nt4,
		    const Category& category,
		    rule_pair_compact_type& rule_pair,
		    const bool sentential=false)
  {
    typedef std::pair<span_type, symbol_type> span_category_type;
    
    // sort by source-side span...
    
    boost::array<span_pair_type, 4> spans_nt;
    spans_nt[0] = __spans_nt1;
    spans_nt[1] = __spans_nt2;
    spans_nt[2] = __spans_nt3;
    spans_nt[3] = __spans_nt4;
    
    std::sort(spans_nt.begin(), spans_nt.end(), less_source());
    
    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    
    boost::array<symbol_type, 4> cat;
    cat[0] = symbol_type(category(spans_nt[0])).non_terminal(1);
    cat[1] = symbol_type(category(spans_nt[1])).non_terminal(2);
    cat[2] = symbol_type(category(spans_nt[2])).non_terminal(3);
    cat[3] = symbol_type(category(spans_nt[3])).non_terminal(4);
    
    builder.clear();
    builder << lhs;
    for (int src = spans.source.first; src != spans_nt[0].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[0];
    for (int src = spans_nt[0].source.second; src != spans_nt[1].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[1];
    for (int src = spans_nt[1].source.second; src != spans_nt[2].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[2];
    for (int src = spans_nt[2].source.second; src != spans_nt[3].source.first; ++ src)
      builder << ' '  << source[src];
    builder << ' ' << cat[3];
    for (int src = spans_nt[3].source.second; src != spans.source.second; ++ src)
      builder << ' '  << source[src];
    rule_pair.source = rules_source.insert(builder);
    
    // sort by target-side span with category,...
    boost::array<span_category_type, 4> spans_cat;
    spans_cat[0] = std::make_pair(spans_nt[0].target, cat[0]);
    spans_cat[1] = std::make_pair(spans_nt[1].target, cat[1]);
    spans_cat[2] = std::make_pair(spans_nt[2].target, cat[2]);
    spans_cat[3] = std::make_pair(spans_nt[3].target, cat[3]);
    
    std::sort(spans_cat.begin(), spans_cat.end(), less_first());
    
    builder.clear();
    builder << lhs;
    for (int trg = spans.target.first; trg != spans_cat[0].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[0].second;
    for (int trg = spans_cat[0].first.second; trg != spans_cat[1].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[1].second;
    for (int trg = spans_cat[1].first.second; trg != spans_cat[2].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[2].second;
    for (int trg = spans_cat[2].first.second; trg != spans_cat[3].first.first; ++ trg)
      builder << ' ' << target[trg];
    builder << ' ' << spans_cat[3].second;
    for (int trg = spans_cat[3].first.second; trg != spans.target.second; ++ trg)
      builder << ' ' << target[trg];
    rule_pair.target = rules_target.insert(builder);
    
    const int nt1_source_size = spans_nt[0].source.second - spans_nt[0].source.first;
    const int nt2_source_size = spans_nt[1].source.second - spans_nt[1].source.first;
    const int nt3_source_size = spans_nt[2].source.second - spans_nt[2].source.first;
    const int nt4_source_size = spans_nt[3].source.second - spans_nt[3].source.first;
    
    const int nt1_target_size = spans_nt[0].target.second - spans_nt[0].target.first;
    const int nt2_target_size = spans_nt[1].target.second - spans_nt[1].target.first;
    const int nt3_target_size = spans_nt[2].target.second - spans_nt[2].target.first;
    const int nt4_target_size = spans_nt[3].target.second - spans_nt[3].target.first;
    
    rule_pair_type::alignment_type alignment;
    for (int src = spans.source.first; src != spans.source.second; ++ src)
      if (is_out_of_span(spans_nt[0].source, src)
	  && is_out_of_span(spans_nt[1].source, src)
	  && is_out_of_span(spans_nt[2].source, src)
	  && is_out_of_span(spans_nt[3].source, src)) {
	point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	
	for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (is_out_of_span(spans_nt[0].target, *aiter)
	      && is_out_of_span(spans_nt[1].target, *aiter)
	      && is_out_of_span(spans_nt[2].target, *aiter)
	      && is_out_of_span(spans_nt[3].target, *aiter)) {
	    const int mask_source1 = - (src >= spans_nt[0].source.second);
	    const int mask_source2 = - (src >= spans_nt[1].source.second);
	    const int mask_source3 = - (src >= spans_nt[2].source.second);
	    const int mask_source4 = - (src >= spans_nt[3].source.second);
	    
	    const int mask_target1 = - (*aiter >= spans_nt[0].target.second);
	    const int mask_target2 = - (*aiter >= spans_nt[1].target.second);
	    const int mask_target3 = - (*aiter >= spans_nt[2].target.second);
	    const int mask_target4 = - (*aiter >= spans_nt[3].target.second);
	    
	    // we shift -1, since we have to take into account the <x1> token... (also for <x2>)
	    const int shift_source = (spans.source.first
				      + (mask_source1 & (nt1_source_size - 1))
				      + (mask_source2 & (nt2_source_size - 1))
				      + (mask_source3 & (nt3_source_size - 1))
				      + (mask_source4 & (nt4_source_size - 1)));
	    const int shift_target = (spans.target.first
				      + (mask_target1 & (nt1_target_size - 1))
				      + (mask_target2 & (nt2_target_size - 1))
				      + (mask_target3 & (nt3_target_size - 1))
				      + (mask_target4 & (nt4_target_size - 1)));
	    
	    alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }
    
    rule_pair.alignment = rules_alignment.insert(alignment);
  }
  

  bool is_out_of_span(const span_type& span, const int pos) const
  {
    return pos < span.first || span.second <= pos;
  }
  
  bool is_parent(const span_type& parent, const span_type& child) const
  {
    return parent.first <= child.first && child.second <= parent.second;
  }
  
  bool is_ordered(const span_type& span1, const span_type& span2) const 
  {
    return span1.second <= span2.first;
  }

  bool is_adjacent(const span_type& span1, const span_type& span2) const 
  {
    return span1.second == span2.first || span2.second == span1.first;
  }
  
  bool is_disjoint(const span_type& span1, const span_type& span2) const 
  {
    return span1.second <= span2.first || span2.second <= span1.first;
  }
  
  double fertility(int length1, int length2)
  {
    return double(utils::bithack::max(length1, length2)) / utils::bithack::min(length1, length2);
  }
  
  void compute_spans(const sentence_type& source,
		     const sentence_type& target,
		     const alignment_type& alignment)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    alignment_source_target.clear();
    alignment_target_source.clear();
    alignment_source_target.reserve(source_size);
    alignment_target_source.reserve(target_size);
    alignment_source_target.resize(source_size);
    alignment_target_source.resize(target_size);
    
    alignment_count_source.clear();
    alignment_count_target.clear();
    alignment_count_source.reserve(source_size + 1);
    alignment_count_target.reserve(target_size + 1);
    alignment_count_source.resize(source_size + 1, 0);
    alignment_count_target.resize(target_size + 1, 0);
    
    spans.clear();
    spans_unique.clear();
    
    if (inverse) {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
	if (aiter->target >= static_cast<int>(alignment_source_target.size()))
	  throw std::runtime_error("invalid alignment");
	if (aiter->source >= static_cast<int>(alignment_target_source.size()))
	  throw std::runtime_error("invalid alignment");
	
	alignment_source_target[aiter->target].push_back(aiter->source);
	alignment_target_source[aiter->source].push_back(aiter->target);
      }
    } else {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
	if (aiter->source >= static_cast<int>(alignment_source_target.size()))
	  throw std::runtime_error("invalid alignment");
	if (aiter->target >= static_cast<int>(alignment_target_source.size()))
	  throw std::runtime_error("invalid alignment");
	
	alignment_source_target[aiter->source].push_back(aiter->target);
	alignment_target_source[aiter->target].push_back(aiter->source);
      }
    }
    
    for (size_type src = 0; src != source_size; ++ src) {
      std::sort(alignment_source_target[src].begin(), alignment_source_target[src].end());
      alignment_count_source[src + 1] = alignment_count_source[src] + alignment_source_target[src].size();
    }
    
    for (size_type trg = 0; trg != target_size; ++ trg) {
      std::sort(alignment_target_source[trg].begin(), alignment_target_source[trg].end());
      alignment_count_target[trg + 1] = alignment_count_target[trg] + alignment_target_source[trg].size();
    }
    
    span_target_chart.clear();
    span_target_chart.reserve(target_size + 1);
    span_target_chart.resize(target_size + 1, span_type(source_size, 0));
    
    // first, collect alignment count only for the target side
    for (int target_first = 0; target_first < static_cast<int>(target_size); ++ target_first) {
      span_type span_source(source_size, 0);
	
      for (int target_last = target_first + 1; target_last <= static_cast<int>(target_size); ++ target_last) {
	
	if (! alignment_target_source[target_last - 1].empty()) {
	  span_source.first  = utils::bithack::min(int(span_source.first),  int(alignment_target_source[target_last - 1].front()));
	  span_source.second = utils::bithack::max(int(span_source.second), int(alignment_target_source[target_last - 1].back()) + 1);
	}
	  
	span_target_chart(target_first, target_last) = span_source;
      }
    }

    // then, collect for source-target
    
    for (int source_first = 0; source_first < static_cast<int>(source_size); ++ source_first) {
      span_type span_target(target_size, 0);
	
      for (int source_last = source_first + 1; source_last <= static_cast<int>(source_size); ++ source_last) {
	  
	if (! alignment_source_target[source_last - 1].empty()) {
	  span_target.first  = utils::bithack::min(int(span_target.first),  int(alignment_source_target[source_last - 1].front()));
	  span_target.second = utils::bithack::max(int(span_target.second), int(alignment_source_target[source_last - 1].back()) + 1);
	}
	  
	const int span_count_source = alignment_count_source[source_last] - alignment_count_source[source_first];
	
	if (span_count_source > 0 && span_target.second - span_target.first > 0) {
	  
	  const int span_count_target = alignment_count_target[span_target.second] - alignment_count_target[span_target.first];
	  
	  if (span_count_source == span_count_target) {
	    
	    // unique span-pair
	    const span_type& span_source = span_target_chart(span_target.first, span_target.second);
	    if (span_source.first == source_first && span_source.second == source_last)
	      spans_unique.push_back(span_pair_type(span_source, span_target));
	    
	    // enlarge the target-span...
	    for (int target_first = span_target.first; target_first >= 0; -- target_first)
	      for (int target_last = span_target.second; target_last <= static_cast<int>(target_size); ++ target_last) {
		
		const int span_count_target = alignment_count_target[target_last] - alignment_count_target[target_first];
		
		if (span_count_source != span_count_target)
		  break;

		spans.push_back(span_pair_type(span_type(source_first, source_last), span_type(target_first, target_last)));
	      }
	  }
	}
      }
    }
  }

};


struct Task
{
  typedef boost::filesystem::path path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef Bitext bitext_type;
  
  typedef ExtractSCFG::rule_pair_type     rule_pair_type;
  typedef ExtractSCFG::rule_pair_set_type rule_pair_set_type;
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> > queue_type;
  
  
  Task(queue_type& __queue,
       const path_type& __output,
       const int max_length,
       const int max_fertility,
       const int max_span_source,
       const int max_span_target,
       const int min_hole_source,
       const int min_hole_target,
       const int max_rank,
       const int max_scope,
       const bool exhaustive,
       const bool constrained,
       const bool exclude,
       const bool sentential,
       const bool inverse,
       const double __max_malloc)
    : queue(__queue),
      output(__output),
      extractor(max_length, max_fertility,
		max_span_source, max_span_target, 
		min_hole_source, min_hole_target,
		max_rank, max_scope,
		exhaustive, constrained, exclude, sentential, inverse),
      max_malloc(__max_malloc) {}
  
  queue_type&   queue;
  path_type     output;
  path_set_type paths;
  
  ExtractSCFG extractor;
  
  double max_malloc;

  struct Dumper
  {
    void operator()(rule_pair_set_type& rule_pairs) const
    {
      if (rule_pairs.size() < 1024 * 1024) return;
      
      const size_t allocated = utils::malloc_stats::used() + rule_pairs.size() * sizeof(void*);
      
      if ((min_counts_size && rule_pairs.size() < min_counts_size && allocated <= (malloc_threshold + (malloc_threshold >> 1)))
	  || allocated <= malloc_threshold) return;
      
      if (! min_counts_size)
	const_cast<size_t&>(min_counts_size) = rule_pairs.size() >> 4;
      
      dump(rule_pairs);
      
      rule_pairs.clear();
      rule_pair_set_type(rule_pairs).swap(rule_pairs);
    }

    typedef std::vector<const rule_pair_type*, std::allocator<const rule_pair_type*> > sorted_type;
    
    template <typename Tp>
    struct less_ptr
    {
      bool operator()(const Tp* x, const Tp* y) const
      {
	return *x < *y;
      }
    };
        
    void dump(const rule_pair_set_type& rule_pairs) const
    {
      if (rule_pairs.empty()) return;
      
      // sorting...
      sorted_type sorted(rule_pairs.size());
      {
	sorted_type::iterator siter = sorted.begin();
	rule_pair_set_type::const_iterator citer_end = rule_pairs.end();
	for (rule_pair_set_type::const_iterator citer = rule_pairs.begin(); citer != citer_end; ++ citer, ++ siter)
	  *siter = &(*citer);
      }
      std::sort(sorted.begin(), sorted.end(), less_ptr<rule_pair_type>());
      
      const path_type path_tmp = utils::tempfile::file_name(output / "counts-XXXXXX");
      utils::tempfile::insert(path_tmp);
      const path_type path = path_tmp.string() + ".gz";
      utils::tempfile::insert(path);
      
      utils::compress_ostream os(path, 1024 * 1024);
      sorted_type::const_iterator siter_end = sorted.end();
      for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
	generator(os, *(*siter)) << '\n';
      
      const_cast<path_set_type&>(paths).push_back(path);
    }

    Dumper(const path_type& __output,
	   path_set_type& __paths,
	   const size_t __malloc_threshold)
      : output(__output),
	paths(__paths),
	malloc_threshold(__malloc_threshold),
	min_counts_size(0) {}
    
    const path_type& output;
    path_set_type& paths;
    const size_t malloc_threshold;
    size_t min_counts_size;

    RulePairGenerator generator;
  };

  
  void operator()()
  {
    bitext_type bitext;
    rule_pair_set_type rule_pairs;
    
    Dumper dumper(output, paths, max_malloc * 1024 * 1024 * 1024);
    
    for (int iter = 0;/**/; ++ iter) {
      bitext.clear();
      queue.pop_swap(bitext);
      
      if (bitext.source.empty()) break;
      
      extractor(bitext.source, bitext.target, bitext.alignment, bitext.spans_source, bitext.spans_target, rule_pairs, dumper);
    }
    
    dumper.dump(rule_pairs);
  }
};

#endif
