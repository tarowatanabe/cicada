//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXTRACT_SCFG_IMPL__HPP__
#define __CICADA__EXTRACT_SCFG_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>

#include <boost/array.hpp>

#include "cicada/sentence.hpp"
#include "cicada/alignment.hpp"
#include "cicada/span_vector.hpp"
#include "cicada/vocab.hpp"

#include "utils/sgi_hash_map.hpp"
#include "utils/sgi_hash_set.hpp"
#include "utils/chart.hpp"

#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/chunk_vector.hpp>

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
    if (std::getline(is, line)) {
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
  
  phrase_type    source;
  phrase_type    target;
  alignment_type alignment;
  count_type     count;

  RulePair() : source(), target(), alignment(), count(0) {}

  friend
  size_t hash_value(RulePair const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    return hasher_type()(x.source.begin(), x.source.end(),
			 hasher_type()(x.target.begin(), x.target.end(),
				       hasher_type()(x.alignment.begin(), x.alignment.end(), 0)));
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
    return (x.source < y.source
	    || (!(y.source < x.source)
		&& (x.target < y.target
		    || (!(y.target < x.target)
			&& x.alignment < y.alignment))));
  }
  
  friend
  bool operator>(const RulePair& x, const RulePair& y)
  {
    return y < x;
  }
};

BOOST_FUSION_ADAPT_STRUCT(RulePair::alignment_type::point_type,
			  (RulePair::alignment_type::index_type, source)
			  (RulePair::alignment_type::index_type, target)
			  )
BOOST_FUSION_ADAPT_STRUCT(
			  RulePair,
			  (RulePair::phrase_type, source)
			  (RulePair::phrase_type, target)
			  (RulePair::alignment_type, alignment)
			  (RulePair::count_type, count)
			  )

struct RulePairGenerator
{
  typedef RulePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::count_type    count_type;
  
  RulePairGenerator() : grammar() {}
  RulePairGenerator(const RulePairGenerator& x) : grammar() {}

  
  template <typename Iterator>
  struct phrase_pair_generator : boost::spirit::karma::grammar<Iterator, phrase_pair_type()>
  {
    phrase_pair_generator() : phrase_pair_generator::base_type(phrase_pair)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      alignment %= -((karma::int_ << '-' << karma::int_) % ' ');
      phrase_pair %= standard::string << " ||| " << standard::string << " ||| " << alignment << " ||| " << double20;
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, alignment_type()> alignment;
    boost::spirit::karma::rule<Iterator, phrase_pair_type()> phrase_pair;
  };

  typedef std::ostream_iterator<char> iterator_type;
  
  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair) const
  {
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, phrase_pair))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  phrase_pair_generator<iterator_type> grammar;
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
      return utils::hashmurmur<size_t>()(x, 0);
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
  

#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
				  std::allocator<rule_pair_type> > rule_pair_set_type;
#else
  typedef sgi::hash_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
			std::allocator<rule_pair_type> > rule_pair_set_type;
#endif
  
  typedef utils::chart<span_type, std::allocator<span_type> >          span_chart_type;
  typedef std::vector<int, std::allocator<int> >                       alignment_count_set_type;
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  
  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > alignment_multiple_type;

  typedef utils::chart<phrase_type, std::allocator<phrase_type> >      phrase_chart_type;
  
  typedef std::pair<span_type, span_type> span_parent_type;


  ExtractSCFG(const int __max_length,
	      const int __max_fertility,
	      const int __max_span,
	      const int __min_hole,
	      const bool __ternary,
	      const bool __sentential,
	      const bool __inverse)
    : max_length(__max_length),
      max_fertility(__max_fertility),
      max_span(__max_span),
      min_hole(__min_hole),
      ternary(__ternary),
      sentential(__sentential),
      inverse(__inverse) {}
		
  int max_length;
  int max_fertility;
  int max_span;
  int min_hole;
  bool ternary;
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
  };
  
  struct ExtractCategoryBase
  {
    typedef utils::chart<symbol_type, std::allocator<std::string> > label_chart_type;
    
#ifdef HAVE_TR1_UNORDERED_MAP
    typedef std::tr1::unordered_map<span_type, symbol_type, utils::hashmurmur<size_t>, std::equal_to<span_type>,
				    std::allocator<std::pair<const span_type, symbol_type> > > label_map_type;
#else
    typedef sgi::hash_map<span_type, symbol_type, utils::hashmurmur<size_t>, std::equal_to<span_type>,
			  std::allocator<std::pair<const span_type, symbol_type> > > label_map_type;
  
#endif
    
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
	  label = '[' + liter->second.non_terminal_strip() + '\\' + siter->second.non_terminal_strip() + ']';
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
  };
  
  struct ExtractCategoryTarget : public ExtractCategoryBase
  {
    ExtractCategoryTarget(const sentence_type& sentence, const span_set_type& spans)
      : ExtractCategoryBase(sentence, spans) {}
    
    symbol_type operator()(const span_pair_type& spans) const
    {
      return extract(spans.target);
    }
  };
  
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
    
    if (! spans_source.empty())
      extract_rules(source, target, alignment, rule_pairs, ExtractCategorySource(source, spans_source), dumper);
    else if (! spans_target.empty())
      extract_rules(source, target, alignment, rule_pairs, ExtractCategoryTarget(target, spans_target), dumper);
    else
      extract_rules(source, target, alignment, rule_pairs, ExtractCategory(), dumper);
  }
  
  template <typename Category, typename Dumper>
  void extract_rules(const sentence_type& source,
		     const sentence_type& target,
		     const alignment_type& alignment,
		     rule_pair_set_type& rule_pairs,
		     const Category& category,
		     const Dumper& dumper)
  {
    typedef utils::chunk_vector<rule_pair_type, 4096 / sizeof(rule_pair_type), std::allocator<rule_pair_type> > rule_pair_list_type;

    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    rule_pair_list_type rule_pair_list;
    rule_pair_list_type sentential_pair_list;
    rule_pair_type rule_pair;
    
    span_pair_set_type::const_iterator iter_end = spans.end();
    for (span_pair_set_type::const_iterator iter = spans.begin(); iter != iter_end; ++ iter) {
      const int source_length = iter->source.second - iter->source.first;
      const int target_length = iter->target.second - iter->target.first;
      
      const int source_count = alignment_count_source[iter->source.second] - alignment_count_source[iter->source.first];
      //const int target_count = alignment_count_target[iter->target.second] - alignment_count_target[iter->target.first];
      
      const bool sentential_mode = sentential && source_length == static_cast<int>(source_size) && target_length == static_cast<int>(target_size);
      const bool rule_mode = (! sentential) && (max_span <= 0 || source_length <= max_span);
      
      if (max_length <= 0 || source_length <= max_length || target_length <= max_length)
	if (max_fertility <= 0 || fertility(source_length, target_length) < max_fertility) {
	  // extract rule...

	  if (! sentential) {
	    extract_rule(source, target, *iter, category, rule_pair);
	    const_cast<rule_pair_type&>(*(rule_pairs.insert(rule_pair).first)).count += 1;
	  } else if (sentential_mode) {
	    extract_rule(source, target, *iter, category, rule_pair, true);
	    const_cast<rule_pair_type&>(*(rule_pairs.insert(rule_pair).first)).count += 1;
	  }
	}
      
      // consider hole!
      if (rule_mode || sentential_mode) {
	rule_pair_list.clear();
	sentential_pair_list.clear();
	
	span_pair_set_type::const_iterator niter_end = spans_unique.end();
	
	// first non-terminal...
	for (span_pair_set_type::const_iterator niter1 = spans_unique.begin(); niter1 != niter_end; ++ niter1) 
	  if (*iter != *niter1
	      && (min_hole <= 1 || (niter1->source.second - niter1->source.first) >= min_hole)
	      && is_parent(iter->source, niter1->source)
	      && is_parent(iter->target, niter1->target)) {
	    
	    const int source_count1 = alignment_count_source[niter1->source.second] - alignment_count_source[niter1->source.first];
	    //const int target_count1 = alignment_count_target[niter1->target.second] - alignment_count_target[niter1->target.first];
	    
	    if (source_count1 == source_count) continue;
	    
	    const int source_length1 = source_length - (niter1->source.second - niter1->source.first);
	    const int target_length1 = target_length - (niter1->target.second - niter1->target.first);
	    
	    if (max_length <= 0 || source_length1 <= max_length || target_length1 <= max_length)
	      if (max_fertility <= 0 || fertility(source_length1, target_length1) < max_fertility) {
		// extract rule...
		
		if (rule_mode) {
		  extract_rule(source, target, *iter, *niter1, category, rule_pair);
		  rule_pair_list.push_back(rule_pair);
		} else if (sentential_mode) {
		  extract_rule(source, target, *iter, *niter1, category, rule_pair, true);
		  sentential_pair_list.push_back(rule_pair);
		}
	      }
	    
	    // second non-terminal...
	    for (span_pair_set_type::const_iterator niter2 = niter1 + 1; niter2 != niter_end; ++ niter2) 
	      if (*iter != *niter2
		  && (min_hole <= 1 || (niter2->source.second - niter2->source.first) >= min_hole)
		  && is_parent(iter->source, niter2->source)
		  && is_parent(iter->target, niter2->target)
		  && is_disjoint(niter1->source, niter2->source)
		  && is_disjoint(niter1->target, niter2->target)
		  && ! is_adjacent(niter1->source, niter2->source)) {
		
		const int source_count2 = alignment_count_source[niter2->source.second] - alignment_count_source[niter2->source.first];
		//const int target_count2 = alignment_count_target[niter2->target.second] - alignment_count_target[niter2->target.first];
		
		if (source_count == source_count1 + source_count2) continue;
		
		const int source_length2 = source_length1 - (niter2->source.second - niter2->source.first);
		const int target_length2 = target_length1 - (niter2->target.second - niter2->target.first);
		
		if (max_length <= 0 || source_length2 <= max_length || target_length2 <= max_length)
		  if (max_fertility <= 0 || fertility(source_length2, target_length2) < max_fertility) {
		    // extract rule...
		    
		    if (rule_mode) {
		      extract_rule(source, target, *iter, *niter1, *niter2, category, rule_pair);
		      rule_pair_list.push_back(rule_pair);
		    } else if (sentential_mode) {
		      extract_rule(source, target, *iter, *niter1, *niter2, category, rule_pair, true);
		      sentential_pair_list.push_back(rule_pair);
		    }
		  }

		
		if (ternary)
		  for (span_pair_set_type::const_iterator niter3 = niter2 + 1; niter3 != niter_end; ++ niter3) 
		    if (*iter != *niter3
			&& (min_hole <= 1 || (niter3->source.second - niter3->source.first) >= min_hole)
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
		      
		      if (source_count == source_count1 + source_count2 + source_count3) continue;
		      
		      const int source_length3 = source_length2 - (niter3->source.second - niter3->source.first);
		      const int target_length3 = target_length2 - (niter3->target.second - niter3->target.first);
		      
		      if (max_length <= 0 || source_length3 <= max_length || target_length3 <= max_length)
			if (max_fertility <= 0 || fertility(source_length3, target_length3) < max_fertility) {
			  // extract rule...
			  
			  if (rule_mode) {
			    extract_rule(source, target, *iter, *niter1, *niter2, *niter3, category, rule_pair);
			    rule_pair_list.push_back(rule_pair);
			  } else if (sentential_mode) {
			    extract_rule(source, target, *iter, *niter1, *niter2, *niter3, category, rule_pair, true);
			    sentential_pair_list.push_back(rule_pair);
			  }
			}
		    }
	      }
	  }
	
	// add into rule-pairs...
	if (! rule_pair_list.empty()) {
	  const double count = 1.0 / rule_pair_list.size();
	  
	  rule_pair_list_type::const_iterator riter_end = rule_pair_list.end();
	  for (rule_pair_list_type::const_iterator riter = rule_pair_list.begin(); riter != riter_end; ++ riter)
	    const_cast<rule_pair_type&>(*(rule_pairs.insert(*riter).first)).count += count;

	  rule_pair_list.clear();
	}
	
	// add into rule-pairs...
	if (! sentential_pair_list.empty()) {
	  const double count = 1.0 / sentential_pair_list.size();
	  
	  rule_pair_list_type::const_iterator riter_end = sentential_pair_list.end();
	  for (rule_pair_list_type::const_iterator riter = sentential_pair_list.begin(); riter != riter_end; ++ riter)
	    const_cast<rule_pair_type&>(*(rule_pairs.insert(*riter).first)).count += count;
	  
	  sentential_pair_list.clear();
	}
      }
      
      dumper(rule_pairs);
    }
  }
  
  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const Category& category,
		    rule_pair_type& rule_pair,
		    const bool sentential=false)
  {
    const symbol_type& lhs = (sentential ? vocab_type::S : category(spans));

    rule_pair.source = static_cast<const std::string&>(lhs);
    for (int src = spans.source.first; src != spans.source.second; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    
    rule_pair.target = static_cast<const std::string&>(lhs);
    for (int trg = spans.target.first; trg != spans.target.second; ++ trg)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
    
    rule_pair.alignment.clear();
    for (int src = spans.source.first; src != spans.source.second; ++ src) {
      point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
      point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
      
      for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	rule_pair.alignment.push_back(std::make_pair(src - spans.source.first, *aiter - spans.target.first));
    }
  }

  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const span_pair_type& spans_nt1,
		    const Category& category,
		    rule_pair_type& rule_pair,
		    const bool sentential=false)
  {
    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    const symbol_type nt1 = category(spans_nt1);

    rule_pair.source = static_cast<const std::string&>(lhs);
    for (int src = spans.source.first; src != spans_nt1.source.first; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    rule_pair.source += ' ' + static_cast<const std::string&>(nt1.non_terminal(1));
    for (int src = spans_nt1.source.second; src != spans.source.second; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    
    rule_pair.target = static_cast<const std::string&>(lhs);
    for (int trg = spans.target.first; trg != spans_nt1.target.first; ++ trg)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
    rule_pair.target += ' ' + static_cast<const std::string&>(nt1.non_terminal(1));
    for (int trg = spans_nt1.target.second; trg != spans.target.second; ++ trg)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);

    const int nt1_source_size = spans_nt1.source.second - spans_nt1.source.first;
    const int nt1_target_size = spans_nt1.target.second - spans_nt1.target.first;
    
    rule_pair.alignment.clear();
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
	    
	    rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }
  }
  
  template <typename Category>
  void extract_rule(const sentence_type& source,
		    const sentence_type& target,
		    const span_pair_type& spans,
		    const span_pair_type& __spans_nt1,
		    const span_pair_type& __spans_nt2,
		    const Category& category,
		    rule_pair_type& rule_pair,
		    const bool sentential=false)
  {
    // we will reorder wrt source side spans...
    const bool __is_ordered = is_ordered(__spans_nt1.source, __spans_nt2.source);
    const span_pair_type& spans_nt1 = (__is_ordered ? __spans_nt1 : __spans_nt2);
    const span_pair_type& spans_nt2 = (__is_ordered ? __spans_nt2 : __spans_nt1);
    
    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    const symbol_type nt1 = category(spans_nt1);
    const symbol_type nt2 = category(spans_nt2);
    
    rule_pair.source = static_cast<const std::string&>(lhs);
    for (int src = spans.source.first; src != spans_nt1.source.first; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    rule_pair.source += ' ' + static_cast<const std::string&>(nt1.non_terminal(1));
    for (int src = spans_nt1.source.second; src != spans_nt2.source.first; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    rule_pair.source += ' ' + static_cast<const std::string&>(nt2.non_terminal(2));
    for (int src = spans_nt2.source.second; src != spans.source.second; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);

    if (is_ordered(spans_nt1.target, spans_nt2.target)) {
      rule_pair.target = static_cast<const std::string&>(lhs);
      for (int trg = spans.target.first; trg != spans_nt1.target.first; ++ trg)
	rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
      rule_pair.target += ' ' + static_cast<const std::string&>(nt1.non_terminal(1));
      for (int trg = spans_nt1.target.second; trg != spans_nt2.target.first; ++ trg)
	rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
      rule_pair.target += ' ' + static_cast<const std::string&>(nt2.non_terminal(2));
      for (int trg = spans_nt2.target.second; trg != spans.target.second; ++ trg)
	rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
    } else {
      rule_pair.target = static_cast<const std::string&>(lhs);
      for (int trg = spans.target.first; trg != spans_nt2.target.first; ++ trg)
	rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
      rule_pair.target += ' ' + static_cast<const std::string&>(nt2.non_terminal(2));
      for (int trg = spans_nt2.target.second; trg != spans_nt1.target.first; ++ trg)
	rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
      rule_pair.target += ' ' + static_cast<const std::string&>(nt1.non_terminal(1));
      for (int trg = spans_nt1.target.second; trg != spans.target.second; ++ trg)
	rule_pair.target += ' ' + static_cast<const std::string&>(target[trg]);
    }
    
    const int nt1_source_size = spans_nt1.source.second - spans_nt1.source.first;
    const int nt2_source_size = spans_nt2.source.second - spans_nt2.source.first;
    const int nt1_target_size = spans_nt1.target.second - spans_nt1.target.first;
    const int nt2_target_size = spans_nt2.target.second - spans_nt2.target.first;

    rule_pair.alignment.clear();
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
	    
	    rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }
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
		    rule_pair_type& rule_pair,
		    const bool sentential=false)
  {
    typedef std::pair<span_type, symbol_type> span_category_type;
    
    // sort by source-side span...
    std::vector<span_pair_type, std::allocator<span_pair_type> > spans_nt(3);
    spans_nt[0] = __spans_nt1;
    spans_nt[1] = __spans_nt2;
    spans_nt[2] = __spans_nt3;
    
    std::sort(spans_nt.begin(), spans_nt.end(), less_source());
    
    const span_pair_type& spans_nt1 = spans_nt[0];
    const span_pair_type& spans_nt2 = spans_nt[1];
    const span_pair_type& spans_nt3 = spans_nt[2];

    const symbol_type lhs = (sentential ? vocab_type::S : category(spans));
    const symbol_type nt1 = category(spans_nt1);
    const symbol_type nt2 = category(spans_nt2);
    const symbol_type nt3 = category(spans_nt3);
    
    rule_pair.source = static_cast<const std::string&>(lhs);
    for (int src = spans.source.first; src != spans_nt1.source.first; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    rule_pair.source += ' ' + static_cast<const std::string&>(nt1.non_terminal(1));
    for (int src = spans_nt1.source.second; src != spans_nt2.source.first; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    rule_pair.source += ' ' + static_cast<const std::string&>(nt2.non_terminal(2));
    for (int src = spans_nt2.source.second; src != spans_nt3.source.first; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    rule_pair.source += ' ' + static_cast<const std::string&>(nt3.non_terminal(3));
    for (int src = spans_nt3.source.second; src != spans.source.second; ++ src)
      rule_pair.source += ' ' + static_cast<const std::string&>(source[src]);
    
    // sort by target-side span with category,...
    std::vector<span_category_type, std::allocator<span_category_type> > spans_cat(3);
    spans_cat[0] = std::make_pair(spans_nt1.target, nt1.non_terminal(1));
    spans_cat[1] = std::make_pair(spans_nt2.target, nt2.non_terminal(2));
    spans_cat[2] = std::make_pair(spans_nt3.target, nt3.non_terminal(3));
    
    std::sort(spans_cat.begin(), spans_cat.end(), less_first());
    
    rule_pair.target = static_cast<const std::string&>(lhs);
    for (int src = spans.target.first; src != spans_cat[0].first.first; ++ src)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[src]);
    rule_pair.target += ' ' + static_cast<const std::string&>(spans_cat[0].second);
    for (int src = spans_cat[0].first.second; src != spans_cat[1].first.first; ++ src)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[src]);
    rule_pair.target += ' ' + static_cast<const std::string&>(spans_cat[1].second);
    for (int src = spans_cat[1].first.second; src != spans_cat[2].first.first; ++ src)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[src]);
    rule_pair.target += ' ' + static_cast<const std::string&>(spans_cat[2].second);
    for (int src = spans_cat[2].first.second; src != spans.target.second; ++ src)
      rule_pair.target += ' ' + static_cast<const std::string&>(target[src]);
    
    const int nt1_source_size = spans_nt1.source.second - spans_nt1.source.first;
    const int nt2_source_size = spans_nt2.source.second - spans_nt2.source.first;
    const int nt3_source_size = spans_nt3.source.second - spans_nt3.source.first;
    const int nt1_target_size = spans_nt1.target.second - spans_nt1.target.first;
    const int nt2_target_size = spans_nt2.target.second - spans_nt2.target.first;
    const int nt3_target_size = spans_nt3.target.second - spans_nt3.target.first;
    
    rule_pair.alignment.clear();
    for (int src = spans.source.first; src != spans.source.second; ++ src)
      if (is_out_of_span(spans_nt1.source, src) && is_out_of_span(spans_nt2.source, src) && is_out_of_span(spans_nt3.source, src)) {
	point_set_type::const_iterator aiter_begin = alignment_source_target[src].begin();
	point_set_type::const_iterator aiter_end   = alignment_source_target[src].end();
	
	for (point_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  if (is_out_of_span(spans_nt1.target, *aiter) && is_out_of_span(spans_nt2.target, *aiter) && is_out_of_span(spans_nt3.target, *aiter)) {
	    const int mask_source1 = - (src >= spans_nt1.source.second);
	    const int mask_source2 = - (src >= spans_nt2.source.second);
	    const int mask_source3 = - (src >= spans_nt3.source.second);
	    
	    const int mask_target1 = - (*aiter >= spans_nt1.target.second);
	    const int mask_target2 = - (*aiter >= spans_nt2.target.second);
	    const int mask_target3 = - (*aiter >= spans_nt3.target.second);
	    
	    // we shift -1, since we have to take into account the <x1> token... (also for <x2>)
	    const int shift_source = (spans.source.first
				      + (mask_source1 & (nt1_source_size - 1))
				      + (mask_source2 & (nt2_source_size - 1))
				      + (mask_source3 & (nt3_source_size - 1)));
	    const int shift_target = (spans.target.first
				      + (mask_target1 & (nt1_target_size - 1))
				      + (mask_target2 & (nt2_target_size - 1))
				      + (mask_target3 & (nt3_target_size - 1)));
	    
	    rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	  }
      }
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
	alignment_source_target[aiter->target].push_back(aiter->source);
	alignment_target_source[aiter->source].push_back(aiter->target);
      }
    } else {
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
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
       const int max_span,
       const int min_hole,
       const bool ternary,
       const bool sentential,
       const bool inverse,
       const double __max_malloc)
    : queue(__queue),
      output(__output),
      extractor(max_length, max_fertility, max_span, min_hole, ternary, sentential, inverse),
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
      if (rule_pairs.size() >= 1024 || utils::malloc_stats::used() <= malloc_threshold) return;
      
      dump(rule_pairs);
      rule_pairs.clear();
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
	malloc_threshold(__malloc_threshold) {}
    
    const path_type& output;
    path_set_type& paths;
    const size_t malloc_threshold;

    RulePairGenerator generator;
  };

  
  void operator()()
  {
    bitext_type bitext;
    rule_pair_set_type rule_pairs;
    
    Dumper dumper(output, paths, max_malloc * 1024 * 1024 * 1024);
    
    for (;;) {
      queue.pop_swap(bitext);
      
      if (bitext.source.empty()) break;
      
      extractor(bitext.source, bitext.target, bitext.alignment, bitext.spans_source, bitext.spans_target, rule_pairs, dumper);
    }
    
    dumper.dump(rule_pairs);
  }
};

#endif
