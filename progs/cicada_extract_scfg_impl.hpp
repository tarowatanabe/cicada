#ifndef __CICADA__EXTRACT_SCFG_IMPL__HPP__
#define __CICADA__EXTRACT_SCFG_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>

#include <string>
#include <vector>

#include <boost/array.hpp>

#include "cicada/sentence.hpp"
#include "cicada/alignment.hpp"
#include "cicada/span_vector.hpp"

#include "utils/sgi_hash_set.hpp"
#include "utils/chart.hpp"

#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>

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
    
    using qi::phrase_parse;
    using qi::lit;
    using standard::space;

    bitext.clear();
    std::string line;
    if (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end  = line.end();
      
      if((!bitext.source.assign(iter, end))
	 || (!phrase_parse(iter, end, lit("|||"), space))
	 || (!bitext.target.assign(iter, end))
	 || (!phrase_parse(iter, end, lit("|||"), space))
	 || (!bitext.alignment.assign(iter, end))
	 || (!phrase_parse(iter, end, lit("|||"), space))
	 || (!bitext.spans_source.assign(iter, end))
	 || (!phrase_parse(iter, end, lit("|||"), space))
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
  typedef boost::array<double, 5> counts_type;
  
  phrase_type    source;
  phrase_type    target;
  alignment_type alignment;
  counts_type    counts;

  RulePair() : source(), target(), alignment(), counts() {}

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
			  (RulePair::counts_type, counts)
			  )

struct RulePairGenerator
{
  typedef RulePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  RulePairGenerator() : grammar() {}
  RulePairGenerator(const RulePairGenerator& x) : grammar() {}

  
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
      using karma::double_;
      using karma::int_;
      using standard::space;
      
      phrase %= +char_;
      alignment %= -((int_ << '-' << int_) % ' ');
      counts %= double20 % ' ';
      phrase_pair %= phrase << " ||| " << phrase << " ||| " << alignment << " ||| " << counts;
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
    boost::spirit::karma::rule<Iterator, alignment_type()> alignment;
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

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  
  typedef std::string phrase_type;
  typedef boost::array<double, 5> counts_type;

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
      return x.source == y.source && x.taret == y.target;
    }
    
    friend
    bool operator!=(const span_pair_type& x, const span_pair_type& y)
    {
      return x.source != y.source || x.taret != y.target;
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
  

#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
				  std::allocator<rule_pair_type> > rule_pair_set_type;
#else
  typedef sgi::hash_set<rule_pair_type, boost::hash<rule_pair_type>, std::equal_to<rule_pair_type>,
			std::allocator<rule_pair_type> > rule_pair_set_type;
#endif
  
  typedef utils::chart<phrase_type, std::allocator<phrase_type> >      phrase_chart_type;
  typedef utils::chart<span_type, std::allocator<span_type> >          span_chart_type;
  typedef std::vector<int, std::allocator<int> >                       alignment_count_set_type;
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  
  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > alignment_multiple_type;

  ExtractSCFG(const int __max_length,
	      const int __max_fertility,
	      const int __max_span,
	      const int __min_hole)
    : max_length(__max_length),
      max_fertility(__max_fertility),
      max_span(__max_span),
      min_hole(__min_hole) {}
		
  int max_length;
  int max_fertility;
  int max_span;
  int min_hole;

  phrase_chart_type phrases_source;
  phrase_chart_type phrases_target;
  
  alignment_multiple_type alignment_source_target;
  alignment_multiple_type alignment_target_source;
  
  span_chart_type span_source_chart;
  span_chart_type span_target_chart;
  
  alignment_count_set_type alignment_count_source;
  alignment_count_set_type alignment_count_target;
  
  span_pair_set_type        spans;
  span_pair_set_type        spans_unique;
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  const alignment_type& alignment,
		  rule_pair_set_type& rule_pairs)
  {
    // first, extract spans...
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    phrases_source.clear();
    phrases_target.clear();
    
    phrases_source.resize(source_size + 1);
    phrases_target.resize(target_size + 1);
    
    alignment_source_target.clear();
    alignment_target_source.clear();
    alignment_source_target.resize(source_size);
    alignment_target_source.resize(target_size);
    
    alignment_count_source.clear();
    alignment_count_target.clear();
    alignment_count_source.resize(source_size + 1, 0);
    alignment_count_target.resize(target_size + 1, 0);
    
    spans.clear();
    spans_unique.clear();
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      alignment_source_target[aiter->source].push_back(aiter->target);
      alignment_target_source[aiter->target].push_back(aiter->source);
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
    
    // now, iterate again...
    span_pair_set_type::const_iterator iter_end = spans.end();
    for (span_pair_set_type::const_iterator iter = spans.begin(); iter != iter_end; ++ iter) {
      
      const int source_length = iter->source.second - iter->source.first;
      const int target_length = iter->target.second - iter->target.first;

      const int source_count = alignment_count_source[iter->source.second] - alignment_count_source[iter->source.first];
      const int target_count = alignment_count_target[iter->target.second] - alignment_count_target[iter->target.first];
      
      if (max_length <= 0 || source_length <= max_length)
	if (max_fertility <= 0 || fertility(source_length, target_length) < max_fertility) {
	  // extract rule...
	  
	}
      
      // consider hole!
      if (max_span <= 0 || source_length <= max_span) {
	span_pair_set_type::const_iterator niter_end = spans_unique.end();
	
	// first non-terminal...
	for (span_pair_set_type::const_iterator niter1 = spans_unique.begin(); niter1 != niter_end; ++ niter1) 
	  if (*iter != *niter1
	      && is_parent(iter->source, niter1->source)
	      && is_parent(iter->target, niter1->target)) {
	    
	    const int source_count1 = alignment_count_source[niter1->source.second] - alignment_count_source[niter1->source.first];
	    const int target_count1 = alignment_count_target[niter1->target.second] - alignment_count_target[niter1->target.first];
	    
	    if (source_count1 == source_count) continue;
	    
	    const int source_length1 = source_length - (niter1->source.second - niter1->source.first);
	    const int target_length1 = target_length - (niter1->target.second - niter1->target.first);
	    
	    if (max_length <= 0 || source_length1 <= max_length)
	      if (max_fertility <= 0 || fertility(source_length1, target_length1) < max_fertility) {
		// extract rule...
		
		
	      }
	    
	    // second non-terminal...
	    for (span_pair_set_type::const_iterator niter2 = niter1 + 1; niter2 != niter_end; ++ niter2) 
	      if (*iter != *niter2
		  && is_parent(iter->source, niter2->source)
		  && is_parent(iter->target, niter2->target)
		  && is_disjoint(niter1->source, niter2->source)
		  && is_disjoint(niter1->taregt, niter2->target)) {
		
		const int source_count2 = alignment_count_source[niter2->source.second] - alignment_count_source[niter2->source.first];
		const int target_count2 = alignment_count_target[niter2->target.second] - alignment_count_target[niter2->target.first];
		
		if (source_count == source_count1 + source_count2) continue;
		
		const int source_length2 = source_length1 - (niter2->source.second - niter2->source.first);
		const int target_length2 = target_length1 - (niter2->target.second - niter2->target.first);
		
		if (max_length <= 0 || source_length2 <= max_length)
		  if (max_fertility <= 0 || fertility(source_length2, target_length2) < max_fertility) {
		    // extract rule...
		    
		    
		  }
	      }
	  }
      }
    }
  }
  
  bool is_parent(const span_type& parent, const span_type& child)
  {
    return parent.first <= child.first && child.second <= parent.second;
  }
  
  bool is_disjoint(const span_type& child1, const span_type& child2)
  {
    return child1.second <= child2.first || child2.second <= child1.first;
  }
  
  double fertility(int length1, int length2)
  {
    return double(utils::bithack::max(length1, length2)) / utils::bithacl::min(length1, length2);
  }
  
  bool is_aligned(const int source, const int target) const
  {
    const int source_size = alignment_source_target.size();
    const int target_size = alignment_target_source.size();
    
    if (source == -1 && target == -1) return true; // aligned at BOS
    if (source <= -1 || target <= -1) return false;
    if (source == source_size && target == target_size) return true; // aligned at EOS
    if (source >= source_size || target >= target_size) return false;
    
    point_set_type::const_iterator aiter_begin = alignment_source_target[source].begin();
    point_set_type::const_iterator aiter_end   = alignment_source_target[source].end();
    
    // check if there exists alignment point!
    return std::find(aiter_begin, aiter_end, target) != aiter_end;
  }
  
};


struct Task
{
  typedef boost::filesystem::path path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef Bitext bitext_type;
  
  typedef ExtractPhrase::rule_pair_type     rule_pair_type;
  typedef ExtractPhrase::rule_pair_set_type rule_pair_set_type;
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> > queue_type;
  
  Task(queue_type& __queue,
       const path_type& __output,
       const int max_length,
       const int max_fertility,
       const double __max_malloc)
    : queue(__queue),
      output(__output),
      extractor(max_length, max_fertility),
      max_malloc(__max_malloc) {}
  
  queue_type&   queue;
  path_type     output;
  path_set_type paths;
  
  ExtractPhrase extractor;
  RulePairGenerator generator;
  
  double max_malloc;
  
  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };
  
  void dump(const rule_pair_set_type& rule_pairs)
  {
    typedef std::vector<const rule_pair_type*, std::allocator<const rule_pair_type*> > sorted_type;
    
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
    const path_type path = path_tmp.file_string() + ".gz";
    utils::tempfile::insert(path);
    
    utils::compress_ostream os(path, 1024 * 1024);
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      generator(os, *(*siter)) << '\n';
    
    paths.push_back(path);
  }
  
  void operator()()
  {
    bitext_type bitext;
    rule_pair_set_type rule_pairs;
    
    const int iteration_mask = (1 << 10) - 1;
    
    for (int iter = 0;/**/; ++ iter) {
      queue.pop_swap(bitext);
      if (bitext.source.empty()) break;
      
      extractor(bitext.source, bitext.target, bitext.alignment, rule_pairs);
      
      if (((iter & iteration_mask) == iteration_mask) && (utils::malloc_stats::used() > size_t(max_malloc * 1024 * 1024 * 1024))) {
	dump(rule_pairs);
	rule_pairs.clear();
      }
    }
    
    if (! rule_pairs.empty()) {
      dump(rule_pairs);
      rule_pairs.clear();
    }
  }
};

#endif
