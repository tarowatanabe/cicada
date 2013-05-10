//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EXTRACT_PHRASE_IMPL__HPP__
#define __CICADA__EXTRACT_PHRASE_IMPL__HPP__ 1

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

#include "utils/unordered_set.hpp"
#include "utils/unordered_map.hpp"
#include "utils/chart.hpp"
#include "utils/compact_set.hpp"
#include "utils/indexed_set.hpp"

#include <utils/lockfree_list_queue.hpp>
#include <utils/bithack.hpp>
#include <utils/compress_stream.hpp>
#include <utils/tempfile.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/vector2.hpp>
#include <utils/hashmurmur3.hpp>
#include <utils/getline.hpp>

struct Bitext
{
  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;

  sentence_type  source;
  sentence_type  target;
  alignment_type alignment;
  
  Bitext() : source(), target(), alignment() {}
  Bitext(const sentence_type& __source,
	 const sentence_type& __target,
	 const alignment_type& __alignment)
    : source(__source), target(__target), alignment(__alignment) {}
  
  void swap(Bitext& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    alignment.swap(x.alignment);
  }

  void clear()
  {
    source.clear();
    target.clear();
    alignment.clear();
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
	 || iter != end)
	bitext.clear();
    }
    return is;
  }
  
  friend
  std::ostream& operator<<(std::ostream& os, const Bitext& bitext)
  {
    os << bitext.source << " ||| " << bitext.target << " ||| " << bitext.alignment;
    return os;
  }
  
};


struct PhrasePair
{
  typedef std::string phrase_type;
  typedef cicada::Alignment alignment_type;
  typedef int64_t count_type;
  typedef boost::array<count_type, 5+3> counts_type;
  
  const phrase_type*    source;
  const phrase_type*    target;
  const alignment_type* alignment;
  counts_type    counts;

  PhrasePair() : source(0), target(0), alignment(0), counts() {}
  PhrasePair(const phrase_type& __source,
	     const phrase_type& __target,
	     const alignment_type& __alignment,
	     const counts_type& __counts)
    : source(&__source), target(&__target), alignment(&__alignment), counts(__counts) {}
  
  friend
  size_t hash_value(PhrasePair const& x)
  {
    typedef utils::hashmurmur3<size_t> hasher_type;
    
    return hasher_type()(x.source, hasher_type()(x.target, (uintptr_t) x.alignment));
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
    return (PhrasePair::compare(x.source, y.source)
	    || (!PhrasePair::compare(y.source, x.source)
		&& (PhrasePair::compare(x.target, y.target)
		    || (!PhrasePair::compare(y.target, x.target)
			&& PhrasePair::compare(x.alignment, y.alignment)))));
  }
  
  friend
  bool operator>(const PhrasePair& x, const PhrasePair& y)
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

BOOST_FUSION_ADAPT_STRUCT(PhrasePair::alignment_type::point_type,
			  (PhrasePair::alignment_type::index_type, source)
			  (PhrasePair::alignment_type::index_type, target)
			  )

struct PhrasePairGenerator
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
  typedef phrase_pair_type::count_type     count_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  PhrasePairGenerator() {}

  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair)
  {
    typedef std::ostream_iterator<char> iterator_type;
  
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    iterator_type iter(os);
    
    if (! karma::generate(iter,
			  standard::string << " ||| " << standard::string
			  << " ||| " << -((karma::int_ << '-' << karma::int_) % ' ')
			  << " ||| " << (count % ' '),
			  *phrase_pair.source,
			  *phrase_pair.target,
			  *phrase_pair.alignment,
			  phrase_pair.counts))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  boost::spirit::karma::int_generator<count_type> count;
};


namespace std
{
  inline
  void swap(Bitext& x, Bitext& y)
  {
    x.swap(y);
  }
};

struct ExtractPhrase
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  
  typedef std::string phrase_type;
  typedef int64_t count_type;
  typedef boost::array<count_type, 5+3> counts_type;

  typedef std::pair<int, int> span_type;
  struct span_pair_type
  {
    span_type source;
    span_type target;
    
    span_pair_type() : source(), target() {}
    span_pair_type(const span_type& __source, const span_type& __target) : source(__source), target(__target) {}
  };
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  
  struct corner_type
  {
    char left_top;
    char right_top;
    char left_bottom;
    char right_bottom;
    
    corner_type() : left_top(0), right_top(0), left_bottom(0), right_bottom(0) {}
  };
  
  typedef utils::vector2<corner_type, std::allocator<corner_type> > corner_set_type;
  
  typedef PhrasePair phrase_pair_type;

  struct phrase_pair_set_type
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef phrase_pair_type::alignment_type alignment_type;
    
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

    struct phrase_pair_unassigned
    {
      phrase_pair_type operator()() const
      {
	return phrase_pair_type();
      }
    };
    
    typedef utils::compact_set<phrase_pair_type,
			       phrase_pair_unassigned, phrase_pair_unassigned,
			       boost::hash<phrase_pair_type>, std::equal_to<phrase_pair_type>,
			       std::allocator<phrase_pair_type> > pair_set_type;

    typedef pair_set_type::const_iterator const_iterator;
    typedef pair_set_type::iterator       iterator;
    
    bool empty() const { return phrases.empty(); }
    size_type size() const { return phrases.size(); }

    const_iterator begin() const { return phrases.begin(); }
    const_iterator end() const { return phrases.end(); }

    std::pair<iterator, bool> insert(const phrase_pair_type& x)
    {
      return phrases.insert(phrase_pair_type(*sources.insert(*x.source).first,
					     *targets.insert(*x.target).first,
					     *alignments.insert(*x.alignment).first,
					     x.counts));
    }
    
    void clear()
    {
      phrases.clear();
      sources.clear();
      targets.clear();
      alignments.clear();
    }

    void swap(phrase_pair_set_type& x)
    {
      phrases.swap(x.phrases);
      sources.swap(x.sources);
      targets.swap(x.targets);
      alignments.swap(x.alignments);
    }
    
    pair_set_type      phrases;
    phrase_set_type    sources;
    phrase_set_type    targets;
    alignment_set_type alignments;
  };  
  
  struct string_hash : public utils::hashmurmur3<size_t>
  {
    typedef utils::hashmurmur3<size_t> hasher_type;
    size_t operator()(const std::string& x) const
    {
      return hasher_type::operator()(x.begin(), x.end(), 0);
    }
  };
  
  typedef std::pair<const phrase_type, bool> phrase_compact_type;
  
  typedef utils::unordered_map<phrase_type,
			       bool,
			       string_hash,
			       std::equal_to<phrase_type>,
			       std::allocator<std::pair<const phrase_type, bool> > >::type phrase_compact_set_type;
  
  struct PhrasePairCompact
  {
    typedef phrase_pair_type::alignment_type alignment_type;
    
    const phrase_compact_type* source;
    const phrase_compact_type* target;
    alignment_type             alignment;
    counts_type                counts;
    
    PhrasePairCompact() : source(), target(), alignment(), counts() {}
    
    friend
    size_t hash_value(PhrasePairCompact const& x)
    {
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      return hasher_type()(x.target, hasher_type()(x.alignment.begin(), x.alignment.end(), (uintptr_t) x.source));
    }
    
    friend
    bool operator==(const PhrasePairCompact& x, const PhrasePairCompact& y) 
    {
      return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
    }
    
    friend
    bool operator!=(const PhrasePairCompact& x, const PhrasePairCompact& y)
    {
      return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
    }
  };
  
  typedef PhrasePairCompact phrase_pair_compact_type;

  typedef utils::unordered_set<phrase_pair_compact_type, boost::hash<phrase_pair_compact_type>, std::equal_to<phrase_pair_compact_type>,
			       std::allocator<phrase_pair_compact_type> >::type phrase_pair_compact_set_type;


  typedef std::pair<const phrase_compact_type*, const phrase_compact_type*> unique_pair_type;
  
  struct unique_pair_unassigned
  {
    unique_pair_type operator()() const { return unique_pair_type(0, 0); }
  };
  
  typedef utils::compact_set<unique_pair_type,
			     unique_pair_unassigned, unique_pair_unassigned,
			     utils::hashmurmur3<size_t>, std::equal_to<unique_pair_type>,
			     std::allocator<unique_pair_type> > unique_pair_set_type;  
  
  typedef utils::chart<const phrase_compact_type*, std::allocator<const phrase_compact_type*> > phrase_chart_type;
  typedef utils::chart<span_type, std::allocator<span_type> >                                   span_chart_type;
  typedef std::vector<int, std::allocator<int> >                       alignment_count_set_type;
  
  typedef std::vector<int, std::allocator<int> > point_set_type;
  typedef std::vector<point_set_type, std::allocator<point_set_type> > alignment_multiple_type;

  ExtractPhrase(const int __max_length,
		const int __max_fertility,
		const bool __inverse)
    : max_length(__max_length),
      max_fertility(__max_fertility),
      inverse(__inverse) {}
		
  int max_length;
  int max_fertility;
  bool inverse;

  phrase_chart_type phrases_source;
  phrase_chart_type phrases_target;
  
  alignment_multiple_type alignment_source_target;
  alignment_multiple_type alignment_target_source;
  
  span_chart_type span_source_chart;
  span_chart_type span_target_chart;
  
  alignment_count_set_type alignment_count_source;
  alignment_count_set_type alignment_count_target;

  span_pair_set_type spans;
  corner_set_type    corners;
  
  // local structs
  phrase_pair_compact_set_type phrase_pairs_local;
  
  phrase_compact_set_type uniques_source;
  phrase_compact_set_type uniques_target;
  unique_pair_set_type    uniques_pair;

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
    Builder& operator<<(Builder& builder, const sentence_type::value_type& x)
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
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  const alignment_type& alignment,
		  phrase_pair_set_type& phrase_pairs)
  {
    // first, extract spans...
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    phrase_pairs_local.clear();
    uniques_source.clear();
    uniques_target.clear();
    
    phrases_source.clear();
    phrases_target.clear();

    phrases_source.reserve(source_size + 1);
    phrases_target.reserve(target_size + 1);
    
    phrases_source.resize(source_size + 1, 0);
    phrases_target.resize(target_size + 1, 0);
    
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
    
    phrase_pair_compact_type phrase_pair;
    
    // clear span..
    spans.clear();
    
    // clear corners.
    corners.clear();
    corners.resize(source_size + 2, target_size + 2);

    corners(0, 0).right_bottom = true; //BOS
    corners(source_size + 1, target_size + 1).left_top = true; // EOS
    
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
	    
#if 0
	    // unique span-pair
	    const span_type& span_source = span_target_chart(span_target.first, span_target.second);
	    if (span_source.first == source_first && span_source.second == source_last) {
	      // we will assign four-corners of a phrase in the alignment matrix...
	      corners(span_source.first + 1, span_target.first + 1).left_top  = true;
	      corners(span_source.second,    span_target.first + 1).right_top = true;
	      corners(span_source.first + 1, span_target.second).left_bottom  = true;
	      corners(span_source.second,    span_target.second).right_bottom = true;
	    }
#endif
	    
	    // enlarge the target-span...
	    for (int target_first = span_target.first; target_first >= 0; -- target_first)
	      for (int target_last = span_target.second; target_last <= static_cast<int>(target_size); ++ target_last) {
		  
		const int span_count_target = alignment_count_target[target_last] - alignment_count_target[target_first];
		
		if (span_count_source != span_count_target)
		  break;
		
		if (max_length > 0 && source_last - source_first > max_length) continue;
		if (max_length > 0 && target_last - target_first > max_length) continue;
		if (max_fertility > 0
		    && (double(utils::bithack::max(source_last - source_first, target_last - target_first))
			/ double(utils::bithack::min(source_last - source_first, target_last - target_first))) >= max_fertility) continue;
		
		corners(source_first + 1, target_first + 1).left_top  = true;
		corners(source_last,      target_first + 1).right_top = true;
		corners(source_first + 1, target_last).left_bottom  = true;
		corners(source_last,      target_last).right_bottom = true;
		
		spans.push_back(span_pair_type(span_type(source_first, source_last),
					       span_type(target_first, target_last)));
	      }
	  }
	}
      }
    }
    
    span_pair_set_type::const_iterator siter_end = spans.end();
    for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter) {
      const int& source_first = siter->source.first;
      const int& source_last  = siter->source.second;
      const int& target_first = siter->target.first;
      const int& target_last  = siter->target.second;
      
      if (! phrases_source(source_first, source_last)) {

	builder.clear();
	for (int i = source_first; i != source_last - 1; ++ i)
	  builder << source[i] << ' ';
	builder << source[source_last - 1];
	
	phrases_source(source_first, source_last) = &(*uniques_source.insert(phrase_compact_type(builder, false)).first);
      }
      
      if (! phrases_target(target_first, target_last)) {
	
	builder.clear();
	for (int i = target_first; i != target_last - 1; ++ i)
	  builder << target[i] << ' ';
	builder << target[target_last - 1];

	phrases_target(target_first, target_last) = &(*uniques_target.insert(phrase_compact_type(builder, false)).first);
      }
      
      // work with this span!
      phrase_pair.source = phrases_source(source_first, source_last);
      phrase_pair.target = phrases_target(target_first, target_last);
      
      phrase_pair.alignment.clear();
      for (int src = source_first; src != source_last; ++ src) {
	point_set_type::const_iterator titer_end = alignment_source_target[src].end();
	for (point_set_type::const_iterator titer = alignment_source_target[src].begin(); titer != titer_end; ++ titer)
	  phrase_pair.alignment.push_back(std::make_pair(src - source_first, *titer - target_first));
      }
      
      // connected implies we have a "phrase"...corners is one-based!
      const bool connected_left_top     = corners(source_first,    target_first).right_bottom;
      const bool connected_right_top    = corners(source_last + 1, target_first).left_bottom;
      const bool connected_left_bottom  = corners(source_first,    target_last + 1).right_top;
      const bool connected_right_bottom = corners(source_last + 1, target_last + 1).left_top;
      
      phrase_pair_compact_set_type::iterator iter = phrase_pairs_local.find(phrase_pair);
      if (iter == phrase_pairs_local.end())
	iter = phrase_pairs_local.insert(phrase_pair).first;
      
      counts_type& counts = const_cast<counts_type&>(iter->counts);
      
      // monotone/swap counts are computed by "unique" monotone and "unique" swap
      
      counts[0] += 1;
      counts[1] += (  connected_left_top && ! connected_right_top);
      counts[2] += (! connected_left_top &&   connected_right_top);
      counts[3] += (  connected_left_bottom && ! connected_right_bottom);
      counts[4] += (! connected_left_bottom &&   connected_right_bottom);
    }
    
    uniques_pair.clear();
    
    phrase_pair_compact_set_type::const_iterator piter_end = phrase_pairs_local.end();
    for (phrase_pair_compact_set_type::const_iterator piter = phrase_pairs_local.begin(); piter != piter_end; /**/) {
      const bool unique_source = ! piter->source->second;
      const bool unique_target = ! piter->target->second;
      
      const_cast<bool&>(piter->source->second) = true;
      const_cast<bool&>(piter->target->second) = true;
      
      std::pair<phrase_pair_set_type::iterator, bool> result = phrase_pairs.insert(phrase_pair_type(piter->source->first,
												    piter->target->first,
												    piter->alignment,
												    piter->counts));
      
      counts_type& counts = const_cast<counts_type&>(result.first->counts);
      
      if (! result.second)
	std::transform(piter->counts.begin(), piter->counts.end(), counts.begin(), counts.begin(), std::plus<count_type>());
      
      counts[5] += uniques_pair.insert(std::make_pair(piter->source, piter->target)).second;
      counts[6] += unique_source;
      counts[7] += unique_target;
      
      phrase_pairs_local.erase(piter ++);
    }
    
    phrase_pairs_local.clear();
    uniques_source.clear();
    uniques_target.clear();
    uniques_pair.clear();
  }
};


struct Task
{
  typedef boost::filesystem::path path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef Bitext bitext_type;
  
  typedef ExtractPhrase::phrase_pair_type     phrase_pair_type;
  typedef ExtractPhrase::phrase_pair_set_type phrase_pair_set_type;
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> > queue_type;
  
  Task(queue_type& __queue,
       const path_type& __output,
       const int max_length,
       const int max_fertility,
       const bool inverse,
       const double __max_malloc)
    : queue(__queue),
      output(__output),
      extractor(max_length, max_fertility, inverse),
      max_malloc(__max_malloc) {}
  
  queue_type&   queue;
  path_type     output;
  path_set_type paths;
  
  ExtractPhrase extractor;
  PhrasePairGenerator generator;
  
  double max_malloc;
  
  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };
  
  void dump(const phrase_pair_set_type& phrase_pairs)
  {
    typedef std::vector<const phrase_pair_type*, std::allocator<const phrase_pair_type*> > sorted_type;
    
    // sorting...
    sorted_type sorted(phrase_pairs.size());
    {
      sorted_type::iterator siter = sorted.begin();
      phrase_pair_set_type::const_iterator citer_end = phrase_pairs.end();
      for (phrase_pair_set_type::const_iterator citer = phrase_pairs.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    std::sort(sorted.begin(), sorted.end(), less_ptr<phrase_pair_type>());
    
    const path_type path_tmp = utils::tempfile::file_name(output / "counts-XXXXXX");
    utils::tempfile::insert(path_tmp);
    const path_type path = path_tmp.string() + ".gz";
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
    phrase_pair_set_type phrase_pairs;
    
    const int iteration_mask = (1 << 4) - 1;
    const size_t malloc_threshold = size_t(max_malloc * 1024 * 1024 * 1024);
    size_t min_counts_size = 0;
    
    for (int iter = 0;/**/; ++ iter) {
      bitext.clear();
      queue.pop_swap(bitext);
      
      if (bitext.source.empty()) break;
      
      extractor(bitext.source, bitext.target, bitext.alignment, phrase_pairs);
      
      if (((iter & iteration_mask) == iteration_mask)
	  && (phrase_pairs.size() > 1024 * 1024)
	  && (! min_counts_size || phrase_pairs.size() > min_counts_size)
	  && (utils::malloc_stats::used() + phrase_pairs.size() * sizeof(void*) > malloc_threshold)) {
	
	if (! min_counts_size)
	  min_counts_size = phrase_pairs.size() >> 3;
	
	dump(phrase_pairs);
	
	phrase_pairs.clear();
	phrase_pair_set_type(phrase_pairs).swap(phrase_pairs);
      }
    }
    
    if (! phrase_pairs.empty()) {
      dump(phrase_pairs);
      phrase_pairs.clear();
    }
  }
};

#endif
