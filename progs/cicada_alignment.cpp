//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// cicada alignment tool:
//
// we can produce intersection/union/grow-{,diag}-{final,final-and}/source/target/itg/max-match from GIZA++ alingment
//        invert alignment
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_object.hpp>

#include <boost/spirit/include/support_istream_iterator.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

//#include <boost/numeric/ublas/matrix.hpp>

#include <stdexcept>
#include <iostream>
#include <string>
#include <vector>
#include <queue>
#include <set>

#include <cicada/alignment.hpp>

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/bithack.hpp"
#include "utils/mathop.hpp"
#include "utils/vector2.hpp"

#include "kuhn_munkres.hpp"
#include "itg_prune.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;

typedef cicada::Alignment alignment_type;
typedef alignment_type::point_type point_type;

struct BitextGiza
{
  typedef std::set<int, std::less<int>, std::allocator<int> > point_set_type;
  typedef std::pair<std::string, point_set_type> word_align_type;
  typedef std::vector<word_align_type, std::allocator<word_align_type> > word_align_set_type;
  typedef std::vector<std::string, std::allocator<std::string> > sent_type;
  
  std::string         comment;
  sent_type           target;
  word_align_set_type source;
  
  BitextGiza() {}

  void clear()
  {
    comment.clear();
    target.clear();
    source.clear();
  }

  void swap(BitextGiza& x)
  {
    comment.swap(x.comment);
    target.swap(x.target);
    source.swap(x.source);
  }
};

namespace std
{
  inline
  void swap(BitextGiza& x, BitextGiza& y)
  {
    x.swap(y);
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  BitextGiza,
			  (std::string,                     comment)
			  (BitextGiza::sent_type,           target)
			  (BitextGiza::word_align_set_type, source)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  cicada::Alignment::point_type,
			  (cicada::Alignment::index_type, source)
			  (cicada::Alignment::index_type, target)
			  )

typedef BitextGiza bitext_giza_type;


path_type source_target_file;
path_type target_source_file;
path_type span_source_file;
path_type span_target_file;
path_type input_file;
path_type output_file = "-";

bool source_target_mode = false;
bool target_source_mode = false;

bool itg_mode = false;
bool max_match_mode = false;

bool intersection_mode = false;
bool union_mode = false;
bool grow_mode = false;
bool final_mode = false;
bool diag_mode = false;
bool final_and_mode = false;
bool invert_mode = false;
bool moses_mode = false;

double prob_null         = 0.01;
double prob_union        = 0.5;
double prob_intersection = 1.0;

double score_null = 0.0;
double score_union = 0.0;
double score_intersection = 0.0;

int threads = 1;
int debug = 0;

void process_giza(std::istream& is_src_trg, std::istream& is_trg_src, std::ostream& os);
void process_giza(std::istream& is, std::ostream& os);
void process_alignment(std::istream& is, std::ostream& os);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);

    score_null         = utils::mathop::log(prob_null);
    score_union        = utils::mathop::log(prob_union);
    score_intersection = utils::mathop::log(prob_intersection);
    
    const bool flush_output = (output_file == "-"
			       || (boost::filesystem::exists(output_file)
				   && ! boost::filesystem::is_regular_file(output_file)));
    
    if (source_target_mode || target_source_mode) {
      if (source_target_mode && target_source_mode)
	throw std::runtime_error("which mode f2e or e2f?");
      
      if (source_target_mode) {
	if (source_target_file != "-" && ! boost::filesystem::exists(source_target_file))
	  throw std::runtime_error("no f2e file?" + source_target_file.file_string());
	
	utils::compress_istream is(source_target_file, 1024 * 1024);
	utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
	
	process_giza(is, os);
      } else {
	if (target_source_file != "-" && ! boost::filesystem::exists(target_source_file))
	  throw std::runtime_error("no e2f file?" + target_source_file.file_string());

	utils::compress_istream is(target_source_file, 1024 * 1024);
	utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
	
	process_giza(is, os);
      }
    } if (input_file == "-" || boost::filesystem::exists(input_file)) {
      utils::compress_istream is(input_file, 1024 * 1024);
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
      
      process_alignment(is, os);
    } else {
      if (source_target_file != "-" && ! boost::filesystem::exists(source_target_file))
	throw std::runtime_error("no f2e file?" + source_target_file.file_string());
      if (target_source_file != "-" && ! boost::filesystem::exists(target_source_file))
	throw std::runtime_error("no e2f file?" + target_source_file.file_string());
      
      utils::compress_istream is_src_trg(source_target_file, 1024 * 1024);
      utils::compress_istream is_trg_src(target_source_file, 1024 * 1024);
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
      
      process_giza(is_src_trg, is_trg_src, os);
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

template <typename Iterator>
struct alignment_parser : boost::spirit::qi::grammar<Iterator, alignment_type(), boost::spirit::standard::blank_type>
{
  alignment_parser() : alignment_parser::base_type(alignment)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    point     %= qi::int_ >> '-' >> qi::int_;
    alignment %= *point >> (qi::eol | qi::eoi);
  };
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, alignment_type::point_type(), blank_type> point;
  boost::spirit::qi::rule<Iterator, alignment_type(), blank_type>             alignment;
};

template <typename Iterator>
struct bitext_giza_parser : boost::spirit::qi::grammar<Iterator, bitext_giza_type(), boost::spirit::standard::blank_type>
{
  bitext_giza_parser() : bitext_giza_parser::base_type(bitext)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    comment %= qi::lexeme[*(standard::char_ - qi::eol)] >> qi::eol;
    target  %= *qi::lexeme[+(standard::char_ - standard::space - qi::eol)] >> qi::eol;
    
    points     %= "({" >> *qi::int_ >> "})";
    word_align %= qi::lexeme[+(standard::char_ - standard::space - qi::eol)] >> points;
    source     %= *word_align >> (qi::eol | qi::eoi);
    
    bitext %= comment >> target >> source;
  }
  
  typedef boost::spirit::standard::blank_type blank_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), blank_type>                           comment;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::sent_type(), blank_type>           target;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::point_set_type(), blank_type>      points;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::word_align_type(), blank_type>     word_align;
  boost::spirit::qi::rule<Iterator, bitext_giza_type::word_align_set_type(), blank_type> source;
  
  boost::spirit::qi::rule<Iterator, bitext_giza_type(), blank_type> bitext;
};

struct AlignmentInserter
{
  AlignmentInserter(alignment_type& __align) : align(__align) {}

  void insert(const point_type& point)
  {
    align.push_back(point);
  }
  
  alignment_type& align;
};

struct Invert
{
  template <typename Alignment>
  void operator()(const alignment_type& align, Alignment& inverted)
  {
    alignment_type::const_iterator aiter_end = align.end();
    for (alignment_type::const_iterator aiter = align.begin(); aiter != aiter_end; ++ aiter)
      inverted.insert(std::make_pair(aiter->target, aiter->source));
  }
};

struct SourceTarget
{
  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext, Alignment& align)
  {
    for (int src = 1; src < static_cast<int>(bitext.source.size()); ++ src) {
      const bitext_giza_type::point_set_type& aligns = bitext.source[src].second;
      
      bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) 
	align.insert(std::make_pair(src - 1, *titer - 1));
    }
  }
};

struct TargetSource
{
  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext, Alignment& align)
  {
    for (int src = 1; src < static_cast<int>(bitext.source.size()); ++ src) {
      const bitext_giza_type::point_set_type& aligns = bitext.source[src].second;
      
      bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) 
	align.insert(std::make_pair(*titer - 1, src - 1));
    }
  }
};

struct Intersect
{
  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    // - 1 for NULL
    const int source_size = utils::bithack::min(bitext_source_target.source.size() - 1, bitext_target_source.target.size());
    const int target_size = utils::bithack::min(bitext_source_target.target.size(),     bitext_target_source.source.size() - 1);
    
    for (int src = 1; src <= source_size; ++ src) {
      const bitext_giza_type::point_set_type& aligns = bitext_source_target.source[src].second;
      
      bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) 
	if (*titer < static_cast<int>(bitext_target_source.source.size())) {
	  const bitext_giza_type::point_set_type& invert = bitext_target_source.source[*titer].second;
	  
	  if (invert.find(src) != invert.end())
	    align.insert(std::make_pair(src - 1, *titer - 1));
	}
    }
  }
};

struct Union
{
  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    for (int src = 1; src < static_cast<int>(bitext_source_target.source.size()); ++ src) {
      const bitext_giza_type::point_set_type& aligns = bitext_source_target.source[src].second;
      
      bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) 
	align.insert(std::make_pair(src - 1, *titer - 1));
    }
    
    for (int trg = 1; trg < static_cast<int>(bitext_target_source.source.size()); ++ trg) {
      const bitext_giza_type::point_set_type& aligns = bitext_target_source.source[trg].second;
      
      bitext_giza_type::point_set_type::const_iterator siter_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator siter = aligns.begin(); siter != siter_end; ++ siter) 
	align.insert(std::make_pair(*siter - 1, trg - 1));
    }
  }
};

struct __Grow
{
  
  std::vector<bool, std::allocator<bool> > aligned_source;
  std::vector<bool, std::allocator<bool> > aligned_target;
  
  template <typename Alignment, size_t N>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align,
		  const point_type (&neighbours)[N])
  {
    const int source_size = utils::bithack::max(bitext_source_target.source.size() - 1, bitext_target_source.target.size());
    const int target_size = utils::bithack::max(bitext_source_target.target.size(),     bitext_target_source.source.size() - 1);

    aligned_source.clear();
    aligned_target.clear();
    
    aligned_source.resize(source_size, false);
    aligned_target.resize(target_size, false);
    
    typename Alignment::const_iterator aiter_end = align.end();
    for (typename Alignment::const_iterator aiter = align.begin(); aiter != aiter_end; ++ aiter) {
      aligned_source[aiter->source] = true;
      aligned_target[aiter->target] = true;
    }
    
    Alignment align_new;
    bool added = true;
    do {
      added = false;
      align_new.clear();
      
      typename Alignment::const_iterator aiter_end = align.end();
      for (typename Alignment::const_iterator aiter = align.begin(); aiter != aiter_end; ++ aiter) {
	for (size_t n = 0; n != N; ++ n) {
	  point_type point = *aiter;
	  point.source += neighbours[n].source;
	  point.target += neighbours[n].target;
	  
	  if (0 <= point.source && point.source < source_size && 0 <= point.target && point.target < target_size) 
	    if ((! aligned_source[point.source]) || (! aligned_target[point.target]))
	      if (has_alignment(bitext_source_target, point.source, point.target)
		  || has_alignment(bitext_target_source, point.target, point.source)) {
		align_new.insert(point);
		
		aligned_source[point.source] = true;
		aligned_target[point.target] = true;
		
		added = true;
	      }
	}
      }
      
      align.insert(align_new.begin(), align_new.end());
      
    } while (added);
  }

  
  bool has_alignment(const bitext_giza_type& bitext, const int source, const int target)
  {
    if (source + 1 >= static_cast<int>(bitext.source.size())) return false;
    
    const bitext_giza_type::point_set_type& aligns = bitext.source[source + 1].second;
    
    return aligns.find(target + 1) != aligns.end();
  }
};

struct Grow
{
  static const point_type points[4];
  
  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    grow(bitext_source_target, bitext_target_source, align, points);
  }
  
  __Grow grow;
};

struct GrowDiag
{
  static const point_type points[8];
  
  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    grow(bitext_source_target, bitext_target_source, align, points);
  }
  
  __Grow grow;
};

const point_type Grow::points[4] = {point_type(-1, 0), point_type(0, -1), point_type(1, 0), point_type(0, 1)};
const point_type GrowDiag::points[8] = {point_type(-1, 0), point_type(0, -1), point_type(1, 0), point_type(0, 1),
					point_type(-1, -1), point_type(-1, 1), point_type(1, -1), point_type(1, 1)};

struct __Final
{
  std::vector<bool, std::allocator<bool> > aligned_source;
  std::vector<bool, std::allocator<bool> > aligned_target;
  
  template <typename Alignment, typename Filter>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align,
		  Filter filter)
  {
    const int source_size = utils::bithack::max(bitext_source_target.source.size() - 1, bitext_target_source.target.size());
    const int target_size = utils::bithack::max(bitext_source_target.target.size(),     bitext_target_source.source.size() - 1);

    aligned_source.clear();
    aligned_target.clear();
    
    aligned_source.resize(source_size, false);
    aligned_target.resize(target_size, false);
    
    typename Alignment::const_iterator aiter_end = align.end();
    for (typename Alignment::const_iterator aiter = align.begin(); aiter != aiter_end; ++ aiter) {
      aligned_source[aiter->source] = true;
      aligned_target[aiter->target] = true;
    }
    
    // grow finally, but not extending from previously infered alignments...
    for (int src = 1; src < static_cast<int>(bitext_source_target.source.size()); ++ src) {
      const bitext_giza_type::point_set_type& aligns =  bitext_source_target.source[src].second;
      
      bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) {
	const int trg = *titer;
	
	if (filter(aligned_source[src - 1], aligned_target[trg - 1])) continue;
	
	align.insert(point_type(src - 1, trg - 1));
	aligned_source[src - 1] = true;
	aligned_target[trg - 1] = true;
      }
    }
    
    for (int trg = 1; trg < static_cast<int>(bitext_target_source.source.size()); ++ trg) {
      const bitext_giza_type::point_set_type& aligns =  bitext_target_source.source[trg].second;
      
      bitext_giza_type::point_set_type::const_iterator siter_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator siter = aligns.begin(); siter != siter_end; ++ siter) {
	const int src = *siter;
	
	if (filter(aligned_source[src - 1], aligned_target[trg - 1])) continue;
	
	align.insert(point_type(src - 1, trg - 1));
	aligned_source[src - 1] = true;
	aligned_target[trg - 1] = true;
      }
    }
  }
};

struct Final
{
  struct Filter
  {
    bool operator()(const bool source, const bool target) const
    {
      return source && target;
    }
  };

  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    final(bitext_source_target, bitext_target_source, align, Filter());
  }

  __Final final;
};

struct FinalAnd
{

  struct Filter
  {
    bool operator()(const bool source, const bool target) const
    {
      return source || target;
    }
  };

  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    final(bitext_source_target, bitext_target_source, align, Filter());
  }
  
  __Final final;
};


struct MaxMatch
{
  //typedef boost::numeric::ublas::matrix<double> matrix_type;
  
  typedef utils::vector2<double, std::allocator<double> > matrix_type;
  typedef utils::vector2<char, std::allocator<char> > assigned_type;

  template <typename Alignment>
  class insert_align
  {
    typedef Alignment align_type;

    int source_size;
    int target_size;
      
    align_type*   align;
      
  public:
    insert_align(const int& _source_size,
		 const int& _target_size,
		 align_type& __align)
      : source_size(_source_size), target_size(_target_size),
	align(&__align) {}
    
    template <typename Edge>
    insert_align& operator=(const Edge& edge)
    {	
      if (edge.first < source_size && edge.second < target_size)
	align->insert(edge);
	
      return *this;
    }
      
    insert_align& operator*() { return *this; }
    insert_align& operator++() { return *this; }
    insert_align operator++(int) { return *this; }
  };


  template <typename Alignment>
  void operator()(const bitext_giza_type& bitext_source_target,
		  const bitext_giza_type& bitext_target_source,
		  Alignment& align)
  {
    const int source_size = utils::bithack::max(bitext_source_target.source.size() - 1, bitext_target_source.target.size());
    const int target_size = utils::bithack::max(bitext_source_target.target.size(),     bitext_target_source.source.size() - 1);
    
    costs.clear();
    costs.resize(source_size + target_size, source_size + target_size, 0.0);
    
    assigned.clear();
    assigned.resize(source_size, target_size, false);
    
    for (int src = 1; src < static_cast<int>(bitext_source_target.source.size()); ++ src) {
      const bitext_giza_type::point_set_type& aligns = bitext_source_target.source[src].second;
      
      bitext_giza_type::point_set_type::const_iterator titer_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator titer = aligns.begin(); titer != titer_end; ++ titer) {
	costs(src - 1, *titer - 1) = score_union;
	assigned(src - 1, *titer - 1) = true;
      }
    }
    
    for (int trg = 1; trg < static_cast<int>(bitext_target_source.source.size()); ++ trg) {
      const bitext_giza_type::point_set_type& aligns = bitext_target_source.source[trg].second;
      
      bitext_giza_type::point_set_type::const_iterator siter_end = aligns.end();
      for (bitext_giza_type::point_set_type::const_iterator siter = aligns.begin(); siter != siter_end; ++ siter)
	costs(*siter - 1, trg - 1) = (assigned(*siter - 1, trg - 1) ? score_intersection : score_union);
    }
    
    for (int trg = 0; trg < target_size; ++ trg)
      for (int src = 0; src < source_size; ++ src) {
	// NULL network...
	costs(src, trg + source_size) = score_null;
	costs(src + target_size, trg) = score_null;
      }
    
    kuhn_munkres_assignment(costs, insert_align<Alignment>(source_size, target_size, align));
  }
  
  matrix_type   costs;
  assigned_type assigned;
};


void process_alignment(std::istream& is, std::ostream& os)
{
  alignment_type align;
  alignment_type inverted;
  AlignmentInserter inserter(inverted);
  Invert invert;
  
  if (invert_mode) {
    while (is >> align) {
      inverted.clear();
      invert(align, inserter);
      std::sort(inverted.begin(), inverted.end());
      os << inverted << '\n';
    }
  } else {
    while (is >> align) {
      std::sort(align.begin(), align.end());
      os << align << '\n';
    }
  }
}

inline
void alignment_to_giza(const alignment_type& alignment, bitext_giza_type& bitext)
{
  bitext.clear();
  
  alignment_type::const_iterator aiter_end = alignment.end();
  for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
    bitext.target.resize(utils::bithack::max(static_cast<size_t>(aiter->target + 1), bitext.target.size()));
    bitext.source.resize(utils::bithack::max(static_cast<size_t>(aiter->source + 2), bitext.source.size()));
    bitext.source[aiter->source + 1].second.insert(aiter->target + 1);
  }
}

void process_giza(std::istream& is, std::ostream& os)
{
  typedef boost::spirit::istream_iterator iter_type;
  
  alignment_parser<iter_type> parser_alignment;
  bitext_giza_parser<iter_type> parser_giza;
  
  bitext_giza_type bitext_giza;
  
  is.unsetf(std::ios::skipws);
  
  iter_type iter(is);
  iter_type iter_end;
  
  alignment_type    alignment;
  alignment_type    inverted;
  
  AlignmentInserter inserter(alignment);
  AlignmentInserter inverted_inserter(inverted);
  
  if (source_target_mode) {
    SourceTarget process;
    Invert       invert;
    
    while (iter != iter_end) {
      bitext_giza.clear();
      alignment.clear();
      
      if (moses_mode) {
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser_alignment, boost::spirit::standard::blank, alignment))
	  if (iter != iter_end)
	    throw std::runtime_error("parsing failed");
	
	alignment_to_giza(alignment, bitext_giza);
	alignment.clear();
      } else {
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser_giza, boost::spirit::standard::blank, bitext_giza))
	  if (iter != iter_end)
	    throw std::runtime_error("parsing failed");
      }
      
      process(bitext_giza, inserter);

      if (invert_mode) {
	inverted.clear();
	invert(alignment, inverted_inserter);
	std::sort(inverted.begin(), inverted.end());
	alignment.swap(inverted);
      }
    
      os << alignment << '\n';
    }
  } else {
    TargetSource process;
    Invert       invert;
  
    while (iter != iter_end) {
      bitext_giza.clear();
      alignment.clear();
      
      if (moses_mode) {
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser_alignment, boost::spirit::standard::blank, alignment))
	  if (iter != iter_end)
	    throw std::runtime_error("parsing failed");
	
	alignment_to_giza(alignment, bitext_giza);
	alignment.clear();
      } else {
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser_giza, boost::spirit::standard::blank, bitext_giza))
	  if (iter != iter_end)
	    throw std::runtime_error("parsing failed");
      }

      process(bitext_giza, inserter);

      if (invert_mode) {
	inverted.clear();
	invert(alignment, inverted_inserter);
	std::sort(inverted.begin(), inverted.end());
	alignment.swap(inverted);
      }
      
      os << alignment << '\n';
    }
  }
}

struct Task
{
  
  void operator()()
  {
    
    
  }
  
};

void process_giza(std::istream& is_src_trg, std::istream& is_trg_src, std::ostream& os)
{
  typedef std::set<point_type, std::less<point_type>, std::allocator<point_type> > align_set_type;

  typedef boost::spirit::istream_iterator iter_type;
  
  alignment_parser<iter_type> parser_alignment;
  bitext_giza_parser<iter_type> parser_giza;
  
  bitext_giza_type bitext_source_target;
  bitext_giza_type bitext_target_source;
  
  is_src_trg.unsetf(std::ios::skipws);
  is_trg_src.unsetf(std::ios::skipws);
  
  iter_type siter(is_src_trg);
  iter_type siter_end;

  iter_type titer(is_trg_src);
  iter_type titer_end;
  
  align_set_type aligns;
  alignment_type alignment;
  alignment_type inverted;
  AlignmentInserter inserter(alignment);
  AlignmentInserter inverted_inserter(inverted);
  
  MaxMatch  __max_match;
  Intersect __intersect;
  Union     __union;
  Grow      __grow;
  GrowDiag  __grow_diag;
  Final     __final;
  FinalAnd  __final_and;
  
  Invert    __invert;

  while (siter != siter_end && titer != titer_end) {
    bitext_source_target.clear();
    bitext_target_source.clear();
    alignment.clear();

    if (moses_mode) {
      if (! boost::spirit::qi::phrase_parse(siter, siter_end, parser_alignment, boost::spirit::standard::blank, alignment))
	if (siter != siter_end)
	  throw std::runtime_error("source-target parsing failed");

      alignment_to_giza(alignment, bitext_source_target);
      alignment.clear();
      
      if (! boost::spirit::qi::phrase_parse(titer, titer_end, parser_alignment, boost::spirit::standard::blank, alignment))
	if (titer != titer_end)
	  throw std::runtime_error("target-source parsing failed");

      alignment_to_giza(alignment, bitext_target_source);
      alignment.clear();
    } else {
      if (! boost::spirit::qi::phrase_parse(siter, siter_end, parser_giza, boost::spirit::standard::blank, bitext_source_target))
	if (siter != siter_end)
	  throw std::runtime_error("source-target parsing failed");
      
      if (! boost::spirit::qi::phrase_parse(titer, titer_end, parser_giza, boost::spirit::standard::blank, bitext_target_source))
	if (titer != titer_end)
	  throw std::runtime_error("target-source parsing failed");
    }
    
    aligns.clear();
    
    if (max_match_mode) {
      __max_match(bitext_source_target, bitext_target_source, aligns);
      
      alignment.insert(alignment.end(), aligns.begin(), aligns.end());
    } else if (intersection_mode) 
      __intersect(bitext_source_target, bitext_target_source, inserter);
    else if (union_mode) {
      __union(bitext_source_target, bitext_target_source, aligns);
      
      alignment.insert(alignment.end(), aligns.begin(), aligns.end());
    } else {
      // first, compute intersection
      __intersect(bitext_source_target, bitext_target_source, aligns);
      
      // grow...
      if (grow_mode) {
	if (diag_mode)
	  __grow_diag(bitext_source_target, bitext_target_source, aligns);
	else
	  __grow(bitext_source_target, bitext_target_source, aligns);
      }
      
      // final...
      if (final_mode)
	__final(bitext_source_target, bitext_target_source, aligns);
      else if (final_and_mode)
	__final_and(bitext_source_target, bitext_target_source, aligns);
      
      alignment.insert(alignment.end(), aligns.begin(), aligns.end());
    }
    
    // invert this alignment...
    if (invert_mode) {
      inverted.clear();
      __invert(alignment, inverted_inserter);
      std::sort(inverted.begin(), inverted.end());
      alignment.swap(inverted);
    }
    
    os << alignment << '\n';
  }
  
  if (siter != siter_end || titer != titer_end)
    throw std::runtime_error("# of samples do not match");
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source-target", po::value<path_type>(&source_target_file), "P(target | source) viterbi output")
    ("target-source", po::value<path_type>(&target_source_file), "P(source | target) viterbi output")
    ("span-source",   po::value<path_type>(&span_source_file),   "source span data")
    ("span-target",   po::value<path_type>(&span_target_file),   "target span data")
    ("input",         po::value<path_type>(&input_file),                      "input alignement")
    ("output",        po::value<path_type>(&output_file)->default_value("-"), "output alignment")
    
    ("f2e", po::bool_switch(&source_target_mode), "source target")
    ("e2f", po::bool_switch(&target_source_mode), "target source")
    
    ("itg",          po::bool_switch(&itg_mode),          "itg")
    ("max-match",    po::bool_switch(&max_match_mode),    "max-match")
    ("intersection", po::bool_switch(&intersection_mode), "intersection")
    ("union",        po::bool_switch(&union_mode),        "union")
    ("grow",         po::bool_switch(&grow_mode),         "grow")
    ("diag",         po::bool_switch(&diag_mode),         "diag")
    ("final",        po::bool_switch(&final_mode),        "final")
    ("final-and",    po::bool_switch(&final_and_mode),    "final-and")
    ("invert",       po::bool_switch(&invert_mode),       "invert alignment")

    ("moses", po::bool_switch(&moses_mode), "moses alignment (not GIZA++ alignment)")

    ("prob-null",         po::value<double>(&prob_null),         "NULL probability")
    ("prob-union",        po::value<double>(&prob_union),        "union probability")
    ("prob-intersection", po::value<double>(&prob_intersection), "intersection probability")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
