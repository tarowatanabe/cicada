//
//  Copyright(C) 2011 Graham Neubig <neubig@ar.media.kyoto-u.ac.jp>
//                    Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// Extract hiero rule from pialign derivation
//
// Given ITG, we have only to traverse children and compute holes.
// We will always insert terminal in source-side, but not target-side..
// How to handle counts? Do we use original-hiero-style fractional counts?
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>
#include <set>
#include <stdexcept>
#include <memory>

#include <cicada/symbol.hpp>
#include <cicada/sentence.hpp>

#include <utils/compress_stream.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>


typedef boost::filesystem::path path_type;
typedef std::pair<int, int> span_type;
struct span_pair_type
{
  span_type source;
  span_type target;
  
  span_pair_type()
    : source(), target() {}
  span_pair_type(const span_type& __source, const span_type& __target)
    : source(__source), target(__target) {}
};

typedef cicada::Symbol   word_type;
typedef cicada::Sentence sentence_type;

struct itg_type
{
  typedef std::vector<itg_type, std::allocator<itg_type> > itg_pair_type;
  
  std::string source;
  std::string target;
  
  itg_pair_type antecedent;
  
  span_pair_type spans;
  
  bool inverse;
  bool block;

  itg_type() : source(), target(), antecedent(), inverse(false), block(false) {}

  void clear()
  {
    source.clear();
    target.clear();
    antecedent.clear();
    
    spans = span_pair_type();

    inverse = false;
    block = false;
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  itg_type,
			  (std::string, source)
			  (std::string, target)
			  (itg_type::itg_pair_type, antecedent)
			  (bool, inverse)
			  (bool, block)
			  )

template <typename Iterator>
struct derivation_parser : boost::spirit::qi::grammar<Iterator, itg_type(), boost::spirit::standard::space_type>
{
  derivation_parser() : derivation_parser::base_type(derivation)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    terminal %= qi::lexeme[+(standard::char_ - standard::space) - "(((" - ")))" - "|||"] | qi::attr("");
    
    derivation_terminal %= "(((" >> terminal >> "|||" >> terminal >> ")))";
    derivation_regular  %= '[' >> qi::attr("") >> qi::attr("") >> derivation_pair >> qi::attr(false) >> ']';
    derivation_inverse  %= '<' >> qi::attr("") >> qi::attr("") >> derivation_pair >> qi::attr(true) >> '>';
    derivation_pair     %= qi::repeat(2)[derivation];
    
    derivation_block %= ('{' >> (qi::hold[derivation_regular] | qi::hold[derivation_inverse] | derivation_terminal) >> '}')[phoenix::at_c<4>(qi::_val) = true];
    derivation_break  %= (qi::hold[derivation_regular] | qi::hold[derivation_inverse] | derivation_terminal)[phoenix::at_c<4>(qi::_val) = false];
    
    derivation %= qi::hold[derivation_block] | derivation_break;
    
    //qi::on_error<qi::fail>(derivation, std::cerr << qi::_4 << phoenix::val(": ") << phoenix::construct<std::string>(qi::_3, qi::_2) << std::endl);
  }
  
  typedef boost::spirit::standard::space_type space_type;
  
  boost::spirit::qi::rule<Iterator, std::string(), space_type>  terminal;

  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_block;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_break;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_terminal;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_regular;
  boost::spirit::qi::rule<Iterator, itg_type(), space_type>     derivation_inverse;
  
  boost::spirit::qi::rule<Iterator, itg_type::itg_pair_type(), space_type> derivation_pair;
};

std::ostream& print_tree(std::ostream& os, const itg_type& itg)
{
  if (itg.block)
    os << '{';

#if 0
  os << ' ' << itg.spans.source.first << ".." << itg.spans.source.second
     << ' ' << itg.spans.target.first << ".." << itg.spans.target.second
     << ' ';
#endif

  if (itg.antecedent.empty())
    os << "((( " << itg.source << " ||| " << itg.target << " )))";
  else {
    os << (itg.inverse ? '<' : '[');
    print_tree(os, itg.antecedent.front());
    print_tree(os, itg.antecedent.back());
    os << (itg.inverse ? '>' : ']');
  }
  if (itg.block)
    os << '}';
  
  return os;
}


void span_derivation_source(itg_type& itg, sentence_type& sentence)
{
  itg.spans.source.first = sentence.size();
  if (itg.antecedent.empty()) {
    if (! itg.source.empty())
      sentence.push_back(itg.source);
  } else {
    span_derivation_source(itg.antecedent.front(), sentence);
    span_derivation_source(itg.antecedent.back(),  sentence);
  }
  itg.spans.source.second = sentence.size();
}

void span_derivation_target(itg_type& itg, sentence_type& sentence)
{
  itg.spans.target.first = sentence.size();
  if (itg.antecedent.empty()) {
    if (! itg.target.empty())
      sentence.push_back(itg.target);
  } else {
    if (! itg.inverse) {
      span_derivation_target(itg.antecedent.front(), sentence);
      span_derivation_target(itg.antecedent.back(),  sentence);
    } else {
      span_derivation_target(itg.antecedent.back(),  sentence);
      span_derivation_target(itg.antecedent.front(), sentence);
    }
  }
  itg.spans.target.second = sentence.size();
}

struct HieroGrammar
{
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  
  typedef std::set<int, std::less<int>, std::allocator<int> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > alignment_type;

  typedef std::pair<int, int> point_type;
  typedef std::vector<point_type, std::allocator<point_type> > point_set_type;

  typedef sentence_type phrase_type;
  
  struct rule_pair_type
  {
    phrase_type source;
    phrase_type target;
    point_set_type alignment;
    
    rule_pair_type() : source(), target(), alignment() {}
    
    void clear()
    {
      source.clear();
      target.clear();
      alignment.clear();
    }

    friend
    std::ostream& operator<<(std::ostream& os, const rule_pair_type& x)
    {
      std::copy(x.source.begin(), x.source.end(), std::ostream_iterator<word_type>(os, " "));
      os << "||| ";
      std::copy(x.target.begin(), x.target.end(), std::ostream_iterator<word_type>(os, " "));
      os << "|||";
      point_set_type::const_iterator piter_end = x.alignment.end();
      for (point_set_type::const_iterator piter = x.alignment.begin(); piter != piter_end; ++ piter)
	os << ' ' << piter->first << '-' << piter->second;
      return os;
    }
  };
  
  typedef std::vector<rule_pair_type, std::allocator<rule_pair_type> > rule_pair_set_type;
  
  HieroGrammar(std::ostream& __os,
	       const int __max_span,
	       const int __max_length)
    : os(__os),
      max_span(__max_span),
      max_length(__max_length) {}
  
  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    span_pair_set_type spans;
    
    alignment.clear();
    alignment.resize(source.size());
    
    operator()(itg, source, target, spans, alignment, blocker);
  }

  template <typename Blocker>
  void operator()(const itg_type& itg,
		  const sentence_type& source,
		  const sentence_type& target,
		  span_pair_set_type& spans,
		  alignment_type& alignment,
		  Blocker blocker)
  {
    // post-traversal order to collect paired-span
  
    if (itg.spans.source.first == itg.spans.source.second
	|| itg.spans.target.first == itg.spans.target.second)
      return;
    
    const int length_source = itg.spans.source.second - itg.spans.source.first;
    const int length_target = itg.spans.target.second - itg.spans.target.first;
    
    if (blocker(itg) || itg.antecedent.empty()) {
      spans.push_back(itg.spans);
      
      for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src)
	for (int trg = itg.spans.target.first; trg != itg.spans.target.second; ++ trg)
	  alignment[src].insert(trg);

    } else if (max_span <= 0 || length_source <= max_span) {
      span_pair_set_type spans1;
      span_pair_set_type spans2;
      
      operator()(itg.antecedent.front(), source, target, spans1, alignment, blocker);
      operator()(itg.antecedent.back(),  source, target, spans2, alignment, blocker);
      
      // we will create hiero rule from spans1 and spans2!

      rule_pair_set_type rule_pairs;
      
      // first, single non-terminal
      if (! spans1.empty())
	hiero_rule(source, target, itg.spans, spans1, alignment, rule_pairs);
      
      if (! spans2.empty())
	hiero_rule(source, target, itg.spans, spans2, alignment, rule_pairs);
      
      // second, two non-terminals
      if (! spans1.empty() && ! spans2.empty())
	hiero_rule(source, target, itg.spans, spans1, spans2, alignment, rule_pairs);
    
      spans.swap(spans1);
      spans.insert(spans.end(), spans2.begin(), spans2.end());
      
      spans.push_back(itg.spans);
      
      if (! rule_pairs.empty()) {
	const double weight = 1.0 / rule_pairs.size();
	
	rule_pair_set_type::const_iterator riter_end = rule_pairs.end();
	for (rule_pair_set_type::const_iterator riter = rule_pairs.begin(); riter != riter_end; ++ riter)
	  os << *riter << " ||| " << weight << '\n';
      }
    }
    
    if (max_length <= 0 || (length_source <= max_length && length_target <= max_length)) {
      rule_pair_type rule_pair;
      rule_pair.source.push_back("[x]");
      rule_pair.target.push_back("[x]");
      rule_pair.source.insert(rule_pair.source.end(), source.begin() + itg.spans.source.first, source.begin() + itg.spans.source.second);
      rule_pair.target.insert(rule_pair.target.end(), target.begin() + itg.spans.target.first, target.begin() + itg.spans.target.second);
      
      for (int src = itg.spans.source.first; src != itg.spans.source.second; ++ src) {
	index_set_type::const_iterator aiter_begin = alignment[src].begin();
	index_set_type::const_iterator aiter_end   = alignment[src].end();
	
	for (index_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	  rule_pair.alignment.push_back(std::make_pair(src - itg.spans.source.first, *aiter - itg.spans.target.first));
      }
      
      os << rule_pair << " ||| " << 1.0 << '\n';
    }
  }
  
  void hiero_rule(const sentence_type& source,
		  const sentence_type& target,
		  const span_pair_type& span,
		  const span_pair_set_type& spans,
		  const alignment_type& alignment,
		  rule_pair_set_type& rule_pairs)
  {
    const int length_source = span.source.second - span.source.first;
    const int length_target = span.target.second - span.target.first;
    
    // we need at least one terminal...
    span_pair_set_type::const_iterator siter_end = spans.end();
    for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)
      if (siter->source != span.source) {
	
	const int length_source1 = siter->source.second - siter->source.first;
	const int length_target1 = siter->target.second - siter->target.first;
	
	if (max_length > 0 && length_source - length_source1 > max_length) continue;
	if (max_length > 0 && length_target - length_target1 > max_length) continue;

	rule_pairs.resize(rule_pairs.size() + 1);
	rule_pair_type& rule_pair = rule_pairs.back();
	
	rule_pair.source.clear();
      
	rule_pair.source.push_back("[x]");
	rule_pair.source.insert(rule_pair.source.end(), source.begin() + span.source.first, source.begin() + siter->source.first);
	rule_pair.source.push_back("[x,1]");
	rule_pair.source.insert(rule_pair.source.end(), source.begin() + siter->source.second, source.begin() + span.source.second);
      
	rule_pair.target.clear();
	rule_pair.target.push_back("[x]");
	rule_pair.target.insert(rule_pair.target.end(), target.begin() + span.target.first, target.begin() + siter->target.first);
	rule_pair.target.push_back("[x,1]");
	rule_pair.target.insert(rule_pair.target.end(), target.begin() + siter->target.second, target.begin() + span.target.second);
	
	for (int src = span.source.first; src != span.source.second; ++ src) 
	  if (is_out_of_span(siter->source, src)) {
	    index_set_type::const_iterator aiter_begin = alignment[src].begin();
	    index_set_type::const_iterator aiter_end   = alignment[src].end();
	    
	    for (index_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
	      if (is_out_of_span(siter->target, *aiter)) {
		const int mask_source = - (src    >= siter->source.second);
		const int mask_target = - (*aiter >= siter->target.second);
		
		// -1 for non-terminal label, [x,1]
		const int shift_source = (mask_source & (length_source1 - 1)) + span.source.first;
		const int shift_target = (mask_target & (length_target1 - 1)) + span.target.first;
		
		rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
	      }
	  }
      }
  }

  void hiero_rule(const sentence_type& source,
		  const sentence_type& target,
		  const span_pair_type& span,
		  const span_pair_set_type& spans1,
		  const span_pair_set_type& spans2,
		  const alignment_type& alignment,
		  rule_pair_set_type& rule_pairs)
  {
    const int length_source = span.source.second - span.source.first;
    const int length_target = span.target.second - span.target.first;
    
    span_pair_set_type::const_iterator siter1_begin = spans1.begin();
    span_pair_set_type::const_iterator siter1_end = spans1.end();
    span_pair_set_type::const_iterator siter2_begin = spans2.begin();
    span_pair_set_type::const_iterator siter2_end = spans2.end();
    for (span_pair_set_type::const_iterator siter1 = siter1_begin; siter1 != siter1_end; ++ siter1) {
      const span_pair_type& span1 = *siter1;

      const int length_source1 = span1.source.second - span1.source.first;
      const int length_target1 = span1.target.second - span1.target.first;
      
      for (span_pair_set_type::const_iterator siter2 = siter2_begin; siter2 != siter2_end; ++ siter2) 
	if (span1.source.second != siter2->source.first) {
	  const span_pair_type& span2 = *siter2;
	  
	  const int length_source2 = span2.source.second - span2.source.first;
	  const int length_target2 = span2.target.second - span2.target.first;
	  
	  if (max_length > 0 && length_source - length_source1 - length_source2 > max_length) continue;
	  if (max_length > 0 && length_target - length_target1 - length_target2 > max_length) continue;

	  rule_pairs.resize(rule_pairs.size() + 1);
	  rule_pair_type& rule_pair = rule_pairs.back();

	  const bool inverse = ! (span1.target.second <= span2.target.first);
	
	  const span_type& span_source1 = span1.source;
	  const span_type& span_source2 = span2.source;
	
	  const span_type& span_target1 = (inverse ? span2.target : span1.target);
	  const span_type& span_target2 = (inverse ? span1.target : span2.target);
	
	  rule_pair.source.clear();
	  rule_pair.source.push_back("[x]");
	  rule_pair.source.insert(rule_pair.source.end(), source.begin() + span.source.first, source.begin() + span_source1.first);
	  rule_pair.source.push_back("[x,1]");
	  rule_pair.source.insert(rule_pair.source.end(), source.begin() + span_source1.second, source.begin() + span_source2.first);
	  rule_pair.source.push_back("[x,2]");
	  rule_pair.source.insert(rule_pair.source.end(), source.begin() + span_source2.second, source.begin() + span.source.second);
	
	  rule_pair.target.clear();
	  rule_pair.target.push_back("[x]");
	  rule_pair.target.insert(rule_pair.target.end(), target.begin() + span.target.first, target.begin() + span_target1.first);
	  rule_pair.target.push_back(inverse ? "[x,2]" : "[x,1]");
	  rule_pair.target.insert(rule_pair.target.end(), target.begin() + span_target1.second, target.begin() + span_target2.first);
	  rule_pair.target.push_back(inverse ? "[x,1]" : "[x,2]");
	  rule_pair.target.insert(rule_pair.target.end(), target.begin() + span_target2.second, target.begin() + span.target.second);
	  
	  for (int src = span.source.first; src != span.source.second; ++ src) 
	    if (is_out_of_span(span1.source, src) && is_out_of_span(span2.source, src)) {
	      index_set_type::const_iterator aiter_begin = alignment[src].begin();
	      index_set_type::const_iterator aiter_end   = alignment[src].end();
	      
	      for (index_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter)
		if (is_out_of_span(span1.target, *aiter) && is_out_of_span(span2.target, *aiter)) {
		  const int mask_source1 = - (src >= span1.source.second);
		  const int mask_source2 = - (src >= span2.source.second);
		  
		  const int mask_target1 = - (*aiter >= span1.target.second);
		  const int mask_target2 = - (*aiter >= span2.target.second);
		  
		  // -1 for non-terminal label, [x,1] or [x,2]
		  const int shift_source = (span.source.first
					    + (mask_source1 & (length_source1 - 1))
					    + (mask_source2 & (length_source2 - 1)));
		  const int shift_target = (span.target.first
					    + (mask_target1 & (length_target1 - 1))
					    + (mask_target2 & (length_target2 - 1)));
		  
		  rule_pair.alignment.push_back(std::make_pair(src - shift_source, *aiter - shift_target));
		}
	    }
	}
    }
  }

  bool is_out_of_span(const span_type& span, const int pos) const
  {
    return pos < span.first || span.second <= pos;
  }

  std::ostream& os;
  const int max_span;
  const int max_length;
};


struct BlockerModel
{
  bool operator()(const itg_type& itg) const
  {
    return itg.block;
  }
};

struct BlockerBlock
{
  bool operator()(const itg_type& itg) const
  {
    return ((itg.spans.source.second - itg.spans.source.first) == 1) || ((itg.spans.target.second - itg.spans.target.first) == 1);
  }
};

struct BlockerTerminal
{
  bool operator()(const itg_type& itg) const
  {
    return itg.antecedent.empty();
  }
};

int main(int argc, char** argv)
{
  try {
    namespace po = boost::program_options;
    
    path_type input_file = "-";
    path_type output_file = "-";
    path_type output_source_file;
    path_type output_target_file;
    path_type output_alignment_file;
    int max_length = 7;
    int max_span = 15;
    
    int max_nodes   = 15;
    int max_height  = 4;
    int max_compose = 0;
    int max_scope   = 0;
    
    bool scfg_mode = false;
    bool ghkm_mode = false;

    bool frontier_source_mode = false;
    bool frontier_target_mode = false;

    bool phrase_mode = false;
    bool block_mode = false;
    bool exhaustive_mode = false;

    int debug = 0;
  
    po::options_description desc("options");
    desc.add_options()
      ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
      ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")
      
      ("source",    po::value<path_type>(&output_source_file),    "output source yield")
      ("target",    po::value<path_type>(&output_target_file),    "output target yield")
      ("alignment", po::value<path_type>(&output_alignment_file), "output word-for-word alignment")
      
      ("max-length", po::value<int>(&max_length)->default_value(max_length), "max terminal length")
      ("max-span",   po::value<int>(&max_span)->default_value(max_span),     "max span")
      
      ("max-nodes",   po::value<int>(&max_nodes)->default_value(max_nodes),     "max nodes")
      ("max-height",  po::value<int>(&max_height)->default_value(max_height),   "max height")
      ("max-compose", po::value<int>(&max_compose)->default_value(max_compose), "max compose")
      ("max-scope",   po::value<int>(&max_scope)->default_value(max_scope),     "max scope")
      ("frontier-source", po::bool_switch(&frontier_source_mode),               "take frontier of source side (string-to-* model)")
      ("frontier-target", po::bool_switch(&frontier_target_mode),               "take frontier of target side (*-to-string model)")
      
      ("scfg", po::bool_switch(&scfg_mode), "extract SCFG rules")
      ("ghkm", po::bool_switch(&ghkm_mode), "extract GHKM rules")
    
      ("phrase",     po::bool_switch(&phrase_mode),     "phrase-wise model alignment (many-to-many)")
      ("block",      po::bool_switch(&block_mode),      "block-wise alignment (one-to-many)")
      ("exhaustive", po::bool_switch(&exhaustive_mode), "exhaustive alignment (one-to-one)")
      
      ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
      ("help", "help message");
  
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
    po::notify(vm);
  
    if (vm.count("help")) {
      std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
      return 0;
    }
  
    if (int(phrase_mode) + block_mode + exhaustive_mode == 0)
      phrase_mode = true;
    
    if (int(phrase_mode) + block_mode + exhaustive_mode > 1)
      throw std::runtime_error("either phrase|block|exhaustive");

    if (scfg_mode && ghkm_mode)
      throw std::runtime_error("either scfg|ghkm");
    
    if (int(scfg_mode) + ghkm_mode == 0)
      scfg_mode = true;
    
    typedef boost::spirit::istream_iterator iter_type;
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    std::auto_ptr<std::ostream> os_src(! output_source_file.empty() ? new utils::compress_ostream(output_source_file, 1024 * 1024) : 0);
    std::auto_ptr<std::ostream> os_trg(! output_target_file.empty() ? new utils::compress_ostream(output_target_file, 1024 * 1024) : 0);
    std::auto_ptr<std::ostream> os_align(! output_alignment_file.empty() ? new utils::compress_ostream(output_alignment_file, 1024 * 1024) : 0);

    is.unsetf(std::ios::skipws);
    os.precision(20);
    
    derivation_parser<iter_type> parser;
    
    iter_type iter(is);
    iter_type end;
  
    itg_type itg;
    sentence_type source;
    sentence_type target;

    HieroGrammar::alignment_type alignment;

    HieroGrammar scfg_grammar(os, max_span, max_length);
  
    while (iter != end) {
      itg.clear();
    
      if (! boost::spirit::qi::phrase_parse(iter, end, parser, boost::spirit::standard::space, itg))
	throw std::runtime_error("parsing failed");

      source.clear();
      target.clear();

      span_derivation_source(itg, source);
      span_derivation_target(itg, target);

      if (os_src.get())
	*os_src << source << '\n';
      if (os_trg.get())
	*os_trg << target << '\n';
      
      //print_tree(std::cout, itg) << std::endl;
    
      if (ghkm_mode) {
	if (phrase_mode)
	  scfg_grammar(itg, source, target, alignment, BlockerModel());
	else if (block_mode)
	  scfg_grammar(itg, source, target, alignment, BlockerBlock());
	else
	  scfg_grammar(itg, source, target, alignment, BlockerTerminal());
      } else {
	
	
	
      }
      
      // dump alignment...
      if (os_align.get()) {
	bool initial = true;
	for (size_t src = 0; src != alignment.size(); ++ src) 
	  if (! alignment[src].empty()) {
	    const HieroGrammar::alignment_type::value_type& align = alignment[src];
	    
	    HieroGrammar::alignment_type::value_type::const_iterator aiter_end = align.end();
	    for (HieroGrammar::alignment_type::value_type::const_iterator aiter = align.begin(); aiter != aiter_end; ++ aiter) {
	      if (! initial)
		*os_align << ' ';
	      *os_align << src << '-' << *aiter;
	      
	      initial = false;
	    }
	  }
	*os_align << '\n';
      }
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}
