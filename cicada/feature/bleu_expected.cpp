
#include <map>


#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"

#include "cicada/feature/bleu_expected.hpp"
#include "cicada/parameter.hpp"
#include "cicada/eval/bleu.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/semiring.hpp"


#include <boost/numeric/conversion/bounds.hpp>

#include <unicode/uchar.h>
#include <unicode/unistr.h>
#include <unicode/schriter.h>
#include <unicode/bytestream.h>

namespace cicada
{
  namespace feature
  {

    class BleuExpectedImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;

      typedef cicada::eval::Score score_type;
      typedef score_type::score_ptr_type score_ptr_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;

      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
            
      // this implementation specific...
      typedef uint32_t id_type;
      typedef uint32_t count_type;

      typedef symbol_type word_type;
      
      typedef std::allocator<std::pair<const word_type, count_type> >  ngram_allocator_type;
      typedef utils::compact_trie_dense<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>, ngram_allocator_type> ngram_set_type;
      
      struct Node
      {
	Node() : word(), parent(id_type(-1)), order(0) {}
	
	word_type word;
	id_type   parent;
	int       order;
      };
      typedef Node node_type;
      typedef std::vector<node_type, std::allocator<node_type> > node_set_type;
      typedef std::vector<int, std::allocator<int> > size_set_type;
    
      typedef utils::simple_vector<count_type, std::allocator<count_type> > count_set_type;
      
      struct count_set_hash : public utils::hashmurmur<size_t>
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	size_t operator()(const count_set_type& x) const
	{
	  return hasher_type::operator()(x.begin(), x.end(), 0);
	}
      };
      
      typedef utils::indexed_set<count_set_type, count_set_hash, std::equal_to<count_set_type>, std::allocator<count_set_type> > states_count_set_type;


    public:
      BleuExpectedImpl(const int __order,
		       const bool __exact,
		       const bool __split,
		       const bool __yield_source)
	: ngrams(word_type()), nodes(), sizes(),
	  order(__order), exact(__exact), split(__split), yield_source(__yield_source) {}
      
      template <typename Sentence>
      void split_non_ascii_characters(const Sentence& phrase, Sentence& result) const
      {
	std::vector<word_type, std::allocator<word_type> > tokens;
	std::string buffer;
	
	typename Sentence::const_iterator piter_end = phrase.end();
	for (typename Sentence::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) {
	  
	  UnicodeString uword = UnicodeString::fromUTF8(static_cast<const std::string&>(*piter));
	  
	  StringCharacterIterator iter(uword);
	  for (iter.setToStart(); iter.hasNext(); /**/) {
	    const UChar32 c = iter.next32PostInc();
	    
	    if (c < 128)
	      buffer.push_back(c);
	    else {
	      // we will split...
	      if (! buffer.empty())
		tokens.push_back(word_type(buffer.begin(), buffer.end()));
	      buffer.clear();
	      
	      StringByteSink<std::string> __sink(&buffer);
	      UnicodeString(c).toUTF8(__sink);
	      
	      tokens.push_back(word_type(buffer.begin(), buffer.end()));
	      buffer.clear();
	    }
	  }
	  
	  if (! buffer.empty())
	    tokens.push_back(word_type(buffer.begin(), buffer.end()));
	  buffer.clear();
	}
	
	result.assign(tokens.begin(), tokens.end());
      }


      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			const bool final=false) const
      {
	if (ngrams.empty()) {
	  char* buf = reinterpret_cast<char*>(state);
	  std::fill(buf, buf + sizeof(symbol_type) * order * 2 + sizeof(int) * 2 + sizeof(id_type) * 2, 0);
	  return 0.0;
	}

	const rule_type& rule = *edge.rule;
	
	phrase_type target_split;
	if (split) {
	  if (yield_source)
	    split_non_ascii_characters(rule.source, target_split);
	  else
	    split_non_ascii_characters(rule.target, target_split);
	}
	
	const phrase_type& target = (split ? target_split : (yield_source ? rule.source : rule.target));
	const phrase_type& source = rule.source;
	
	count_set_type counts;

	symbol_type* context_first = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_last  = context_first + order * 2;

	std::fill(context_first, context_last, vocab_type::EMPTY);
	
	int*     context_parsed     = reinterpret_cast<int*>(context_last);
	int*     context_hypothesis = context_parsed + 1;
	id_type* context_count       = reinterpret_cast<id_type*>(context_hypothesis + 1);

	const int context_size = order - 1;
	
	if (states.empty()) {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  
	  phrase_type::const_iterator titer_end = target.end();
	  for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	    if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  
	  collect_counts(buffer.begin(), buffer.end(), counts);
	  
	  if (buffer.size() <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context_first);
	  else {
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram_prefix(biter_begin, biter_begin + context_size);
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix = ngram_suffix(biter_end - context_size, biter_end);
	    
	    std::copy(prefix.first, prefix.second, context_first);
	    context_first[prefix.second - prefix.first] = vocab_type::STAR;
	    std::copy(suffix.first, suffix.second, context_first + (prefix.second - prefix.first) + 1);
	  }
	  
	  *context_parsed = 0;
	  phrase_type::const_iterator siter_end = source.end();
	  for (phrase_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
	    *context_parsed += (*siter != vocab_type::EPSILON);
	  
	  *context_hypothesis = buffer.size();
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();

	  const double bleu = bleu_score(counts, *context_hypothesis, *context_parsed, source_size, ! final);
	  
	  //std::cerr << "bleu: " << bleu << ' ';
	  //std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
	  //std::cerr << std::endl;
	  
	  return bleu;

	} else {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  buffer.reserve(target.size() + (order * 2) * states.size());
	  
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  phrase_spans.clear();
	  target.terminals(std::back_inserter(phrase_spans));
	  
	  int star_first = -1;
	  int star_last  = -1;
	  
	  *context_parsed = 0;
	  for (phrase_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter)
	    *context_parsed += (*siter != vocab_type::EPSILON && siter->is_terminal());
	  
	  for (phrase_type::const_iterator iter = phrase_spans.front().first; iter != phrase_spans.front().second; ++ iter)
	    if (*iter != vocab_type::EPSILON)
	      buffer.push_back(*iter);
	  
	  *context_hypothesis = buffer.size();
	  
	  collect_counts(buffer.begin(), buffer.end(), counts);
	  
	  double bleu_antecedent = 0.0;

	  buffer_type::const_iterator biter_first = buffer.begin();
	  
	  phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = siter - (siter_begin + 1);
	    
	    const symbol_type* antecedent_first = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    const symbol_type* antecedent_last  = antecedent_first + order * 2;
	    
	    const symbol_type* antecedent_end  = std::find(antecedent_first, antecedent_last, vocab_type::EMPTY);
	    const symbol_type* antecedent_star = std::find(antecedent_first, antecedent_end, vocab_type::STAR);
	    
	    const int*     antecedent_parsed     = reinterpret_cast<const int*>(antecedent_last);
	    const int*     antecedent_hypothesis = antecedent_parsed + 1;
	    const id_type* antecedent_count       = reinterpret_cast<const id_type*>(antecedent_hypothesis + 1);

	    const count_set_type& counts_antecedent = states_counts[*antecedent_count];
	    
	    const double bleu_ant = bleu_score(counts_antecedent, *antecedent_hypothesis, *antecedent_parsed, source_size);
	    
	    //std::cerr << "bleu antecedent: " << bleu_ant << ' ';
	    //std::copy(antecedent_first, antecedent_end, std::ostream_iterator<symbol_type>(std::cerr, " "));
	    //std::cerr << std::endl;
	    
	    bleu_antecedent += bleu_ant;
	    
	    // merge statistics...
	    counts.resize(utils::bithack::max(counts.size(), counts_antecedent.size()), count_type(0));
	    std::transform(counts_antecedent.begin(), counts_antecedent.end(), counts.begin(), counts.begin(), std::plus<count_type>());
	    
	    *context_parsed     += *antecedent_parsed;
	    *context_hypothesis += *antecedent_hypothesis;
	    
	    buffer_type::const_iterator biter = buffer.end();
	    
	    buffer.insert(buffer.end(), antecedent_first, antecedent_star);
	    
	    buffer_type::const_iterator biter_end = buffer.end();
	    
	    collect_counts(biter_first, biter, biter_end, counts);
	    
	    // insert context after star
	    if (antecedent_star != antecedent_end) {
	      biter_first = buffer.end() + 1;
	      
	      star_last = buffer.size() + 1;
	      if (star_first < 0)
		star_first = buffer.size() + 1;
	      
	      buffer.insert(buffer.end(), antecedent_star, antecedent_end);
	    }
	    
	    // collect counts from edge's terminals
	    {
	      buffer_type::const_iterator biter = buffer.end();
	      
	      for (phrase_type::const_iterator iter = span.first; iter != span.second; ++ iter)
		if (*iter != vocab_type::EPSILON)
		  buffer.push_back(*iter);
	      
	      buffer_type::const_iterator biter_end = buffer.end();
	      
	      collect_counts(biter_first, biter, biter_end, counts);
	      
	      *context_hypothesis += biter_end - biter;
	    }
	  }
	  
	  if (star_first >= 0) {
	    const int prefix_size = utils::bithack::min(star_first, context_size);
	    const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	    
	    buffer_type::const_iterator biter_begin = buffer.begin();
	    buffer_type::const_iterator biter_end   = buffer.end();
	    
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram_prefix(biter_begin, biter_begin + prefix_size);
	    std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix = ngram_suffix(biter_end - suffix_size, biter_end);
	    
	    std::copy(prefix.first, prefix.second, context_first);
	    context_first[prefix.second - prefix.first] = vocab_type::STAR;
	    std::copy(suffix.first, suffix.second, context_first + (prefix.second - prefix.first) + 1);
	  } else {
	    if (buffer.size() <= context_size)
	      std::copy(buffer.begin(), buffer.end(), context_first);
	    else {
	      buffer_type::const_iterator biter_begin = buffer.begin();
	      buffer_type::const_iterator biter_end   = buffer.end();
	      
	      std::pair<buffer_type::const_iterator, buffer_type::const_iterator> prefix = ngram_prefix(biter_begin, biter_begin + context_size);
	      std::pair<buffer_type::const_iterator, buffer_type::const_iterator> suffix = ngram_suffix(biter_end - context_size, biter_end);
	      
	      std::copy(prefix.first, prefix.second, context_first);
	      context_first[prefix.second - prefix.first] = vocab_type::STAR;
	      std::copy(suffix.first, suffix.second, context_first + (prefix.second - prefix.first) + 1);
	    }
	  }
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();

	  const double bleu = bleu_score(counts, *context_hypothesis, *context_parsed, source_size, ! final);
	  
	  //std::cerr << "bleu: " << bleu;
	  //std::cerr << " antecedent: " << bleu_antecedent << ' ';
	  //std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
	  //std::cerr << std::endl;
	  
	  return  bleu - bleu_antecedent;
	}
      }

      void initialize()
      {
	states_counts.clear();
      }
      
      void clear()
      {
	// for ngram handling
	ngrams.clear();
	nodes.clear();
	sizes.clear();
	
	// for count set representation
	states_counts.clear();
	
	source_size = 0;

	score.reset();
      }

      void insert(const score_ptr_type& __score)
      {
	score = __score;
	
	if (score && ! dynamic_cast<const cicada::eval::Bleu*>(score.get()))
	  throw std::runtime_error("this is not a bleu-score!");
      }
      
      void insert(const sentence_type& __sentence)
      {
	typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
	sentence_type sentence_split;
	if (split)
	  split_non_ascii_characters(__sentence, sentence_split);
	const sentence_type& sentence = (split ? sentence_split : __sentence);
	
	counts_type counts;
	sentence_type::const_iterator siter_end = sentence.end();
	for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  ngram_set_type::id_type id = ngrams.root();
	  
	  int n = 1;
	  for (sentence_type::const_iterator iter = siter; iter != std::min(siter + order, siter_end); ++ iter, ++ n) {
	    const ngram_set_type::id_type id_next = ngrams.insert(id, *iter);
	    
	    ++ counts[id_next];
	    
	    if (id_next >= nodes.size())
	      nodes.resize(id_next + 1, node_type());
	    
	    nodes[id_next].word = *iter;
	    nodes[id_next].parent = id;
	    nodes[id_next].order  = n;
	    
	    id = id_next;
	  }
	}
	
	// collect clipped ngram counts
	counts_type::const_iterator citer_end = counts.end();
	for (counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
	  ngrams[citer->first] = utils::bithack::max(ngrams[citer->first], citer->second);
	
	// keep sizes...
	sizes.push_back(sentence.size());
	std::sort(sizes.begin(), sizes.end());
      }

    private:
      template <typename Iterator>
      std::pair<Iterator, Iterator> ngram_prefix(Iterator first, Iterator last) const
      {
	if (first == last || first + 1 == last) return std::make_pair(first, last);
	
	Iterator iter = first;
	ngram_set_type::id_type id = ngrams.root();
	for (/**/; iter != last; ++ iter) {
	  const ngram_set_type::id_type id_next = ngrams.find(id, *iter);
	  if (ngrams.is_root(id)) break;
	}
	return std::make_pair(first, std::min(std::max(first + 1, iter), last));
      }
      
      template <typename Iterator>
      std::pair<Iterator, Iterator> ngram_suffix(Iterator first, Iterator last) const
      {
	if (first == last || first + 1 == last) return std::make_pair(first, last);
	
	first = std::max(first, last - order);
	
	for (/**/; first != last - 1; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter = first; iter != last; ++ iter) {
	    id = ngrams.find(id, *iter);
	    if (ngrams.is_root(id)) break;
	  }
	  if (! ngrams.is_root(id))
	    return std::make_pair(first, last);
	}
	
	return std::make_pair(first, last);
      }
      
      template <typename Iterator>
      void collect_counts(Iterator first, Iterator iter, Iterator last, count_set_type& counts) const
      {
	const int context_size = order - 1;
	
	first = std::max(first, iter - context_size);
	
	if (exact)
	  counts.resize(nodes.size(), count_type(0));
	else
	  counts.resize(order, count_type(0));
	
	// we will collect counts at [iter, last) with context from [first, iter)
	for (/**/; first != iter; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	    id = ngrams.find(id, *iter2);
	    
	    if (ngrams.is_root(id)) break;
	    if (iter2 < iter) continue;
	    
	    if (exact)
	      ++ counts[id];
	    else
	      ++ counts[nodes[id].order - 1];
	  }
	}
      }
      
      template <typename Iterator>
      void collect_counts(Iterator first, Iterator last, count_set_type& counts) const
      {
	if (exact)
	  counts.resize(nodes.size(), count_type(0));
	else
	  counts.resize(order, count_type(0));
	
	// we will collect counts at [first, last)
	for (/**/; first != last; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter = first; iter != std::min(first + order, last); ++ iter) {
	    id = ngrams.find(id, *iter);
	    
	    if (ngrams.is_root(id)) break;
	    
	    if (exact)
	      ++ counts[id];
	    else
	      ++ counts[nodes[id].order - 1];
	  }
	}
      }
      

      double brevity_penalty(const double hypothesis_size, const double reference_size) const
      {
	if (hypothesis_size == 0 || reference_size == 0)
	  return 0.0;
	else
	  return std::min(1.0 - reference_size / hypothesis_size, 0.0);
      }
      
      double bleu_score(const count_set_type& __counts, const int hypothesis_size, const int parsed_size, const int source_size, const bool scaling=true) const
      {
	count_set_type counts_bleu(order, count_type(0));
	if (exact)
	  for (ngram_set_type::id_type id = 0; id < __counts.size();++ id)
	    counts_bleu[nodes[id].order - 1] += utils::bithack::min(int(__counts[id]), int(ngrams[id]));
	
	const count_set_type& counts = (exact ? counts_bleu : __counts);

	const cicada::eval::Bleu* __bleu = (score ? dynamic_cast<const cicada::eval::Bleu*>(score.get()) : 0);
	
	if (__bleu && __bleu->length_reference > 0) {
	  
	  const double hypothesis_length = tst_size(hypothesis_size, parsed_size, source_size, scaling);
	  const double reference_length  = ref_size(hypothesis_length);
	  
	  double smooth = 0.5;
	  double bleu = brevity_penalty(hypothesis_length + __bleu->length_hypothesis, reference_length + __bleu->length_reference);
	  
	  const int ngram_size = utils::bithack::min(int(counts.size()), hypothesis_size);
	  const int bleu_order = utils::bithack::max(counts.size(), __bleu->ngrams_hypothesis.size());
	  
	  const double factor = 1.0 / order;
	  for (int n = 1; n <= bleu_order; ++ n) {
	    const double count = (double(n <= ngram_size ? double(counts[n - 1]) : 0.0)
				  + (n <= __bleu->ngrams_hypothesis.size() ? double(__bleu->ngrams_hypothesis[n - 1]) : 0.0));
	    const double norm  = (double(n <= ngram_size ? double(hypothesis_size + 1 - n) : 0.0)
				  + (n <= __bleu->ngrams_reference.size() ? double(__bleu->ngrams_reference[n - 1]) : 0.0));
	    const double p = (norm > 0.0 ? ((count > 0.0 ? count : smooth) / norm) : 0.0);
	    
	    bleu += (p > 0.0 ? std::log(p) : 0.0) * factor;
	    smooth *= 0.5;
	  }
	  
	  return std::exp(bleu);
	} else {
	  if (hypothesis_size == 0 || counts.empty()) return 0.0;
	  
	  const double hypothesis_length = tst_size(hypothesis_size, parsed_size, source_size, scaling);
	  const double reference_length  = ref_size(hypothesis_length);
	  
	  double smooth = 0.5;
	  double bleu = brevity_penalty(hypothesis_length, reference_length);
	  
	  const int ngram_size = utils::bithack::min(int(counts.size()), hypothesis_size);
	  
	  const double factor = 1.0 / order;
	  for (int n = 1; n <= ngram_size; ++ n) {
	    const int count = counts[n - 1];
	    
	    bleu += std::log((count ? double(count) : smooth) / (hypothesis_size + 1 - n)) * factor;
	    smooth *= 0.5;
	  }
	  
	  return std::exp(bleu);
	}
      }

      double tst_size(int length, int parsed, int source_size, const bool scaling=true) const
      {
	if (length == 0 || parsed == 0) return 0.0;
	
	// we will scale hypothesis length by the # of parsed words
	
	return (parsed < source_size && scaling
		? (double(source_size) / parsed) * length
		: double(length));
      }
      
      double ref_size(const double hypothesis_size) const
      {
	if (hypothesis_size == 0) return 0;
	
	int reference_size = 0;
	double min_diff = boost::numeric::bounds<double>::highest();
	for (size_set_type::const_iterator siter = sizes.begin(); siter != sizes.end(); ++ siter) {
	  const double diff = std::fabs(hypothesis_size - *siter);

	  if (diff < min_diff) {
	    min_diff = diff;
	    reference_size = *siter;
	  } else if (diff == min_diff)
	    reference_size = utils::bithack::min(reference_size, *siter);
	}
	
	return reference_size;
      }
      
    public:
      
      buffer_type          buffer_impl;
      phrase_span_set_type phrase_spans_impl;

      ngram_set_type ngrams;
      node_set_type  nodes;
      size_set_type  sizes;
      
      states_count_set_type states_counts;
      
      int source_size;

      score_ptr_type score;

      int order;
      bool exact;
      bool split;
      bool yield_source;
    };
    
    
    BleuExpected::BleuExpected(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "bleu")
	throw std::runtime_error("this is not BleuExpected feature: " + parameter);

      int order = 4;
      bool exact = false;
      bool split = false;
      
      bool yield_source = false;
      bool yield_target = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "order") == 0)
	  order = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "exact") == 0)
	  exact = utils::lexical_cast<bool>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "split") == 0)
	  split = utils::lexical_cast<bool>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    yield_source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    yield_target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	} else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (yield_source && yield_target)
	throw std::runtime_error("you cannot specify both source/target yield");
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, exact, split, yield_source));
      
      // two-side context + length (hypothesis/reference) + counts-id (hypothesis/reference)
      base_type::__state_size = sizeof(symbol_type) * order * 2 + sizeof(int) * 2 + sizeof(impl_type::id_type) * 2;
      base_type::__feature_name = "bleu";
      
      pimpl = bleu_impl.release();
    }
    
    BleuExpected::~BleuExpected() { std::auto_ptr<impl_type> tmp(pimpl); }

    BleuExpected::BleuExpected(const BleuExpected& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    BleuExpected& BleuExpected::operator=(const BleuExpected& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void BleuExpected::apply(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {
      double score = pimpl->bleu_score(state, states, edge, final);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
    }

    void BleuExpected::apply_coarse(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      
    }
    
    void BleuExpected::initialize()
    {
      pimpl->initialize();
    }
    
    void BleuExpected::clear()
    {
      pimpl->clear();
    }
    
    struct source_length_function
    {
      typedef cicada::Vocab vocab_type;
      typedef cicada::Rule rule_type;
      typedef cicada::semiring::Tropical<int> value_type;
      
      template <typename Edge>
      value_type operator()(const Edge& edge) const
      {
	int length = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->source.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->source.begin(); siter != siter_end; ++ siter)
	  length += (*siter != vocab_type::EPSILON && siter->is_terminal());
    
	// since we will "max" at operator+, we will collect negative length
	return cicada::semiring::traits<value_type>::log(- length);
      }
    };


    void BleuExpected::assign(const hypergraph_type& hypergraph,
			      const lattice_type& lattice,
			      const span_set_type& spans,
			      const sentence_set_type& targets,
			      const ngram_count_set_type& ngram_counts)
    {
      int source_length = lattice.shortest_distance();
      if (hypergraph.is_valid()) {
	// we will enumerate forest structure... and collect min-size...
	std::vector<source_length_function::value_type, std::allocator<source_length_function::value_type> > lengths(hypergraph.nodes.size());
	
	cicada::inside(hypergraph, lengths, source_length_function());
	
	source_length = - log(lengths.back());
      }
      
      pimpl->clear();
      sentence_set_type::const_iterator titer_end = targets.end();
      for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
	pimpl->insert(*titer);

      pimpl->source_size = source_length;
    }
    

    void BleuExpected::assign(const score_ptr_type& score)
    {
      pimpl->insert(score);
    }

  };
};
