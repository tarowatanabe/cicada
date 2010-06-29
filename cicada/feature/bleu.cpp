
#include <map>

#include "cicada/feature/bleu.hpp"
#include "cicada/parameter.hpp"

#include "utils/hashmurmur.hpp"
#include "utils/compact_trie.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"

#include "cicada/eval/bleu.hpp"

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace feature
  {

    class BleuImpl
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
      typedef utils::compact_trie<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>, ngram_allocator_type> ngram_set_type;
      
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
      
      enum {
	NON_TERMINAL = 0,
	PHRASE_LEFT = 1,
	PHRASE_LEFT_RIGHT = 2,
      };
      

    public:
      BleuImpl(const int __order,
	       const bool __exact)
	: order(__order), exact(__exact) {}

      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge) const
      {
	if (ngrams.empty()) {
	  char* buf = reinterpret_cast<char*>(state);
	  std::fill(buf, buf + sizeof(int) * 2 + sizeof(id_type) * 4, 0);
	  return 0.0;
	}

	const rule_type& rule = *edge.rule;
	const phrase_type& target = rule.target;
	const phrase_type& source = rule.source;
	
	count_set_type counts;
	
	const int context_size = order - 1;

	int*     context_parsed     = reinterpret_cast<int*>(state);
        int*     context_hypothesis = context_parsed + 1;
        id_type* context_count      = reinterpret_cast<id_type*>(context_hypothesis + 1);
	id_type& context_prefix     = *(context_count + 1);
	id_type& context_suffix     = *(context_count + 2);
	id_type& context_flag       = *(context_count + 3);
	
	if (states.empty()) {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  
	  phrase_type::const_iterator titer_end = target.end();
	  for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	    if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  
	  collect_counts(buffer.begin(), buffer.end(), counts);
	  state_context(buffer.begin(), buffer.end(), context_prefix, context_suffix, context_flag);
	  
	  *context_parsed = 0;
	  phrase_type::const_iterator siter_end = source.end();
	  for (phrase_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
	    *context_parsed += (*siter != vocab_type::EPSILON);
	  
	  *context_hypothesis = buffer.size();
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();

	  const double bleu = bleu_score(counts, *context_hypothesis, *context_parsed, source_size);
	  
	  return bleu;

	} else {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  buffer.reserve(target.size() + (order * 2) * states.size());
	  
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  phrase_spans.clear();
	  target.terminals(std::back_inserter(phrase_spans));
	  
	  for (phrase_type::const_iterator iter = phrase_spans.front().first; iter != phrase_spans.front().second; ++ iter)
	    if (*iter != vocab_type::EPSILON)
	      buffer.push_back(*iter);
	  
	  *context_parsed = 0;
	  for (phrase_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter)
	    *context_parsed += (*siter != vocab_type::EPSILON && siter->is_terminal());
	  
	  *context_hypothesis = buffer.size();
	  
	  collect_counts(buffer.begin(), buffer.end(), counts);
	  state_context(buffer.begin(), buffer.end(), context_prefix, context_suffix, context_flag);
	  
	  double bleu_antecedent = 0.0;
	  
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = phrase_spans.begin() + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    const int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      throw std::runtime_error("this is a non-terminal, but no index!");
	    
	    const int*     antecedent_parsed     = reinterpret_cast<const int*>(states[antecedent_index]);
	    const int*     antecedent_hypothesis = antecedent_parsed + 1;
	    const id_type* antecedent_count      = reinterpret_cast<const id_type*>(antecedent_hypothesis + 1);
	    const id_type& antecedent_prefix     = *(antecedent_count + 1);
	    const id_type& antecedent_suffix     = *(antecedent_count + 2);
	    const id_type& antecedent_flag       = *(antecedent_count + 3);
	    
	    *context_parsed     += *antecedent_parsed;
	    *context_hypothesis += *antecedent_hypothesis + (span.second - span.first);
	    
	    const count_set_type& counts_antecedent = states_counts[*antecedent_count];
	    
	    collect_counts(span.first, span.second, counts);
	    counts.resize(utils::bithack::max(counts.size(), counts_antecedent.size()), count_type(0));
	    std::transform(counts_antecedent.begin(), counts_antecedent.end(), counts.begin(), counts.begin(), std::plus<count_type>());
	    
	    const double bleu_ant = bleu_score(counts_antecedent, *antecedent_hypothesis, *antecedent_parsed, source_size);
	    
	    //std::cerr << "bleu antecedent: " << bleu_ant << ' ';
	    //std::copy(antecedent_first, antecedent_end, std::ostream_iterator<symbol_type>(std::cerr, " "));
	    //std::cerr << std::endl;
	    
	    bleu_antecedent += bleu_ant;
	    
	    if (antecedent_flag != PHRASE_LEFT) {
	      const id_type suffix = (context_flag == PHRASE_LEFT ? context_prefix : context_suffix);
	      
	      collect_counts(suffix, antecedent_prefix, counts);
	      if (! ngrams.is_root(antecedent_suffix) && span.first != span.second)
		collect_counts(antecedent_suffix, longest_prefix(span.first, std::min(span.first + context_size, span.second)).first, counts);
	      
	      if (context_flag == PHRASE_LEFT) {
		buffer.clear();
		
		extract_contexts(antecedent_prefix, std::back_inserter(buffer));
		extract_contexts(suffix, std::back_inserter(buffer));
		std::reverse(buffer.begin(), buffer.end());
		
		if (! buffer.empty())
		  context_prefix = longest_prefix(buffer.begin(), std::min(buffer.begin() + context_size, buffer.end())).first;
		else
		  context_prefix = id_type(-1);
	      }
	      
	      buffer.clear();
	      extract_contexts(antecedent_suffix, std::back_inserter(buffer));
	      std::reverse(buffer.begin(), buffer.end());
	      buffer.insert(buffer.end(), span.first, span.second);
	      
	      if (! buffer.empty())
		context_suffix = longest_suffix(std::max(buffer.begin(), buffer.end() - context_size), buffer.end()).first;
	      else
		context_suffix = id_type(-1);
	      
	      context_flag = PHRASE_LEFT_RIGHT;
	    } else {
	      const id_type suffix = (context_flag == PHRASE_LEFT ? context_prefix : context_suffix);
	      
	      collect_counts(suffix, antecedent_prefix, counts);
	      
	      buffer.clear();
	      extract_contexts(antecedent_prefix, std::back_inserter(buffer));
	      extract_contexts(suffix, std::back_inserter(buffer));
	      std::reverse(buffer.begin(), buffer.end());
	      
	      if (span.first != span.second && ! buffer.empty())
		collect_counts(longest_suffix(std::max(buffer.end() - context_size, buffer.begin()), buffer.end()).first,
			       longest_prefix(span.first, std::min(span.first + context_size, span.second)).first,
			       counts);
	      
	      buffer.insert(buffer.end(), span.first, span.second);
	      
	      if (context_flag == PHRASE_LEFT)
		state_context(buffer.begin(), buffer.end(), context_prefix, context_suffix, context_flag);
	      else {
		if (! buffer.empty())
		  context_suffix = longest_suffix(std::max(buffer.begin(), buffer.end() - context_size), buffer.end()).first;
		else
		  context_suffix = id_type(-1);
	      }
	    }
	  }
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();
	  
	  const double bleu = bleu_score(counts, *context_hypothesis, *context_parsed, source_size);
	  	  
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
      
      void insert(const int __source_size, const sentence_type& sentence)
      {
	typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
	source_size = __source_size;
	
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
      void state_context(Iterator first, Iterator last, id_type& prefix, id_type& suffix, id_type& flag) const
      {
	const int context_size = order - 1;
	
	if (first == last) {
	  prefix = id_type(-1);
	  suffix = id_type(-1);
	  flag   = PHRASE_LEFT;
	} else if (last - first > context_size) {
	  prefix = longest_prefix(first, first + context_size).first;
	  suffix = longest_suffix(last - context_size, last).first;
	  flag   = PHRASE_LEFT_RIGHT;
	} else {
	  const std::pair<id_type, bool> result = longest_prefix(first, last);
	  if (result.second) {
	    prefix = result.first;
	    suffix = id_type(-1);
	    flag = PHRASE_LEFT;
	  } else {
	    prefix = result.first;
	    suffix = longest_suffix(first, last).first;
	    flag = PHRASE_LEFT_RIGHT;
	  }
	}
      }
      
      template <typename Iterator>
      void extract_contexts(id_type id, Iterator result) const
      {
	while (! ngrams.is_root(id)) {
	  *result = nodes[id].word;
	  ++ result;
	  id = nodes[id].parent;
	}
      }

      
      // compute longest suffix/prefix...
      template <typename Iterator>
      std::pair<id_type, bool> longest_suffix(Iterator first, Iterator last) const
      {
	// longest suffix of [first, last) where ngram matches..
	// we use this to supply additional conetxt for the future ngrams.
	bool complete = true;
	for (/**/; first != last; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter = first; iter != last; ++ iter) {
	    id = ngrams.find(id, *iter);
	    if (ngrams.is_root(id)) break;
	  }
	  if (! ngrams.is_root(id))
	    return std::make_pair(id, complete);
	  
	  complete = false;
	}
	return std::make_pair(ngrams.root(), complete);
      }
      
      template <typename Iterator>
      std::pair<id_type, bool> longest_prefix(Iterator first, Iterator last) const
      {
	// longest prefix of [first, last) where ngram matches...
	// we use this to compute the ngram ranges where we need to rescore...
	ngram_set_type::id_type id = ngrams.root();
	for (/**/; first != last; ++ first) {
	  const ngram_set_type::id_type id_next = ngrams.find(id, *first);
	  if (ngrams.is_root(id_next))
	    return std::make_pair(id, false);
	  id = id_next;
	}
	return std::make_pair(id, true);
      }
      
      template <typename Iterator>
      void collect_counts(Iterator first, Iterator last, count_set_type& counts) const
      {
	typedef typename boost::is_integral<Iterator>::type __integral;
	
	__collect_counts_dispatch(first, last, counts, __integral());
      }

      void __collect_counts_dispatch(id_type id_prefix, id_type id_suffix, count_set_type& counts, boost::true_type) const
      {
	typedef std::vector<word_type, std::allocator<word_type> > context_type;
      
	// we do not consider additonal contexts for id
	if (ngrams.is_root(id_prefix) || ngrams.is_root(id_suffix)) return;
      
	if (exact)
	  counts.resize(nodes.size(), count_type(0));
	else
	  counts.resize(order, count_type(0));
      
	// traverse back id_prefix and id_suffix do we construct a phrase?
      
	context_type context;
      
	// suffix to rescore with additional prefix...
	while (! ngrams.is_root(id_suffix)) {
	  context.push_back(nodes[id_suffix].word);
	  id_suffix = nodes[id_suffix].parent;
	}
	const int size_suffix = context.size();
      
	// now try collect count using additional prefix
	while (! ngrams.is_root(id_prefix)) {
	  context.push_back(nodes[id_prefix].word);
	  id_prefix = nodes[id_prefix].parent;
	
	  // we start from context.back(), and enumerate ngrams from context.size() - size_suffix to order...
	  const int order_prefix = context.size() - size_suffix;
	  int n = 1;
	  ngram_set_type::id_type id = ngrams.root();
	  context_type::const_reverse_iterator citer_begin = context.rbegin();
	  context_type::const_reverse_iterator citer_end = context.rend();
	  context_type::const_reverse_iterator citer_last = std::min(citer_begin + order, citer_end);
	  for (context_type::const_reverse_iterator citer = citer_begin; citer != citer_last; ++ citer, ++ n) {
	    id = ngrams.find(id, *citer);
	    if (ngrams.is_root(id)) break;
	  
	    if (n > order_prefix) {

	      if (exact)
		++ counts[id];
	      else
		++ counts[nodes[id].order - 1];
	    }
	  }
	}
      }
    
      template <typename Iterator>
      void __collect_counts_dispatch(Iterator first, Iterator last, count_set_type& counts, boost::false_type) const
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
      
      double bleu_score(const count_set_type& __counts, const int hypothesis_size, const int parsed_size, const int source_size) const
      {
	count_set_type counts_bleu(order, count_type(0));
	if (exact)
	  for (ngram_set_type::id_type id = 0; id < __counts.size();++ id)
	    counts_bleu[nodes[id].order - 1] += utils::bithack::min(int(__counts[id]), int(ngrams[id]));
	
	const count_set_type& counts = (exact ? counts_bleu : __counts);

	const cicada::eval::Bleu* __bleu = (score ? dynamic_cast<const cicada::eval::Bleu*>(score.get()) : 0);
	
	if (__bleu && __bleu->length_reference > 0) {
	  
	  const double hypothesis_length = tst_size(hypothesis_size, parsed_size, source_size);
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
	  
	  const double hypothesis_length = tst_size(hypothesis_size, parsed_size, source_size);
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

      double tst_size(int length, int parsed, int source_size) const
      {
	if (length == 0 || parsed == 0) return 0.0;
	
	// we will scale hypothesis length by the # of parsed words
	
	return (parsed < source_size
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
      
    private:
      
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
    };
    
    inline
    bool true_false(const std::string& token)
    {
      if (strcasecmp(token.c_str(), "true") == 0)
	return true;
      if (strcasecmp(token.c_str(), "yes") == 0)
	return true;
      if (atoi(token.c_str()) > 0)
	return true;
      return false;
    }
    
    Bleu::Bleu(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "bleu")
	throw std::runtime_error("this is not Bleu feature: " + parameter);
      
      const int order = (param.find("order") != param.end() ? boost::lexical_cast<int>(param.find("order")->second) : 4);
      const bool exact = (param.find("exact") != param.end() ? true_false(param.find("exact")->second) : false);
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, exact));
      
      // two-side context + counts-id (prefix/suffix/hypothesis) + flag
      base_type::__state_size = sizeof(int) * 2 + sizeof(impl_type::id_type) * 4;
      base_type::__feature_name = "bleu";
      
      pimpl = bleu_impl.release();
    }
    
    Bleu::~Bleu() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    void Bleu::operator()(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates) const
    {
      const double score = pimpl->bleu_score(state, states, edge);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
    }
    
    void Bleu::operator()(const state_ptr_type& state,
			  feature_set_type& features) const
    {
      // we do nothing...
    }

    void Bleu::initialize()
    {
      pimpl->initialize();
    }
    
    void Bleu::clear()
    {
      pimpl->clear();
    }
    
    void Bleu::insert(const int source_size, const sentence_type& sentence)
    {
      pimpl->insert(source_size, sentence);
    }

    void Bleu::insert(const score_ptr_type& score)
    {
      pimpl->insert(score);
    }

  };
};
