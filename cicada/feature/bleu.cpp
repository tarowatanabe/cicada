
#include <map>

#include "cicada/feature/bleu.hpp"
#include "cicada/parameter.hpp"

#include "utils/hashmurmur.hpp"
#include "utils/compact_trie.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"

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


    public:
      BleuImpl(const int __order,
	       const bool __exact)
	: order(__order), exact(__exact) {}

      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge) const
      {
	const rule_type& rule = *edge.rule;
	const phrase_type& target = rule.target;
	const phrase_type& source = rule.source;
	
	count_set_type counts;

	symbol_type* context_first = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_last  = context_first + order * 2;
	
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
	  
	  std::fill(context_first, context_last, vocab_type::EMPTY);
	  
	  if (buffer.size() <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context_first);
	  else {
	    std::copy(buffer.begin(), buffer.begin() + context_size, context_first);
	    context_first[context_size] = vocab_type::STAR;
	    std::copy(buffer.end() - context_size, buffer.end(), context_first + order);
	  }
	  
	  *context_parsed = 0;
	  phrase_type::const_iterator siter_end = source.end();
	  for (phrase_type::const_iterator siter = source.begin(); siter != siter_end; ++ siter)
	    *context_parsed += (*siter != vocab_type::EPSILON);
	  
	  *context_hypothesis = buffer.size();
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();
	  
	  return bleu_score(counts, *context_hypothesis, *context_parsed, source_size);

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
	  
	  double score = 1.0;

	  buffer_type::const_iterator biter_first = buffer.begin();
	  
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = phrase_spans.begin() + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    const int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      throw std::runtime_error("this is a non-terminal, but no index!");
	    
	    const symbol_type* antecedent_first = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    const symbol_type* antecedent_last  = antecedent_first + order * 2;
	    
	    const symbol_type* antecedent_end  = std::find(antecedent_first, antecedent_last, vocab_type::EMPTY);
	    const symbol_type* antecedent_star = std::find(antecedent_first, antecedent_end, vocab_type::STAR);
	    
	    const int*     antecedent_parsed     = reinterpret_cast<const int*>(antecedent_last);
	    const int*     antecedent_hypothesis = antecedent_parsed + 1;
	    const id_type* antecedent_count       = reinterpret_cast<const id_type*>(antecedent_hypothesis + 1);
	    
	    *context_parsed     += *antecedent_parsed;
	    *context_hypothesis += *antecedent_hypothesis;

	    score /= bleu_score(states_counts[*antecedent_count], *antecedent_hypothesis, *antecedent_parsed, source_size);
	    
	    buffer_type::const_iterator biter = buffer.end();
	    
	    buffer.insert(buffer.end(), antecedent_first, antecedent_star);
	    
	    buffer_type::const_iterator biter_end = buffer.end();
	    
	    collect_counts(std::max(biter_first, biter - context_size), biter, biter_end, counts);
	    
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
	      
	      collect_counts(std::max(biter_first, biter - context_size), biter, biter_end, counts);
	      
	      *context_hypothesis += biter_end - biter;
	    }
	  }
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();
	  
	  return bleu_score(counts, *context_hypothesis, *context_parsed, source_size);
	}
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
      void collect_counts(Iterator first, Iterator iter, Iterator last, count_set_type& counts) const
      {
	if (exact)
	  counts.resize(nodes.size(), count_type(0));
	else
	  counts.resize(order, count_type(0));
	
	// we will collect counts at [first, last)
	for (/**/; first != last; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	    id = ngrams.find(id, *iter2);
	    
	    if (ngrams.is_root(id)) break;
	    
	    // additional constraint: iter2 must be greater or equal to iter for count collection
	    if (iter2 >= iter) {
	      if (exact)
		++ counts[id];
	      else
		++ counts[nodes[id].order - 1];
	    }
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
      
      double bleu_score(const count_set_type& __counts, const int hypothesis_size, const int parsed_size, const int source_size) const
      {
	count_set_type counts_bleu(order, count_type(0));
	if (exact)
	  for (ngram_set_type::id_type id = 0; id < __counts.size();++ id)
	    counts_bleu[nodes[id].order - 1] += utils::bithack::min(int(__counts[id]), int(ngrams[id]));
	
	const count_set_type& counts = (exact ? counts_bleu : __counts);
	
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
      
      // two-side context + length (hypothesis/reference) + counts-id (hypothesis/reference)
      base_type::__state_size = sizeof(symbol_type) * order * 2 + sizeof(int) * 2 + sizeof(impl_type::id_type) * 2;
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
    
    void Bleu::operator()(state_ptr_type& state,
			  feature_set_type& features) const
    {
      // we do nothing...
    }
    
    void Bleu::clear()
    {
      pimpl->clear();
    }
    
    void Bleu::insert(const int source_size, const sentence_type& sentence)
    {
      pimpl->insert(source_size, sentence);
    }

  };
};
