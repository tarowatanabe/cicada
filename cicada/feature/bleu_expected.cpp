//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"
#include "utils/bithack.hpp"

#include "cicada/feature/bleu_expected.hpp"
#include "cicada/parameter.hpp"
#include "cicada/eval/bleu.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/semiring.hpp"

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

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
      
      // this implementation specific...
      typedef uint32_t id_type;
      typedef uint16_t count_type;
      typedef double   expected_type;

      typedef symbol_type word_type;
      
      typedef std::allocator<std::pair<const word_type, expected_type> >  ngram_allocator_type;
      typedef utils::compact_trie_dense<word_type, expected_type, boost::hash<word_type>, std::equal_to<word_type>, ngram_allocator_type> ngram_set_type;
      
      struct Node
      {
	Node() : word(), parent(id_type(-1)), order(0) {}
	
	word_type word;
	id_type   parent;
	int       order;
      };
      typedef Node node_type;
      typedef std::vector<node_type, std::allocator<node_type> > node_set_type;
      
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
      BleuExpectedImpl(const int __order)
	: ngrams(word_type()), nodes(), unigram(0),
	  order(__order) {}

      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			const bool final=false) const
      {
	if (ngrams.empty()) {
	  char* buf = reinterpret_cast<char*>(state);
	  std::fill(buf, buf + sizeof(symbol_type) * order * 2 + sizeof(int) + sizeof(id_type), 0);
	  return 0.0;
	}
	
	const rule_type& rule = *edge.rule;
	
	const phrase_type& target = rule.rhs;
	
	count_set_type counts;
	
	symbol_type* context_first = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_last  = context_first + order * 2;
	
	std::fill(context_first, context_last, vocab_type::EMPTY);
	
	int*     context_hypothesis = reinterpret_cast<int*>(context_last);
	id_type* context_count      = reinterpret_cast<id_type*>(context_hypothesis + 1);

	const int context_size = order - 1;
	
	if (states.empty()) {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  
	  phrase_type::const_iterator titer_end = target.end();
	  for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	    if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  
	  collect_counts(buffer.begin(), buffer.end(), counts);
	  
	  if (static_cast<int>(buffer.size()) <= context_size)
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
	  
	  *context_hypothesis = buffer.size();
	  
	  states_count_set_type::iterator citer = const_cast<states_count_set_type&>(states_counts).insert(counts).first;
	  *context_count = citer - const_cast<states_count_set_type&>(states_counts).begin();

	  return bleu_score(counts, *context_hypothesis, ! final);

	} else {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  buffer.reserve(target.size() + (order * 2) * states.size());
	  
	  int star_first = -1;
	  int star_last  = -1;
	  
	  buffer_type::iterator biter_first = buffer.begin();
	  buffer_type::iterator biter       = buffer.begin();

	  double bleu_antecedent = 0.0;
	  
	  int non_terminal_pos = 0;
	  phrase_type::const_iterator titer_end = target.end();
	  for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer) {
	    if (titer->is_non_terminal()) {
	      const int __non_terminal_index = titer->non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;

	      const symbol_type* antecedent_first = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	      const symbol_type* antecedent_last  = antecedent_first + order * 2;
	      
	      const symbol_type* antecedent_end  = std::find(antecedent_first, antecedent_last, vocab_type::EMPTY);
	      const symbol_type* antecedent_star = std::find(antecedent_first, antecedent_end, vocab_type::STAR);
	      
	      const int*     antecedent_hypothesis = reinterpret_cast<const int*>(antecedent_last);
	      const id_type* antecedent_count      = reinterpret_cast<const id_type*>(antecedent_hypothesis + 1);

	      const count_set_type& counts_antecedent = states_counts[*antecedent_count];
	      
	      bleu_antecedent += bleu_score(counts_antecedent, *antecedent_hypothesis, true);

	      // merge statistics...
	      counts.resize(utils::bithack::max(counts.size(), counts_antecedent.size()), count_type(0));
	      std::transform(counts_antecedent.begin(), counts_antecedent.end(), counts.begin(), counts.begin(), std::plus<count_type>());

	      *context_hypothesis += *antecedent_hypothesis;
	      
	      if (biter != buffer.end()) {
		if (biter_first != biter)
		  collect_counts(biter_first, biter, buffer.end(), counts);
		collect_counts(biter, buffer.end(), counts);
		biter = buffer.end();
	      }
	      
	      buffer.insert(buffer.end(), antecedent_first, antecedent_star);
	      if (biter_first != biter && biter != buffer.end())
		collect_counts(biter_first, biter, buffer.end(), counts);
	      biter = buffer.end();
	      
	      if (antecedent_star != antecedent_end) {
		star_last = buffer.size() + 1;
		if (star_first < 0)
		  star_first = buffer.size() + 1;
		
		biter_first = buffer.end() + 1;
		buffer.insert(buffer.end(), antecedent_star, antecedent_end);
		biter = buffer.end();
	      }
	      
	    } else if (*titer != vocab_type::EPSILON) {
	      buffer.push_back(*titer);
	      *context_hypothesis += 1;
	    }
	  }

	  if (biter != buffer.end()) {
	    if (biter_first != biter)
	      collect_counts(biter_first, biter, buffer.end(), counts);
	    collect_counts(biter, buffer.end(), counts);
	    biter = buffer.end();
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
	    if (static_cast<int>(buffer.size()) <= context_size)
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

	  return bleu_score(counts, *context_hypothesis, ! final) - bleu_antecedent;
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
	unigram = 0.0;
	
	// for count set representation
	states_counts.clear();
	
	score.reset();
      }

      void insert(const score_ptr_type& __score)
      {
	score = __score;
	
	if (score && ! dynamic_cast<const cicada::eval::Bleu*>(score.get()))
	  throw std::runtime_error("this is not a bleu-score!");
      }
      
      template <typename Iterator>
      void insert(Iterator first, Iterator last, const double& count)
      {
	if (last - first == 1)
	  unigram += count;

	ngram_set_type::id_type id = ngrams.root();
	int n = 1;
	
	for (/**/; first != last; ++ first, ++ n) {
	  const ngram_set_type::id_type id_next = ngrams.insert(id, *first);
	  
	  if (id_next >= nodes.size())
	    nodes.resize(id_next + 1, node_type());
	  
	  nodes[id_next].word = *first;
	  nodes[id_next].parent = id;
	  nodes[id_next].order  = n;
	  
	  id = id_next;
	}
	
	ngrams[id] = count;
      }

    private:
      template <typename Iterator>
      std::pair<Iterator, Iterator> ngram_prefix(Iterator first, Iterator last) const
      {
	if (first == last || first + 1 == last) return std::make_pair(first, last);
	
	Iterator iter = first;
	ngram_set_type::id_type id = ngrams.root();
	for (/**/; iter != last; ++ iter) {
	  id = ngrams.find(id, *iter);
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
	
	counts.resize(nodes.size(), count_type(0));
	
	// we will collect counts at [iter, last) with context from [first, iter)
	for (/**/; first != iter; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	    id = ngrams.find(id, *iter2);
	    
	    if (ngrams.is_root(id)) break;
	    if (iter2 < iter) continue;

	    counts[id] = utils::bithack::min(1 + counts[id], int(ngrams[id]) + 1);
	  }
	}
      }
      
      template <typename Iterator>
      void collect_counts(Iterator first, Iterator last, count_set_type& counts) const
      {
	counts.resize(nodes.size(), count_type(0));
	
	// we will collect counts at [first, last)
	for (/**/; first != last; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter = first; iter != std::min(first + order, last); ++ iter) {
	    id = ngrams.find(id, *iter);
	    
	    if (ngrams.is_root(id)) break;
	    
	    counts[id] = utils::bithack::min(1 + counts[id], int(ngrams[id]) + 1);
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
      
      double bleu_score(const count_set_type& __counts, const int hypothesis_size, const bool scaling=true) const
      {
	std::vector<double, std::allocator<double> > counts(order);
	for (ngram_set_type::id_type id = 0; id < __counts.size();++ id)
	  counts[nodes[id].order - 1] += std::min(double(__counts[id]), ngrams[id]);
	
	const cicada::eval::Bleu* __bleu = (score ? dynamic_cast<const cicada::eval::Bleu*>(score.get()) : 0);
	
	if (__bleu && __bleu->length_reference > 0) {
	  const double hypothesis_length = tst_size(hypothesis_size, scaling);
	  const double reference_length  = ref_size(hypothesis_length);
	  
	  double smooth = 0.5;
	  double bleu = brevity_penalty(hypothesis_length + __bleu->length_hypothesis, reference_length + __bleu->length_reference);
	  
	  const size_t ngram_size = utils::bithack::min(int(counts.size()), hypothesis_size);
	  const size_t bleu_order = utils::bithack::max(counts.size(), __bleu->ngrams_hypothesis.size());
	  
	  const double factor = 1.0 / order;
	  for (size_t n = 1; n <= bleu_order; ++ n) {
	    const double count = (double(n <= ngram_size ? counts[n - 1] : 0.0)
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
	  
	  const double hypothesis_length = tst_size(hypothesis_size, scaling);
	  const double reference_length  = ref_size(hypothesis_length);
	  
	  double smooth = 0.5;
	  double bleu = brevity_penalty(hypothesis_length, reference_length);
	  
	  const int ngram_size = utils::bithack::min(int(counts.size()), hypothesis_size);
	  
	  const double factor = 1.0 / order;
	  for (int n = 1; n <= ngram_size; ++ n) {
	    const double& count = counts[n - 1];
	    
	    bleu += std::log((count > 0.0 ? count : smooth) / (hypothesis_size + 1 - n)) * factor;
	    smooth *= 0.5;
	  }
	  
	  return std::exp(bleu);
	}
      }

      double tst_size(int length, const bool scaling=true) const
      {
	if (length == 0) return 0.0;
	
	return (length < unigram && scaling ? unigram : double(length));
      }
      
      double ref_size(const double hypothesis_size) const
      {
	return unigram;
      }
      
    public:
      
      buffer_type          buffer_impl;
      
      ngram_set_type ngrams;
      node_set_type  nodes;
      count_type     unigram;
      
      states_count_set_type states_counts;
      
      score_ptr_type score;

      int order;
    };
    
    
    BleuExpected::BleuExpected(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "bleu-expected")
	throw std::runtime_error("this is not BleuExpected feature: " + parameter);

      int order = 4;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order));
      
      // two-side context + length + counts-id
      base_type::__state_size = sizeof(symbol_type) * order * 2 + sizeof(int) + sizeof(impl_type::id_type);
      base_type::__feature_name = "bleu-expected";
      
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
    void BleuExpected::apply_predict(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     feature_set_type& estimates,
				     const bool final) const
    {}
    void BleuExpected::apply_scan(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  const int dot,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {}
    void BleuExpected::apply_complete(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      feature_set_type& estimates,
				      const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void BleuExpected::initialize()
    {
      pimpl->initialize();
    }
    
    void BleuExpected::clear()
    {
      pimpl->clear();
    }

    void BleuExpected::assign(const size_type& id,
			      const hypergraph_type& hypergraph,
			      const lattice_type& lattice,
			      const span_set_type& spans,
			      const sentence_set_type& targets,
			      const ngram_count_set_type& ngram_counts)
    {
      pimpl->clear();

      ngram_count_set_type::const_iterator niter_end = ngram_counts.end();
      for (ngram_count_set_type::const_iterator niter = ngram_counts.begin(); niter != niter_end; ++ niter) {
	//std::cerr << "ngram: " << niter->first << " count: " << niter->second << std::endl;
	
	pimpl->insert(niter->first.begin(), niter->first.end(), niter->second);
      }

      //std::cerr << "unigram: " << pimpl->unigram << std::endl;
    }
    

    void BleuExpected::assign(const score_ptr_type& score)
    {
      pimpl->insert(score);
    }

  };
};
