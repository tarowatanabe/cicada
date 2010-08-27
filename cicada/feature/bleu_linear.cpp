
#include <map>

#include "cicada/feature/bleu_linear.hpp"
#include "cicada/parameter.hpp"

#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace feature
  {

    class BleuLinearImpl
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
      
      typedef std::vector<double, std::allocator<double> > ngram_factor_type;

    public:
      BleuLinearImpl(const int __order,
		     const double __precision,
		     const double __ratio)
	: ngrams(word_type()), nodes(), sizes(), order(__order), precision(__precision), ratio(__ratio)
      {
	factors.clear();
	factors.resize(order + 1, 0.0);
	
	factors[0] = - 1.0 / 1.0;
	
	double ratio_factor = 1.0;
	for (int n = 1; n <= order; ++ n) {
	  factors[n] = 1.0 / (4.0 * precision * ratio_factor);
	  ratio_factor *= ratio;
	}
      }

      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge) const
      {
	const rule_type& rule = *edge.rule;
	const phrase_type& target = rule.target;
	const phrase_type& source = rule.source;

	symbol_type* context_first = reinterpret_cast<symbol_type*>(state);
	symbol_type* context_last  = context_first + order * 2;
	
	const int context_size = order - 1;
	
	if (states.empty()) {
	  buffer_type& buffer = const_cast<buffer_type&>(buffer_impl);
	  buffer.clear();
	  
	  phrase_type::const_iterator titer_end = target.end();
	  for (phrase_type::const_iterator titer = target.begin(); titer != titer_end; ++ titer)
	    if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  
	  const double bleu = bleu_score(buffer.begin(), buffer.end());
	  
	  std::fill(context_first, context_last, vocab_type::EMPTY);
	  
	  if (buffer.size() <= context_size)
	    std::copy(buffer.begin(), buffer.end(), context_first);
	  else {
	    std::copy(buffer.begin(), buffer.begin() + context_size, context_first);
	    context_first[context_size] = vocab_type::STAR;
	    std::copy(buffer.end() - context_size, buffer.end(), context_first + order);
	  }
	  
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
	  
	  for (phrase_type::const_iterator iter = phrase_spans.front().first; iter != phrase_spans.front().second; ++ iter)
	    if (*iter != vocab_type::EPSILON)
	      buffer.push_back(*iter);
	  
	  double bleu = bleu_score(buffer.begin(), buffer.end());

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
	    
	    buffer_type::const_iterator biter = buffer.end();
	    
	    buffer.insert(buffer.end(), antecedent_first, antecedent_star);
	    
	    buffer_type::const_iterator biter_end = buffer.end();
	    
	    bleu += bleu_score(biter_first, biter, biter_end);
	    
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
	      
	      bleu += bleu_score(biter_first, biter, biter_end);
	    }
	  }
	  
	  if (star_first >= 0) {
	    const int prefix_size = utils::bithack::min(star_first, context_size);
	    const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	    
	    std::copy(buffer.begin(), buffer.begin() + prefix_size, context_first);
	    context_first[prefix_size] = vocab_type::STAR;
	    std::copy(buffer.end() - suffix_size, buffer.end(), context_first + prefix_size + 1);
	  } else {
	    if (buffer.size() <= context_size)
	      std::copy(buffer.begin(), buffer.end(), context_first);
	    else {
	      std::copy(buffer.begin(), buffer.begin() + context_size, context_first);
	      context_first[context_size] = vocab_type::STAR;
	      std::copy(buffer.end() - context_size, buffer.end(), context_first + order);
	    }
	  }
	  
	  return  bleu;
	}
      }

      void initialize()
      {
	
      }
      
      void clear()
      {
	// for ngram handling
	ngrams.clear();
	nodes.clear();
	sizes.clear();
	
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
      double bleu_score(Iterator first, Iterator iter, Iterator last) const
      {
	if (ngrams.empty()) return 0.0;

	double bleu = factors[0] * (last - iter);
	
	const int context_size = order - 1;
	
	first = std::max(first, iter - context_size);
		
	// we will collect counts at [iter, last) with context from [first, iter)
	for (/**/; first != iter; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	    id = ngrams.find(id, *iter2);
	    
	    if (ngrams.is_root(id)) break;
	    if (iter2 < iter) continue;
	    
	    bleu += factors[nodes[id].order];
	  }
	}
	
	return bleu;
      }
      
      template <typename Iterator>
      double bleu_score(Iterator first, Iterator last) const
      {
	if (ngrams.empty()) return 0.0;

	double bleu = factors[0] * (last - first);
	
	// we will collect counts at [first, last)
	for (/**/; first != last; ++ first) {
	  ngram_set_type::id_type id = ngrams.root();
	  for (Iterator iter = first; iter != std::min(first + order, last); ++ iter) {
	    id = ngrams.find(id, *iter);
	    
	    if (ngrams.is_root(id)) break;
	    
	    bleu += factors[nodes[id].order];
	  }
	}
	return bleu;
      }
      
      
    private:
      
      buffer_type          buffer_impl;
      phrase_span_set_type phrase_spans_impl;

      ngram_set_type ngrams;
      node_set_type  nodes;
      size_set_type  sizes;
      
      int source_size;

      ngram_factor_type factors;

      int order;
      double precision;
      double ratio;
    };
    
    BleuLinear::BleuLinear(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "bleu-linear")
	throw std::runtime_error("this is not BleuLinear feature: " + parameter);

      int order        = 4;
      double precision = 0.8;
      double ratio     = 0.6;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "order") == 0)
	  order = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "precision") == 0)
	  precision = boost::lexical_cast<double>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "ratio") == 0)
	  ratio = boost::lexical_cast<double>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu-linear: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, precision, ratio));
      
      // two-side context + length (hypothesis/reference) + counts-id (hypothesis/reference)
      base_type::__state_size = sizeof(symbol_type) * order * 2;
      base_type::__feature_name = "bleu-linear";
      
      pimpl = bleu_impl.release();
    }
    
    BleuLinear::~BleuLinear() { std::auto_ptr<impl_type> tmp(pimpl); }

    BleuLinear::BleuLinear(const BleuLinear& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    BleuLinear& BleuLinear::operator=(const BleuLinear& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void BleuLinear::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      const double score = pimpl->bleu_score(state, states, edge);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
    }

    void BleuLinear::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
    }
    
    void BleuLinear::initialize()
    {
      pimpl->initialize();
    }
    
    void BleuLinear::clear()
    {
      pimpl->clear();
    }
    
    void BleuLinear::insert(const int source_size, const sentence_type& sentence)
    {
      pimpl->insert(source_size, sentence);
    }
    
    void BleuLinear::insert(const score_ptr_type& score)
    {
      // do nothing...
    }

  };
};
