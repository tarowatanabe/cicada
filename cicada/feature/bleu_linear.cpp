//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <map>

#include "cicada/feature/bleu_linear.hpp"
#include "cicada/parameter.hpp"
#include "cicada/semiring.hpp"
#include "cicada/tokenizer.hpp"

#include "utils/space_separator.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/tokenizer.hpp>

namespace cicada
{
  namespace feature
  {

    class BleuLinearImpl
    {
    public:
      typedef cicada::Symbol         symbol_type;
      typedef cicada::Vocab          vocab_type;
      typedef cicada::Sentence       sentence_type;
      typedef cicada::SentenceVector sentence_set_type;

      typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;
      
      typedef cicada::FeatureFunction feature_function_type;

      typedef cicada::Tokenizer tokenizer_type;
      
      struct tokenizer_wrapper_type
      {
	tokenizer_wrapper_type(const tokenizer_type* __tokenizer)
	  : tokenizer(__tokenizer) {}

	tokenizer_wrapper_type()
	  : tokenizer(0) {}
	tokenizer_wrapper_type(const tokenizer_wrapper_type& x)
	  : tokenizer(0)
	{
	  if (x.tokenizer)
	    tokenizer = &tokenizer_type::create(x.tokenizer->algorithm());
	}
	tokenizer_wrapper_type& operator=(const tokenizer_wrapper_type& x)
	{
	  tokenizer = 0;
	  if (x.tokenizer)
	    tokenizer = &tokenizer_type::create(x.tokenizer->algorithm());
	  return *this;
	}
	
	template <typename Sent>
	void operator()(const Sent& source, Sent& tokenized) const
	{
	  if (tokenizer)
	    tokenizer->operator()(source, tokenized);
	  else
	    tokenized = source;
	}
	
	operator bool() const { return tokenizer; }
	
	const tokenizer_type* tokenizer;
      };
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
            
      // this implementation specific...
      typedef uint32_t id_type;
      typedef uint16_t count_type;

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
		     const double __ratio,
		     const tokenizer_type* __tokenizer,
		     const bool __skip_sgml_tag)
	: ngrams(word_type()), nodes(), sizes(), order(__order), precision(__precision), ratio(__ratio),
	  tokenizer(__tokenizer),
	  skip_sgml_tag(__skip_sgml_tag)
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

      
      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON;
	}
      };

      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || (word != vocab_type::BOS && word != vocab_type::EOS && word.is_sgml_tag());
	}
      };

      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge) const
      {
	if (skip_sgml_tag)
	  return bleu_score(state, states, edge, skipper_sgml());
	else
	  return bleu_score(state, states, edge, skipper_epsilon());
      }
      
      template <typename Skipper>
      double bleu_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			Skipper skipper) const
      {
	const rule_type& rule = *edge.rule;
	
	const phrase_type& __target = rule.rhs;
	
	phrase_type __target_tokenized;
	if (tokenizer)
	  tokenizer(__target, __target_tokenized);
	const phrase_type& target = (tokenizer ? __target_tokenized : __target);
	
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
	  
	  if (static_cast<int>(buffer.size()) <= context_size)
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

	  int star_first = -1;
	  int star_last  = -1;

	  double bleu = 0;
	  
	  buffer_type::iterator biter_first = buffer.begin();
	  buffer_type::iterator biter       = buffer.begin();
	  
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
	      
	      if (biter != buffer.end()) {
		if (biter_first != biter)
		  bleu += bleu_score(biter_first, biter, buffer.end());
		bleu += bleu_score(biter, buffer.end());
		biter = buffer.end();
	      }
	      
	      buffer.insert(buffer.end(), antecedent_first, antecedent_star);
	      if (biter_first != biter && biter != buffer.end())
		bleu += bleu_score(biter_first, biter, buffer.end());
	      biter = buffer.end();
	      
	      if (antecedent_star != antecedent_end) {
		star_first = utils::bithack::branch(star_first < 0, static_cast<int>(buffer.size()) + 1, star_first);
		star_last  = buffer.size() + 1;
		
		biter_first = buffer.end() + 1;
		buffer.insert(buffer.end(), antecedent_star, antecedent_end);
		biter = buffer.end();
	      }
	      
	    } else if (*titer != vocab_type::EPSILON)
	      buffer.push_back(*titer);
	  }
	  
	  if (biter != buffer.end()) {
	    if (biter_first != biter)
	      bleu += bleu_score(biter_first, biter, buffer.end());
	    bleu += bleu_score(biter, buffer.end());
	    biter = buffer.end();
	  }
	  
	  if (star_first >= 0) {
	    const int prefix_size = utils::bithack::min(star_first, context_size);
	    const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	    
	    std::copy(buffer.begin(), buffer.begin() + prefix_size, context_first);
	    context_first[prefix_size] = vocab_type::STAR;
	    std::copy(buffer.end() - suffix_size, buffer.end(), context_first + prefix_size + 1);
	  } else {
	    if (static_cast<int>(buffer.size()) <= context_size)
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
      }
      
      void insert(const sentence_type& __sentence)
      {
	typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
	sentence_type __sentence_tokenized;
	if (tokenizer)
	  tokenizer(__sentence, __sentence_tokenized);
	const sentence_type& sentence = (tokenizer ? __sentence_tokenized : __sentence);
	
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
      
      
    public:
      sentence_document_type refset;
      
      buffer_type          buffer_impl;

      ngram_set_type ngrams;
      node_set_type  nodes;
      size_set_type  sizes;

      ngram_factor_type factors;

      int order;
      double precision;
      double ratio;
      
      tokenizer_wrapper_type tokenizer;

      bool skip_sgml_tag;
    };
    
    BleuLinear::BleuLinear(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "bleu-linear")
	throw std::runtime_error("this is not BleuLinear feature: " + parameter);

      int order        = 4;
      double precision = 0.8;
      double ratio     = 0.6;
      bool skip_sgml_tag = false;
      
      const cicada::Tokenizer* tokenizer = 0;
      
      std::string name;
      path_type   refset_file;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "precision")
	  precision = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "ratio")
	  ratio = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "tokenizer")
	  tokenizer = &cicada::Tokenizer::create(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "refset")
	  refset_file = piter->second;
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu-linear: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! refset_file.empty() && ! boost::filesystem::exists(refset_file))
	throw std::runtime_error("no refset file?: " + refset_file.string());
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, precision, ratio, tokenizer, skip_sgml_tag));
      
      // two-side context + length (hypothesis/reference) + counts-id (hypothesis/reference)
      base_type::__state_size = sizeof(symbol_type) * order * 2;
      base_type::__feature_name = (name.empty() ? std::string("bleu-linear") : name);
      
      pimpl = bleu_impl.release();

      pimpl->refset.clear();
      
      if (! refset_file.empty()) {
	typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	
	utils::compress_istream is(refset_file, 1024 * 1024);
	std::string line;
	
	while (std::getline(is, line)) {
	  utils::piece line_piece(line);
	  tokenizer_type tokenizer(line_piece);
	  
	  tokenizer_type::iterator iter = tokenizer.begin();
	  if (iter == tokenizer.end()) continue;
	  
	  const int id = utils::lexical_cast<int>(*iter);
	  ++ iter;
	  
	  if (iter == tokenizer.end()) continue;
	  if (*iter != "|||") continue;
	  ++ iter;
	  
	  if (id >= static_cast<int>(pimpl->refset.size()))
	    pimpl->refset.resize(id + 1);
	  
	  pimpl->refset[id].push_back(sentence_type(iter, tokenizer.end()));
	}
      }
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
    void BleuLinear::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
    void BleuLinear::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    void BleuLinear::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void BleuLinear::initialize()
    {
      pimpl->initialize();
    }
    
    void BleuLinear::clear()
    {
      pimpl->clear();
    }
    
    void BleuLinear::assign(const size_type& id,
			    const hypergraph_type& hypergraph,
			    const lattice_type& lattice,
			    const span_set_type& spans,
			    const sentence_set_type& targets,
			    const ngram_count_set_type& ngram_counts)
    {
      pimpl->clear();
      
      if (! targets.empty()) {
	sentence_set_type::const_iterator titer_end = targets.end();
	for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
	  pimpl->insert(*titer);
      } else if (! pimpl->refset.empty()) {
	if (id < pimpl->refset.size()) {
	  sentence_set_type::const_iterator titer_end = pimpl->refset[id].end();
	  for (sentence_set_type::const_iterator titer = pimpl->refset[id].begin(); titer != titer_end; ++ titer)
	    pimpl->insert(*titer);
	}
      } else
	throw std::runtime_error("no reference set?");
    }
    
    void BleuLinear::assign(const score_ptr_type& score)
    {
      // do nothing...
    }

  };
};
