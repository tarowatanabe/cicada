//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <map>

#include "utils/space_separator.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"

#include "cicada/feature/bleu.hpp"
#include "cicada/parameter.hpp"
#include "cicada/eval/bleu.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/semiring.hpp"
#include "cicada/sentence_vector.hpp"

#include "cicada/tokenizer.hpp"

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/tokenizer.hpp>


namespace cicada
{
  namespace feature
  {

    class BleuImpl
    {
    public:
      typedef cicada::Symbol         symbol_type;
      typedef cicada::Vocab          vocab_type;
      typedef cicada::Sentence       sentence_type;
      typedef cicada::SentenceVector sentence_set_type;

      typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;

      typedef cicada::eval::Score score_type;
      typedef score_type::score_ptr_type score_ptr_type;
      
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
	       const bool __exact,
	       const tokenizer_type* __tokenizer)
	: ngrams(word_type()), nodes(), sizes(),
	  order(__order), exact(__exact), tokenizer(__tokenizer)
      {
	
      }

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
	
	const phrase_type& __target = rule.rhs;
	
	phrase_type __target_tokenized;
	if (tokenizer)
	  tokenizer(__target, __target_tokenized);
	const phrase_type& target = (tokenizer ? __target_tokenized : __target);
	
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

	  const double bleu = bleu_score(counts, *context_hypothesis, ! final);
	  
	  //std::cerr << "bleu: " << bleu << ' ';
	  //std::copy(buffer.begin(), buffer.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
	  //std::cerr << std::endl;
	  
	  return bleu;

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
	      int antecedent_index = titer->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = non_terminal_pos;
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
	  
	  const double bleu =  bleu_score(counts, *context_hypothesis, ! final) - bleu_antecedent;

	  return bleu;
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
	
	if (exact) {
	  counts.resize(nodes.size(), count_type(0));
	  
	  // we will collect counts at [iter, last) with context from [first, iter)
	  for (/**/; first != iter; ++ first) {
	    ngram_set_type::id_type id = ngrams.root();
	    for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	      id = ngrams.find(id, *iter2);
	      
	      if (ngrams.is_root(id)) break;
	      if (iter2 < iter) continue;

	      counts[id] = utils::bithack::min(count_type(counts[id] + 1), ngrams[id]);
	    }
	  }
	} else {
	  counts.resize(order, count_type(0));

	  // we will collect counts at [iter, last) with context from [first, iter)
	  for (/**/; first != iter; ++ first) {
	    ngram_set_type::id_type id = ngrams.root();
	    for (Iterator iter2 = first; iter2 != std::min(first + order, last); ++ iter2) {
	      id = ngrams.find(id, *iter2);
	      
	      if (ngrams.is_root(id)) break;
	      if (iter2 < iter) continue;

	      ++ counts[nodes[id].order - 1];
	    }
	  }
	}
      }
      
      template <typename Iterator>
      void collect_counts(Iterator first, Iterator last, count_set_type& counts) const
      {
	if (exact) {
	  counts.resize(nodes.size(), count_type(0));
	  
	  // we will collect counts at [first, last)
	  for (/**/; first != last; ++ first) {
	    ngram_set_type::id_type id = ngrams.root();
	    for (Iterator iter = first; iter != std::min(first + order, last); ++ iter) {
	      id = ngrams.find(id, *iter);
	    
	      if (ngrams.is_root(id)) break;
	      
	      counts[id] = utils::bithack::min(count_type(counts[id] + 1), ngrams[id]);
	    }
	  }
	} else {
	  counts.resize(order, count_type(0));

	  // we will collect counts at [first, last)
	  for (/**/; first != last; ++ first) {
	    ngram_set_type::id_type id = ngrams.root();
	    for (Iterator iter = first; iter != std::min(first + order, last); ++ iter) {
	      id = ngrams.find(id, *iter);
	      
	      if (ngrams.is_root(id)) break;
	      
	      ++ counts[nodes[id].order - 1];
	    }
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
	count_set_type counts_bleu(order, count_type(0));
	if (exact)
	  for (ngram_set_type::id_type id = 0; id < __counts.size();++ id)
	    counts_bleu[nodes[id].order - 1] += utils::bithack::min(int(__counts[id]), int(ngrams[id]));
	
	const count_set_type& counts = (exact ? counts_bleu : __counts);

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
	  
	  const double hypothesis_length = tst_size(hypothesis_size, scaling);
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

      double tst_size(int length, const bool scaling=true) const
      {
	if (length == 0 || sizes.empty()) return length;
	
	const int mask = int(scaling) - 1;
	
	return ((~mask) & utils::bithack::max(length, sizes.front())) | (mask & length);
	//return (scaling ? utils::bithack::max(length, sizes.front()) : length);
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
      sentence_document_type refset;
      
      buffer_type          buffer_impl;

      ngram_set_type ngrams;
      node_set_type  nodes;
      size_set_type  sizes;
      
      states_count_set_type states_counts;
      
      score_ptr_type score;

      int order;
      bool exact;

      tokenizer_wrapper_type tokenizer;
    };
    
    
    Bleu::Bleu(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "bleu")
	throw std::runtime_error("this is not Bleu feature: " + parameter);

      int order = 4;
      bool exact = false;

      const cicada::Tokenizer* tokenizer = 0;
      
      std::string name;
      path_type   refset_file;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = boost::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "exact")
	  exact = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "tokenizer")
	  tokenizer = &cicada::Tokenizer::create(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "refset")
	  refset_file = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! refset_file.empty() && ! boost::filesystem::exists(refset_file))
	throw std::runtime_error("no refset file?: " + refset_file.file_string());
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, exact, tokenizer));
      
      // two-side context + length + counts-id 
      base_type::__state_size = sizeof(symbol_type) * order * 2 + sizeof(int) + sizeof(impl_type::id_type);
      base_type::__feature_name = (name.empty() ? std::string("bleu") : name);
      
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
	  
	  const int id = boost::lexical_cast<int>(*iter);
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
    
    Bleu::~Bleu() { std::auto_ptr<impl_type> tmp(pimpl); }

    Bleu::Bleu(const Bleu& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    Bleu& Bleu::operator=(const Bleu& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Bleu::apply(state_ptr_type& state,
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

    void Bleu::apply_coarse(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
    {
      
    }
    void Bleu::apply_predict(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Bleu::apply_scan(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  const int dot,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {}
    void Bleu::apply_complete(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void Bleu::initialize()
    {
      pimpl->initialize();
    }
    
    void Bleu::clear()
    {
      pimpl->clear();
    }

    void Bleu::assign(const size_type& id,
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
    

    void Bleu::assign(const score_ptr_type& score)
    {
      pimpl->insert(score);
    }

  };
};
