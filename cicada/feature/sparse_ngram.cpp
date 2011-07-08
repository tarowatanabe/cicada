//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "sparse_ngram.hpp"
#include "cicada/parameter.hpp"

#include "utils/piece.hpp"
#include "utils/indexed_trie.hpp"
#include "utils/lexical_cast.hpp"

namespace cicada
{
  namespace feature
  {

    class SparseNGramImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef symbol_type word_type;

      typedef utils::indexed_trie<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> > trie_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > cache_type;
      typedef std::vector<bool, std::allocator<bool> > checked_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > buffer_type;
      
      SparseNGramImpl() : trie(), cache(), checked(), prefix("sparse-ngram"), forced_feature(false) {}
      
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
      
      void ngram_final_score(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features)
      {
	if (no_bos_eos) return;
	
	const trie_type::id_type* context = reinterpret_cast<const trie_type::id_type*>(state);

	buffer.clear();
	buffer.reserve(order * 2 + 2);
	
	buffer.push_back(vocab_type::BOS);
	unpack_context(context[0], buffer);
	
	ngram_feature(buffer.begin(), buffer.begin() + 1, buffer.end(), features);
	
	if (context[1] != trie.root()) {
	  buffer.clear();
	  unpack_context(context[1], buffer);
	}
	buffer.push_back(vocab_type::EOS);
	
	ngram_feature(buffer.begin(), buffer.end() - 1, buffer.end(), features);
      }

      void ngram_score(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       feature_set_type& features)
      {
	if (skip_sgml_tag)
	  ngram_score(state, states, edge, features, skipper_sgml());
	else
	  ngram_score(state, states, edge, features, skipper_epsilon());
      }
      
      template <typename Skipper>
      void ngram_score(state_ptr_type& state,
		       const state_ptr_set_type& states,
		       const edge_type& edge,
		       feature_set_type& features,
		       Skipper skipper)
      {
	const int context_size = order - 1;
	const rule_type& rule = *(edge.rule);
	const phrase_type& phrase = rule.rhs;
	
	phrase_type::const_iterator titer_begin = phrase.begin();
	phrase_type::const_iterator titer_end   = phrase.end();

	buffer.clear();
	buffer.reserve((titer_end - titer_begin) + (order * 2) * states.size());

	if (states.empty()) {
	  // we will copy to buffer...
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	    if (! skipper(*titer))
	      buffer.push_back(*titer);
	  
	  ngram_feature(buffer.begin(), buffer.end(), features);
	  
	  trie_type::id_type* context = reinterpret_cast<trie_type::id_type*>(state);
	  
	  if (static_cast<int>(buffer.size()) <= context_size) {
	    context[0] = ngram_context(buffer.begin(), buffer.end());
	    context[1] = trie.root();
	  } else {
	    context[0] = ngram_context(buffer.begin(), buffer.begin() + context_size);
	    context[1] = ngram_context(buffer.end() - context_size, buffer.end());
	  }
	} else {
	  int star_first = -1;
	  int star_last  = -1;

	  buffer_type::iterator biter_first = buffer.begin();
	  buffer_type::iterator biter       = buffer.begin();
	  
	  int non_terminal_pos = 0;
	  for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	    if (titer->is_non_terminal()) {
	      const int __non_terminal_index = titer->non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      const trie_type::id_type* context_antecedent = reinterpret_cast<const trie_type::id_type*>(states[antecedent_index]);
	      
	      if (biter != buffer.end()) {
		if (biter_first != biter)
		  ngram_feature(biter_first, biter, buffer.end(), features);
		ngram_feature(biter, buffer.end(), features);
		biter = buffer.end();
	      }
	      
	      // uncover antecedent_context into buffer.
	      unpack_context(context_antecedent[0], buffer);
	      
	      if (biter_first != biter && biter != buffer.end())
		ngram_feature(biter_first, biter, buffer.end(), features);
	      biter = buffer.end();
	      
	      if (context_antecedent[1] != trie.root()) {
		star_first = utils::bithack::branch(star_first < 0, static_cast<int>(buffer.size()), star_first);
		star_last  = buffer.size();
		
		biter_first = buffer.end();
		
		// uncover suffix context
		unpack_context(context_antecedent[1], buffer);
		
		biter = buffer.end();
	      }
	    } else if (! skipper(*titer))
	      buffer.push_back(*titer);
	  }
	  
	  if (biter != buffer.end()) {
	    if (biter_first != biter)
	      ngram_feature(biter_first, biter, buffer.end(), features);
	    ngram_feature(biter, buffer.end(), features);
	    biter = buffer.end();
	  }
	  
	  trie_type::id_type* context = reinterpret_cast<trie_type::id_type*>(state);

	  if (star_first >= 0) {
	    const int prefix_size = utils::bithack::min(star_first, context_size);
	    const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
	    
	    context[0] = ngram_context(buffer.begin(), buffer.begin() + prefix_size);
	    context[1] = ngram_context(buffer.end() - suffix_size, buffer.end());
	  } else if (buffer.size() <= context_size) {
	    context[0] = ngram_context(buffer.begin(), buffer.end());
	    context[1] = trie.root();
	  } else {
	    context[0] = ngram_context(buffer.begin(), buffer.begin() + context_size);
	    context[1] = ngram_context(buffer.end() - context_size, buffer.end());
	  }
	  
	}
      }

      template <typename Iterator>
      const feature_type& ngram_feature(Iterator first, Iterator iter, Iterator last, feature_set_type& features)
      {
	const int context_size = order - 1;
	
	first = std::max(iter - context_size, first);
	for (/**/; first != iter; ++ first) {
	  trie_type::id_type id = trie.root();
	  
	  Iterator end = std::min(first + order, last);
	  for (Iterator iter2 = first; iter2 != end; ++ iter2) {
	    id = traverse(id, first, iter2);
	    
	    if (iter2 < iter) continue;
	    
	    if (! cache[id].empty())
	      features[cache[id]] += 1.0;
	  }
	}
      }
      
      template <typename Iterator>
      const feature_type& ngram_feature(Iterator first, Iterator last, feature_set_type& features)
      {
	for (/**/; first != last; ++ first) {
	  trie_type::id_type id = trie.root();
	  
	  Iterator end = std::min(first + order, last);
	  for (Iterator iter = first; iter != end; ++ iter) {
	    id = traverse(id, first, iter);
	    
	    if (! cache[id].empty())
	      features[cache[id]] += 1.0;
	  }
	}
      }

      template <typename Iterator>
      trie_type::id_type ngram_context(Iterator first, Iterator last) {
	trie_type::id_type id = trie.root();
	
	for (Iterator iter = first; iter != last; ++ iter)
	  id = traverse(id, first, iter);
	
	return id;
      }
      
      
      void unpack_context(trie_type::id_type id, buffer_type& buffer)
      {
	// we assume that buffer is allocated with enough space
	
	buffer_type::iterator biter = buffer.end();
	
	while (id != trie.root()) {
	  buffer.push_back(trie[id]);
	  id = trie.pop(id);
	}
	std::reverse(biter, buffer.end());
      }
      
      template <typename Iterator>
      trie_type::id_type traverse(trie_type::id_type id, Iterator first, Iterator iter) {
	id = trie.push(id, *iter);
	
	if (id >= cache.size())
	  cache.resize(id + 1, feature_type());
	if (id >= checked.size())
	  checked.resize(id + 1, false);
	
	if (! checked[id]) {
	  // ngram at [first, iter + 1)
	  
	  std::string name = prefix + ":";
	  for (Iterator fiter = first; fiter != iter + 1; ++ fiter)
	    name += "_" + static_cast<const std::string&>(*fiter);
	  	  
	  if (forced_feature || feature_type::exists(name))
	    cache[id] = name;
	  
	  checked[id] = true;
	}
	
	return id;
      }
      
      void clear()
      {
	trie.clear();
	cache.clear();
	checked.clear();
      }

      buffer_type buffer;
      
      trie_type    trie;
      cache_type   cache;
      checked_type checked;
      
      int order;
      bool no_bos_eos;
      bool skip_sgml_tag;
      
      std::string prefix;
      bool forced_feature;
    };
    
    SparseNGram::SparseNGram(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "sparse-ngram")
	throw std::runtime_error("this is not sparse ngram feature: " + parameter);
      
      int order = 3;
      bool skip_sgml_tag = false;
      bool no_bos_eos = false;
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for sparse ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> ngram_impl(new impl_type());

      ngram_impl->order = order;
      ngram_impl->no_bos_eos = no_bos_eos;
      ngram_impl->skip_sgml_tag = skip_sgml_tag;
      ngram_impl->prefix = (name.empty() ? std::string("sparse-ngram") : name);
      
      // two-side context + length + counts-id 
      base_type::__state_size = sizeof(impl_type::trie_type::id_type) * 2;
      base_type::__feature_name = (name.empty() ? std::string("sparse-ngram") : name);
      base_type::__sparse_feature = true;
      
      pimpl = ngram_impl.release();
    }
    
    SparseNGram::~SparseNGram() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    SparseNGram::SparseNGram(const SparseNGram& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    SparseNGram& SparseNGram::operator=(const SparseNGram& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void SparseNGram::apply(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));

      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->ngram_score(state, states, edge, features);
      if (final)
	pimpl->ngram_final_score(state, states, edge, features);
    }

    void SparseNGram::apply_coarse(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {
      
    }
    
    void SparseNGram::apply_predict(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {}
    void SparseNGram::apply_scan(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 const int dot,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {}
    void SparseNGram::apply_complete(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     feature_set_type& estimates,
				     const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    void SparseNGram::initialize()
    {
      // initialize internal data structs
      pimpl->clear();
    }

  };
};
