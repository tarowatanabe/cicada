//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <memory>

#include "cicada/feature/tree_rnn.hpp"
#include "cicada/tree_rnn.hpp"
#include "cicada/parameter.hpp"
#include "cicada/symbol_vector.hpp"

#include "utils/array_power2.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/small_vector.hpp"

namespace cicada
{
  namespace feature
  {
    class TreeRNNImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef cicada::Symbol symbol_type;
      typedef cicada::Vocab  vocab_type;
      
      typedef cicada::TreeRNN      rnn_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;

      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > name_set_type;

      typedef rnn_type::parameter_type parameter_type;
      typedef rnn_type::matrix_type    matrix_type;
      typedef rnn_type::tensor_type    tensor_type;

      typedef rnn_type::path_type path_type;

      typedef rule_type::word_type word_type;
      
      typedef std::vector<parameter_type, std::allocator<parameter_type> > buffer_type;
      typedef std::vector<word_type, std::allocator<word_type> > word_set_type;
      
    public:
      TreeRNNImpl(const path_type& path, const std::string& name)
	: rnn(path), no_bos_eos(false), skip_sgml_tag(false)
      {
	initialize(name);
      }
      
      TreeRNNImpl(const size_type& hidden, const size_type& embedding, const path_type& path, const std::string& name)
	: rnn(hidden, embedding, path), no_bos_eos(false), skip_sgml_tag(false)
      {
	initialize(name);
      }
      
      TreeRNNImpl(const TreeRNNImpl& x)
	: rnn(x.rnn),
	  no_bos_eos(x.no_bos_eos),
	  skip_sgml_tag(x.skip_sgml_tag),
	  feature_names(x.feature_names)
	  
      {
      }

      TreeRNNImpl& operator=(const TreeRNNImpl& x)
      {
	rnn = x.rnn;
	
	no_bos_eos    = x.no_bos_eos;
	skip_sgml_tag = x.skip_sgml_tag;
	feature_names = x.feature_names;
	
	return *this;
      }
      
      void initialize(const std::string& name)
      {
	feature_names.clear();
	for (size_type i = 0; i != rnn.hidden_; ++ i)
	  feature_names.push_back(name + ':' + utils::lexical_cast<std::string>(i));
      }
      
      struct extract_word
      {
	const symbol_type& operator()(const symbol_type& word) const
	{
	  return word;
	}
      };

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
      
      void rnn_score(state_ptr_type& state,
		     const state_ptr_set_type& states,
		     const edge_type& edge,
		     feature_set_type& features,
		     const bool final) const
      {
	if (final && ! no_bos_eos) {
	  word_set_type& words = const_cast<word_set_type&>(words_tmp);
	  
	  words.clear();
	  words.push_back(vocab_type::BOS);
	  words.insert(words.end(), edge.rule->rhs.begin(), edge.rule->rhs.end());
	  words.push_back(vocab_type::EOS);
	  
	  if (skip_sgml_tag)
	    rnn_score(state, states, edge, words.begin(), words.end(),  features, extract_word(), skipper_sgml());
	  else
	    rnn_score(state, states, edge, words.begin(), words.end(), features, extract_word(), skipper_epsilon());
	} else {
	  if (skip_sgml_tag)
	    rnn_score(state, states, edge, edge.rule->rhs.begin(), edge.rule->rhs.end(), features, extract_word(), skipper_sgml());
	  else
	    rnn_score(state, states, edge, edge.rule->rhs.begin(), edge.rule->rhs.end(), features, extract_word(), skipper_epsilon());
	}
      }

      
      template <typename Iterator, typename Extract, typename Skipper>
      void rnn_score(state_ptr_type& state,
		     const state_ptr_set_type& states,
		     const edge_type& edge,
		     Iterator first, Iterator last, 
		     feature_set_type& features,
		     Extract extract,
		     Skipper skipper) const
      {
	const size_type offset1 = 0;
	const size_type offset2 = rnn.hidden_;
	
	const_cast<buffer_type&>(buffer_tmp).resize(rnn.hidden_);
	
	parameter_type* pointer_curr = const_cast<parameter_type*>(&(*buffer_tmp.begin()));
	parameter_type* pointer_next = reinterpret_cast<parameter_type*>(state);
	
	if (states.empty()) {
	  bool is_initial = true;
	  for (/**/; first != last; ++ first)
	    if (! skipper(*first)) {
	      matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
	      matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
	      
	      const word_type word = extract(*first);

	      if (is_initial)
		buffer_next = (rnn.Wp_ * rnn.input_->operator()(word) + rnn.Bp_).array().unaryExpr(rnn_type::shtanh());
	      else
		buffer_next = (rnn.Bt_
			       + rnn.Wt_.block(0, offset1, rnn.hidden_, rnn.hidden_) * buffer_curr
			       + rnn.Wt_.block(0, offset2, rnn.hidden_, rnn.embedding_) * rnn.input_->operator()(word)
			       ).array().unaryExpr(rnn_type::shtanh());
	      
	      std::swap(pointer_curr, pointer_next);
	      is_initial = false;
	    }
	} else {
	  if (states.size() == 1 && std::distance(first, last) == 1) {
	    // special handling of unary rules

	    matrix_type buffer_prev(const_cast<parameter_type*>(reinterpret_cast<const parameter_type*>(states[0])), rnn.hidden_, 1);
	    matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
	    
	    buffer_next = (rnn.Bu_ + rnn.Wu_ * buffer_prev).array().unaryExpr(rnn_type::shtanh());
	    
	    std::swap(pointer_curr, pointer_next);
	  } else {
	    bool is_initial = true;
	    int non_terminal_pos = 0;
	    for (/**/; first != last; ++ first) {
	      if (first->is_non_terminal()) {
		matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
		matrix_type buffer_next(pointer_next, rnn.hidden_, 1);

		const int __non_terminal_index = first->non_terminal_index();
		const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0,
								    non_terminal_pos,
								    __non_terminal_index - 1);
		++ non_terminal_pos;

		matrix_type buffer_prev(const_cast<parameter_type*>(reinterpret_cast<const parameter_type*>(states[antecedent_index])),
					rnn.hidden_, 1);
		
		if (is_initial)
		  buffer_next = buffer_prev;
		else 
		  buffer_next = (rnn.Bn_
				 + rnn.Wn_.block(0, offset1, rnn.hidden_, rnn.hidden_) * buffer_curr
				 + rnn.Wn_.block(0, offset2, rnn.hidden_, rnn.hidden_) * buffer_prev
				 ).array().unaryExpr(rnn_type::shtanh());
		
		std::swap(pointer_curr, pointer_next);
		is_initial = false;
	      } else if (! skipper(*first)) {
		matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
		matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
		
		const word_type word = extract(*first);

		if (is_initial)
		  buffer_next = (rnn.Wp_ * rnn.input_->operator()(word) + rnn.Bp_).array().unaryExpr(rnn_type::shtanh());
		else
		  buffer_next = (rnn.Bt_
				 + rnn.Wt_.block(0, offset1, rnn.hidden_, rnn.hidden_) * buffer_curr
				 + rnn.Wt_.block(0, offset2, rnn.hidden_, rnn.embedding_) * rnn.input_->operator()(word)
				 ).array().unaryExpr(rnn_type::shtanh());
		
		std::swap(pointer_curr, pointer_next);
		is_initial = false;
	      }
	    }
	  }
	}
	
	// copy into state buffer when necessary..
	if (pointer_curr != reinterpret_cast<parameter_type*>(state)) {
	  matrix_type buffer_curr(pointer_curr, rnn.hidden_, 1);
	  matrix_type buffer_next(pointer_next, rnn.hidden_, 1);
	  
	  buffer_next = buffer_curr;
	}
	
	// add features...
	features.reserve(features.size() + rnn.hidden_);
	
	matrix_type buffer(reinterpret_cast<parameter_type*>(state), rnn.hidden_, 1);
	for (size_type i = 0; i != rnn.hidden_; ++ i)
	  if (buffer(i, 0) != parameter_type(0)) 
	    features[feature_names[i]] = buffer(i, 0);
      }
      
    public:
      rnn_type rnn;

      buffer_type buffer_tmp;
      word_set_type words_tmp;
      
      bool no_bos_eos;
      bool skip_sgml_tag;
      
      // names...
      name_set_type feature_names;
    };
    
    TreeRNN::TreeRNN(const std::string& parameter)
      : pimpl(0)
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "tree-rnn")
	throw std::runtime_error("is this really tree-rnn feature function? " + parameter);
      
      path_type   path;
      path_type   path_embedding;
      size_type   hidden = 0;
      size_type   embedding = 0;
      bool        skip_sgml_tag = false;
      bool        no_bos_eos = false;
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "embedding-file")
	  path_embedding = piter->second;
	else if (utils::ipiece(piter->first) == "dimension-hidden")
	  hidden = utils::lexical_cast<size_type>(piter->second);
	else if (utils::ipiece(piter->first) == "dimension-embedding")
	  embedding = utils::lexical_cast<size_type>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-bos-eos")
	  no_bos_eos = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for ngram: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! path.empty() && ! path_embedding.empty())
	throw std::runtime_error("either one of model file or embedding file");
      
      if (! path.empty()) {
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error("no model file? " + path.string());
	
      } else if (! path_embedding.empty()) {
	if (! boost::filesystem::exists(path_embedding))
	  throw std::runtime_error("no embedding file? " + path_embedding.string());
	
	if (hidden == 0 || embedding == 0)
	  throw std::runtime_error("invalid dimension");
	
      } else
	throw std::runtime_error("no model file nor embedding file?");
      
      if (name.empty())
	name = "tree-rnn";
      
      std::auto_ptr<impl_type> rnn_impl(! path.empty()
					? new impl_type(path, name)
					: new impl_type(hidden, embedding, path_embedding, name));
      
      rnn_impl->no_bos_eos    = no_bos_eos;
      rnn_impl->skip_sgml_tag = skip_sgml_tag;
      
      base_type::__state_size = sizeof(tree_rnn_type::parameter_type) * rnn_impl->rnn.hidden_;
      base_type::__feature_name = name;
      
      pimpl = rnn_impl.release();
    }
    
    TreeRNN::~TreeRNN()
    {
      std::auto_ptr<impl_type> tmp(pimpl);
    }
    
    TreeRNN::TreeRNN(const TreeRNN& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    TreeRNN::tree_rnn_type& TreeRNN::model() const
    {
      return const_cast<tree_rnn_type&>(pimpl->rnn);
    }
    
    TreeRNN& TreeRNN::operator=(const TreeRNN& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    
    void TreeRNN::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			const bool final) const
    {
      pimpl->rnn_score(state, states, edge, features, final);
    }

    void TreeRNN::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       const bool final) const
    { }
    
    void TreeRNN::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const
    { }
    
    void TreeRNN::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     const bool final) const
    { }
    
    void TreeRNN::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const
    {
      apply(state, states, edge, features, final);
    }
  };
};
