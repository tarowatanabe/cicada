//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>
#include <vector>

#include "cicada/feature/vocabulary.hpp"
#include "cicada/parameter.hpp"

#include "utils/piece.hpp"
#include "utils/bithack.hpp"
#include "utils/compress_stream.hpp"

#include <boost/filesystem.hpp>

namespace cicada
{
  namespace feature
  {
    class VocabularyImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Sentence sentence_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef std::vector<bool, std::allocator<bool> > vocab_type;
      
      VocabularyImpl() {}
      
      void vocabulary_score(const edge_type& edge,
			    feature_set_type& features) const
      {
	int count = 0;
	if (oov_penalty) {
	  rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	    if (riter->is_terminal() && *riter != cicada::Vocab::BOS && *riter != cicada::Vocab::EOS) {
	      const symbol_type::id_type id = riter->id();
	      
	      if (id >= vocab.size() || ! vocab[id])
		++ count;
	    }
	  
	  if (count)
	    features[feature_name] = - count;
	  else
	    features.erase(feature_name);
	} else {
	  rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	    if (riter->is_terminal() && *riter != cicada::Vocab::BOS && *riter != cicada::Vocab::EOS) {
	      const symbol_type::id_type id = riter->id();
	      
	      if (id < vocab.size() && vocab[id])
		++ count;
	    }
	  
	  if (count)
	    features[feature_name] = count;
	  else
	    features.erase(feature_name);
	}
      }
      
      vocab_type vocab;
      feature_type feature_name;
      bool oov_penalty;
    };
    
    
    Vocabulary::Vocabulary(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "vocabulary")
	throw std::runtime_error("is this really vocabulary feature function? " + parameter);
      
      path_type   path;
      std::string name;
      bool oov_penalty = false;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for vocabulary: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (path.empty() || (path != "-" && ! boost::filesystem::exists(path)))
	throw std::runtime_error("no vocabulary file: " + path.string());
      
      std::auto_ptr<impl_type> vocabulary_impl(new impl_type());

      {
	utils::compress_istream is(path, 1024 * 1024);
	symbol_type word;
	while (is >> word) {
	  const symbol_type::id_type id = word.id();
	  
	  if (id >= vocabulary_impl->vocab.size())
	    vocabulary_impl->vocab.resize(id + 1, false);
	  vocabulary_impl->vocab[id] = true;
	}
	
	impl_type::vocab_type(vocabulary_impl->vocab).swap(vocabulary_impl->vocab);
      }
      
      vocabulary_impl->feature_name = (name.empty() ? std::string("vocabulary") : name);
      vocabulary_impl->oov_penalty = oov_penalty;
      
      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("vocabulary") : name);
      
      pimpl = vocabulary_impl.release();
    }
    
    Vocabulary::~Vocabulary() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    Vocabulary::Vocabulary(const Vocabulary& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    Vocabulary& Vocabulary::operator=(const Vocabulary& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Vocabulary::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      pimpl->vocabulary_score(edge, features);
    }

    void Vocabulary::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void Vocabulary::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void Vocabulary::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    
    void Vocabulary::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {}
  };
};
