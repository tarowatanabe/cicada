//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "utils/space_separator.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"

#include "cicada/feature/depeval.hpp"
#include "cicada/parameter.hpp"
#include "cicada/eval/depeval.hpp"
#include "cicada/semiring.hpp"
#include "cicada/sentence_vector.hpp"

#include "cicada/tokenizer.hpp"

#include <boost/numeric/conversion/bounds.hpp>
#include <boost/tokenizer.hpp>

namespace cicada
{
  namespace feature
  {
    class DepevalImpl
    {
      
    };
    
    
    Depeval::Depeval(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "depeval")
	throw std::runtime_error("this is not Depeval feature: " + parameter);

      bool skip_sgml_tag = false;
      const cicada::Tokenizer* tokenizer = 0;
      
      std::string name;
      path_type   refset_file;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "tokenizer")
	  tokenizer = &cicada::Tokenizer::create(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "refset")
	  refset_file = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for depeval: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! refset_file.empty() && ! boost::filesystem::exists(refset_file))
	throw std::runtime_error("no refset file?: " + refset_file.string());
      
      std::auto_ptr<impl_type> depeval_impl(new impl_type(skip_sgml_tag, tokenizer));
      
      // matched count + total count
      base_type::__state_size = sizeof(int) * 2;
      base_type::__feature_name = (name.empty() ? std::string("depeval") : name);
      
      pimpl = depeval_impl.release();
      
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
    
    Depeval::~Depeval() { std::auto_ptr<impl_type> tmp(pimpl); }

    Depeval::Depeval(const Depeval& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    Depeval& Depeval::operator=(const Depeval& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Depeval::apply(state_ptr_type& state,
		     const state_ptr_set_type& states,
		     const edge_type& edge,
		     feature_set_type& features,
		     feature_set_type& estimates,
		     const bool final) const
    {
      double score = pimpl->depeval_score(state, states, edge, final);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
    }

    void Depeval::apply_coarse(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
    {
      
    }
    void Depeval::apply_predict(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Depeval::apply_scan(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  const int dot,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {}
    void Depeval::apply_complete(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void Depeval::initialize()
    {
      pimpl->initialize();
    }
    
    void Depeval::clear()
    {
      pimpl->clear();
    }

    void Depeval::assign(const size_type& id,
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
    
    
    void Depeval::assign(const score_ptr_type& score)
    {
      pimpl->insert(score);
    }

  };
};
