//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/math/special_functions/expm1.hpp>
#include <boost/filesystem.hpp>

#include "lexicon.hpp"

#include "cicada/parameter.hpp"
#include "cicada/lexicon.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/mathop.hpp"
#include "utils/compact_set.hpp"

namespace cicada
{
  namespace feature
  {
    class LexiconImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef boost::filesystem::path path_type;

      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

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
      
      typedef utils::compact_set<word_type,
				 utils::unassigned<word_type>, utils::unassigned<word_type>,
				 boost::hash<word_type>, std::equal_to<word_type>,
				 std::allocator<word_type> > word_set_type;
      struct score_set_type
      {
	double model1;
	double viterbi;
	double noisy_or;
	
	score_set_type() 
	  : model1(- std::numeric_limits<double>::infinity()),
	    viterbi(- std::numeric_limits<double>::infinity()),
	    noisy_or(- std::numeric_limits<double>::infinity()) {}
      };
      typedef std::vector<score_set_type, std::allocator<score_set_type> > cache_set_type;
      
      typedef cicada::Lexicon lexicon_type;
      
      LexiconImpl()
	: lexicon(0), uniques(), words(), caches(), skip_sgml_tag(false), unique_source(false), feature_model1(), feature_viterbi(), feature_noisy_or()
      {  }

      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS;
	}
      };
      
      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS || word.is_sgml_tag();
	}
      };
      
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features)
      {
	if (skip_sgml_tag)
	  lexicon_score(edge, features, skipper_sgml());
	else
	  lexicon_score(edge, features, skipper_epsilon());
      }
      
      
      template <typename Skipper>
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features,
			 Skipper skipper)
      {
	const phrase_type& phrase = edge.rule->rhs;
	const double inf = std::numeric_limits<double>::infinity();

	double score_model1 = 0.0;
	double score_viterbi = 0.0;
	double score_noisy_or = 0.0;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) 
	  if (piter->is_terminal() && ! skipper(*piter)) {
	    const symbol_type& target = *piter;

	    if (target.id() >= caches.size())
	      caches.resize(target.id() + 1);
	    
	    if (caches[target.id()].model1 == - inf) {
	      //
	      // 1.0 - exp(score) == - expm1(score)
	      //
	      
	      double score_model1   = lexicon->operator()(vocab_type::EPSILON, target);
	      double score_viterbi  = score_model1;
	      double score_noisy_or = 0.0;
	      
	      sentence_type::const_iterator siter_end = words.end();
	      for (sentence_type::const_iterator siter = words.begin(); siter != siter_end; ++ siter) {
		const double score = lexicon->operator()(*siter, target);
		
		score_model1   += score;
		score_viterbi   = std::max(score_viterbi, score);
		score_noisy_or += utils::mathop::log(1.0 - score);
	      }
	      
	      caches[target.id()].model1   = utils::mathop::log(score_model1);
	      caches[target.id()].viterbi  = utils::mathop::log(score_viterbi);
	      caches[target.id()].noisy_or = utils::mathop::log(- boost::math::expm1(score_noisy_or));
	    }
	    
	    score_model1   += caches[target.id()].model1;
	    score_viterbi  += caches[target.id()].viterbi;
	    score_noisy_or += caches[target.id()].noisy_or;
	  }
	
	if (score_model1 != 0.0)
	  features[feature_model1] = score_model1;
	else
	  features.erase(feature_model1);

	if (score_viterbi != 0.0)
	  features[feature_viterbi] = score_viterbi;
	else
	  features.erase(feature_viterbi);
	
	if (score_noisy_or != 0.0)
	  features[feature_noisy_or] = score_noisy_or;
	else
	  features.erase(feature_noisy_or);
      }

      void assign(const lattice_type& lattice)
      {
	if (skip_sgml_tag)
	  assign(lattice, skipper_sgml());
	else
	  assign(lattice, skipper_epsilon());
      }
      
      template <typename Skipper>
      void assign(const lattice_type& lattice, Skipper skipper)
      {
	clear();
	
	if (unique_source) {
	  lattice_type::const_iterator liter_end = lattice.end();
	  for (lattice_type::const_iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	    lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	    for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	      if (! skipper(aiter->label))
		uniques.insert(aiter->label);
	  }
	  
	  word_set_type::const_iterator uiter_end = uniques.end();
	  for (word_set_type::const_iterator uiter = uniques.begin(); uiter != uiter_end; ++ uiter)
	    words.push_back(*uiter);
	  
	} else {
	  lattice_type::const_iterator liter_end = lattice.end();
	  for (lattice_type::const_iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	    lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	    for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	      if (! skipper(aiter->label))
		words.push_back(aiter->label);
	  }
	}
	
	std::sort(words.begin(), words.end());
      }
      
      void clear()
      {
	uniques.clear();
	words.clear();
	caches.clear();
      }
      
      lexicon_type* lexicon;
      
      word_set_type  uniques;
      sentence_type  words;
      cache_set_type caches;
      
      bool skip_sgml_tag;
      bool unique_source;
      
      feature_type feature_model1;
      feature_type feature_viterbi;
      feature_type feature_noisy_or;
    };
    
    Lexicon::Lexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "lexicon")
	throw std::runtime_error("this is not a lexicon feature: " + parameter);

      std::string path;
      
      bool populate = false;
      bool skip_sgml_tag = false;
      bool unique_source = false;
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "populate")
	  populate = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "unique-source")
	  unique_source = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (path.empty())
	throw std::runtime_error("no lexicon file? " + path);
      
      std::unique_ptr<impl_type> lexicon_impl(new impl_type());
      
      lexicon_impl->skip_sgml_tag = skip_sgml_tag;
      lexicon_impl->unique_source = unique_source;
      lexicon_impl->feature_model1   = (name.empty() ? std::string("lexicon") : name) + ":model1";
      lexicon_impl->feature_viterbi  = (name.empty() ? std::string("lexicon") : name) + ":viterbi";
      lexicon_impl->feature_noisy_or = (name.empty() ? std::string("lexicon") : name) + ":noisy-or";
      
      // open!
      lexicon_impl->lexicon = &cicada::Lexicon::create(path);

      if (! lexicon_impl->lexicon)
	throw std::runtime_error("no lexicon");

      if (populate)
	lexicon_impl->lexicon->populate();
      
      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("lexicon") : name);
      
      pimpl = lexicon_impl.release();
    }
    
    Lexicon::~Lexicon() { std::unique_ptr<impl_type> tmp(pimpl); }
    
    Lexicon::Lexicon(const Lexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {
      pimpl->lexicon = 0;
      if (x.pimpl->lexicon)
	pimpl->lexicon = &cicada::Lexicon::create(x.pimpl->lexicon->path().string());
    }

    Lexicon& Lexicon::operator=(const Lexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      pimpl->lexicon = 0;
      if (x.pimpl->lexicon)
	pimpl->lexicon = &cicada::Lexicon::create(x.pimpl->lexicon->path().string());
      
      return *this;
    }
    
    void Lexicon::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			const bool final) const
    {
      pimpl->lexicon_score(edge, features);
    }

    void Lexicon::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void Lexicon::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void Lexicon::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     const bool final) const
    {}
    
    void Lexicon::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const
    {}
    
    void Lexicon::assign(const size_type& id,
			 const hypergraph_type& hypergraph,
			 const lattice_type& lattice,
			 const span_set_type& spans,
			 const sentence_set_type& targets,
			 const ngram_count_set_type& ngram_counts)
    {
      //
      // how do we assign lexion feature from hypergraph...???
      // we assume that the lattice is always filled with the source-word...!
      //
      
      pimpl->clear();
      if (! lattice.empty())
	pimpl->assign(lattice);
    }
    
  };
};
