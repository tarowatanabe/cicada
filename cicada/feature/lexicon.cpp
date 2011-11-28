//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/filesystem.hpp>

#include "lexicon.hpp"

#include "cicada/parameter.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/mathop.hpp"

#include <google/dense_hash_set>
#include <google/dense_hash_map>

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
      
      typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;
      typedef std::vector<double, std::allocator<double> > cache_set_type;
      
      struct ttable_type
      {
	struct table_type
	{
	  typedef google::dense_hash_map<word_type, double, boost::hash<word_type> , std::equal_to<word_type> > __table_type;
	  
	  table_type() : table() { table.set_empty_key(word_type()); }
	  
	  __table_type table;
	};
	
	typedef utils::alloc_vector<table_type, std::allocator<table_type> > table_set_type;
	
	ttable_type() : smooth(1e-30), tables() {}
	ttable_type(const path_type& path) : smooth(), tables() { open(path); }
	
	void open(const path_type& path)
	{
	  namespace qi = boost::spirit::qi;
	  namespace standard = boost::spirit::standard;
	  
	  typedef boost::fusion::tuple<std::string, std::string, double> parsed_type;
	  
	  qi::rule<std::string::const_iterator, std::string(), standard::space_type> word;
	  qi::rule<std::string::const_iterator, parsed_type(), standard::space_type> lexicon;
	  
	  word %= qi::lexeme[+(standard::char_ - standard::space)];
	  lexicon %= word >> word >> qi::double_;
	  
	  smooth = std::numeric_limits<double>::infinity();
	  tables.clear();
	  
	  utils::compress_istream is(path, 1024 * 1024);
	  std::string line;
	  parsed_type parsed;
	  
	  while (std::getline(is, line)) {
	    std::string::const_iterator iter = line.begin();
	    std::string::const_iterator end = line.end();
	    
	    boost::fusion::get<0>(parsed).clear();
	    boost::fusion::get<1>(parsed).clear();
	    
	    const bool result = qi::phrase_parse(iter, end, lexicon, standard::space, parsed);
	    if (! result || iter != end) continue;
	    
	    tables[word_type(boost::fusion::get<1>(parsed)).id()].table[boost::fusion::get<0>(parsed)] = boost::fusion::get<2>(parsed);
	    
	    smooth = std::min(smooth, boost::fusion::get<2>(parsed));
	  }
	}
	
	double operator()(const word_type& source, const word_type& target) const
	{
	  if (tables.exists(source.id())) {
	    const table_type& table = tables[source.id()];
	    
	    table_type::__table_type::const_iterator iter = table.table.find(target);
	    return (iter != table.table.end() ? iter->second : smooth);
	  } else 
	    return smooth;
	}
	
	
	double smooth;
	table_set_type tables;
      };

      
      LexiconImpl()
	: uniques(), words(), caches_model1(), caches_noisy_or(), skip_sgml_tag(false), unique_source(false), feature_model1(), feature_noisy_or()
      { uniques.set_empty_key(word_type());  }

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
	double score_noisy_or = 0.0;
	
	phrase_type::const_iterator piter_end = phrase.end();
	for (phrase_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) 
	  if (piter->is_terminal() && ! skipper(*piter)) {
	    const symbol_type& target = *piter;

	    if (target.id() >= caches_model1.size())
	      caches_model1.resize(target.id() + 1, - inf);
	    if (target.id() >= caches_noisy_or.size())
	      caches_noisy_or.resize(target.id() + 1, - inf);
	    
	    if (caches_model1[target.id()] == - inf) {
	      double score = ttable(vocab_type::EPSILON, target);
	      
	      sentence_type::const_iterator siter_end = words.end();
	      for (sentence_type::const_iterator siter = words.begin(); siter != siter_end; ++ siter)
		score += ttable(*siter, target);
	      
	      caches_model1[target.id()] = utils::mathop::log(score);
	    }
	    
	    if (caches_noisy_or[target.id()] == - inf) {
	      //
	      // 1.0 - exp(score) == - expm1(score)
	      //
	      
	      double score = 0.0;
	      sentence_type::const_iterator siter_end = words.end();
	      for (sentence_type::const_iterator siter = words.begin(); siter != siter_end; ++ siter)
		score += utils::mathop::log(1.0 - ttable(*siter, target));
	      
	      caches_noisy_or[target.id()] = utils::mathop::log(- boost::math::expm1(score));
	    }
	    
	    score_model1   += caches_model1[target.id()];
	    score_noisy_or += caches_noisy_or[target.id()];
	  }
	
	if (score_model1 != 0.0)
	  features[feature_model1] = score_model1;
	else
	  features.erase(feature_model1);
	
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
	caches_model1.clear();
	caches_noisy_or.clear();
      }

      ttable_type ttable;
      
      word_set_type  uniques;
      sentence_type  words;
      cache_set_type caches_model1;
      cache_set_type caches_noisy_or;
      
      bool skip_sgml_tag;
      bool unique_source;
      
      feature_type feature_model1;
      feature_type feature_noisy_or;
    };
    
    Lexicon::Lexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "lexicon")
	throw std::runtime_error("this is not sparse lexicon feature: " + parameter);

      path_type path;
      
      bool skip_sgml_tag = false;
      bool unique_source = false;
      
      std::string name;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "file")
	  path = piter->second;
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "unique-source")
	  unique_source = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (path.empty() || ! boost::filesystem::exists(path))
	throw std::runtime_error("no lexicon file? " + path.string());
      
      std::auto_ptr<impl_type> lexicon_impl(new impl_type());
      
      lexicon_impl->skip_sgml_tag = skip_sgml_tag;
      lexicon_impl->unique_source = unique_source;
      lexicon_impl->feature_model1   = (name.empty() ? std::string("lexicon") : name) + ":model1";
      lexicon_impl->feature_noisy_or = (name.empty() ? std::string("lexicon") : name) + ":noisy-or";
      
      // open!
      lexicon_impl->ttable.open(path);
      
      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("lexicon") : name);
      
      pimpl = lexicon_impl.release();
    }
    
    Lexicon::~Lexicon() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    Lexicon::Lexicon(const Lexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    Lexicon& Lexicon::operator=(const Lexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
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
