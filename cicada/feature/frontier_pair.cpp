//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "frontier_pair.hpp"
#include "feature_builder.hpp"

#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/array_power2.hpp"
#include "utils/small_vector.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/space_separator.hpp"
#include "utils/unordered_map.hpp"

#include <boost/tokenizer.hpp>

namespace cicada
{
  namespace feature
  {


    class FrontierPairImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;


      typedef utils::hashmurmur3<size_t> hasher_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef symbol_type word_type;
      
      struct cache_phrase_type
      {
	std::string frontier;
	std::string phrase;
	
	cache_phrase_type() : frontier(), phrase() {}
      };
      typedef utils::array_power2<cache_phrase_type, 1024 * 4, std::allocator<cache_phrase_type> > cache_phrase_set_type;
      
      typedef std::pair<std::string, std::string> phrase_pair_type;
      
      struct phrase_pair_hash : public utils::hashmurmur3<size_t>
      {
	typedef utils::hashmurmur3<size_t> hasher_type;
	
	size_t operator()(const phrase_pair_type& x) const
	{
	  return hasher_type::operator()(x.first.begin(), x.first.end(), hasher_type::operator()(x.second.begin(), x.second.end(), 0));
	}
      };

      typedef utils::unordered_map<phrase_pair_type, feature_type, phrase_pair_hash, std::equal_to<phrase_pair_type>,
				   std::allocator<std::pair<const phrase_pair_type, feature_type> > >::type cache_feature_set_type;
      
      
      typedef FeatureBuilder feature_builder_type;

      FrontierPairImpl()
	: skip_sgml_tag(false), prefix("frontier-pair"), forced_feature(false),
	  attr_frontier_source("frontier-source"),
	  attr_frontier_target("frontier-target")
      { }
      
      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word.is_terminal() && (word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS);
	}
      };

      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word.is_terminal() && (word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS || word.is_sgml_tag());
	}
      };
      
      void pair_score(const edge_type& edge,
		      feature_set_type& features)
      {
	if (skip_sgml_tag)
	  pair_score(edge, features, skipper_sgml());
	else
	  pair_score(edge, features, skipper_epsilon());
      }
      
      struct __attribute_string : public boost::static_visitor<const cicada::AttributeVector::string_type&>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	static 
	const std::string& empty()
	{
	  static std::string __empty;
	  return __empty;
	}
	
	const attribute_set_type::string_type& operator()(const attribute_set_type::int_type& x) const { return empty(); }
	const attribute_set_type::string_type& operator()(const attribute_set_type::float_type& x) const { return empty(); }
	const attribute_set_type::string_type& operator()(const attribute_set_type::string_type& x) const { return x; }
      };
      
      template <typename Skipper>
      void pair_score(const edge_type& edge,
		      feature_set_type& features,
		      Skipper skipper)
      {
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_frontier_source);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_frontier_target);
	
	if (siter == edge.attributes.end() || titer == edge.attributes.end()) return;
	
	const std::string& frontier_source = boost::apply_visitor(__attribute_string(), siter->second);
	const std::string& frontier_target = boost::apply_visitor(__attribute_string(), titer->second);
	
	if (frontier_source.empty() || frontier_target.empty()) return;
	
	cache_feature_set_type::iterator fiter = cache_features.find(phrase_pair_type(frontier_source, frontier_target));
	
	if (fiter == cache_features.end()) {
	  fiter = cache_features.insert(std::make_pair(phrase_pair_type(frontier_source, frontier_target), feature_type())).first;
	  
	  const std::string& phrase_source = phrase(frontier_source, cache_phrase_source, skipper);
	  const std::string& phrase_target = phrase(frontier_target, cache_phrase_target, skipper);
	  
	  feature_builder.clear();
	  feature_builder << prefix << ":" << phrase_source << "|" << phrase_target;
	  
	  if (forced_feature || feature_builder.exists())
	    fiter->second = feature_builder;
	}
	
	if (! fiter->second.empty())
	  features[fiter->second] += 1.0;
      }
      
      template <typename Skipper>
      const std::string& phrase(const std::string& frontier,
				cache_phrase_set_type& caches,
				Skipper skipper)
      {
	cache_phrase_type& cache = caches[hasher_type::operator()(frontier.begin(), frontier.end(), 0) & (caches.size() - 1)];
	
	if (cache.frontier != frontier) {
	  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	  
	  utils::piece frontier_piece(frontier);
	  tokenizer_type tokenizer(frontier_piece);
	  
	  cache.frontier = frontier;
	  
	  feature_builder.clear();
	  bool initial = true;
	  
	  tokenizer_type::iterator titer_end = tokenizer.end();
	  for (tokenizer_type::iterator titer = tokenizer.begin(); titer != titer_end; ++ titer) {
	    const symbol_type word = *titer;
	    
	    if (! skipper(word)) {
	      if (! initial)
		feature_builder << "_";
	      
	      feature_builder << word;
	      initial = false;
	    }
	  }
	  
	  cache.phrase = feature_builder;
	}
	
	return cache.phrase;
      }
            
      void clear()
      {
	cache_features.clear();
      }
      
      void clear_cache()
      {
	cache_phrase_source.clear();
	cache_phrase_target.clear();

	cache_features.clear();
      }

      cache_phrase_set_type cache_phrase_source;
      cache_phrase_set_type cache_phrase_target;

      cache_feature_set_type cache_features;

      feature_builder_type feature_builder;
      
      bool skip_sgml_tag;
      
      std::string prefix;
      bool forced_feature;

      attribute_type attr_frontier_source;
      attribute_type attr_frontier_target;
    };
    
    FrontierPair::FrontierPair(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "frontier-pair")
	throw std::runtime_error("this is not frontier pair feature: " + parameter);

      bool skip_sgml_tag = false;
      
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for frontier pair: " << piter->first << "=" << piter->second << std::endl;
      }
            
      std::unique_ptr<impl_type> pair_impl(new impl_type());

      pair_impl->skip_sgml_tag = skip_sgml_tag;
      pair_impl->prefix = (name.empty() ? std::string("frontier-pair") : name);

      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("frontier-pair") : name);
      base_type::__sparse_feature = true;
      
      pimpl = pair_impl.release();
    }
    
    FrontierPair::~FrontierPair() { std::unique_ptr<impl_type> tmp(pimpl); }
    
    FrontierPair::FrontierPair(const FrontierPair& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {
      pimpl->clear_cache();
    }

    FrontierPair& FrontierPair::operator=(const FrontierPair& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      pimpl->clear_cache();
            
      return *this;
    }
    
    void FrontierPair::apply(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     const bool final) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      feature_set_type feats;
      
      pimpl->pair_score(edge, feats);

      features.update(feats, static_cast<const std::string&>(base_type::feature_name()));
    }

    void FrontierPair::apply_coarse(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void FrontierPair::apply_predict(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void FrontierPair::apply_scan(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  const int dot,
				  feature_set_type& features,
				  const bool final) const
    {}
    void FrontierPair::apply_complete(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      const bool final) const
    {}


    void FrontierPair::initialize()
    {
      pimpl->clear();
    }
    
  };
};
