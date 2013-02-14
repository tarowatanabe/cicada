//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "frontier_bigram.hpp"
#include "feature_builder.hpp"

#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"
#include "cicada/feature_vector_linear.hpp"

#include "utils/indexed_set.hpp"
#include "utils/chunk_vector.hpp"
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


    class FrontierBigramImpl : public utils::hashmurmur3<size_t>
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

      typedef FeatureVectorLinear<feature_set_type::mapped_type> feature_linear_set_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef symbol_type word_type;
      
      template <typename Tp>
      struct hash_sequence : public utils::hashmurmur3<size_t>
      {
	typedef utils::hashmurmur3<size_t> hasher_type;
	
	size_t operator()(const Tp& x) const
	{
	  return hasher_type::operator()(x.begin(), x.end(), 0);
	}
      };
      
      typedef utils::indexed_set<std::string, hash_sequence<std::string>, std::equal_to<std::string>,
				 std::allocator<std::string> > frontier_map_type;

      typedef utils::chunk_vector<std::string, 4096 / sizeof(std::string), std::allocator<std::string> > frontier_string_type;
      
      typedef utils::simple_vector<frontier_map_type::index_type, std::allocator<frontier_map_type::index_type> > node_set_type;
      
      typedef utils::indexed_set<node_set_type, hash_sequence<node_set_type>, std::equal_to<node_set_type>,
				 std::allocator<node_set_type> > node_map_type;
      
      typedef node_map_type::index_type id_type;

      typedef std::pair<id_type, id_type> node_pair_type;

      struct node_pair_hash : public utils::hashmurmur3<size_t>
      {
	typedef utils::hashmurmur3<size_t> hasher_type;
	
	size_t operator()(const node_pair_type& x) const
	{
	  return hasher_type::operator()(x.first, x.second);
	}
      };

      typedef utils::unordered_map<node_pair_type, feature_type, node_pair_hash, std::equal_to<node_pair_type>,
				   std::allocator<std::pair<const node_pair_type, feature_type> > >::type cache_feature_set_type;
      
      
      typedef FeatureBuilder feature_builder_type;

      FrontierBigramImpl()
	: skip_sgml_tag(false), prefix("frontier-bigram"), forced_feature(false),
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
      
      void bigram_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features)
      {
	if (skip_sgml_tag)
	  bigram_score(state, states, edge, features, skipper_sgml());
	else
	  bigram_score(state, states, edge, features, skipper_epsilon());
      }
      
      struct __attribute_string : public boost::static_visitor<cicada::AttributeVector::string_type>
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
      void bigram_score(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			Skipper skipper)
      {
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_frontier_source);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_frontier_target);
	
	if (siter == edge.attributes.end() || titer == edge.attributes.end()) return;
	
	const std::string& frontier_source = boost::apply_visitor(__attribute_string(), siter->second);
	const std::string& frontier_target = boost::apply_visitor(__attribute_string(), titer->second);

	static const frontier_map_type::index_type root = frontier_map_type::index_type(-1);
	
	frontier_map_type::index_type node(root);
	
	if (! frontier_source.empty() && ! frontier_target.empty()) {
	  frontier_map_type::iterator iter = const_cast<frontier_map_type&>(frontier_map).insert(frontier_source).first;
	  node = iter - frontier_map.begin();
	  
	  if (frontier_string.size() <= node)
	    frontier_string.resize(node + 1);
	  
	  if (frontier_string[node].empty()) {
	    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	    
	    utils::piece frontier_piece(frontier_source);
	    tokenizer_type tokenizer(frontier_piece);

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
	    
	    frontier_string[node] = feature_builder;
	  }
	}
	
	if (node != root) {
	  if (! states.empty()) {
	    state_ptr_set_type::const_iterator siter_end = states.end();
	    for (state_ptr_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
	      const node_set_type& antecedents = node_map[*reinterpret_cast<const id_type*>(*siter)];
	      
	      // we use the pair of node and antecedents as a feature!
	      
	      node_set_type::const_iterator aiter_end = antecedents.end();
	      for (node_set_type::const_iterator aiter = antecedents.begin(); aiter != aiter_end; ++ aiter) {

		cache_feature_set_type::iterator fiter = cache_features.find(node_pair_type(node, *aiter));
		if (fiter == cache_features.end()) {
		  fiter = cache_features.insert(std::make_pair(node_pair_type(node, *aiter), feature_type())).first;
		  
		  feature_builder.clear();
		  feature_builder << prefix << ":" << frontier_string[node] << "|" << frontier_string[*aiter];
		  
		  if (forced_feature || feature_builder.exists())
		    fiter->second = feature_builder;
		}
		
		if (! fiter->second.empty())
		  features[fiter->second] += 1.0;
	      }
	    }
	  }
	  
	  // new state...
	  node_map_type::iterator niter = const_cast<node_map_type&>(node_map).insert(node_set_type(1, node)).first;
	  *reinterpret_cast<id_type*>(state) = niter - node_map.begin();
	} else {
	  if (states.empty()) {
	    // no states from antecedents
	    
	    node_map_type::iterator niter = const_cast<node_map_type&>(node_map).insert(node_set_type()).first;
	    *reinterpret_cast<id_type*>(state) = niter - node_map.begin();
	  } else {
	    // simply copy from states into state...
	    
	    node_set_type nodes;
	    
	    state_ptr_set_type::const_iterator siter_end = states.end();
	    for (state_ptr_set_type::const_iterator siter = states.begin(); siter != siter_end; ++ siter) {
	      const node_set_type& antecedents = node_map[*reinterpret_cast<const id_type*>(*siter)];
	      
	      if (! antecedents.empty())
		nodes.insert(nodes.end(), antecedents.begin(), antecedents.end());
	    }
	    
	    std::sort(nodes.begin(), nodes.end());
	    
	    node_map_type::iterator niter = const_cast<node_map_type&>(node_map).insert(nodes).first;
	    *reinterpret_cast<id_type*>(state) = niter - node_map.begin();
	  }
	}
      }
	
            
      void clear()
      {
	cache_features.clear();
	
	frontier_map.clear();
	frontier_string.clear();
	node_map.clear();
      }
      
      void clear_cache()
      {
	clear();
      }

      frontier_map_type    frontier_map;
      frontier_string_type frontier_string;
      node_map_type        node_map;

      cache_feature_set_type cache_features;

      feature_builder_type feature_builder;
      
      bool skip_sgml_tag;
      
      std::string prefix;
      bool forced_feature;

      attribute_type attr_frontier_source;
      attribute_type attr_frontier_target;
    };
    
    FrontierBigram::FrontierBigram(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "frontier-bigram")
	throw std::runtime_error("this is not frontier bigram feature: " + parameter);

      bool skip_sgml_tag = false;
      
      std::string name;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for frontier bigram: " << piter->first << "=" << piter->second << std::endl;
      }
            
      std::auto_ptr<impl_type> bigram_impl(new impl_type());

      bigram_impl->skip_sgml_tag = skip_sgml_tag;
      bigram_impl->prefix = (name.empty() ? std::string("frontier-bigram") : name);

      base_type::__state_size = sizeof(impl_type::id_type);
      base_type::__feature_name = (name.empty() ? std::string("frontier-bigram") : name);
      base_type::__sparse_feature = true;
      
      pimpl = bigram_impl.release();
    }
    
    FrontierBigram::~FrontierBigram() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    FrontierBigram::FrontierBigram(const FrontierBigram& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {
      pimpl->clear_cache();
    }

    FrontierBigram& FrontierBigram::operator=(const FrontierBigram& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      pimpl->clear_cache();
      
      return *this;
    }
    
    void FrontierBigram::apply(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       const bool final) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      feature_set_type feats;
      
      pimpl->bigram_score(state, states, edge, feats);

      features.update(feats, static_cast<const std::string&>(base_type::feature_name()));
    }

    void FrontierBigram::apply_coarse(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void FrontierBigram::apply_predict(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void FrontierBigram::apply_scan(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    const int dot,
				    feature_set_type& features,
				    const bool final) const
    {}
    void FrontierBigram::apply_complete(state_ptr_type& state,
					const state_ptr_set_type& states,
					const edge_type& edge,
					feature_set_type& features,
					const bool final) const
    {}


    void FrontierBigram::initialize()
    {
      pimpl->clear();
    }
    
  };
};
