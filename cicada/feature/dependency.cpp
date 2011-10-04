
#include <vector>

#include "feature/dependency.hpp"

#include "parameter.hpp"

#include "utils/piece.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/simple_vector.hpp"
#include "utils/indexed_trie.hpp"
#include "utils/chunk_vector.hpp"
#include "utils/sgi_hash_map.hpp"

#include <boost/fusion/tuple.hpp>
#include <boost/array.hpp>

namespace cicada
{
  namespace feature
  {
    class DependencyImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;

      typedef Dependency::feature_function_type feature_function_type;

      typedef feature_function_type::symbol_type symbol_type;
      typedef feature_function_type::vocab_type  vocab_type;

      typedef feature_function_type::hypergraph_type hypergraph_type;
      typedef feature_function_type::lattice_type    lattice_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef std::pair<int, int> lattice_edge_type;
      typedef std::vector<lattice_edge_type, std::allocator<lattice_edge_type> > lattice_edge_set_type;

      typedef std::vector<int, std::allocator<int> > lattice_node_set_type;
      typedef std::vector<lattice_node_set_type, std::allocator<lattice_node_set_type> > lattice_node_map_type;
      
      typedef std::pair<symbol_type, symbol_type> terminal_pos_type;
      typedef std::vector<terminal_pos_type, std::allocator<terminal_pos_type> > terminal_pos_set_type;
      
      typedef std::pair<int, int> dependency_type;
      typedef std::vector<dependency_type, std::allocator<dependency_type> >         dependency_set_type;
      typedef std::vector<dependency_set_type, std::allocator<dependency_set_type> > dependency_map_type;
      
      
      typedef utils::indexed_trie<dependency_type, utils::hashmurmur<size_t>, std::equal_to<dependency_type>, std::allocator<dependency_type> > dependency_index_type;

      typedef utils::simple_vector<feature_type, std::allocator<feature_type> > feature_list_type;
      //typedef std::vector<feature_type, std::allocator<feature_type> > feature_list_type;
      typedef std::vector<feature_list_type, std::allocator<feature_list_type> > feature_map_type;
      typedef utils::chunk_vector<feature_list_type, 4096 / sizeof(feature_list_type), std::allocator<feature_list_type> > feature_pair_map_type;

      
      typedef std::pair<dependency_type, dependency_type> dependency_pair_type;
#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<dependency_pair_type, feature_list_type, utils::hashmurmur<size_t>, std::equal_to<dependency_pair_type>,
				      std::allocator<std::pair<const dependency_pair_type, feature_list_type> > > feature_order_map_type;
#else
      typedef sgi::hash_map<dependency_pair_type, feature_list_type, utils::hashmurmur<size_t>, std::equal_to<dependency_pair_type>,
			    std::allocator<std::pair<const dependency_pair_type, feature_list_type> > > feature_order_map_type;
#endif
      
      // temporary...
      typedef std::vector<std::string, std::allocator<std::string> > feats_type;

      struct __attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
      {
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
      };
      
      DependencyImpl(const int __order)
	: order(__order),
	  lattice(0),
	  forced_feature(false),
	  feat_none("dependency"),
	  feat_root_multiple("dependency:root-multiple"),
	  attr_dependency_pos("dependency-pos"),
	  attr_dependency_head("dependency-head"),
	  attr_dependency_dependent("dependency-dependent") {}
      
      void dependency_score(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features)
      {
	typedef dependency_index_type::id_type id_type;
	
	int pos_head = -1;
	int pos_dep  = -1;
	
	attribute_set_type::const_iterator hiter = edge.attributes.find(attr_dependency_head);
	attribute_set_type::const_iterator diter = edge.attributes.find(attr_dependency_dependent);
	
	if (hiter != edge.attributes.end() && diter != edge.attributes.end()) {
	  pos_head = boost::apply_visitor(__attribute_integer(), hiter->second);
	  pos_dep  = boost::apply_visitor(__attribute_integer(), diter->second);
	}

	std::cerr << "hypergraph head: " << edge.head << std::endl;
	
	const id_type state_dep = (pos_head >= 0 && pos_dep >= 0 
				   ? dependency_index.push(dependency_index.root(), dependency_type(pos_head, pos_dep))
				   : dependency_index.root());
	
	if (pos_head >= 0 && pos_dep >= 0) {
	  std::cerr << "dep: " << pos_head << "->" << pos_dep << std::endl;

	  apply_features(state_dep, pos_head, pos_dep, features);
	}
		
	// root-count...
	int root_count = (pos_head == 0 && pos_dep >= 0);
	bool fired_root_multiple = false;
		
	dependency_antecedents.clear();
	dependency_antecedents.resize(order - 1);
	
	std::cerr << "states sizes: " << states.size() << std::endl;
	
	for (size_t i = 0; i != states.size(); ++ i) {
	  // root count...
	  const int& root_count_antecedent = *reinterpret_cast<const int*>(states[i]);
	  
	  root_count += root_count_antecedent;
	  fired_root_multiple |= root_count_antecedent > 1;
	  
	  // convert state representation into ordered dependencies...
	  const id_type* id_antecedent = reinterpret_cast<const id_type*>(reinterpret_cast<const int*>(states[i]) + 1);
	  const id_type* id_antecedent_last = id_antecedent + order;
	  
	  // TODO: how to represent state space...?
	  for (int k = 0; k != order - 1 && *id_antecedent != dependency_index.root(); ++ k, ++ id_antecedent) {
	    id_type id = *id_antecedent;
	    while (id != dependency_index.root()) {
	      dependency_antecedents[k].push_back(dependency_index[id]);
	      
	      std::cerr << "order: " << k << " dep: " << dependency_index[id].first << "->" << dependency_index[id].second << std::endl;

	      id = dependency_index.parent(id);
	    }
	  }

	  if (id_antecedent > id_antecedent_last)
	    throw std::runtime_error("invalid access for antecedent state");
	}
	
	// we will fire for higher order features...
	// in this implementation, we will simply ignore the boundary of antecedents...
	if (pos_head >= 0 && pos_dep >= 0)
	  apply_features(state_dep, pos_head, pos_dep, dependency_antecedents, features);
	  
	// fire root multiple
	if (! fired_root_multiple && root_count > 1) 
	  features[feat_root_multiple] = 1.0;
	
	// update states...
	
	// root-count
	*reinterpret_cast<int*>(state) = root_count;
	
	id_type* state_id = reinterpret_cast<id_type*>(reinterpret_cast<int*>(state) + 1);
	id_type* state_id_last = state_id + order;
	std::fill(state_id, state_id + order, dependency_index.root());
	
	// antecedents..
	int order_adjusted = order - 1;
	if (pos_head >= 0 && pos_dep >= 0) {
	  *state_id = state_dep;
	  ++ state_id;
	  order_adjusted = order - 2;
	}
	
	for (int k = 0; k < order_adjusted && ! dependency_antecedents[k].empty(); ++ k, ++ state_id) {
	  std::sort(dependency_antecedents[k].begin(), dependency_antecedents[k].end());
	  
	  id_type id = dependency_index.root();
	  dependency_set_type::const_iterator diter_end = dependency_antecedents[k].end();
	  for (dependency_set_type::const_iterator diter = dependency_antecedents[k].begin(); diter != diter_end; ++ diter)
	    id = dependency_index.push(id, *diter);
	  
	  *state_id = id;
	}
	
	if (state_id > state_id_last)
	  throw std::runtime_error("exceed max?");
	
	std::cerr << "finished" << std::endl;
      }	

      template <typename Iterator>
      void apply_features(Iterator first, Iterator last, feature_set_type& features)
      {
	for (/**/; first != last; ++ first)
	  if (*first != feat_none)
	    features[*first] += 1.0;
      }
      
      void apply_features(const dependency_index_type::id_type& state,
			  const int& pos_head,
			  const int& pos_tail,
			  feature_set_type& features)
      {
	// we will do caching for base features....
	
	
	// head...
	if (features_heads[pos_head].empty()) {
	  const std::string& word = terminals[pos_head].first;
	  const std::string& pos = terminals[pos_head].second;
	  
	  const std::string feat_word_pos = "dependency:head-word-pos:" + word + '|' + pos;
	  const std::string feat_word     = "dependency:head-word:" + word;
	  const std::string feat_pos      = "dependency:head-pos:" + pos;
	  
	  if (forced_feature) {
	    features_heads[pos_head].push_back(feat_word_pos);
	    features_heads[pos_head].push_back(feat_word);
	    features_heads[pos_head].push_back(feat_pos);
	  } else {
	    if (feature_type::exists(feat_word_pos))
	      features_heads[pos_head].push_back(feat_word_pos);
	    if (feature_type::exists(feat_word))
	      features_heads[pos_head].push_back(feat_word);
	    if (feature_type::exists(feat_pos))
	      features_heads[pos_head].push_back(feat_pos);
	    
	    // fallback to NONE
	    if (features_heads[pos_head].empty())
	      features_heads[pos_head].push_back(feat_none);
	  }
	}
	
	// apply features...
	apply_features(features_heads[pos_head].begin(), features_heads[pos_head].end(), features);
	
	// dependent...
	if (features_tails[pos_tail].empty()) {
	  const std::string& word = terminals[pos_tail].first;
	  const std::string& pos = terminals[pos_tail].second;
	  
	  const std::string feat_word_pos = "dependency:dep-word-pos:" + word + '|' + pos;
	  const std::string feat_word     = "dependency:dep-word:" + word;
	  const std::string feat_pos      = "dependency:dep-pos:" + pos;
	  
	  if (forced_feature) {
	    features_tails[pos_tail].push_back(feat_word_pos);
	    features_tails[pos_tail].push_back(feat_word);
	    features_tails[pos_tail].push_back(feat_pos);
	  } else {
	    if (feature_type::exists(feat_word_pos))
	      features_tails[pos_tail].push_back(feat_word_pos);
	    if (feature_type::exists(feat_word))
	      features_tails[pos_tail].push_back(feat_word);
	    if (feature_type::exists(feat_pos))
	      features_tails[pos_tail].push_back(feat_pos);
	    
	    // fallback to NONE
	    if (features_tails[pos_tail].empty())
	      features_tails[pos_tail].push_back(feat_none);
	  }
	}
	
	// apply features...
	apply_features(features_tails[pos_tail].begin(), features_tails[pos_tail].end(), features);
	
	// pairs...
	if (state >= features_pairs.size())
	  features_pairs.resize(state + 1);
	
	if (features_pairs[state].empty()) {
	  feats.clear();
	  
	  // bigram features...
	  {
	    const std::string& head_word = terminals[pos_head].first;
	    const std::string& head_pos  = terminals[pos_head].second;
	    const std::string& tail_word = terminals[pos_tail].first;
	    const std::string& tail_pos  = terminals[pos_tail].second;
	    static const std::string empty;
	    
	    feats.push_back("dependency:head-dep:" + head_word + '|' + head_pos + '+' + tail_word + '|' + tail_pos);
	    feats.push_back("dependency:head-dep:" + empty     + '|' + head_pos + '+' + tail_word + '|' + tail_pos);
	    feats.push_back("dependency:head-dep:" + head_word + '|' + empty    + '+' + tail_word + '|' + tail_pos);
	    feats.push_back("dependency:head-dep:" + head_word + '|' + head_pos + '+' + empty     + '|' + tail_pos);
	    feats.push_back("dependency:head-dep:" + head_word + '|' + head_pos + '+' + tail_word + '|' + empty);
	    feats.push_back("dependency:head-dep:" + head_word + '|' + empty    + '+' + tail_word + '|' + empty);
	    feats.push_back("dependency:head-dep:" + empty     + '|' + head_pos + '+' + empty     + '|' + tail_pos);
	    
	    // direction
	    feats.push_back("dependency:head-dep-dir:" + head_pos + '+' + tail_pos + (pos_tail > pos_head ? ":R" : ":L"));
	  }
	  
	  // surrounding POS context features
	  {
	    static const lattice_node_set_type nodes_empty;
	    
	    const lattice_node_set_type& nodes_head_prev = (pos_head == 0 ? nodes_empty : nodes_backward[edges[pos_head].first]);
	    const lattice_node_set_type& nodes_head_next = nodes_forward[edges[pos_head].second];
	    
	    const lattice_node_set_type& nodes_tail_prev = nodes_backward[edges[pos_tail].first];
	    const lattice_node_set_type& nodes_tail_next = nodes_forward[edges[pos_tail].second];
	      
	    const std::string& head_pos = terminals[pos_head].second;
	    const std::string& tail_pos = terminals[pos_tail].second;
	      
	    if (! nodes_head_next.empty() && ! nodes_tail_prev.empty()) {
	      lattice_node_set_type::const_iterator head_niter_end = nodes_head_next.end();
	      for (lattice_node_set_type::const_iterator head_niter = nodes_head_next.begin(); head_niter != head_niter_end; ++ head_niter) {
		lattice_node_set_type::const_iterator tail_piter_end = nodes_tail_prev.end();
		for (lattice_node_set_type::const_iterator tail_piter = nodes_tail_prev.begin(); tail_piter != tail_piter_end; ++ tail_piter) {
		  const std::string& head_pos_next = terminals[*head_niter].second;
		  const std::string& tail_pos_prev = terminals[*tail_piter].second;
		    
		  feats.push_back("dependency:+1-1:" + head_pos + '|' + head_pos_next + '+' + tail_pos_prev + '|' + tail_pos);
		}
	      }
	    }
	      
	    if (! nodes_head_prev.empty() && ! nodes_tail_prev.empty()) {
	      lattice_node_set_type::const_iterator head_piter_end = nodes_head_prev.end();
	      for (lattice_node_set_type::const_iterator head_piter = nodes_head_prev.begin(); head_piter != head_piter_end; ++ head_piter) {
		lattice_node_set_type::const_iterator tail_piter_end = nodes_tail_prev.end();
		for (lattice_node_set_type::const_iterator tail_piter = nodes_tail_prev.begin(); tail_piter != tail_piter_end; ++ tail_piter) {
		  const std::string& head_pos_prev = terminals[*head_piter].second;
		  const std::string& tail_pos_prev = terminals[*tail_piter].second;
		    
		  feats.push_back("dependency:-1-1:" + head_pos_prev + '|' + head_pos + '+' + tail_pos_prev + '|' + tail_pos);
		}
	      }
	    }
	      
	    if (! nodes_head_next.empty() && ! nodes_tail_next.empty()) {
	      lattice_node_set_type::const_iterator head_niter_end = nodes_head_next.end();
	      for (lattice_node_set_type::const_iterator head_niter = nodes_head_next.begin(); head_niter != head_niter_end; ++ head_niter) {
		lattice_node_set_type::const_iterator tail_niter_end = nodes_tail_next.end();
		for (lattice_node_set_type::const_iterator tail_niter = nodes_tail_next.begin(); tail_niter != tail_niter_end; ++ tail_niter) {
		  const std::string& head_pos_next = terminals[*head_niter].second;
		  const std::string& tail_pos_next = terminals[*tail_niter].second;
		    
		  feats.push_back("dependency:+1+1:" + head_pos + '|' + head_pos_next + '+' + tail_pos + '|' + tail_pos_next);
		}
	      }
	    }
	      
	    if (! nodes_head_prev.empty() && ! nodes_tail_next.empty()) {
	      lattice_node_set_type::const_iterator head_piter_end = nodes_head_prev.end();
	      for (lattice_node_set_type::const_iterator head_piter = nodes_head_prev.begin(); head_piter != head_piter_end; ++ head_piter) {
		lattice_node_set_type::const_iterator tail_niter_end = nodes_tail_next.end();
		for (lattice_node_set_type::const_iterator tail_niter = nodes_tail_next.begin(); tail_niter != tail_niter_end; ++ tail_niter) {
		  const std::string& head_pos_prev = terminals[*head_piter].second;
		  const std::string& tail_pos_next = terminals[*tail_niter].second;
		  
		  feats.push_back("dependency:-1+1:" + head_pos_prev + '|' + head_pos + '+' + tail_pos + '|' + tail_pos_next);
		}
	      }
	    }
	  }
	  
	  if (forced_feature)
	    features_pairs[state].insert(features_pairs[state].end(), feats.begin(), feats.end());
	  else {
	    feats_type::const_iterator fiter_end = feats.end();
	    for (feats_type::const_iterator fiter = feats.begin(); fiter != fiter_end; ++ fiter)
	      if (feature_type::exists(*fiter))
		features_pairs[state].push_back(*fiter);
	    
	    // fallback to NONE
	    if (features_pairs[state].empty())
	      features_pairs[state].push_back(feat_none);
	  }
	}
	
	// apply features...
	apply_features(features_pairs[state].begin(), features_pairs[state].end(), features);
      }
      
      void apply_features(const dependency_index_type::id_type& state,
			  const int& pos_head,
			  const int& pos_tail,
			  const dependency_map_type& antecedents,
			  feature_set_type& features)
      {
	static const std::string empty;
	const dependency_type parent(pos_head, pos_tail);
	
	// we will do caching...
	for (size_t k = 0; k != antecedents.size(); ++ k)
	  if (! antecedents[k].empty()) {
	    dependency_set_type::const_iterator aiter_end = antecedents[k].end();
	    for (dependency_set_type::const_iterator aiter = antecedents[k].begin(); aiter != aiter_end; ++ aiter) {
	      const dependency_type& antecedent = *aiter;
	      
	      std::pair<feature_order_map_type::iterator, bool> result = features_order.insert(std::make_pair(dependency_pair_type(parent, antecedent),
													      feature_list_type()));
	      
	      if (result.second) {
		feats.clear();
		
		if (pos_tail == antecedent.first) {
		  const int pos1 = pos_head;
		  const int pos2 = pos_tail;
		  const int pos3 = antecedent.second;

		  const std::string& head_word = terminals[pos_head].first;
		  const std::string& head_pos  = terminals[pos_head].second;
		  const std::string& mid_word  = terminals[pos_tail].first;
		  const std::string& mid_pos   = terminals[pos_tail].second;
		  const std::string& tail_word = terminals[antecedent.second].first;
		  const std::string& tail_pos  = terminals[antecedent.second].second;
		  
		  feats.push_back("dependency:parent:" + head_word + '|' + head_pos
				  + '&' + empty     + '|' + mid_pos
				  + '&' + empty     + '|' + tail_pos);
		  feats.push_back("dependency:parent:" + empty     + '|' + head_pos
				  + '&' + mid_word  + '|' + mid_pos
				  + '&' + empty     + '|' + tail_pos);
		  feats.push_back("dependency:parent:" + empty     + '|' + head_pos
				  + '&' + empty     + '|' + mid_pos
				  + '&' + tail_word + '|' + tail_pos);
		  
		  feats.push_back("dependency:parent:" + head_pos + '&' + mid_pos + '&' + tail_pos
				  + std::string(pos2 > pos1 ? ":R" : ":L")
				  + std::string(pos3 > pos2 ? "&R" : "&L")
				  + std::string(pos3 > pos1 ? "&R" : "&L"));
		} else if (antecedent.second == pos_head) {
		  // actually, this will not hapen...?
		  const int pos1 = antecedent.first;
		  const int pos2 = antecedent.second;
		  const int pos3 = pos_tail;
		  
		  const std::string& head_word = terminals[antecedent.first].first;
		  const std::string& head_pos  = terminals[antecedent.first].second;
		  const std::string& mid_word  = terminals[antecedent.second].first;
		  const std::string& mid_pos   = terminals[antecedent.second].second;
		  const std::string& tail_word = terminals[pos_tail].first;
		  const std::string& tail_pos  = terminals[pos_tail].second;
		  
		  feats.push_back("dependency:child:" + head_word + '|' + head_pos
				  + '&' + empty     + '|' + mid_pos
				  + '&' + empty     + '|' + tail_pos);
		  feats.push_back("dependency:child:" + empty     + '|' + head_pos
				  + '&' + mid_word  + '|' + mid_pos
				  + '&' + empty     + '|' + tail_pos);
		  feats.push_back("dependency:child:" + empty     + '|' + head_pos
				  + '&' + empty     + '|' + mid_pos
				  + '&' + tail_word + '|' + tail_pos);
		  
		  feats.push_back("dependency:child:" + head_pos + '&' + mid_pos + '&' + tail_pos
				  + std::string(pos2 > pos1 ? ":R" : ":L")
				  + std::string(pos3 > pos2 ? "&R" : "&L")
				  + std::string(pos3 > pos1 ? "&R" : "&L"));
		} else if (pos_head == antecedent.first) {
		  const std::string& head_word = terminals[pos_head].first;
		  const std::string& head_pos  = terminals[pos_head].second;
		  const std::string& tail1_word = terminals[pos_tail].first;
		  const std::string& tail1_pos  = terminals[pos_tail].second;
		  const std::string& tail2_word = terminals[antecedent.second].first;
		  const std::string& tail2_pos  = terminals[antecedent.second].second;
		  
		  feats.push_back("dependency:sibling:" + head_word + '|' + head_pos
				  + '+' + empty      + '|' + tail1_pos
				  + '&' + empty      + '|' + tail2_pos);
		  feats.push_back("dependency:sibling:" + empty     + '|' + head_pos
				  + '+' + tail1_word + '|' + tail1_pos
				  + '&' + empty      + '|' + tail2_pos);
		  feats.push_back("dependency:sibling:" + empty     + '|' + head_pos
				  + '+' + empty      + '|' + tail1_pos
				  + '&' + tail2_word + '|' + tail2_pos);
		  feats.push_back("dependency:sibling:" + head_pos + '+' + tail1_pos + '&' + tail2_pos
				  + std::string(pos_tail > pos_head ? ":R" : ":L")
				  + std::string(antecedent.second > pos_head ? "+R" : "+L")
				  + std::string(pos_tail > antecedent.second ? "&R" : "&L"));
		} else {
		  // floating... do we construct a feature set...?
#if 0
		  const std::string& head1_word = terminals[pos_head].first;
		  const std::string& head1_pos  = terminals[pos_head].second;
		  const std::string& tail1_word = terminals[pos_tail].first;
		  const std::string& tail1_pos  = terminals[pos_tail].second;
		  const std::string& tail1_word = terminals[antecedent.first].first;
		  const std::string& tail1_pos  = terminals[antecedent.first].second;
		  const std::string& tail2_word = terminals[antecedent.second].first;
		  const std::string& tail2_pos  = terminals[antecedent.second].second;
#endif
		  
		  
		}
		
		if (forced_feature)
		  result.first->second.insert(result.first->second.end(), feats.begin(), feats.end());
		else {
		  feats_type::const_iterator fiter_end = feats.end();
		  for (feats_type::const_iterator fiter = feats.begin(); fiter != fiter_end; ++ fiter)
		    if (feature_type::exists(*fiter))
		      result.first->second.push_back(*fiter);
		}
	      }
	      
	      apply_features(result.first->second.begin(), result.first->second.end(), features);
	    }
	  }
      }
      
      void clear()
      {
      }
      
      void assign(const lattice_type& __lattice)
      {
	lattice = &__lattice;
	
	edges.clear();
	terminals.clear();
	
	nodes_forward.clear();
	nodes_backward.clear();
	nodes_forward.resize(lattice->size() + 1);
	nodes_backward.resize(lattice->size() + 1);
	
	//
	// we need to compute adjacent nodes and edges...
	//
	
	// ROOT
	nodes_backward.front().push_back(0);
	edges.push_back(std::make_pair(-1, 0));
	terminals.push_back(std::make_pair(vocab_type::BOS, vocab_type::BOS));
	
	// terminals/POSs
	for (size_type pos = 0; pos != lattice->size(); ++ pos) {
	  lattice_type::arc_set_type::const_iterator aiter_end = lattice->operator[](pos).end();
	  for (lattice_type::arc_set_type::const_iterator aiter = lattice->operator[](pos).begin(); aiter != aiter_end; ++ aiter) {
	    // keep forward and backward edges...
	    nodes_forward[pos].push_back(edges.size());
	    nodes_backward[pos + aiter->distance].push_back(edges.size());
	    
	    edges.push_back(std::make_pair(pos, pos + aiter->distance));
	    terminals.push_back(std::make_pair(aiter->label.terminal(), aiter->label.pos()));
	    
	    if (terminals.back().second.empty())
	      terminals.back().second = vocab_type::X;
	  }
	}
	
	// END
	nodes_forward[lattice->size()].push_back(edges.size());
	edges.push_back(std::make_pair(lattice->size(), lattice->size() + 1));
	terminals.push_back(std::make_pair(vocab_type::EOS, vocab_type::EOS));
	
	// clear dependency index
	dependency_index.clear();
	dependency_antecedents.clear();
	
	// caching...
	features_heads.clear();
	features_tails.clear();
	features_pairs.clear();
	features_order.clear();
	
	features_heads.resize(edges.size());
	features_tails.resize(edges.size());
      }
      
      int order;
      
      const lattice_type*   lattice;
      lattice_edge_set_type edges;
      terminal_pos_set_type terminals;
      
      lattice_node_map_type nodes_forward;
      lattice_node_map_type nodes_backward;
      
      bool forced_feature;

      feature_type feat_none;
      feature_type feat_root_multiple;
      
      attribute_type attr_dependency_pos;
      attribute_type attr_dependency_head;
      attribute_type attr_dependency_dependent;
      
      // internal use only...
      dependency_index_type dependency_index;
      dependency_map_type   dependency_antecedents;
      
      feature_map_type       features_heads;
      feature_map_type       features_tails;
      feature_pair_map_type  features_pairs;
      feature_order_map_type features_order;

      feats_type feats;
    };

    Dependency::Dependency(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "dependency")
	throw std::runtime_error("is this really dependency feature function? " + parameter);
      
      int order = 2;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for dependency: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (order <= 0)
	throw std::runtime_error("we do not support zero or negative orders");
      
      pimpl = new impl_type(order);
      
      base_type::__state_size = sizeof(int) + sizeof(impl_type::dependency_index_type::id_type) * order;
      base_type::__feature_name = "dependency";
      base_type::__sparse_feature = true;
    }
    
    Dependency::Dependency(const Dependency& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl)) {}
    
    Dependency::~Dependency() { if (pimpl) delete pimpl; }
    
    Dependency& Dependency::operator=(const Dependency& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    
    void Dependency::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->dependency_score(state, states, edge, features);
    }
    
    void Dependency::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {

    }
    
    void Dependency::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {

    }
    
    void Dependency::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {

    }
    
    void Dependency::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    void Dependency::initialize()
    {
      pimpl->clear();
    }
    
    void Dependency::assign(const size_type& id,
			    const hypergraph_type& hypergraph,
			    const lattice_type& lattice,
			    const span_set_type& spans,
			    const sentence_set_type& targets,
			    const ngram_count_set_type& ngram_counts)
    {
      pimpl->assign(lattice);
    }
    
  };
};
