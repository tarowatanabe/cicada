//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/semiring.hpp>
#include <cicada/kbest.hpp>
#include <cicada/graphviz.hpp>
#include <cicada/treebank.hpp>
#include <cicada/inside_outside.hpp>
#include <cicada/span_node.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/debinarize.hpp>

#include <cicada/operation/output.hpp>
#include <cicada/operation/functional.hpp>
#include <cicada/operation/traversal.hpp>

#include <utils/lexical_cast.hpp>
#include <utils/resource.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/piece.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  namespace operation
  {
    
    template <typename Hypergraph, typename Function, typename Filter>
    inline
    void kbest_derivations(std::ostream& os,
			   const size_t id,
			   const Hypergraph& graph,
			   const int kbest_size,
			   const Function& function,
			   const Filter& filter,
			   const bool no_id,
			   const bool graphviz_mode,
			   const bool treebank_mode,
			   const bool debinarize)
    {
      typedef Hypergraph hypergraph_type;
      typedef typename hypergraph_type::rule_type rule_type;
      
      if (! graph.is_valid()) {
	hypergraph_type graph_empty;
	if (! no_id)
	  os << id << " ||| ";
	if (graphviz_mode)
	  cicada::graphviz(os, graph_empty) << '\n';
	else if (treebank_mode)
	  cicada::treebank(os, graph_empty) << " ||| ||| 0" << '\n';
	else
	  os << graph_empty << " ||| ||| 0" << '\n';
	return;
      }
      
      cicada::KBest<edge_feature_traversal, Function, Filter> derivations(graph, kbest_size, edge_feature_traversal(), function, filter);
      
      typedef edge_feature_traversal::value_type    derivation_type;
      typedef edge_feature_traversal::edge_set_type edge_set_type;
  
      typedef typename hypergraph_type::id_type id_type;

      typedef typename Function::value_type weight_type;

      typedef google::dense_hash_map<id_type, id_type, utils::hashmurmur<size_t>, std::equal_to<id_type> > node_map_type;

      typedef std::vector<id_type, std::allocator<id_type> > head_set_type;
      
      derivation_type derivation;
      weight_type     weight;
      node_map_type   node_maps;
      head_set_type   heads;
      hypergraph_type graph_kbest;
      
      node_maps.set_empty_key(id_type(-1));

      edge_set_type tails;
  
      for (int k = 0; k < kbest_size; ++ k) {
	if (! derivations(k, derivation, weight))
	  break;
    
	const edge_set_type& edges = boost::get<0>(derivation);
	
	heads.clear();
	node_maps.clear();
	graph_kbest.clear();

	heads.reserve(edges.size());
    
	id_type node_id = 0;
	edge_set_type::const_iterator eiter_end = edges.end();
	for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
	  std::pair<typename node_map_type::iterator, bool> result = node_maps.insert(std::make_pair(graph.edges[*eiter].head, node_id));
	  
	  heads.push_back(result.first->second);
	  node_id += result.second;
	}
    
	for (id_type node = 0; node != node_id; ++ node)
	  graph_kbest.add_node();
	
	id_type edge_id = 0;
	for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter, ++ edge_id) {
	  const typename hypergraph_type::edge_type& edge = graph.edges[*eiter];
      
	  tails.clear();
	  typename hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	  for (typename hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
	    typename node_map_type::const_iterator niter = node_maps.find(*titer);
	    if (niter == node_maps.end())
	      throw std::runtime_error("no node?");
	
	    tails.push_back(niter->second);
	  }
      
	  typename hypergraph_type::edge_type& edge_kbest = graph_kbest.add_edge(tails.begin(), tails.end());
	  edge_kbest.rule = edge.rule;
	  edge_kbest.features = edge.features;
	  edge_kbest.attributes = edge.attributes;
      
	  graph_kbest.connect_edge(edge_kbest.id, heads[edge_id]);
	}
    
	typename node_map_type::const_iterator niter = node_maps.find(graph.goal);
	if (niter == node_maps.end())
	  throw std::runtime_error("did not reach goal?");
    
	graph_kbest.goal = niter->second;

	graph_kbest.topologically_sort();


	if (debinarize)
	  cicada::debinarize(graph_kbest);
	
	if (! no_id)
	  os << id << " ||| ";
	
	if (graphviz_mode)
	  os << cicada::graphviz(os, graph_kbest) << '\n';
	else {
	  if (treebank_mode)
	    cicada::treebank(os, graph_kbest);
	  else
	    os << graph_kbest;
	  os << " |||";
	  
	  typename hypergraph_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
	  for (typename hypergraph_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
	    os << ' ' << fiter->first << '=' << fiter->second;
	  os << " ||| ";
	  os << weight;
	  os << '\n';
	}
      }
    }


    template <typename Hypergraph, typename Traversal, typename Function, typename Filter>
    inline
    void kbest_derivations(std::ostream& os,
			   const size_t id,
			   const Hypergraph& graph,
			   const int kbest_size,
			   const Traversal& traversal, 
			   const Function& function,
			   const Filter& filter,
			   const bool no_id)
    {
      typedef Hypergraph hypergraph_type;
      typedef typename hypergraph_type::rule_type rule_type;
      typedef typename hypergraph_type::feature_set_type feature_set_type;
      
      if (! graph.is_valid()) {
	if (! no_id)
	  os << id << " |||";
	os << " ||| ||| 0" << '\n';
	return;
      }
      
      cicada::KBest<Traversal, Function, Filter> derivations(graph, kbest_size, traversal, function, filter);
  
      typename Traversal::value_type derivation;
      typename Function::value_type  weight;
  
      for (int k = 0; k < kbest_size; ++ k) {
	if (! derivations(k, derivation, weight))
	  break;
    
	if (! no_id)
	  os << id << " ||| ";
	os << boost::get<0>(derivation) << " |||";
	typename hypergraph_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
	for (typename hypergraph_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
	  os << ' ' << fiter->first << '=' << fiter->second;
	os << " ||| ";
	os << weight;
	os << '\n';
      }
    }

    Output::Output(const std::string& parameter, output_data_type& __output_data, const int __debug)
      : base_type("output"),
	output_data(__output_data), file(), directory(), weights(0), weights_assigned(0), weights_one(false), weights_fixed(false),
	kbest_size(0), kbest_unique(false),
	insertion_prefix(),
	yield_string(false),
	yield_terminal_pos(false),
	yield_tree(false),
	yield_graphviz(false),
	yield_treebank(false),
	yield_alignment(false),
	yield_dependency(false),
	yield_span(false),
	
	debinarize(false),
	graphviz(false),
	statistics(false),
	lattice_mode(false),
	forest_mode(false),
	no_id(false),
	debug(__debug)
    {
      typedef cicada::Parameter param_type;
    
      param_type param(parameter);
      if (utils::ipiece(param.name()) != "output")
	throw std::runtime_error("this is not a outputter");
      
      for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "kbest")
	  kbest_size = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "unique")
	  kbest_unique = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  weights = &base_type::weights(piter->second);
	else if (utils::ipiece(piter->first) == "weights-one")
	  weights_one = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "debinarize")
	  debinarize = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "graphviz")
	  graphviz = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "statistics")
	  statistics = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "lattice")
	  lattice_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "forest")
	  forest_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-id")
	  no_id = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "file")
	  file = piter->second;
	else if (utils::ipiece(piter->first) == "directory")
	  directory = piter->second;
	else if (utils::ipiece(piter->first) == "insertion-prefix")
	  insertion_prefix = piter->second;
	else if (utils::ipiece(piter->first) == "yield") {
	  const utils::ipiece value = piter->second;
	
	  if (value == "sentence" || value == "string")
	    yield_string = true;
	  else if (value == "sentence-pos" || value == "terminal-pos")
	    yield_terminal_pos = true;
	  else if (value == "derivation" || value == "tree")
	    yield_tree = true;
	  else if (value == "graphviz")
	    yield_graphviz = true;
	  else if (value == "treebank")
	    yield_treebank = true;
	  else if (value == "alignment" || value == "align")
	    yield_alignment = true;
	  else if (value == "dependency" || value == "dep")
	    yield_dependency = true;
	  else if (value == "span")
	    yield_span = true;
	  else
	    throw std::runtime_error("unknown yield: " + piter->second);
	} else
	  std::cerr << "WARNING: unsupported parameter for output: " << piter->first << "=" << piter->second << std::endl;
      }

    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");
      
      if (weights || weights_one)
	weights_fixed = true;
      
      if (! weights)
	weights = &base_type::weights();
      
      // default to stdout
      if (directory.empty() && file.empty())
	file = "-";
    
      if (! directory.empty() && ! file.empty())
	throw std::runtime_error("you cannot output both in directory and file");
	
      if (int(yield_string) + yield_terminal_pos + yield_tree + yield_graphviz + yield_treebank + yield_alignment + yield_dependency + yield_span > 1)
	throw std::runtime_error("only string, tree or alignment yield for kbest");
	
      if (int(yield_string) + yield_terminal_pos + yield_tree + yield_graphviz + yield_treebank + yield_alignment + yield_dependency + yield_span == 0)
	yield_string = true;
	
      if (graphviz && statistics)
	throw std::runtime_error("only one of graphviz or statistics can be specified...");

      if (graphviz || yield_graphviz)
	no_id = true;
      
      if (lattice_mode && forest_mode && graphviz)
	throw std::runtime_error("only one of lattice or forest can be dumped");
      
      if (int(lattice_mode) + forest_mode == 0)
	forest_mode = true;
    }

    void Output::assign(const weight_set_type& __weights)
    {
      if (! weights_fixed)
	weights_assigned = &__weights;
    }

    void Output::clear()
    {
      output_data.buffer.clear();
	
      if (output_data.os) {
	if (! directory.empty())
	  output_data.os.reset();
	else
	  *output_data.os << std::flush;
      }
    }

    void Output::operator()(data_type& data) const
    {
      typedef cicada::semiring::Logprob<double> weight_type;

      const size_type& id = data.id;
      const hypergraph_type& hypergraph = data.hypergraph;
    	
      boost::iostreams::filtering_ostream os_buffer;
      if (output_data.use_buffer)
	os_buffer.push(boost::iostreams::back_inserter(const_cast<std::string&>(output_data.buffer)));
      else if (! output_data.os) {
	const path_type path = (! file.empty() ? file  : directory / (utils::lexical_cast<std::string>(id) + ".gz"));
	const_cast<boost::shared_ptr<std::ostream>&>(output_data.os).reset(new utils::compress_ostream(path, 1024 * 1024));
      }
	
      std::ostream& os = (output_data.use_buffer
			  ? static_cast<std::ostream&>(os_buffer)
			  : *output_data.os);
      
      os.precision(10);

      if (debug)
	std::cerr << name << ": " << data.id << std::endl;

      utils::resource start;

      if (statistics) {
	if (lattice_mode) {
	  const size_t num_nodes = data.lattice.size() + 1;
	  size_t num_edges = 0;
	  size_t num_epsilon = 0;
	    
	  lattice_type::const_iterator liter_end = data.lattice.end();
	  for (lattice_type::const_iterator liter = data.lattice.begin(); liter != liter_end; ++ liter) {
	    num_edges += liter->size();
	    
	    lattice_type::arc_set_type::const_iterator aiter_end = liter->end();
	    for (lattice_type::arc_set_type::const_iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	      num_epsilon += aiter->label == vocab_type::EPSILON;
	  }
	  
	  if (no_id)
	    os << "lattice-num-node: "    << num_nodes << '\n'
	       << "lattice-num-edge: "    << num_edges << '\n'
	       << "lattice-num-epsilon: " << num_epsilon << '\n'
	       << "lattice-shortest-distance: " << data.lattice.shortest_distance() << '\n'
	       << "lattice-longest-distance: "  << data.lattice.longest_distance() << '\n';
	  else
	    os << id << " ||| lattice-num-node: "    << num_nodes << '\n'
	       << id << " ||| lattice-num-edge: "    << num_edges << '\n'
	       << id << " ||| lattice-num-epsilon: " << num_epsilon << '\n'
	       << id << " ||| lattice-shortest-distance: " << data.lattice.shortest_distance() << '\n'
	       << id << " ||| lattice-longest-distance: "  << data.lattice.longest_distance() << '\n';
	}
	
	if (forest_mode) {
	  if (data.hypergraph.is_valid()) {
	    std::vector<shortest_length_function::value_type, std::allocator<shortest_length_function::value_type> > lengths_shortest(data.hypergraph.nodes.size());
	    std::vector<longest_length_function::value_type, std::allocator<longest_length_function::value_type> >   lengths_longest(data.hypergraph.nodes.size());
	    
	    cicada::inside(data.hypergraph, lengths_shortest, shortest_length_function());
	    cicada::inside(data.hypergraph, lengths_longest, longest_length_function());
	    
	    const int length_shortest = - log(lengths_shortest.back());
	    const int length_longest  =   log(lengths_longest.back());
	    
	    if (no_id)
	      os << "hypergraph-num-node: " << data.hypergraph.nodes.size() << '\n'
		 << "hypergraph-num-edge: " << data.hypergraph.edges.size() << '\n'
		 << "hypergraph-shortest-leaf: " << length_shortest << '\n'
		 << "hypergraph-longest-leaf: "  << length_longest << '\n';
	    else
	      os << id << " ||| hypergraph-num-node: " << data.hypergraph.nodes.size() << '\n'
		 << id << " ||| hypergraph-num-edge: " << data.hypergraph.edges.size() << '\n'
		 << id << " ||| hypergraph-shortest-leaf: " << length_shortest << '\n'
		 << id << " ||| hypergraph-longest-leaf: "  << length_longest << '\n';
	  } else {
	    if (no_id)
	      os << "hypergraph-num-node: " << 0 << '\n'
		 << "hypergraph-num-edge: " << 0 << '\n'
		 << "hypergraph-shortest-leaf: " << 0 << '\n'
		 << "hypergraph-longest-leaf: "  << 0 << '\n';
	    else
	      os << id << " ||| hypergraph-num-node: " << 0 << '\n'
		 << id << " ||| hypergraph-num-edge: " << 0 << '\n'
		 << id << " ||| hypergraph-shortest-leaf: " << 0 << '\n'
		 << id << " ||| hypergraph-longest-leaf: "  << 0 << '\n';
	  }
	}
      } else if (graphviz) {
	if (! no_id)
	  os << id << " ||| ";
	if (lattice_mode)
	  cicada::graphviz(os, data.lattice);
	else
	  cicada::graphviz(os, hypergraph);
	os << '\n';
      } else if (kbest_size <= 0) {	
	if (! no_id)
	  os << id << " ||| ";
	if (lattice_mode && forest_mode)
	  os << data.lattice << " ||| " << hypergraph;
	else if (lattice_mode)
	  os << data.lattice;
	else 
	  os << hypergraph;
	
	os << '\n';
      } else {
	const weight_set_type* weights_kbest = (weights_assigned ? weights_assigned : &(weights->weights));
	
	if (weights_one) {
	  if (kbest_unique) {
	    if (yield_alignment)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				alignment_feature_traversal(),
				weight_function_one<weight_type>(),
				kbest_alignment_filter_unique(hypergraph),
				no_id);
	    else if (yield_dependency)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				dependency_feature_traversal(),
				weight_function_one<weight_type>(),
				kbest_dependency_filter_unique(hypergraph),
				no_id);
	    else if (yield_span)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				span_feature_traversal(),
				weight_function_one<weight_type>(),
				kbest_span_filter_unique(hypergraph),
				no_id);
	    else if (yield_string)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_feature_traversal(insertion_prefix),
				weight_function_one<weight_type>(),
				kbest_sentence_filter_unique(hypergraph),
				no_id);
	    else if (yield_terminal_pos)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_pos_feature_traversal(insertion_prefix),
				weight_function_one<weight_type>(),
				kbest_sentence_filter_unique(hypergraph),
				no_id);
	    else
	      kbest_derivations(os, id, hypergraph, kbest_size,
				weight_function_one<weight_type>(),
				kbest_sentence_filter(),
				no_id,
				yield_graphviz,
				yield_treebank,
				debinarize);
	  } else {
	    if (yield_alignment)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				alignment_feature_traversal(),
				weight_function_one<weight_type>(),
				kbest_alignment_filter(),
				no_id);
	    else if (yield_dependency)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				dependency_feature_traversal(),
				weight_function_one<weight_type>(),
				kbest_dependency_filter(),
				no_id);
	    else if (yield_span)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				span_feature_traversal(),
				weight_function_one<weight_type>(),
				kbest_span_filter(),
				no_id);
	    else if (yield_string)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_feature_traversal(insertion_prefix),
				weight_function_one<weight_type>(),
				kbest_sentence_filter(),
				no_id);
	    else if (yield_terminal_pos)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_pos_feature_traversal(insertion_prefix),
				weight_function_one<weight_type>(),
				kbest_sentence_filter(),
				no_id);
	    else
	      kbest_derivations(os, id, hypergraph, kbest_size,
				weight_function_one<weight_type>(),
				kbest_sentence_filter(),
				no_id,
				yield_graphviz,
				yield_treebank,
				debinarize);
	  }
	} else {
	  if (kbest_unique) {
	    if (yield_alignment)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				alignment_feature_traversal(),
				weight_function<weight_type>(*weights_kbest),
				kbest_alignment_filter_unique(hypergraph),
				no_id);
	    else if (yield_dependency)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				dependency_feature_traversal(),
				weight_function<weight_type>(*weights_kbest),
				kbest_dependency_filter_unique(hypergraph),
				no_id);
	    else if (yield_span)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				span_feature_traversal(),
				weight_function<weight_type>(*weights_kbest),
				kbest_span_filter_unique(hypergraph),
				no_id);
	    else if (yield_string)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_feature_traversal(insertion_prefix),
				weight_function<weight_type>(*weights_kbest),
				kbest_sentence_filter_unique(hypergraph),
				no_id);
	    else if (yield_terminal_pos)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_pos_feature_traversal(insertion_prefix),
				weight_function<weight_type>(*weights_kbest),
				kbest_sentence_filter_unique(hypergraph),
				no_id);
	    else
	      kbest_derivations(os, id, hypergraph, kbest_size,
				weight_function<weight_type>(*weights_kbest),
				kbest_sentence_filter(),
				no_id,
				yield_graphviz,
				yield_treebank,
				debinarize);
	  } else {
	    if (yield_alignment)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				alignment_feature_traversal(),
				weight_function<weight_type>(*weights_kbest),
				kbest_alignment_filter(),
				no_id);
	    else if (yield_dependency)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				dependency_feature_traversal(),
				weight_function<weight_type>(*weights_kbest),
				kbest_dependency_filter(),
				no_id);
	    else if (yield_span)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				span_feature_traversal(),
				weight_function<weight_type>(*weights_kbest),
				kbest_span_filter(),
				no_id);
	    else if (yield_string)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_feature_traversal(insertion_prefix),
				weight_function<weight_type>(*weights_kbest),
				kbest_sentence_filter(),
				no_id);
	    else if (yield_terminal_pos)
	      kbest_derivations(os, id, hypergraph, kbest_size,
				sentence_pos_feature_traversal(insertion_prefix),
				weight_function<weight_type>(*weights_kbest),
				kbest_sentence_filter(),
				no_id);
	    else
	      kbest_derivations(os, id, hypergraph, kbest_size,
				weight_function<weight_type>(*weights_kbest),
				kbest_sentence_filter(),
				no_id,
				yield_graphviz,
				yield_treebank,
				debinarize);
	  }
	}
      }
	
      utils::resource end;
	
      if (debug)
	std::cerr << name << ": " << data.id
		  << " cpu time: " << (end.cpu_time() - start.cpu_time())
		  << " user time: " << (end.user_time() - start.user_time())
		  << " thread time: " << (end.thread_time() - start.thread_time())
		  << std::endl;
      
      statistics_type::statistic_type& stat = data.statistics[name];
      
      ++ stat.count;
      stat.user_time += (end.user_time() - start.user_time());
      stat.cpu_time  += (end.cpu_time() - start.cpu_time());
      stat.thread_time  += (end.thread_time() - start.thread_time());
    }
  };
};
