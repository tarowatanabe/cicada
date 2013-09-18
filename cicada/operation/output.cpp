//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <unistd.h>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/functional/hash/hash.hpp>

#include <iostream>
#include <iterator>

#include <cicada/operation.hpp>
#include <cicada/parameter.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sample.hpp>
#include <cicada/sample_uniform.hpp>
#include <cicada/kbest.hpp>
#include <cicada/kbest_diverse.hpp>
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
#include <utils/piece.hpp>
#include <utils/compact_map.hpp>
#include <utils/sampler.hpp>

namespace cicada
{
  namespace operation
  {
    template <typename Tp>
    struct unassigned_id
    {
      Tp operator()() const { return Tp(-1); }
    };
    
    template <typename Iterator>
    struct feature_generator : boost::spirit::karma::grammar<Iterator, HyperGraph::feature_set_type()>
    {
      feature_generator() : feature_generator::base_type(features) 
      {
	namespace karma = boost::spirit::karma;
	namespace standard = boost::spirit::standard;
	
	features %= (standard::string << '=' << double10) % ' ';
      }
      
      struct real_precision : boost::spirit::karma::real_policies<double>
      {
	static unsigned int precision(double) 
	{ 
	  return 10;
	}
      };
      
      boost::spirit::karma::real_generator<double, real_precision> double10;
      boost::spirit::karma::rule<Iterator, HyperGraph::feature_set_type()> features;
    };

    struct KBestParam
    {
      // kbest parameter
      int    kbest_size;
      bool   kbest_sample;
      bool   kbest_uniform;
      double kbest_diversity;
      
      // output modification
      bool no_id;
      bool graphviz;
      bool treebank;
      bool debinarize;
      
      KBestParam() {}
    };

    template <typename Hypergraph, typename Derivations, typename Removes>
    void kbest_derivations_tree(std::ostream& os,
				const Operation::id_type id,
				const Hypergraph& graph,
				Derivations& derivations,
				const Removes& removes,
				const KBestParam& param)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      typedef Hypergraph hypergraph_type;
      typedef typename hypergraph_type::rule_type rule_type;
      typedef typename hypergraph_type::id_type   id_type;
      
      typedef edge_feature_traversal::value_type    derivation_type;
      typedef edge_feature_traversal::edge_set_type edge_set_type;
      
      typedef typename Derivations::function_type::value_type weight_type;
      
      typedef utils::compact_map<id_type, id_type,
				 unassigned_id<id_type>, unassigned_id<id_type>,
				 boost::hash<id_type>, std::equal_to<id_type>,
				 std::allocator<std::pair<const id_type, id_type> > > node_map_type;
      
      typedef std::vector<id_type, std::allocator<id_type> > head_set_type;
      
      typedef std::ostream_iterator<char> iterator_type;
      
      feature_generator<iterator_type> features;
      
      node_map_type   node_maps;
      head_set_type   heads;
      hypergraph_type graph_kbest;
      
      edge_set_type tails;

      typename Derivations::const_iterator diter_end = derivations.end();
      for (typename Derivations::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
	const weight_type& weight     = diter->first;
	derivation_type&   derivation = const_cast<derivation_type&>(diter->second);
    
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

	if (param.debinarize)
	  cicada::debinarize(graph_kbest);
	
	if (! param.no_id)
	  os << id << " ||| ";
	
	if (param.graphviz) {
	  cicada::graphviz(os, graph_kbest);
	  os << '\n';
	} else {
	  if (param.treebank)
	    cicada::treebank(os, graph_kbest);
	  else
	    os << graph_kbest;
	  
	  typename Removes::const_iterator riter_end = removes.end();
	  for (typename Removes::const_iterator riter = removes.begin(); riter != riter_end; ++ riter)
	    boost::get<1>(derivation).erase(*riter);
	  
	  karma::generate(iterator_type(os), " ||| " << features << " ||| ", boost::get<1>(derivation));
	  
	  os << weight << '\n';
	}
      }
    }
    
    template <typename Hypergraph, typename Function, typename Sampler, typename Removes>
    inline
    void kbest_derivations_tree(std::ostream& os,
				const Operation::id_type id,
				const Hypergraph& graph,
				const Function& function,
				const Sampler& sampler,
				const Removes& removes,
				const KBestParam& param)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      typedef Hypergraph hypergraph_type;
      typedef typename hypergraph_type::rule_type rule_type;
      
      if (! graph.is_valid()) {
	hypergraph_type graph_empty;
	if (! param.no_id)
	  os << id << " ||| ";
	if (param.graphviz) {
	  cicada::graphviz(os, graph_empty);
	  os << '\n';
	} else if (param.treebank) {
	  cicada::treebank(os, graph_empty);
	  os << " ||| ||| 0" << '\n';
	} else
	  os << graph_empty << " ||| ||| 0" << '\n';
	return;
      }
      
      if (param.kbest_sample) {
	cicada::Sample<edge_feature_traversal, Function, Sampler> derivations(graph, param.kbest_size, edge_feature_traversal(), function, const_cast<Sampler&>(sampler));
	
	kbest_derivations_tree(os, id, graph, derivations, removes, param);
      } else if (param.kbest_uniform) {
	cicada::SampleUniform<edge_feature_traversal, Function, Sampler> derivations(graph, param.kbest_size, edge_feature_traversal(), function, const_cast<Sampler&>(sampler));
	
	kbest_derivations_tree(os, id, graph, derivations, removes, param);
      } else if (param.kbest_diversity != 0.0) {
	cicada::KBestDiverse<edge_feature_traversal, Function, kbest_sentence_filter> derivations(graph, param.kbest_size, edge_feature_traversal(), function, kbest_sentence_filter(), param.kbest_diversity);

	kbest_derivations_tree(os, id, graph, derivations, removes, param);
      } else {
	cicada::KBest<edge_feature_traversal, Function, kbest_sentence_filter> derivations(graph, param.kbest_size, edge_feature_traversal(), function, kbest_sentence_filter());
	
	kbest_derivations_tree(os, id, graph, derivations, removes, param);
      }
    }
    
    template <typename Hypergraph, typename Derivations, typename Removes>
    inline
    void kbest_derivations(std::ostream& os,
			   const Operation::id_type id,
			   const Hypergraph& graph,
			   Derivations& derivations,
			   const Removes& removes,
			   const KBestParam& param)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      typedef Hypergraph hypergraph_type;
      typedef typename hypergraph_type::rule_type rule_type;
      typedef typename hypergraph_type::feature_set_type feature_set_type;
      
      typedef std::ostream_iterator<char> iterator_type;
	
      feature_generator<iterator_type> features;
      
      typedef typename Derivations::traversal_type::value_type derivation_type;
      typedef typename Derivations::function_type::value_type  weight_type;
      
      typename Derivations::const_iterator diter_end = derivations.end();
      for (typename Derivations::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
	const weight_type&  weight  = diter->first;
	derivation_type& derivation = const_cast<derivation_type&>(diter->second);
	
	if (! param.no_id)
	  os << id << " ||| ";
	os << boost::get<0>(derivation);
	
	typename Removes::const_iterator riter_end = removes.end();
	for (typename Removes::const_iterator riter = removes.begin(); riter != riter_end; ++ riter)
	  boost::get<1>(derivation).erase(*riter);

	karma::generate(iterator_type(os), " ||| " << features << " ||| ", boost::get<1>(derivation));
	
	os << weight << '\n';
      }
    }

    template <typename Hypergraph, typename Traversal, typename Function, typename Filter, typename Sampler, typename Removes>
    inline
    void kbest_derivations(std::ostream& os,
			   const Operation::id_type id,
			   const Hypergraph& graph,
			   const Traversal& traversal, 
			   const Function& function,
			   const Filter& filter,
			   const Sampler& sampler,
			   const Removes& removes,
			   const KBestParam& param)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;

      typedef Hypergraph hypergraph_type;
      typedef typename hypergraph_type::rule_type rule_type;
      typedef typename hypergraph_type::feature_set_type feature_set_type;
      
      if (! graph.is_valid()) {
	if (! param.no_id)
	  os << id << " |||";
	os << " ||| ||| 0" << '\n';
	return;
      }
      
      if (param.kbest_sample) {
	cicada::Sample<Traversal, Function, Sampler> derivations(graph, param.kbest_size, traversal, function, const_cast<Sampler&>(sampler));
	
	kbest_derivations(os, id, graph, derivations, removes, param);
      } else if (param.kbest_uniform) {
	cicada::SampleUniform<Traversal, Function, Sampler> derivations(graph, param.kbest_size, traversal, function, const_cast<Sampler&>(sampler));
	
	kbest_derivations(os, id, graph, derivations, removes, param);
      } else if (param.kbest_diversity != 0.0) {
	cicada::KBestDiverse<Traversal, Function, Filter> derivations(graph, param.kbest_size, traversal, function, filter, param.kbest_diversity);
	
	kbest_derivations(os, id, graph, derivations, removes, param);
      } else {
	cicada::KBest<Traversal, Function, Filter> derivations(graph, param.kbest_size, traversal, function, filter);
	
	kbest_derivations(os, id, graph, derivations, removes, param);
      }
    }

    template <typename Hypergraph, typename Sampler, typename Removes>
    inline
    void kbest_derivations_tree(std::ostream& os,
				const Operation::id_type id,
				const Hypergraph& graph,
				const Operation::weight_set_type* weights,
				const Operation::feature_set_type& weights_extra,
				const Sampler& sampler,
				const Removes& removes,
				const KBestParam& param)
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      
      if (weights) {
	if (! weights_extra.empty())
	  kbest_derivations_tree(os, id, graph, weight_function_extra<weight_type>(*weights, weights_extra.begin(), weights_extra.end()), sampler, removes, param);
	else
	  kbest_derivations_tree(os, id, graph, weight_function<weight_type>(*weights), sampler, removes, param);
      } else
	kbest_derivations_tree(os, id, graph, weight_function_one<weight_type>(), sampler, removes, param);
    }

    template <typename Hypergraph, typename Traversal, typename Filter, typename Sampler, typename Removes>
    inline
    void kbest_derivations(std::ostream& os,
			   const Operation::id_type id,
			   const Hypergraph& graph,
			   const Traversal& traversal, 
			   const Filter& filter,
			   const Operation::weight_set_type* weights,
			   const Operation::feature_set_type& weights_extra,
			   const Sampler& sampler,
			   const Removes& removes,
			   const KBestParam& param)
    {
      typedef cicada::semiring::Logprob<double> weight_type;
      
      if (weights) {
	if (! weights_extra.empty())
	  kbest_derivations(os, id, graph, traversal, weight_function_extra<weight_type>(*weights, weights_extra.begin(), weights_extra.end()), filter, sampler, removes, param);
	else
	  kbest_derivations(os, id, graph, traversal, weight_function<weight_type>(*weights), filter, sampler, removes, param);
      } else
	kbest_derivations(os, id, graph, traversal, weight_function_one<weight_type>(), filter, sampler, removes, param);
    }

    Output::Output(const std::string& parameter, output_data_type& __output_data, const int __debug)
      : base_type("output"),
	output_data(__output_data), file(), directory(), weights(0), weights_assigned(0),
	weights_one(false), weights_fixed(false), weights_extra(),
	kbest_size(0), kbest_unique(false), kbest_sample(false), kbest_uniform(false), diversity(0.0),
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
        span_mode(false),
        alignment_mode(false),
        dependency_mode(false),
        bitext_mode(false),
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
	else if (utils::ipiece(piter->first) == "diversity")
	  diversity = utils::lexical_cast<double>(piter->second);
	else if (utils::ipiece(piter->first) == "unique")
	  kbest_unique = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "sample")
	  kbest_sample = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "uniform")
	  kbest_uniform = utils::lexical_cast<bool>(piter->second);
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
	else if (utils::ipiece(piter->first) == "span")
	  span_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "alignment")
	  alignment_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "dependency")
	  dependency_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "bitext")
	  bitext_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "no-id")
	  no_id = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "file")
	  file = piter->second;
	else if (utils::ipiece(piter->first) == "directory")
	  directory = piter->second;
	else if (utils::ipiece(piter->first) == "insertion-prefix")
	  insertion_prefix = piter->second;
	else if (utils::ipiece(piter->first) == "remove-feature")
	  removes.push_back(piter->second);
	else if (utils::ipiece(piter->first) == "yield") {
	  const utils::ipiece value = piter->second;
	
	  if (value == "sentence" || value == "string")
	    yield_string = true;
	  else if (value == "sentence-pos" || value == "terminal-pos")
	    yield_terminal_pos = true;
	  else if (value == "derivation" || value == "tree" || value == "hypergraph")
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
	} else if (utils::ipiece(piter->first) == "weight") {
	  namespace qi = boost::spirit::qi;
	  namespace standard = boost::spirit::standard;
	  
	  std::string::const_iterator iter = piter->second.begin();
	  std::string::const_iterator iter_end = piter->second.end();

	  std::string name;
	  double      value;
	  
	  if (! qi::phrase_parse(iter, iter_end,
				 qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi))
					      >> (standard::char_ - standard::space))]
				 >> '='
				 >> qi::double_,
				 standard::blank, name, value) || iter != iter_end)
	    throw std::runtime_error("weight parameter parsing failed");
	  
	  weights_extra[name] = value;
	} else
	  std::cerr << "WARNING: unsupported parameter for output: " << piter->first << "=" << piter->second << std::endl;
      }

    
      if (weights && weights_one)
	throw std::runtime_error("you have weights, but specified all-one parameter");

      if (weights_one && ! weights_extra.empty())
	throw std::runtime_error("you have extra weights, but specified all-one parameter");
      
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
	throw std::runtime_error("only one of string, tree, graphviz, treebank, alignment, dependency or span yield for kbest");
	
      if (int(yield_string) + yield_terminal_pos + yield_tree + yield_graphviz + yield_treebank + yield_alignment + yield_dependency + yield_span == 0)
	yield_string = true;

      if (kbest_size > 0)
	if (int(kbest_sample) + kbest_uniform + (diversity != 0.0) > 1)
	  throw std::runtime_error("when kbest, only one of sample or uniform, diversity can be specified");
	
      if (graphviz && statistics)
	throw std::runtime_error("only one of graphviz or statistics can be specified...");

      if (graphviz || yield_graphviz)
	no_id = true;
      
      if (lattice_mode && forest_mode && graphviz)
	throw std::runtime_error("only one of lattice or forest can be dumped");
      
      if (int(lattice_mode) + forest_mode + span_mode + alignment_mode + dependency_mode + bitext_mode == 0)
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

      const id_type& id = data.id;
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
	    std::vector<length_function::value_type, std::allocator<length_function::value_type> >lengths(data.hypergraph.nodes.size());
	    
	    cicada::inside(data.hypergraph, lengths, length_function());
	    
	    const int length_shortest = - cicada::semiring::log(lengths.back().first);
	    const int length_longest  =   cicada::semiring::log(lengths.back().second);
	    
#if 0
	    std::vector<shortest_length_function::value_type, std::allocator<shortest_length_function::value_type> > lengths_shortest(data.hypergraph.nodes.size());
	    std::vector<longest_length_function::value_type, std::allocator<longest_length_function::value_type> >   lengths_longest(data.hypergraph.nodes.size());
	    
	    cicada::inside(data.hypergraph, lengths_shortest, shortest_length_function());
	    cicada::inside(data.hypergraph, lengths_longest, longest_length_function());
	    
	    const int length_shortest = - log(lengths_shortest.back());
	    const int length_longest  =   log(lengths_longest.back());
#endif
	    
	    // # of unaries
	    int unaries = 0;
	    {
	      hypergraph_type::edge_set_type::const_iterator eiter_end = data.hypergraph.edges.end();
	      for (hypergraph_type::edge_set_type::const_iterator eiter = data.hypergraph.edges.begin(); eiter != eiter_end; ++ eiter) {
		const hypergraph_type::edge_type& edge = *eiter;
		
		unaries += (edge.tails.size() == 1 && edge.rule->rhs.size() == 1);
	      }
	    }
	    
	    if (no_id)
	      os << "hypergraph-num-node: " << data.hypergraph.nodes.size() << '\n'
		 << "hypergraph-num-edge: " << data.hypergraph.edges.size() << '\n'
		 << "hypergraph-num-unary: " << unaries << '\n'
		 << "hypergraph-shortest-leaf: " << length_shortest << '\n'
		 << "hypergraph-longest-leaf: "  << length_longest << '\n';
	    else
	      os << id << " ||| hypergraph-num-node: " << data.hypergraph.nodes.size() << '\n'
		 << id << " ||| hypergraph-num-edge: " << data.hypergraph.edges.size() << '\n'
		 << id << " ||| hypergraph-num-unary: " << unaries << '\n'
		 << id << " ||| hypergraph-shortest-leaf: " << length_shortest << '\n'
		 << id << " ||| hypergraph-longest-leaf: "  << length_longest << '\n';
	  } else {
	    if (no_id)
	      os << "hypergraph-num-node: " << 0 << '\n'
		 << "hypergraph-num-edge: " << 0 << '\n'
		 << "hypergraph-num-unary: " << 0 << '\n'
		 << "hypergraph-shortest-leaf: " << 0 << '\n'
		 << "hypergraph-longest-leaf: "  << 0 << '\n';
	    else
	      os << id << " ||| hypergraph-num-node: " << 0 << '\n'
		 << id << " ||| hypergraph-num-edge: " << 0 << '\n'
		 << id << " ||| hypergraph-num-unary: " << 0 << '\n'
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
	
	bool need_separator = false;
	if (! no_id)
	  os << id << " ||| ";
	
	if (lattice_mode && forest_mode) {
	  os << data.lattice << " ||| " << hypergraph;
	  need_separator = true;
	} else if (lattice_mode) {
	  os << data.lattice;
	  need_separator = true;
	} else {
	  os << hypergraph;
	  need_separator = true;
	}
	
	if (span_mode) {
	  if (need_separator)
	    os << " ||| ";
	  os << data.spans;
	  need_separator = true;
	}
	
	if (alignment_mode) {
	  if (need_separator)
	    os << " ||| ";
	  os << data.alignment;
	  need_separator = true;
	}
	
	if (dependency_mode) {
	  if (need_separator)
	    os << " ||| ";
	  os << data.dependency;
	  need_separator = true;
	}
	
	if (bitext_mode) {
	  if (need_separator)
	    os << " ||| ";
	  os << data.targets;
	  need_separator = true;
	}
	
	os << '\n';
      } else {
	const weight_set_type* weights_kbest = (weights_one ? 0 : (weights_assigned ? weights_assigned : &(weights->weights)));
	
	// kbest-parameters...
	KBestParam param;
	param.kbest_size      = kbest_size;
	param.kbest_sample    = kbest_sample;
	param.kbest_uniform   = kbest_uniform;
	param.kbest_diversity = diversity;
	param.no_id      = no_id;
	param.graphviz   = yield_graphviz;
	param.treebank   = yield_treebank;
	param.debinarize = debinarize;
	
	if (kbest_unique) {
	  if (yield_alignment)
	    kbest_derivations(os, id, hypergraph,
			      alignment_feature_traversal(),
			      kbest_alignment_filter_unique(hypergraph),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_dependency)
	    kbest_derivations(os, id, hypergraph,
			      dependency_feature_traversal(),
			      kbest_dependency_filter_unique(hypergraph),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_span)
	    kbest_derivations(os, id, hypergraph,
			      span_feature_traversal(),
			      kbest_span_filter_unique(hypergraph),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_string)
	    kbest_derivations(os, id, hypergraph,
			      sentence_feature_traversal(insertion_prefix),
			      kbest_sentence_filter_unique(hypergraph),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_terminal_pos)
	    kbest_derivations(os, id, hypergraph,
			      sentence_pos_feature_traversal(insertion_prefix),
			      kbest_sentence_filter_unique(hypergraph),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else
	    kbest_derivations_tree(os, id, hypergraph,
				   weights_kbest,
				   weights_extra,
				   sampler,
				   removes,
				   param);
	} else {
	  if (yield_alignment)
	    kbest_derivations(os, id, hypergraph,
			      alignment_feature_traversal(),
			      kbest_alignment_filter(),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_dependency)
	    kbest_derivations(os, id, hypergraph,
			      dependency_feature_traversal(),
			      kbest_dependency_filter(),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_span)
	    kbest_derivations(os, id, hypergraph,
			      span_feature_traversal(),
			      kbest_span_filter(),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_string)
	    kbest_derivations(os, id, hypergraph,
			      sentence_feature_traversal(insertion_prefix),
			      kbest_sentence_filter(),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else if (yield_terminal_pos)
	    kbest_derivations(os, id, hypergraph,
			      sentence_pos_feature_traversal(insertion_prefix),
			      kbest_sentence_filter(),
			      weights_kbest,
			      weights_extra,
			      sampler,
			      removes,
			      param);
	  else
	    kbest_derivations_tree(os, id, hypergraph,
				   weights_kbest,
				   weights_extra,
				   sampler,
				   removes,
				   param);
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
