#include <stdexcept>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>

#include <utils/compress_stream.hpp>
#include <utils/filesystem.hpp>
#include <utils/lexical_cast.hpp>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"

#include "cicada/kbest.hpp"
#include "cicada/parameter.hpp"
#include "cicada/graphviz.hpp"

#include "cicada/model.hpp"
#include "cicada/grammar.hpp"

#include "cicada/apply.hpp"
#include "cicada/compose.hpp"
#include "cicada/inside_outside.hpp"
#include "cicada/intersect.hpp"
#include "cicada/binarize.hpp"
#include "cicada/permute.hpp"
#include "cicada/sort.hpp"
#include "cicada/prune.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"
#include "cicada/span_vector.hpp"

#include "cicada/feature/variational.hpp"
#include "cicada/feature/bleu.hpp"

#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/sgi_hash_set.hpp"

typedef boost::filesystem::path path_type;

typedef std::vector<std::string, std::allocator<std::string> > grammar_file_set_type;
typedef std::vector<std::string, std::allocator<std::string> > feature_parameter_set_type;

typedef cicada::Symbol          symbol_type;
typedef cicada::Vocab           vocab_type;
typedef cicada::Sentence        sentence_type;
typedef cicada::Lattice         lattice_type;
typedef cicada::Rule            rule_type;
typedef cicada::HyperGraph      hypergraph_type;
typedef cicada::Grammar         grammar_type;
typedef cicada::Model           model_type;
typedef cicada::FeatureFunction feature_function_type;

typedef feature_function_type::feature_function_ptr_type feature_function_ptr_type;

typedef rule_type::feature_set_type    feature_set_type;
typedef feature_set_type::feature_type feature_type;
typedef cicada::WeightVector<double>   weight_set_type;

typedef cicada::SpanVector span_set_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;

struct source_length_function
{
  typedef cicada::semiring::Tropical<int> value_type;
  
  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    int length = 0;
    rule_type::symbol_set_type::const_iterator siter_end = edge.rule->source.end();
    for (rule_type::symbol_set_type::const_iterator siter = edge.rule->source.begin(); siter != siter_end; ++ siter)
      length += (*siter != vocab_type::EPSILON && siter->is_terminal());
    
    // since we will "max" at operator+, we will collect negative length
    return cicada::semiring::traits<value_type>::log(- length);
  }
};

template <typename Weight>
struct weight_set_scaled_function
{
  typedef Weight value_type;
  
  weight_set_scaled_function(const weight_set_type& __weights, const double& __scale)
    : weights(__weights), scale(__scale) {}
  
  const weight_set_type& weights;
  const double scale;
  
  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<value_type>::log(edge.features.dot(weights) * scale);
  }
  
};

struct weight_set_function
{
  typedef cicada::semiring::Logprob<double> value_type;

  weight_set_function(const weight_set_type& __weights)
    : weights(__weights) {}

  const weight_set_type& weights;
  
  template <typename FeatureSet>
  value_type operator()(const FeatureSet& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.dot(weights));
  }
};


struct weight_set_function_one
{
  typedef cicada::semiring::Logprob<double> value_type;

  weight_set_function_one(const weight_set_type& __weights) {}
  
  template <typename FeatureSet>
  value_type operator()(const FeatureSet& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.dot());
  }
};


struct kbest_function
{
  typedef rule_type::feature_set_type feature_set_type;

  typedef cicada::semiring::Logprob<double> value_type;

  kbest_function(const weight_set_type& __weights)
    : weights(__weights) {}

  const weight_set_type& weights;
  
  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
  }

  value_type operator()(const feature_set_type& features) const
  {
    return cicada::semiring::traits<value_type>::log(features.dot(weights));
  }

};

struct kbest_function_one
{
  typedef rule_type::feature_set_type feature_set_type;

  typedef cicada::semiring::Logprob<double> value_type;
  
  kbest_function_one(const weight_set_type& __weights) {}

  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<value_type>::log(edge.features.dot());
  }

  value_type operator()(const feature_set_type& features) const
  {
    return cicada::semiring::traits<value_type>::log(features.dot());
  }
};

struct kbest_traversal_edges
{
  typedef rule_type::feature_set_type feature_set_type;
  typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;

  typedef boost::tuple<edge_set_type, feature_set_type> value_type;
  
  template <typename Edge, typename Iterator>
  void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
  {
    boost::get<0>(yield).clear();

    boost::get<0>(yield).push_back(edge.id);
    boost::get<1>(yield) = edge.features;
    
    // collect edge and features
    for (/**/; first != last; ++ first) {
      boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*first).begin(), boost::get<0>(*first).end());
      boost::get<1>(yield) += boost::get<1>(*first);
    }
  }
};

struct kbest_traversal_source
{
  typedef rule_type::feature_set_type feature_set_type;
  
  typedef boost::tuple<sentence_type, feature_set_type> value_type;
  
  template <typename Edge, typename Iterator>
  void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
  {
    // extract source-yield, features

    boost::get<0>(yield).clear();
    boost::get<1>(yield) = edge.features;
    
    int non_terminal_pos = 0;
    rule_type::symbol_set_type::const_iterator siter_end = edge.rule->source.end();
    for (rule_type::symbol_set_type::const_iterator siter = edge.rule->source.begin(); siter != siter_end; ++ siter)
      if (siter->is_non_terminal()) {
	const int pos = siter->non_terminal_index() - 1;
	
	if (pos < 0)
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + non_terminal_pos)).begin(), boost::get<0>(*(first + non_terminal_pos)).end());
	else
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + pos)).begin(), boost::get<0>(*(first + pos)).end());

	++ non_terminal_pos;
      } else if (*siter != vocab_type::EPSILON)
	boost::get<0>(yield).push_back(*siter);
    
    // collect features...
    for (/**/; first != last; ++ first)
      boost::get<1>(yield) += boost::get<1>(*first);
  }
};

struct kbest_traversal
{
  typedef rule_type::feature_set_type feature_set_type;
  
  typedef boost::tuple<sentence_type, feature_set_type> value_type;
  
  template <typename Edge, typename Iterator>
  void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
  {
    // extract target-yield, features

    boost::get<0>(yield).clear();
    boost::get<1>(yield) = edge.features;
    
    int non_terminal_pos = 0;
    rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
    for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
      if (titer->is_non_terminal()) {
	const int pos = titer->non_terminal_index() - 1;
	
	if (pos < 0)
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + non_terminal_pos)).begin(), boost::get<0>(*(first + non_terminal_pos)).end());
	else
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + pos)).begin(), boost::get<0>(*(first + pos)).end());

	++ non_terminal_pos;
      } else if (*titer != vocab_type::EPSILON)
	boost::get<0>(yield).push_back(*titer);
    
    // collect features...
    for (/**/; first != last; ++ first)
      boost::get<1>(yield) += boost::get<1>(*first);
  }
};

struct kbest_filter
{
  kbest_filter(const hypergraph_type& graph) {}
  
  template <typename Node, typename Yield>
  bool operator()(const Node& node, const Yield& yield) const
  {
    return false;
  }
};

struct kbest_filter_unique
{
 #ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<sentence_type, boost::hash<sentence_type>, std::equal_to<sentence_type>, std::allocator<sentence_type> > unique_type;
#else
  typedef sgi::hash_set<sentence_type, boost::hash<sentence_type>, std::equal_to<sentence_type>, std::allocator<sentence_type> > unique_type;
#endif
  typedef std::vector<unique_type, std::allocator<unique_type> > unique_set_type;
 

  kbest_filter_unique(const hypergraph_type& graph) : uniques(graph.nodes.size()) {}
  
  template <typename Node, typename Yield>
  bool operator()(const Node& node, const Yield& yield) const
  {
    unique_set_type& sents = const_cast<unique_set_type&>(uniques);
    unique_type::iterator iter = sents[node.id].find(boost::get<0>(yield));
    if (iter == sents[node.id].end()) {
      sents[node.id].insert(boost::get<0>(yield));
      return false;
    } else
      return true;
  }

  unique_set_type uniques;
};

template <typename Function, typename Filter>
void kbest_derivations(std::ostream& os,
		       const size_t id,
		       const hypergraph_type& graph,
		       const int kbest_size,
		       const Function& function,
		       const Filter& filter)
{
  cicada::KBest<kbest_traversal_edges, Function, Filter> derivations(graph, kbest_size, kbest_traversal_edges(), function, filter);
  
  typedef kbest_traversal_edges::value_type    derivation_type;
  typedef kbest_traversal_edges::edge_set_type edge_set_type;
  
  typedef hypergraph_type::id_type id_type;

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<id_type, id_type, utils::hashmurmur<size_t>, std::equal_to<id_type>,
    std::allocator<std::pair<id_type, id_type> > > node_map_type;
#else
  typedef sgi::hash_map<id_type, id_type, utils::hashmurmur<size_t>, std::equal_to<id_type>,
    std::allocator<std::pair<id_type, id_type> > > node_map_type;
#endif

  derivation_type derivation;
  node_map_type   node_maps;
  hypergraph_type graph_kbest;

  edge_set_type tails;
  
  for (int k = 0; k < kbest_size; ++ k) {
    if (! derivations(k, derivation))
      break;
    
    const edge_set_type& edges = boost::get<0>(derivation);
    node_maps.clear();
    graph_kbest.clear();
    
    id_type node_id = 0;
    edge_set_type::const_iterator eiter_end = edges.end();
    for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter)
      if (node_maps.find(graph.edges[*eiter].head) == node_maps.end()) {
	node_maps[graph.edges[*eiter].head] = node_id;
	++ node_id;
      }
    
    for (id_type node = 0; node != node_id; ++ node)
      graph_kbest.add_node();
    
    for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
      const hypergraph_type::edge_type& edge = graph.edges[*eiter];
      
      tails.clear();
      hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
      for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
	node_map_type::const_iterator niter = node_maps.find(*titer);
	if (niter == node_maps.end())
	  throw std::runtime_error("no node?");
	
	tails.push_back(niter->second);
      }
      
      hypergraph_type::edge_type& edge_kbest = graph_kbest.add_edge(tails.begin(), tails.end());
      edge_kbest.rule = edge.rule;
      edge_kbest.features = edge.features;
      
      graph_kbest.connect_edge(edge_kbest.id, node_maps[edge.head]);
    }
    
    node_map_type::const_iterator niter = node_maps.find(graph.goal);
    if (niter == node_maps.end())
      throw std::runtime_error("did not reach goal?");
    
    graph_kbest.goal = niter->second;

    graph_kbest.topologically_sort();
    
    os << id << " ||| " << graph_kbest << " |||";
    rule_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
    for (rule_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
      os << ' ' << fiter->first << '=' << fiter->second;
    os << " ||| ";
    os << function(boost::get<1>(derivation));
    os << '\n';
  }
}


template <typename Traversal, typename Function, typename Filter>
void kbest_derivations(std::ostream& os,
		       const size_t id,
		       const hypergraph_type& graph,
		       const int kbest_size,
		       const Traversal& traversal, 
		       const Function& function,
		       const Filter& filter)
{
  cicada::KBest<Traversal, Function, Filter> derivations(graph, kbest_size, traversal, function, filter);
  
  typename Traversal::value_type derivation;
  
  for (int k = 0; k < kbest_size; ++ k) {
    if (! derivations(k, derivation))
      break;
    
    os << id << " ||| " << boost::get<0>(derivation) << " |||";
    rule_type::feature_set_type::const_iterator fiter_end = boost::get<1>(derivation).end();
    for (rule_type::feature_set_type::const_iterator fiter = boost::get<1>(derivation).begin(); fiter != fiter_end; ++ fiter)
      os << ' ' << fiter->first << '=' << fiter->second;
    os << " ||| ";
    os << function(boost::get<1>(derivation));
    os << '\n';
  }
}


template <typename Iterator>
inline
bool parse_id(size_t& id, Iterator& iter, Iterator end)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  namespace phoenix = boost::phoenix;
  
  using qi::phrase_parse;
  using qi::_1;
  using qi::ulong_;
  using standard::space;
  
  using phoenix::ref;
  
  return phrase_parse(iter, end, ulong_ [ref(id) = _1] >> "|||", space);
}

template <typename Iterator>
inline
bool parse_separator(Iterator& iter, Iterator end)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  namespace phoenix = boost::phoenix;
  
  using qi::phrase_parse;
  using standard::space;
  
  return phrase_parse(iter, end, "|||", space);
}

template <typename HyperGraph, typename Lattice, typename SentenceSet, typename Sentence>
inline
bool parse_line(const std::string& line,
		size_t& id,
		HyperGraph& hypergraph,
		Lattice& lattice,
		Lattice& target,
		SentenceSet& target_sentences,
		Sentence& sentence,
		const bool input_id,
		const bool input_lattice,
		const bool input_forest, 
		const bool input_bitext)
{
  std::string::const_iterator iter = line.begin();
  std::string::const_iterator end = line.end();
  
  if (input_id)
    if (! parse_id(id, iter, end))
      throw std::runtime_error("invalid id-prefixed format");
  
  if (input_lattice) {
    if (! lattice.assign(iter, end))
      throw std::runtime_error("invalid lattive format");
  } else if (input_forest) {
    if (! hypergraph.assign(iter, end))
	  throw std::runtime_error("invalid hypergraph format");
  } else {
    if (! sentence.assign(iter, end))
      throw std::runtime_error("invalid sentence format");
    
    lattice = Lattice(sentence);
  }
  
  if (input_bitext) {
    target_sentences.clear();
    
    while (parse_separator(iter, end)) {
      target_sentences.push_back(Sentence());
      
      if (! target_sentences.back().assign(iter, end))
	throw std::runtime_error("invalid sentence format");
    }
    
    if (target_sentences.empty())
      throw std::runtime_error("no bitext?");
    
    target = Lattice(target_sentences.front());
  }
  
  return iter == end;
}


class Operation
{
public:
  typedef Operation base_type;

  Operation() {}
  virtual ~Operation() {}
  
  virtual void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& graph) const = 0;

  virtual void assign(const weight_set_type& weights) {}

  virtual void clear() {};

private:
  struct hash_string : public utils::hashmurmur<size_t>
  {
    size_t operator()(const std::string& x) const
    {
      return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
    }
  };

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, weight_set_type, hash_string, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, weight_set_type> > > weight_map_type;
#else
  typedef sgi::hash_map<std::string, weight_set_type, hash_string, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, weight_set_type> > > weight_map_type;

#endif

public:
  static const weight_set_type& weights(const path_type& path) 
  {
#ifdef HAVE_TLS
    static __thread weight_map_type* __weights_tls = 0;
    static boost::thread_specific_ptr<weight_map_type> __weights;
    
    if (! __weights_tls) {
      __weights.reset(new weight_map_type());
      __weights_tls = __weights.get();
    }
    weight_map_type& weights_map = *__weights_tls;
#else
    static boost::thread_specific_ptr<weight_map_type> __weights;
    
    if (! __weights.get())
      __weights.reset(new weight_map_type());
    
    weight_map_type& weights_map = *__weights;
#endif
    
    weight_map_type::iterator iter = weights_map.find(path.file_string());
    if (iter == weights_map.end()) {
      iter = weights_map.insert(std::make_pair(path.file_string(), weight_set_type())).first;
      
      if (! path.empty()) {
	if (path != "-" && ! boost::filesystem::exists(path))
	  throw std::runtime_error("no feture weights? " + path.file_string());
	
	utils::compress_istream is(path);
	is >> iter->second;
      }
    }
    return iter->second;
  }
};

class Binarize : public Operation
{
public:
  Binarize(const std::string& parameter, const int __debug)
    : size(0), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "binarize")
      throw std::runtime_error("this is not a binarizer");

    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "size") == 0)
	size = boost::lexical_cast<int>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "direction") == 0) {
	const std::string& dir = piter->second;
	
	if (strcasecmp(dir.c_str(), "left") == 0)
	  left = true;
	else if (strcasecmp(dir.c_str(), "right") == 0)
	  right = true;
	else
	  throw std::runtime_error("unuspported direction: " + parameter);
      } else
	std::cerr << "WARNING: unsupported parameter for binarize: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (! left && ! right)
      right == true;
  }
  
  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type binarized;
    
    if (debug)
      std::cerr << "binarization" << std::endl;
    
    utils::resource start;
    
    if (left)
      cicada::binarize_left(hypergraph, binarized, size);
    else if (right)
      cicada::binarize_right(hypergraph, binarized, size);
    
    utils::resource end;
    
    if (debug)
      std::cerr << "binarize cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << binarized.nodes.size()
		<< " # of edges: " << binarized.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(binarized.is_valid())
		<< std::endl;
    
    hypergraph.swap(binarized);
  }
  
  int size;
  
  bool left;
  bool right;
  
  int debug;
};

class Permute : public Operation
{
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type>, std::allocator<symbol_type> > exclude_set_type;
#else
  typedef sgi::hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type>, std::allocator<symbol_type> > exclude_set_type;
#endif

public:
  Permute(const std::string& parameter, const int __debug)
    : excludes(), weights(0), size(0), feature(false), collapse(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "permute")
      throw std::runtime_error("this is not a permuter");

    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "size") == 0)
	size = boost::lexical_cast<int>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "feature") == 0)
	feature = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "collapse") == 0)
	collapse = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "exclude") == 0)
	excludes.insert(piter->second);
      else
	std::cerr << "WARNING: unsupported parameter for permute: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (collapse && ! weights)
      throw std::runtime_error("collapsing but no weights...");
  }
  
  struct Filter
  {
    Filter(const exclude_set_type& __excludes)
      : excludes(__excludes) {}
    
    const exclude_set_type& excludes;
    
    template <typename Cat>
    bool operator()(const Cat& x) const
    {
      return ! excludes.empty() && excludes.find(x) != excludes.end();
    }
  };

  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type permuted;
    
    if (debug)
      std::cerr << "permute" << std::endl;
    
    utils::resource start;
    
    if (collapse)
      cicada::permute(hypergraph, permuted, cicada::PermuteFeatureCollapsed<weight_set_type>(*weights), Filter(excludes), size);
    else if (feature)
      cicada::permute(hypergraph, permuted, cicada::PermuteFeature(), Filter(excludes), size);
    else
      cicada::permute(hypergraph, permuted, cicada::PermuteNoFeature(), Filter(excludes), size);
    
    utils::resource end;
    
    if (debug)
      std::cerr << "permute cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << permuted.nodes.size()
		<< " # of edges: " << permuted.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(permuted.is_valid())
		<< std::endl;
    
    hypergraph.swap(permuted);
  }

  exclude_set_type excludes;
  
  const weight_set_type* weights;
  
  int size;
  bool feature;
  bool collapse;
  
  int debug;
};

class ComposeEarley : public Operation
{
public:
  ComposeEarley(const grammar_type& __grammar,
		const std::string& __goal,
		const std::string& __non_terminal,
		const bool __insertion,
		const bool __deletion,
		const int __debug)
    : grammar(__grammar),
      goal(__goal), non_terminal(__non_terminal), 
      insertion(__insertion), deletion(__deletion),
      debug(__debug)
  { }
  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type composed;
    
    if (debug)
      std::cerr << "composition" << std::endl;

    utils::resource start;

    grammar_type grammar_translation(grammar);
    
    if (insertion)
      grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(hypergraph, non_terminal)));
    if (deletion)
      grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(hypergraph, non_terminal)));

    
    cicada::compose_earley(grammar_translation, hypergraph, composed);
    
    utils::resource end;
    
    if (debug)
      std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << composed.nodes.size()
		<< " # of edges: " << composed.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		<< std::endl;
    
    hypergraph.swap(composed);
  }
  
  const grammar_type& grammar;
  
  const std::string goal;
  const std::string non_terminal;
  
  const bool insertion;
  const bool deletion;
  
  int debug;
};

class ComposeCKY : public Operation
{
public:
  ComposeCKY(const grammar_type& __grammar,
	     const std::string& __goal,
	     const std::string& __non_terminal,
	     const bool __insertion,
	     const bool __deletion,
	     const int __debug)
    : grammar(__grammar),
      goal(__goal), non_terminal(__non_terminal), 
      insertion(__insertion), deletion(__deletion),
      debug(__debug)
  { }
  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type composed;
    
    if (debug)
      std::cerr << "composition" << std::endl;

    utils::resource start;

    grammar_type grammar_translation(grammar);
    
    if (insertion)
      grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarInsertion(lattice, non_terminal)));
    if (deletion)
      grammar_translation.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarDeletion(lattice, non_terminal)));

    
    cicada::compose_cky(goal, grammar_translation, lattice, composed);
    
    utils::resource end;
    
    if (debug)
      std::cerr << "compose cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << composed.nodes.size()
		<< " # of edges: " << composed.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(composed.is_valid())
		<< std::endl;
    
    hypergraph.swap(composed);
  }
  
  const grammar_type& grammar;
  
  const std::string goal;
  const std::string non_terminal;
  
  const bool insertion;
  const bool deletion;
  
  int debug;
};


class Apply : public Operation
{
public:
  Apply(const std::string& parameter,
	const model_type& __model,
	const int __debug)
    : model(__model), weights(0), size(200), weights_one(false), exact(false), forced(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "apply")
      throw std::runtime_error("this is not a feature-functin applier");

    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "size") == 0)
	size = boost::lexical_cast<int>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "exact") == 0)
	exact = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "forced") == 0)
	forced = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	weights_one = utils::lexical_cast<bool>(piter->second);
      else
	std::cerr << "WARNING: unsupported parameter for apply: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
  }

  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type applied;
    
    // assignment...
    const_cast<model_type&>(model).assign(hypergraph);
    const_cast<model_type&>(model).assign(lattice);
    const_cast<model_type&>(model).assign(spans);
    
    if (debug)
      std::cerr << "apply features" << std::endl;

    if (forced)
      const_cast<model_type&>(model).apply_feature(true);
    
    weight_set_type weights_zero;
    const weight_set_type* weights_apply = (weights ? weights : &weights_zero);
    
    utils::resource start;
    
    // apply...
    if (model.is_stateless())
      cicada::apply_state_less(model, hypergraph, applied);
    else if (exact)
      cicada::apply_exact(model, hypergraph, applied);
    else {
      if (weights_one)
	cicada::apply_cube_prune(model, hypergraph, applied, weight_set_function_one(*weights_apply), size);
      else
	cicada::apply_cube_prune(model, hypergraph, applied, weight_set_function(*weights_apply), size);
    }
    
    utils::resource end;

    const_cast<model_type&>(model).apply_feature(false);
	
    
    if (debug)
      std::cerr << "apply cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;

    if (debug)
      std::cerr << "# of nodes: " << applied.nodes.size()
		<< " # of edges: " << applied.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(applied.is_valid())
		<< std::endl;
	
    hypergraph.swap(applied);
  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }

  const model_type& model;
  const weight_set_type* weights;
  int size;
  bool weights_one;
  bool exact;
  bool forced;
  
  int debug;
};

class Bleu : public Operation
{
public:
  Bleu(const std::string& parameter, const model_type& model, const int __debug)
    : weights(0), size(200), weights_one(false), exact(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "bleu")
      throw std::runtime_error("this is not a bleu-computer");
    
    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "size") == 0)
	size = boost::lexical_cast<int>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "exact") == 0)
	exact = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	weights_one = utils::lexical_cast<bool>(piter->second);
      else
	std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
    
    for (model_type::iterator iter = model.begin(); iter != model.end(); ++ iter) {
      cicada::feature::Bleu* __bleu = dynamic_cast<cicada::feature::Bleu*>(iter->get());
      if (__bleu)
	feature = *iter;
    }
    
    if (! feature)
      throw std::runtime_error("you have no bleu feature function");
  }

  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    int source_length = lattice.shortest_distance();
    if (hypergraph.is_valid()) {
      // we will enumerate forest structure... and collect min-size...
      std::vector<source_length_function::value_type, std::allocator<source_length_function::value_type> > lengths(hypergraph.nodes.size());
      
      cicada::inside(hypergraph, lengths, source_length_function());
      
      source_length = - log(lengths.back());
    }
    
    if (debug)
      std::cerr << "source length: " << source_length << std::endl;
    
    hypergraph_type applied;
    
    cicada::feature::Bleu* __bleu = dynamic_cast<cicada::feature::Bleu*>(feature.get());

    __bleu->clear();
    sentence_set_type::const_iterator titer_end = targets.end();
    for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
      __bleu->insert(source_length, *titer);
    
    model_type model;
    model.push_back(feature);
        
    weight_set_type weights_zero;
    const weight_set_type* weights_apply = (weights ? weights : &weights_zero);
    
    if (debug)
      std::cerr << "bleu features" << std::endl;
	
    utils::resource start;
    
    if (exact) {
      if (weights_one)
	cicada::apply_exact(model, hypergraph, applied);
      else
	cicada::apply_exact(model, hypergraph, applied);
    } else {
      if (weights_one)
	cicada::apply_cube_prune(model, hypergraph, applied, weight_set_function_one(*weights_apply), size);
      else
	cicada::apply_cube_prune(model, hypergraph, applied, weight_set_function(*weights_apply), size);
    }
    
    utils::resource end;

    __bleu->clear();
	
    if (debug)
      std::cerr << "bleu cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << applied.nodes.size()
		<< " # of edges: " << applied.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(applied.is_valid())
		<< std::endl;
    
    hypergraph.swap(applied);
  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }
  
  void clear()
  {
    dynamic_cast<cicada::feature::Bleu*>(feature.get())->clear();
  }
  
  feature_function_type::feature_function_ptr_type feature;
  
  const weight_set_type* weights;
  int size;
  bool weights_one;
  bool exact;

  int debug;
};

class Variational : public Operation
{
public:
  Variational(const std::string& parameter, const model_type& model, const int __debug)
    : weights(0), weights_one(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "variational")
      throw std::runtime_error("this is not a variational decoder");
    
    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	weights_one = utils::lexical_cast<bool>(piter->second);
      else
	std::cerr << "WARNING: unsupported parameter for variational: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");


    for (model_type::iterator iter = model.begin(); iter != model.end(); ++ iter) {
      cicada::feature::Variational* __variational = dynamic_cast<cicada::feature::Variational*>(iter->get());
      if (__variational)
	feature = *iter;
    }
    
    if (! feature)
      throw std::runtime_error("you have no variational feature function");
  }
    
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type variational;

    // clear weights to one if feature-weights-one...
    
    weight_set_type __weights;
    if (weights_one) {
      __weights.allocate();
      for (weight_set_type::feature_type::id_type id = 0; id != __weights.size(); ++ id)
	if (! weight_set_type::feature_type(id).empty())
	  __weights[weight_set_type::feature_type(id)] = 1.0;
    }
    
    const weight_set_type* weights_variational(weights ? weights : &__weights);
    
    if (debug)
      std::cerr << "variational decoding" << std::endl;
    
    cicada::feature::Variational* __variational = dynamic_cast<cicada::feature::Variational*>(feature.get());
    
    utils::resource start;
    
    // first, compute vatiational model
    __variational->insert(hypergraph, *weights_variational);
    
    model_type model;
    model.push_back(feature);
        
    // second, apply again...
    cicada::apply_exact(model, hypergraph, variational);
    
    utils::resource end;

    __variational->clear();

    if (debug)
      std::cerr << "variational cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
	
    if (debug)
      std::cerr << "# of nodes: " << variational.nodes.size()
		<< " # of edges: " << variational.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(variational.is_valid())
		<< std::endl;
    
    hypergraph.swap(variational);
  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }
  
  void clear()
  {
    dynamic_cast<cicada::feature::Variational*>(feature.get())->clear();
  }
  
  feature_function_type::feature_function_ptr_type feature;
  
  const weight_set_type* weights;
  bool weights_one;
  
  int debug;
};

class Prune : public Operation
{
public:
  Prune(const std::string& parameter, const int __debug)
    : weights(0), beam(0.0), density(0.0), scale(1.0), weights_one(false), 
      semiring_tropical(false), semiring_logprob(false), semiring_log(false),
      debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "prune")
      throw std::runtime_error("this is not a pruner");

    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "beam") == 0)
	beam = boost::lexical_cast<double>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "density") == 0)
	density = boost::lexical_cast<double>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "scale") == 0)
	scale = boost::lexical_cast<double>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	weights_one = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "semiring") == 0) {
	const std::string& name = piter->second;
	
	if (strcasecmp(name.c_str(), "tropical") == 0)
	  semiring_tropical = true;
	else if (strcasecmp(name.c_str(), "logprob") == 0)
	  semiring_logprob = true;
	else if (strcasecmp(name.c_str(), "log") == 0)
	  semiring_log = true;
	else
	  throw std::runtime_error("unknown semiring: " + name);
	
      } else
	std::cerr << "WARNING: unsupported parameter for prune: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (beam > 0.0 && density > 1.0)
      throw std::runtime_error("you cannot specify both beam and density pruning");
    
    if (beam <= 0.0 && density <= 1.0)
      throw std::runtime_error("you may want to specify either beam or density pruning");

    if (int(semiring_tropical) + semiring_logprob + semiring_log == 0)
      semiring_tropical = true;
    
    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
  }

  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type pruned;

    weight_set_type __weights;
    if (weights_one) {
      __weights.allocate();
      for (weight_set_type::feature_type::id_type id = 0; id != __weights.size(); ++ id)
	if (! weight_set_type::feature_type(id).empty())
	  __weights[weight_set_type::feature_type(id)] = 1.0;
    }
    
    const weight_set_type* weights_prune = (weights ? weights : &__weights);

    if (debug)
      std::cerr << "pruning:"
		<< " # of nodes: " << hypergraph.nodes.size()
		<< " # of edges: " << hypergraph.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(hypergraph.is_valid())
		<< std::endl;
    
    utils::resource prune_start;

    if (beam > 0.0) {
      if (semiring_tropical)
	cicada::prune_beam(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), beam);
      else if (semiring_logprob)
	cicada::prune_beam(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), beam);
      else
	cicada::prune_beam(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), beam);
    } else if (density > 1.0) {
      if (semiring_tropical)
	cicada::prune_density(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), density);
      else if (semiring_logprob)
	cicada::prune_density(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), density);
      else
	cicada::prune_density(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), density);
    } else
      throw std::runtime_error("what pruning?");
    
	
    utils::resource prune_end;
    
    if (debug)
      std::cerr << "prune cpu time: " << (prune_end.cpu_time() - prune_start.cpu_time())
		<< " user time: " << (prune_end.user_time() - prune_start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "pruned:"
		<< " # of nodes: " << pruned.nodes.size()
		<< " # of edges: " << pruned.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(pruned.is_valid())
		<< std::endl;
    
    hypergraph.swap(pruned);

  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }

  const weight_set_type* weights;
  
  double beam;
  double density;
  double scale;
  
  bool weights_one;

  bool semiring_tropical;
  bool semiring_logprob;
  bool semiring_log;
  
  int debug;
};

class Intersect : public Operation
{
public:
  Intersect(const int __debug)
    : debug(__debug) {}

  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    if (targets.empty())
      throw std::runtime_error("no target?");
    
    lattice_type target(targets.front());

    hypergraph_type intersected;
    
    utils::resource start;
    
    cicada::intersect(hypergraph, target, intersected);
    
    utils::resource end;
    
    if (debug)
      std::cerr << "intersect cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
	
    if (debug)
      std::cerr << "# of nodes: " << intersected.nodes.size()
		<< " # of edges: " << intersected.edges.size()
		<< " valid? " << utils::lexical_cast<std::string>(intersected.is_valid())
		<< std::endl;
    
    hypergraph.swap(intersected);
  }
  
  int debug;
};

class OutputString : public Operation
{
public:
  OutputString(const std::string& parameter, std::string& __buffer, size_t& __id, const int __debug)
    : buffer(__buffer), id(__id), weights(0), weights_one(false),
      kbest_size(0), kbest_unique(false),
      yield_source(false), yield_target(false), yield_tree(false),
      graphviz(false),
      debug(__debug)
  {
    typedef cicada::Parameter param_type;

    param_type param(parameter);
    if (param.name() != "output")
      throw std::runtime_error("this is not a outputter");

    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "kbest") == 0)
	kbest_size = boost::lexical_cast<int>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "unique") == 0)
	kbest_unique = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	weights_one = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "graphviz") == 0)
	graphviz = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "file") == 0)
	;
      else if (strcasecmp(piter->first.c_str(), "directory") == 0)
	; 
      else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	const std::string& value = piter->second;
	
	if (strcasecmp(value.c_str(), "source") == 0)
	  yield_source = true;
	else if (strcasecmp(value.c_str(), "target") == 0)
	  yield_target = true;
	else if (strcasecmp(value.c_str(), "derivation") == 0 || strcasecmp(value.c_str(), "tree") == 0)
	  yield_tree = true;
	else
	  throw std::runtime_error("unknown yield: " + value);
      } else
	std::cerr << "WARNING: unsupported parameter for output: " << piter->first << "=" << piter->second << std::endl;
    }
    
    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
    
    if (int(yield_source) + yield_target + yield_tree > 1)
      throw std::runtime_error("only source, target or tree yield for kbest");
    
    if (! yield_source && ! yield_target && ! yield_tree)
      yield_target = true;
  }

  void clear()
  {
    buffer.clear();
  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }
  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    if (! hypergraph.is_valid()) return;
    
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::back_inserter(const_cast<std::string&>(buffer)));
    
    if (kbest_size <= 0) {
      if (graphviz)
	cicada::graphviz(os, hypergraph);
      else
	os << id << " ||| " << hypergraph << '\n';
    } else {
      weight_set_type weights_zero;
      const weight_set_type* weights_kbest = (weights ? weights : &weights_zero);

      if (weights_one) {
	if (kbest_unique) {
	  if (yield_source)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function_one(*weights_kbest), kbest_filter_unique(hypergraph));
	  else if (yield_target)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter_unique(hypergraph));
	  else
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	} else {
	  if (yield_source)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	  else if (yield_target)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	  else
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	}
      } else {
	if (kbest_unique) {
	  if (yield_source)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function(*weights_kbest), kbest_filter_unique(hypergraph));
	  else if (yield_target)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter_unique(hypergraph));
	  else
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_function(*weights_kbest), kbest_filter(hypergraph));
	} else {
	  if (yield_source)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function(*weights_kbest), kbest_filter(hypergraph));
	  else if (yield_target)
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter(hypergraph));
	  else
	    kbest_derivations(os, id, hypergraph, kbest_size, kbest_function(*weights_kbest), kbest_filter(hypergraph));
	}
      }
    }    
  }
  
  std::string& buffer;
  size_t& id;
  
  const weight_set_type* weights;
  bool weights_one;
  int  kbest_size;
  bool kbest_unique;

  bool yield_source;
  bool yield_target;
  bool yield_tree;

  bool graphviz;

  int debug;
};

class Output : public Operation
{
public:
  Output(const std::string& parameter, boost::shared_ptr<std::ostream>& __os, size_t& __id, const int __debug)
    : os(__os), id(__id), file(), directory(), weights(0), weights_one(false),
      kbest_size(0), kbest_unique(false),
      yield_source(false), yield_target(false), yield_tree(false),
      graphviz(false),
      debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "output")
      throw std::runtime_error("this is not a outputter");

    
    for (param_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "kbest") == 0)
	kbest_size = boost::lexical_cast<int>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "unique") == 0)
	kbest_unique = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights") == 0)
	weights = &base_type::weights(piter->second);
      else if (strcasecmp(piter->first.c_str(), "weights-one") == 0)
	weights_one = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "graphviz") == 0)
	graphviz = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "file") == 0)
	file = piter->second;
      else if (strcasecmp(piter->first.c_str(), "directory") == 0)
	directory = piter->second;
      else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	const std::string& value = piter->second;
	
	if (strcasecmp(value.c_str(), "source") == 0)
	  yield_source = true;
	else if (strcasecmp(value.c_str(), "target") == 0)
	  yield_target = true;
	else if (strcasecmp(value.c_str(), "derivation") == 0 || strcasecmp(value.c_str(), "tree") == 0)
	  yield_tree = true;
	else
	  throw std::runtime_error("unknown yield: " + value);
      } else
	std::cerr << "WARNING: unsupported parameter for output: " << piter->first << "=" << piter->second << std::endl;
    }

    
    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
    
    // default to stdout
    if (directory.empty() && file.empty())
      file = "-";
    
    if (! directory.empty() && ! file.empty())
      throw std::runtime_error("you cannot output both in directory and file");
    
    if (int(yield_source) + yield_target + yield_tree > 1)
      throw std::runtime_error("only source, target or tree yield for kbest");
    
    if (! yield_source && ! yield_target && ! yield_tree)
      yield_target = true;
  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }
  
  void clear()
  {
    if (! os) return;
    
    if (! directory.empty())
      os.reset();
    else
      *os << std::flush;
  }
  
  void operator()(const lattice_type& lattice, const span_set_type& spans, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    if (! hypergraph.is_valid()) return;

    if (! os) {
      const path_type path = (! file.empty() ? file  : directory / (boost::lexical_cast<std::string>(id) + ".gz"));
      
      const_cast<boost::shared_ptr<std::ostream>&>(os).reset(new utils::compress_ostream(path, 1024 * 1024));
    }
    
    if (kbest_size <= 0) {
      if (graphviz)
	cicada::graphviz(*os, hypergraph);
      else
	*os << id << " ||| " << hypergraph << '\n';
    } else {
      weight_set_type weights_zero;
      const weight_set_type* weights_kbest = (weights ? weights : &weights_zero);
      
      if (weights_one) {
	if (kbest_unique) {
	  if (yield_source)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function_one(*weights_kbest), kbest_filter_unique(hypergraph));
	  else if (yield_target)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter_unique(hypergraph));
	  else
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	} else {
	  if (yield_source)	  
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	  else if (yield_target)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	  else
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
	}
      } else {
	if (kbest_unique) {
	  if (yield_source)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function(*weights_kbest), kbest_filter_unique(hypergraph));
	  else if (yield_target)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter_unique(hypergraph));
	  else
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_function(*weights_kbest), kbest_filter(hypergraph));
	} else {
	  if (yield_source)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal_source(), kbest_function(*weights_kbest), kbest_filter(hypergraph));
	  else if (yield_target)
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter(hypergraph));
	  else
	    kbest_derivations(*os, id, hypergraph, kbest_size, kbest_function(*weights_kbest), kbest_filter(hypergraph));
	}
      }
    }    
  }
  
  boost::shared_ptr<std::ostream>& os;
  size_t& id;
  
  path_type file;
  path_type directory;
  
  const weight_set_type* weights;
  bool weights_one;
  int  kbest_size;
  bool kbest_unique;

  bool yield_source;
  bool yield_target;
  bool yield_tree;

  bool graphviz;
  
  int debug;
};

class OperationSet
{
public:

  typedef Operation operation_type;
  typedef boost::shared_ptr<operation_type> operation_ptr_type;
  typedef std::vector<operation_ptr_type, std::allocator<operation_ptr_type> > operation_ptr_set_type;
  
  template <typename Iterator>
  OperationSet(Iterator first, Iterator last,
	       const grammar_type& grammar,
	       const model_type& model,
	       const std::string& goal,
	       const std::string& non_terminal,
	       const bool insertion,
	       const bool deletion,
	       const bool __input_id,
	       const bool __input_lattice,
	       const bool __input_forest,
	       const bool __input_span,
	       const bool __input_bitext,
	       const bool __input_mpi,
	       const int debug)
    : id(size_t(-1)), 
      input_id(__input_id),
      input_lattice(__input_lattice),
      input_forest(__input_forest),
      input_span(__input_span),
      input_bitext(__input_bitext),
      input_mpi(__input_mpi)
  {
    typedef cicada::Parameter param_type;

    bool checked = false;
    
    for (/**/; first != last; ++ first) {
      param_type param(*first);
      
      if (param.name() == "binarize")
	operations.push_back(operation_ptr_type(new Binarize(*first, debug)));
      else if (param.name() == "permute")
	operations.push_back(operation_ptr_type(new Permute(*first, debug)));
      else if (param.name() == "compose-earley")
	operations.push_back(operation_ptr_type(new ComposeEarley(grammar, goal, non_terminal, insertion, deletion, debug)));
      else if (param.name() == "compose-cky")
	operations.push_back(operation_ptr_type(new ComposeCKY(grammar, goal, non_terminal, insertion, deletion, debug)));
      else if (param.name() == "apply")
	operations.push_back(operation_ptr_type(new Apply(*first, model, debug)));
      else if (param.name() == "bleu")
	operations.push_back(operation_ptr_type(new Bleu(*first, model, debug)));
      else if (param.name() == "variational")
	operations.push_back(operation_ptr_type(new Variational(*first, model, debug)));
      else if (param.name() == "prune")
	operations.push_back(operation_ptr_type(new Prune(*first, debug)));
      else if (param.name() == "intersect")
	operations.push_back(operation_ptr_type(new Intersect(debug)));
      else if (param.name() == "output") {
	// we do extra checking so that all the output directed to either the same directory or output-file
	boost::shared_ptr<Output> output(new Output(*first, os, id, debug));
	
	if (! checked) {
	  file      = output->file;
	  directory = output->directory;
	} else {
	  output->file      = file;
	  output->directory = directory;
	}
	
	checked = true;
	
	if (input_mpi && ! file.empty())
	  operations.push_back(operation_ptr_type(new OutputString(*first, buffer, id, debug)));
	else
	  operations.push_back(output);
	
      } else
	throw std::runtime_error("unsupport op: " + std::string(*first));
    }
  }
  
  static std::string lists()
  {
    static const char* desc = "\
binarize: perform binarization (monolingual tree)\n\
\tdirection=[left|right] binarization direction\n\
\tsize=binarization size\n\
permute: permute tree (monolingual tree only)\n\
\tfeature=[true|false] apply feature\n\
\tweights=file weight file for composed feature\n\
\tsize=permute size\n\
\tcollapse=[true|false] collapse sparse features\n\
\texclude=[a non-terminal] to prohibit permutation. You can supply multiple\n\
compose-earley: composition from tree with grammar\n\
compose-cky: composition from lattice (or sentence) with grammar\n\
apply: feature application\n\
\tsize=<cube size>\n\
\texact=[true|false] no pruning feature application\n\
\tforced=[true|false] forced feature application\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialized weight\n\
bleu: BLEU computation\n\
\tsize=<cube size>\n\
\texact=[true|false] no pruning feature application\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialized weight\n\
variational: variational decoding\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialized weight\n\
prune: pruning\n\
\tbeam=beam pruning threshold in threshold > 0.0\n\
\tdensity=density pruning threshold in threshold > 1.0\n\
\tscale=scaling for score\n\
\tsemiring=[tropical|logprob|log] semiring to perform score computation\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialzied weight\n\
intersect: compute intersection\n\
output: kbest or hypergraph output\n\
\tkbest=<kbest size> zero for hypergraph output (default)\n\
\tunique=[true|false] unique translation\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialize weight\n\
\tyield=[source|target|derivation|tree] yield for kbest\n\
\tgraphviz=[true|false] dump in graphviz format\n\
\tdirectory=directory for output\n\
\tfile=file for output\n\
";
    return desc;
  }

  void assign(const weight_set_type& weights)
  {
    operation_ptr_set_type::iterator oiter_end = operations.end();
    for (operation_ptr_set_type::iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
      (*oiter)->assign(weights);
  }
  
  void operator()(const std::string& line) {
    
    // clear...
    {
      operation_ptr_set_type::iterator oiter_end = operations.end();
      for (operation_ptr_set_type::iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
	(*oiter)->clear();
    }

    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    if (input_id) {
      if (! parse_id(id, iter, end))
	throw std::runtime_error("invalid id-prefixed format");
    } else
      ++ id;
    
    if (input_lattice) {
      if (! lattice.assign(iter, end))
	throw std::runtime_error("invalid lattive format");
    } else if (input_forest) {
      if (! hypergraph.assign(iter, end))
	throw std::runtime_error("invalid hypergraph format");
    } else {
      if (! sentence.assign(iter, end))
	throw std::runtime_error("invalid sentence format");
      
      lattice = lattice_type(sentence);
    }
    
    if (input_span) {
      spans.clear();
      
      if (! parse_separator(iter, end))
	throw std::runtime_error("invalid span format (separator)");
      
      if (! spans.assign(iter, end))
	throw std::runtime_error("invalid span format");
    }
    
    if (input_bitext) {
      targets.clear();
      
      while (parse_separator(iter, end)) {
	targets.push_back(sentence_type());
	
	if (! targets.back().assign(iter, end))
	  throw std::runtime_error("invalid sentence format");
      }
      
      if (targets.empty())
	throw std::runtime_error("no bitext?");
    }
    
    if (iter != end)
      throw std::runtime_error("invalid input format");
    
    if (lattice.empty() && ! hypergraph.is_valid()) {
      
    } else {
      operation_ptr_set_type::const_iterator oiter_end = operations.end();
      for (operation_ptr_set_type::const_iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
	(*oiter)->operator()(lattice, spans, targets, hypergraph);
    }
  }
  
  bool input_id;
  bool input_lattice;
  bool input_forest;
  bool input_span;
  bool input_bitext;
  bool input_mpi;
  
  // output related...
  boost::shared_ptr<std::ostream> os;
  std::string                     buffer;

  path_type file;
  path_type directory;
  
  size_t id;
  
  lattice_type      lattice;
  span_set_type     spans;
  sentence_set_type targets;
  hypergraph_type   hypergraph;

  sentence_type sentence;
  
  operation_ptr_set_type operations;
};
