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

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"

#include "cicada/kbest.hpp"
#include "cicada/parameter.hpp"

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

#include "cicada/feature/variational.hpp"
#include "cicada/feature/bleu.hpp"

#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/sgi_hash_map.hpp"

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
    
    rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
    for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
      if (titer->is_non_terminal()) {
	const int pos = titer->non_terminal_index() - 1;
	boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + pos)).begin(), boost::get<0>(*(first + pos)).end());
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


inline
bool true_false(const std::string& token)
{
  if (strcasecmp(token.c_str(), "true") == 0)
    return true;
  if (strcasecmp(token.c_str(), "yes") == 0)
    return true;
  if (atoi(token.c_str()) > 0)
    return true;
  return false;
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
  
  virtual void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& graph) const = 0;

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
    
    if (param.find("size") != param.end())
      size = boost::lexical_cast<int>(param.find("size")->second);

    if (param.find("direction") != param.end()) {
      const std::string& dir = param.find("direction")->second;
      
      if (strcasecmp(dir.c_str(), "left") == 0)
	left = true;
      else if (strcasecmp(dir.c_str(), "right") == 0)
	right = true;
      else
	throw std::runtime_error("unuspported direction: " + parameter);
    }

    if (! left && ! right)
      right == true;
  }
  
  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
		<< " valid? " << (binarized.is_valid() ? "true" : "false")
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
public:
  Permute(const std::string& parameter, const int __debug)
    : weights(0), size(0), feature(false), collapse(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "permute")
      throw std::runtime_error("this is not a permuter");
    
    if (param.find("size") != param.end())
      size = boost::lexical_cast<int>(param.find("size")->second);
    
    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);

    if (param.find("feature") != param.end())
      feature = true_false(param.find("feature")->second);

    if (param.find("collapse") != param.end())
      collapse = true_false(param.find("collapse")->second);

    if (collapse && ! weights)
      throw std::runtime_error("collapsing but no weights...");
  }
  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type permuted;
    
    if (debug)
      std::cerr << "permute" << std::endl;
    
    utils::resource start;
    
    if (collapse)
      cicada::permute(hypergraph, permuted, cicada::PermuteFeatureCollapsed<weight_set_type>(*weights), size);
    else if (feature)
      cicada::permute(hypergraph, permuted, cicada::PermuteFeature(), size);
    else
      cicada::permute(hypergraph, permuted, cicada::PermuteNoFeature(), size);
    
    utils::resource end;
	
    if (debug)
      std::cerr << "permute cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << permuted.nodes.size()
		<< " # of edges: " << permuted.edges.size()
		<< " valid? " << (permuted.is_valid() ? "true" : "false")
		<< std::endl;
    
    hypergraph.swap(permuted);
  }
  
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
  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
		<< " valid? " << (composed.is_valid() ? "true" : "false")
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
  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
		<< " valid? " << (composed.is_valid() ? "true" : "false")
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
    : model(__model), weights(0), size(200), weights_one(false), exact(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "apply")
      throw std::runtime_error("this is not a feature-functin applier");

    if (param.find("size") != param.end())
      size = boost::lexical_cast<int>(param.find("size")->second);

    if (param.find("exact") != param.end())
      exact = true_false(param.find("exact")->second);
    
    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);
    
    if (param.find("weights-one") != param.end())
      weights_one = true_false(param.find("weights-one")->second);

    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
  }

  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    hypergraph_type applied;

    if (debug)
      std::cerr << "apply features" << std::endl;
    
    utils::resource start;

    weight_set_type weights_zero;
    const weight_set_type* weights_apply = (weights ? weights : &weights_zero);
    
    if (exact || model.is_stateless()) {
      if (weights_one)
	cicada::apply_exact<weight_set_function_one>(model, hypergraph, applied, weight_set_function_one(*weights_apply), size);
      else
	cicada::apply_exact<weight_set_function>(model, hypergraph, applied, weight_set_function(*weights_apply), size);
    } else {
      if (weights_one)
	cicada::apply_cube_prune<weight_set_function_one>(model, hypergraph, applied, weight_set_function_one(*weights_apply), size);
      else
	cicada::apply_cube_prune<weight_set_function>(model, hypergraph, applied, weight_set_function(*weights_apply), size);
    }
    
    utils::resource end;
    
    if (debug)
      std::cerr << "apply cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;

    if (debug)
      std::cerr << "# of nodes: " << applied.nodes.size()
		<< " # of edges: " << applied.edges.size()
		<< " valid? " << (applied.is_valid() ? "true" : "false")
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
    
    if (param.find("size") != param.end())
      size = boost::lexical_cast<int>(param.find("size")->second);

    if (param.find("exact") != param.end())
      exact = true_false(param.find("exact")->second);

    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);
    
    if (param.find("weights-one") != param.end())
      weights_one = true_false(param.find("weights-one")->second);

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

  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
	cicada::apply_exact<weight_set_function_one>(model, hypergraph, applied, weight_set_function_one(*weights_apply), size);
      else
	cicada::apply_exact<weight_set_function>(model, hypergraph, applied, weight_set_function(*weights_apply), size);
    } else {
      if (weights_one)
	cicada::apply_cube_prune<weight_set_function_one>(model, hypergraph, applied, weight_set_function_one(*weights_apply), size);
      else
	cicada::apply_cube_prune<weight_set_function>(model, hypergraph, applied, weight_set_function(*weights_apply), size);
    }
    
    utils::resource end;
	
    if (debug)
      std::cerr << "bleu cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << applied.nodes.size()
		<< " # of edges: " << applied.edges.size()
		<< " valid? " << (applied.is_valid() ? "true" : "false")
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
    : weights(0), weights_variational(0), size(200), weights_one(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "variational")
      throw std::runtime_error("this is not a variational decoder");
    
    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);

    if (param.find("weights-variational") != param.end())
      weights_variational = &base_type::weights(param.find("weights-variational")->second);
    
    if (param.find("weights-one") != param.end())
      weights_one = true_false(param.find("weights-one")->second);

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
    
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
    
    const weight_set_type* weights_orig(weights ? weights : &__weights);
    
    if (debug)
      std::cerr << "variational decoding" << std::endl;

    cicada::feature::Variational* __variational = dynamic_cast<cicada::feature::Variational*>(feature.get());
    
    utils::resource start;
    
    __variational->insert(hypergraph, *weights_orig);
    
    model_type model;
    model.push_back(feature);

    const weight_set_type* weights_apply(weights_variational ? weights_variational : weights_orig);
    
    // second, apply again...
    
    cicada::apply_exact<weight_set_function>(model, hypergraph, variational, weight_set_function(*weights_apply), size);
	
    utils::resource end;

    if (debug)
      std::cerr << "variational cpu time: " << (end.cpu_time() - start.cpu_time())
		<< " user time: " << (end.user_time() - start.user_time())
		<< std::endl;
	
    if (debug)
      std::cerr << "# of nodes: " << variational.nodes.size()
		<< " # of edges: " << variational.edges.size()
		<< " valid? " << (variational.is_valid() ? "true" : "false")
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
  const weight_set_type* weights_variational;
  int size;
  bool weights_one;
  
  int debug;
};

class Prune : public Operation
{
public:
  Prune(const std::string& parameter, const int __debug)
    : weights(0), beam(0.0), scale(1.0), weights_one(false), 
      semiring_tropical(false), semiring_logprob(false), semiring_log(false),
      debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "prune")
      throw std::runtime_error("this is not a pruner");
    
    if (param.find("beam") != param.end())
      beam = boost::lexical_cast<double>(param.find("beam")->second);
    
    if (beam <= 0.0)
      throw std::runtime_error("beam threshold must be 0.0 < threshold");

    if (param.find("scale") != param.end())
      scale = boost::lexical_cast<double>(param.find("scale")->second);
    
    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);
    
    if (param.find("weights-one") != param.end())
      weights_one = true_false(param.find("weights-one")->second);
    
    if (param.find("semiring") != param.end()) {
      const std::string& name = param.find("semiring")->second;

      if (strcasecmp(name.c_str(), "tropical") == 0)
	semiring_tropical = true;
      else if (strcasecmp(name.c_str(), "logprob") == 0)
	semiring_logprob = true;
      else if (strcasecmp(name.c_str(), "log") == 0)
	semiring_log = true;
      else
	throw std::runtime_error("unknown semiring: " + name);
    }

    
    if (int(semiring_tropical) + semiring_logprob + semiring_log == 0)
      semiring_tropical = true;
    

    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
  }

  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
    
    utils::resource prune_start;
    
    if (semiring_tropical)
      cicada::beam_prune(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Tropical<double> >(*weights_prune, scale), beam);
    else if (semiring_logprob)
      cicada::beam_prune(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Logprob<double> >(*weights_prune, scale), beam);
    else
      cicada::beam_prune(hypergraph, pruned, weight_set_scaled_function<cicada::semiring::Log<double> >(*weights_prune, scale), beam);
    
	
    utils::resource prune_end;
    
    if (debug)
      std::cerr << "prune cpu time: " << (prune_end.cpu_time() - prune_start.cpu_time())
		<< " user time: " << (prune_end.user_time() - prune_start.user_time())
		<< std::endl;
    
    if (debug)
      std::cerr << "# of nodes: " << pruned.nodes.size()
		<< " # of edges: " << pruned.edges.size()
		<< " valid? " << (pruned.is_valid() ? "true" : "false")
		<< std::endl;
    
    hypergraph.swap(pruned);

  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }

  const weight_set_type* weights;
  
  double beam;
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

  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
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
		<< " valid? " << (intersected.is_valid() ? "true" : "false")
		<< std::endl;
    
    hypergraph.swap(intersected);
  }
  
  int debug;
};

class OutputString : public Operation
{
public:
  OutputString(const std::string& parameter, std::string& __buffer, size_t& __id, const int __debug)
    : buffer(__buffer), id(__id), weights(0), weights_one(false), kbest_size(0), kbest_unique(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;

    param_type param(parameter);
    if (param.name() != "output")
      throw std::runtime_error("this is not a outputter");

    if (param.find("kbest") != param.end())
      kbest_size = boost::lexical_cast<int>(param.find("kbest")->second);
    
    if (param.find("unique") != param.end())
      kbest_unique = true_false(param.find("unique")->second);
    
    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);
   
    if (param.find("weights-one") != param.end())
      weights_one = true_false(param.find("weights-one")->second);

    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
  }

  void clear()
  {
    buffer.clear();
  }

  void assign(const weight_set_type& __weights)
  {
    weights = &__weights;
  }
  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    if (! hypergraph.is_valid()) return;
    
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::back_inserter(const_cast<std::string&>(buffer)));
    
    if (kbest_size <= 0)
      os << id << " ||| " << hypergraph << '\n';
    else {
      weight_set_type weights_zero;
      const weight_set_type* weights_kbest = (weights ? weights : &weights_zero);
      
      if (weights_one) {
	if (kbest_unique)
	  kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter_unique(hypergraph));
	else
	  kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
      } else {
	if (kbest_unique)
	  kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter_unique(hypergraph));
	else
	  kbest_derivations(os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter(hypergraph));
      }
    }    
  }
  
  std::string& buffer;
  size_t& id;
  
  const weight_set_type* weights;
  bool weights_one;
  int  kbest_size;
  bool kbest_unique;

  int debug;
};

class Output : public Operation
{
public:
  Output(const std::string& parameter, boost::shared_ptr<std::ostream>& __os, size_t& __id, const int __debug)
    : os(__os), id(__id), file(), directory(), weights(0), weights_one(false), kbest_size(0), kbest_unique(false), debug(__debug)
  {
    typedef cicada::Parameter param_type;
    
    param_type param(parameter);
    if (param.name() != "output")
      throw std::runtime_error("this is not a outputter");

    if (param.find("kbest") != param.end())
      kbest_size = boost::lexical_cast<int>(param.find("kbest")->second);
    
    if (param.find("unique") != param.end())
      kbest_unique = true_false(param.find("unique")->second);
    
    if (param.find("weights") != param.end())
      weights = &base_type::weights(param.find("weights")->second);
   
    if (param.find("weights-one") != param.end())
      weights_one = true_false(param.find("weights-one")->second);

    if (weights && weights_one)
      throw std::runtime_error("you have weights, but specified all-one parameter");
    
    if (param.find("directory") != param.end())
      directory = param.find("directory")->second;

    if (param.find("file") != param.end())
      file = param.find("file")->second;
    
    // default to stdout
    if (directory.empty() && file.empty())
      file = "-";
    
    if (! directory.empty() && ! file.empty())
      throw std::runtime_error("you cannot output both in directory and file");
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
  
  void operator()(const lattice_type& lattice, const sentence_set_type& targets, hypergraph_type& hypergraph) const
  {
    if (! hypergraph.is_valid()) return;

    if (! os) {
      const path_type path = (! file.empty() ? file  : directory / (boost::lexical_cast<std::string>(id) + ".gz"));
      
      const_cast<boost::shared_ptr<std::ostream>&>(os).reset(new utils::compress_ostream(path, 1024 * 1024));
    }
    
    if (kbest_size <= 0)
      *os << id << " ||| " << hypergraph << '\n';
    else {
      weight_set_type weights_zero;
      const weight_set_type* weights_kbest = (weights ? weights : &weights_zero);
      
      if (weights_one) {
	if (kbest_unique)
	  kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter_unique(hypergraph));
	else
	  kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function_one(*weights_kbest), kbest_filter(hypergraph));
      } else {
	if (kbest_unique)
	  kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter_unique(hypergraph));
	else
	  kbest_derivations(*os, id, hypergraph, kbest_size, kbest_traversal(), kbest_function(*weights_kbest), kbest_filter(hypergraph));
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
	       const bool __input_bitext,
	       const bool __input_mpi,
	       const int debug)
    : id(0), 
      input_id(__input_id),
      input_lattice(__input_lattice),
      input_forest(__input_forest),
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
\tdirection=[left|right] binarization direction,\n\
\tsize=binarization size\n\
permute: permute tree (monolingual tree only)\n\
\tfeature=[true|false] apply feature,\n\
\tweights=file weight file for composed feature,\n\
\tsize=permute size\n\
\tcollapse=[true|false] collapse sparse features\n\
compose-earley: composition from tree with grammar\n\
compose-cky: composition from lattice (or sentence) with grammar\n\
apply: feature application\n\
\tsize=<cube size>,\n\
\texact=[true|false] no pruning feature application,\n\
\tweights=weight file for feature,\n\
\tweights-one=[true|false] one initialized weight\n\
bleu: BLEU computation\n\
\tsize=<cube size>\n\
\texact=[true|false] no pruning feature application,\n\
\tweights=weight file for feature,\n\
\tweights-one=[true|false] one initialized weight\n\
variational: variational decoding\n\
\tweights=weight file for feature,\n\
\tweights-variational=weighs for variational decoding feature,\n\
\tweights-one=[true|false] one initialized weight\n\
prune: beam pruning\n\
\tbeam=beam threshold in 0.0 < threshold,\n\
\tscale=scaling for score,\n\
\tsemiring=[tropical|logprob|log] semiring to perform score computation,\n\
\tweights=weight file for feature,\n\
\tweights-one=[true|false] one initialzied weight\n\
\tintersect\n\
output: kbest or hypergraph output\n\
\tkbest=<kbest size> zero for hypergraph output (default),\n\
\tunique=[true|false] unique translation,\n\
\tweights=weight file for feature,\n\
\tweights-one=[true|false] one initialize weight,\n\
\tdirectory=directory for output,\n\
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
      
      lattice = lattice_type(sentence);
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
    
    if (lattice.empty() && ! hypergraph.is_valid())
      ++ id;
    else {
      operation_ptr_set_type::const_iterator oiter_end = operations.end();
      for (operation_ptr_set_type::const_iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
	(*oiter)->operator()(lattice, targets, hypergraph);
      
      ++ id;
    }
  }
  
  bool input_id;
  bool input_lattice;
  bool input_forest;
  bool input_bitext;
  bool input_mpi;
  
  // output related...
  boost::shared_ptr<std::ostream> os;
  std::string                     buffer;

  path_type file;
  path_type directory;
  
  size_t id;
  
  lattice_type      lattice;
  sentence_set_type targets;
  hypergraph_type   hypergraph;

  sentence_type sentence;
  
  operation_ptr_set_type operations;
};
