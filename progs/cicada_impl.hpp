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

#include "cicada/model.hpp"
#include "cicada/grammar.hpp"

#include "cicada/operation_set.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"
#include "cicada/span_vector.hpp"
#include "cicada/ngram_count_set.hpp"

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


typedef cicada::OperationSet operation_set_type;

typedef feature_function_type::feature_function_ptr_type feature_function_ptr_type;

typedef rule_type::feature_set_type    feature_set_type;
typedef feature_set_type::feature_type feature_type;
typedef cicada::WeightVector<double>   weight_set_type;

typedef cicada::SpanVector span_set_type;
typedef cicada::NGramCountSet ngram_count_set_type;
typedef cicada::SentenceVector sentence_set_type;

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

  value_type operator()(const hypergraph_type::edge_type& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.features.dot(weights));
  }
  
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

  value_type operator()(const hypergraph_type::edge_type& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.features.dot());
  }
  
  template <typename FeatureSet>
  value_type operator()(const FeatureSet& x) const
  {
    return cicada::semiring::traits<value_type>::log(x.dot());
  }
};



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

