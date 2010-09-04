#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include "operation_set.hpp"

#include "utils/hashmurmur.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/compress_stream.hpp"

#include <boost/filesystem.hpp>

namespace cicada
{
  std::string OperationSet::lists()
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
generate-earley: re-generation from tree\n\
apply: feature application\n\
\tsize=<cube size>\n\
\texact=[true|false]  no pruning feature application\n\
\tprune=[true|false]  cube-pruning for feature application\n\
\tgrow=[true|false]   cube-growing for feature application\n\
\tforced=[true|false] forced feature application\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialized weight\n\
\tfeature=feature function\n\
expected-ngram: expected ngram computation\n\
\torder=<ngram order>\n\
\tbos-eos=[true|false] include <s> and </s> in ngrams\n\
\tweights=weight file for feature\n\
\tweights-one=[true|false] one initialized weight\n\
\tyield=[source|target]\n\
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


  OperationSet::OperationSet(const parameter_set_type& parameters,
			     const model_type& model,
			     const grammar_type& grammar,
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
    : input_id(__input_id),
      input_lattice(__input_lattice),
      input_forest(__input_forest),
      input_span(__input_span),
      input_bitext(__input_bitext),
      input_mpi(__input_mpi)
  {
    // operation object creation...

    
    
  }
  
  void OperationSet::assign(const weight_set_type& weights)
  {
    operation_ptr_set_type::iterator oiter_end = operations.end();
    for (operation_ptr_set_type::iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
      (*oiter)->assign(weights);
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
  
  void OperationSet::operator()(const std::string& line)
  {
    // clear...
    {
      operation_ptr_set_type::iterator oiter_end = operations.end();
      for (operation_ptr_set_type::iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
	(*oiter)->clear();
    }
    
    // perform parsing...
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    
    if (input_id) {
      if (! parse_id(data.id, iter, end))
	throw std::runtime_error("invalid id-prefixed format");
    } else
      ++ data.id;
    
    if (input_lattice) {
      if (! data.lattice.assign(iter, end))
	throw std::runtime_error("invalid lattive format");
    } else if (input_forest) {
      if (! data.hypergraph.assign(iter, end))
	throw std::runtime_error("invalid hypergraph format");
    } else {
      sentence_type sentence;
      if (! sentence.assign(iter, end))
	throw std::runtime_error("invalid sentence format");
      
      data.lattice = lattice_type(sentence);
    }
    
    if (input_span) {
      data.spans.clear();
      
      if (! parse_separator(iter, end))
	throw std::runtime_error("invalid span format (separator)");
      
      if (! data.spans.assign(iter, end))
	throw std::runtime_error("invalid span format");
    }
    
    if (input_bitext) {
      data.targets.clear();

      if (! parse_separator(iter, end))
	throw std::runtime_error("invalid bitext format (separator)");

      if (! data.targets.assign(iter, end))
	throw std::runtime_error("invalid sentence set format");
    }
    
    if (iter != end)
      throw std::runtime_error("invalid input format");
    
    // processing...
    
    if (data.lattice.empty() && ! data.hypergraph.is_valid()) {
      
    } else {
      operation_ptr_set_type::const_iterator oiter_end = operations.end();
      for (operation_ptr_set_type::const_iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
	(*oiter)->operator()(data);
    }
  }
};
