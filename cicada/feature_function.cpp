
#include "feature_function.hpp"

#include "parameter.hpp"

#include "feature/antecedent.hpp"
#include "feature/bleu.hpp"
#include "feature/bleu_expected.hpp"
#include "feature/bleu_linear.hpp"
#include "feature/boundary.hpp"
#include "feature/global_lexicon.hpp"
#include "feature/neighbours.hpp"
#include "feature/ngram.hpp"
#include "feature/ngram_tree.hpp"
#include "feature/parent.hpp"
#include "feature/penalty.hpp"
#include "feature/span.hpp"
#include "feature/variational.hpp"


namespace cicada
{
  
  std::string FeatureFunction::lists()
  {
    static const char* desc = "\
antecedent: antecedent feature\n\
\tyield=[source|target]\n\
\tcluster=[word class file]\n\
\tprefix=[prefix stemming size]\n\
\tsuffix=[suffix stemming size]\n\
\tdigits=[perform digits stemming]\n\
bleu: BLEU\n\
\torder=<order>,\n\
\texact=[true|false] clipped ngram computation\n\
\tsplit=[true|false] split non-ascii chars\n\
\tlower=[true|false] perform lower casing\n\
\tname=feature-name(default: bleu)\n\
\trefset=reference set file\n\
\tyield=[source|target]\n\
bleu-expected: expected-BLEU\n\
\torder=<order>,\n\
\tyield=[source|target]\n\
bleu-linear: linear corpus-BLEU\n\
\torder=<order>,\n\
\tprecision=<default 0.8>,\n\
\tratio=<default 0.6>\n\
\tsplit=[true|false] split non-ascii chars\n\
\tlower=[true|false] perform lower casing\n\
\tname=feature-name(default: bleu-linear)\n\
\trefset=reference set file\n\
\tyield=[source|target]\n\
boundary: boundary bigram feature\n\
\tcluster-source=[word class file for source]\n\
\tcluster-target=[word class file for target]\n\
\tprefix-source=[source prefix stemming size]\n\
\tprefix-target=[target prefix stemming size]\n\
\tsuffix-source=[source suffix stemming size]\n\
\tsuffix-target=[target suffix stemming size]\n\
\tdigits-source=[source digits stemming]\n\
\tdigits-target=[target digits stemming]\n\
global-lexicon: global lexicon feature\n\
\tyield=[source|target]\n\
neighbours: neighbour words feature\n\
\tyield=[source|target]\n\
\tcluster=[word class file]\n\
\tprefix=[prefix stemming size]\n\
\tsuffix=[suffix stemming size]\n\
\tdigits=[perform digits stemming]\n\
ngram: ngram language model\n\
\tfile=<file>\n\
\torder=<order>\n\
\tcluster=<word class>\n\
\tyiled=[source|target yield (default target)]\n\
\tname=feature-name(default: ngram)\n\
\tcoarse-order=<order> ngram order for coarse heuristic\n\
\tcoarse-file=<file>   ngram for coarrse heuristic\n\
\tcoarse-cluster=<word class> word class for coarse heuristics\n\
ngram-tree: ngram tree feature\n\
\tyield=[source|target]\n\
\tcluster=[word class file]\n\
\tprefix=[prefix stemming size]\n\
\tsuffix=[suffix stemming size]\n\
\tdigits=[perform digits stemming]\n\
parent: parent feature\n\
\tyield=[source|target]\n\
\tcluster=[word class file]\n\
\tprefix=[prefix stemming size]\n\
\tsuffix=[suffix stemming size]\n\
\tdigits=[perform digits stemming]\n\
\texclude-terminal=[true|false] exclude terminal symbol\n\
span: lexical span feature\n\
variational: variational feature for variational decoding\n\
\torder=<order>\n\
target-word-penalty: target word penalty feature\n\
soruce-word-penalty: source word penalty feature\n\
rule-penalty: rule penalty feature\n\
";
    return desc;
  }

  FeatureFunction::feature_function_ptr_type FeatureFunction::create(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    
    if (param.name() == "ngram")
      return feature_function_ptr_type(new feature::NGram(parameter));
    else if (param.name() == "neighbours" || param.name() == "neighbors")
      return feature_function_ptr_type(new feature::Neighbours(parameter));
    else if (param.name() == "ngram-tree")
      return feature_function_ptr_type(new feature::NGramTree(parameter));
    else if (param.name() == "antecedent")
      return feature_function_ptr_type(new feature::Antecedent(parameter));
    else if (param.name() == "parent")
      return feature_function_ptr_type(new feature::Parent(parameter));
    else if (param.name() == "boundary")
      return feature_function_ptr_type(new feature::Boundary(parameter));
    else if (param.name() == "bleu")
      return feature_function_ptr_type(new feature::Bleu(parameter));
    else if (param.name() == "bleu-expected")
      return feature_function_ptr_type(new feature::BleuExpected(parameter));
    else if (param.name() == "bleu-linear")
      return feature_function_ptr_type(new feature::BleuLinear(parameter));
    else if (param.name() == "global-lexicon")
      return feature_function_ptr_type(new feature::GlobalLexicon(parameter));
    else if (param.name() == "span")
      return feature_function_ptr_type(new feature::Span(parameter));
    else if (param.name() == "variational")
      return feature_function_ptr_type(new feature::Variational(parameter));
    else if (param.name() == "target-word-penalty")
      return feature_function_ptr_type(new feature::TargetWordPenalty());
    else if (param.name() == "source-word-penalty")
      return feature_function_ptr_type(new feature::SourceWordPenalty());
    else if (param.name() == "rule-penalty")
      return feature_function_ptr_type(new feature::RulePenalty());
    else
      throw std::runtime_error("unknown featuer: " + parameter);
    
  }
  
};
