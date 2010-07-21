
#include "feature_function.hpp"

#include "parameter.hpp"

#include "feature/antecedent.hpp"
#include "feature/bleu.hpp"
#include "feature/bleu_linear.hpp"
#include "feature/boundary.hpp"
#include "feature/neighbours.hpp"
#include "feature/ngram.hpp"
#include "feature/ngram_tree.hpp"
#include "feature/penalty.hpp"
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
bleu: BLEU\n\
\torder=<order>,\n\
\texact=[true|false] clipped ngram computation\n\
bleu-linear: linear corpus-BLEU\n\
\torder=<order>,\n\
\tprecision=<default 0.8>,\n\
\tratio=<default 0.6>\n\
boundary: boundary bigram feature\n\
\tcluster-source=[word class file for source]\n\
\tcluster-target=[word class file for target]\n\
neighbours: neighbour words feature\n\
\tyield=[source|target]\n\
\tcluster=[word class file]\n\
\tprefix=[prefix stemming size]\n\
\tsuffix=[suffix stemming size]\n\
ngram: ngram language model\n\
\tfile=<file>,\n\
\torder=<order>,\n\
\tname=feature-name(default: ngram)\n\
ngram-tree: ngram tree feature\n\
\tyield=[source|target]\n\
\tcluster=[word class file]\n\
\tprefix=[prefix stemming size]\n\
\tsuffix=[suffix stemming size]\n\
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
    else if (param.name() == "boundary")
      return feature_function_ptr_type(new feature::Boundary(parameter));
    else if (param.name() == "bleu")
      return feature_function_ptr_type(new feature::Bleu(parameter));
    else if (param.name() == "bleu-linear")
      return feature_function_ptr_type(new feature::BleuLinear(parameter));
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
