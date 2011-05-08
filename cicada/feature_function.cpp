//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "feature_function.hpp"

#include "parameter.hpp"

#include "feature/alignment.hpp"
#include "feature/antecedent.hpp"
#include "feature/bleu.hpp"
#include "feature/bleu_expected.hpp"
#include "feature/bleu_linear.hpp"
#include "feature/bleu_multi.hpp"
#include "feature/distortion.hpp"
#include "feature/global_lexicon.hpp"
#include "feature/lexicalized_reordering.hpp"
#include "feature/neighbours.hpp"
#include "feature/ngram.hpp"
#include "feature/ngram_tree.hpp"
#include "feature/parent.hpp"
#include "feature/penalty.hpp"
#include "feature/permute.hpp"
#include "feature/rule_shape.hpp"
#include "feature/span.hpp"
#include "feature/variational.hpp"

#include "utils/piece.hpp"

namespace cicada
{
  
  const char* FeatureFunction::lists()
  {
    static const char* desc = "\
antecedent: antecedent feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
\talignment=[true|false] alignment forest mode\n\
\tsource-root=[true|false] source root mode (for tree composition)\n\
\tname=feature-name-prefix(default: antecedent)\n\
bleu: BLEU\n\
\torder=<order>\n\
\texact=[true|false] clipped ngram computation\n\
\ttokenizer=[tokenizer spec]\n\
\tname=feature-name(default: bleu)\n\
\trefset=reference set file\n\
bleu-expected: expected-BLEU\n\
\torder=<order>\n\
bleu-linear: linear corpus-BLEU\n\
\torder=<order>\n\
\tprecision=<default 0.8>\n\
\tratio=<default 0.6>\n\
\ttokenizer=[tokenizer spec]\n\
\tname=feature-name(default: bleu-linear)\n\
\trefset=reference set file\n\
bleu-multi: multiple BLEU\n\
\torder=<order>\n\
\texact=[true|false] clipped ngram computation\n\
\ttokenizer=[tokenizer spec]\n\
\tname=feature-name(default: bleu-multi:integer)\n\
\tsize=# of BLEU features\n\
distortion: phrase-based distortion\n\
global-lexicon: global lexicon feature\n\
\tfile=global lexicon file\n\
lexicalized-reordering: lexicalized reordering for phrase-based\n\
\tbidirectional=[true|false]\n\
\tmonotonicity=[true|false]\n\
\tfeature=attribute name mapping\n\
neighbours: neighbour words feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
\talignment=[true|false] alignment forest mode\n\
\tsource-root=[true|false] source root mode (for tree composition)\n\
\tname=feature-name-prefix(default: neighbours)\n\
ngram: ngram language model\n\
\tfile=<file>\n\
\torder=<order>\n\
\tcluster=<word class>\n\
\tname=feature-name(default: ngram)\n\
\tno-bos-eos=[true|false] do not add bos/eos\n\
\tcoarse-order=<order> ngram order for coarse heuristic\n\
\tcoarse-file=<file>   ngram for coarrse heuristic\n\
\tcoarse-cluster=<word class> word class for coarse heuristics\n\
ngram-tree: ngram tree feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
\talignment=[true|false] alignment forest mode\n\
\tsource-root=[true|false] source root mode (for tree composition)\n\
\tname=feature-name-prefix(default: ngram-tree)\n\
parent: parent feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
\texclude-terminal=[true|false] exclude terminal symbol\n\
\talignment=[true|false] alignment forest mode\n\
\tsource-root=[true|false] source root mode (for tree composition)\n\
\tname=feature-name-prefix(default: parent)\n\
permute: permutation feature\n\
\tweights=weight file for collapsed feature\n\
\tcollapse=[true|false] collapsed feature\n\
span: lexical span feature\n\
variational: variational feature for variational decoding\n\
\torder=<order>\n\
\tno-bos-eos=[true|false] do not add bos/eos\n\
word-penalty: word penalty feature\n\
rule-penalty: rule penalty feature\n\
arity-penalty: rule arity penalty feature\n\
glue-tree-penalty: glue tree penalty feature\n\
non-latin-penalty: non-latin word penalty feature\n\
rule-shape: rule shape feature\n\
relative-position: relative alignment feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
null-path: path involving epsilon\n\
target-bigram: target bigram feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
word-pair: word pair feature\n\
\tcluster-source=[word class file]\n\
\tcluster-target=[word class file]\n\
\tstemmer-source=[stemmer spec]\n\
\tstemmer-target=[stemmer spec]\n\
";
    return desc;
  }

  FeatureFunction::feature_function_ptr_type FeatureFunction::create(const utils::piece& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const utils::ipiece param_name(param.name());
    
    if (param_name == "ngram")
      return feature_function_ptr_type(new feature::NGram(parameter));
    else if (param_name == "neighbours" || param_name == "neighbors")
      return feature_function_ptr_type(new feature::Neighbours(parameter));
    else if (param_name == "ngram-tree")
      return feature_function_ptr_type(new feature::NGramTree(parameter));
    else if (param_name == "antecedent")
      return feature_function_ptr_type(new feature::Antecedent(parameter));
    else if (param_name == "parent")
      return feature_function_ptr_type(new feature::Parent(parameter));
    else if (param_name == "permute")
      return feature_function_ptr_type(new feature::Permute(parameter));
    else if (param_name == "bleu")
      return feature_function_ptr_type(new feature::Bleu(parameter));
    else if (param_name == "bleu-expected")
      return feature_function_ptr_type(new feature::BleuExpected(parameter));
    else if (param_name == "bleu-linear")
      return feature_function_ptr_type(new feature::BleuLinear(parameter));
    else if (param_name == "bleu-multi" || param_name == "bleu-multiple")
      return feature_function_ptr_type(new feature::BleuMulti(parameter));
    else if (param_name == "distortion")
      return feature_function_ptr_type(new feature::Distortion(parameter));
    else if (param_name == "global-lexicon")
      return feature_function_ptr_type(new feature::GlobalLexicon(parameter));
    else if (param_name == "lexicalized-reordering"
	     || param_name == "lexicalized-reorder"
	     || param_name == "lexical-reordering"
	     || param_name == "lexical-reorder")
      return feature_function_ptr_type(new feature::LexicalizedReordering(parameter));
    else if (param_name == "span")
      return feature_function_ptr_type(new feature::Span(parameter));
    else if (param_name == "variational")
      return feature_function_ptr_type(new feature::Variational(parameter));
    else if (param_name == "word-penalty")
      return feature_function_ptr_type(new feature::WordPenalty());
    else if (param_name == "rule-penalty")
      return feature_function_ptr_type(new feature::RulePenalty());
    else if (param_name == "arity-penalty")
      return feature_function_ptr_type(new feature::ArityPenalty());
    else if (param_name == "glue-tree-penalty")
      return feature_function_ptr_type(new feature::GlueTreePenalty());
    else if (param_name == "non-latin-penalty")
      return feature_function_ptr_type(new feature::NonLatinPenalty());
    else if (param_name == "relative-position")
      return feature_function_ptr_type(new feature::RelativePosition(parameter));
    else if (param_name == "null-path")
      return feature_function_ptr_type(new feature::NullPath(parameter));
    else if (param_name == "target-bigram")
      return feature_function_ptr_type(new feature::TargetBigram(parameter));
    else if (param_name == "word-pair")
      return feature_function_ptr_type(new feature::WordPair(parameter));
    else if (param_name == "rule-shape")
      return feature_function_ptr_type(new feature::RuleShape(parameter));
    else
      throw std::runtime_error("unknown featuer: " + parameter);
    
  }
  
};
