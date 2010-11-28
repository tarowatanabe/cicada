//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "feature_function.hpp"

#include "parameter.hpp"

#include "feature/antecedent.hpp"
#include "feature/bleu.hpp"
#include "feature/bleu_expected.hpp"
#include "feature/bleu_linear.hpp"
#include "feature/distortion.hpp"
#include "feature/global_lexicon.hpp"
#include "feature/lexicalized_reordering.hpp"
#include "feature/neighbours.hpp"
#include "feature/ngram.hpp"
#include "feature/ngram_tree.hpp"
#include "feature/parent.hpp"
#include "feature/penalty.hpp"
#include "feature/span.hpp"
#include "feature/variational.hpp"


namespace cicada
{
  
  const char* FeatureFunction::lists()
  {
    static const char* desc = "\
antecedent: antecedent feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
bleu: BLEU\n\
\torder=<order>\n\
\texact=[true|false] clipped ngram computation\n\
\tokenizer=[tokenizer spec]\n\
\tname=feature-name(default: bleu)\n\
\trefset=reference set file\n\
bleu-expected: expected-BLEU\n\
\torder=<order>\n\
bleu-linear: linear corpus-BLEU\n\
\torder=<order>\n\
\tprecision=<default 0.8>\n\
\tratio=<default 0.6>\n\
\tokenizer=[tokenizer spec]\n\
\tname=feature-name(default: bleu-linear)\n\
\trefset=reference set file\n\
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
ngram: ngram language model\n\
\tfile=<file>\n\
\torder=<order>\n\
\tcluster=<word class>\n\
\tname=feature-name(default: ngram)\n\
\tcoarse-order=<order> ngram order for coarse heuristic\n\
\tcoarse-file=<file>   ngram for coarrse heuristic\n\
\tcoarse-cluster=<word class> word class for coarse heuristics\n\
ngram-tree: ngram tree feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
parent: parent feature\n\
\tcluster=[word class file]\n\
\tstemmer=[stemmer spec]\n\
\texclude-terminal=[true|false] exclude terminal symbol\n\
span: lexical span feature\n\
variational: variational feature for variational decoding\n\
\torder=<order>\n\
word-penalty: word penalty feature\n\
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
    else if (param.name() == "bleu")
      return feature_function_ptr_type(new feature::Bleu(parameter));
    else if (param.name() == "bleu-expected")
      return feature_function_ptr_type(new feature::BleuExpected(parameter));
    else if (param.name() == "bleu-linear")
      return feature_function_ptr_type(new feature::BleuLinear(parameter));
    else if (param.name() == "distortion")
      return feature_function_ptr_type(new feature::Distortion(parameter));
    else if (param.name() == "global-lexicon")
      return feature_function_ptr_type(new feature::GlobalLexicon(parameter));
    else if (param.name() == "lexicalized-reordering"
	     || param.name() == "lexicalized-reorder"
	     || param.name() == "lexical-reordering"
	     || param.name() == "lexical-reorder")
      return feature_function_ptr_type(new feature::LexicalizedReordering(parameter));
    else if (param.name() == "span")
      return feature_function_ptr_type(new feature::Span(parameter));
    else if (param.name() == "variational")
      return feature_function_ptr_type(new feature::Variational(parameter));
    else if (param.name() == "word-penalty")
      return feature_function_ptr_type(new feature::WordPenalty());
    else if (param.name() == "rule-penalty")
      return feature_function_ptr_type(new feature::RulePenalty());
    else
      throw std::runtime_error("unknown featuer: " + parameter);
    
  }
  
};
