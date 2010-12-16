//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/assign/list_of.hpp>

#include "collins.hpp"

namespace cicada
{
  namespace head
  {
    static const char* __collins_ADJP[]   = {"NNS", "QP", "NN", "$", "ADVP", "JJ", "VBN", "VBG",
					    "ADJP", "JJR", "NP", "JJS", "DT", "FW", "RBR", "RBS",
					    "SBAR", "RB"};
    static const char* __collins_ADVP[]   = {"RB", "RBR", "RBS", "FW", "ADVP", "TO", "CD", "JJR",
					    "JJ", "IN", "NP", "JJS", "NN"};
    static const char* __collins_CONJP[]  = {"CC", "RB", "IN"};
    static const char* __collins_FRAG[]   = {};
    static const char* __collins_INTJ[]   = {};
    static const char* __collins_LST[]    = {"LS", ":", "COLON"};
    static const char* __collins_NAC[]    = {"NN", "NNS", "NNP", "NNPS", "NP", "NAC", "EX", "$",
					    "CD", "QP", "PRP", "VBG", "JJ", "JJS", "JJR", "ADJP",
					    "FW"};
    static const char* __collins_NX[]     = {};
    static const char* __collins_PP[]     = {"right", "IN", "TO", "VBG", "VBN", "RP", "FW"};
    static const char* __collins_PRN[]    = {};
    static const char* __collins_PRT[]    = {"RP"};
    static const char* __collins_QP[]     = {"$", "IN", "NNS", "NN", "JJ", "RB", "DT", "CD",
					    "NCD", "QP", "JJR", "JJS"};
    static const char* __collins_RRC[]    = {"VP", "NP", "ADVP", "ADJP", "PP"};
    static const char* __collins_S[]      = {"TO", "IN", "VP", "S", "SBAR", "ADJP", "UCP", "NP"};
    static const char* __collins_SBAR[]   = {"WHNP", "WHPP", "WHADVP", "WHADJP", "IN", "DT", "S", "SQ",
					    "SINV", "SBAR", "FRAG"};
    static const char* __collins_SBARQ[]  = {"SQ", "S", "SINV", "SBARQ", "FRAG"};
    static const char* __collins_SINV[]   = {"VBZ", "VBD", "VBP", "VB", "MD", "VP", "S", "SINV",
					    "ADJP", "NP"};
    static const char* __collins_SQ[]     = {"VBZ", "VBD", "VBP", "VB", "MD", "VP", "SQ"};
    static const char* __collins_UCP[]    = {};
    static const char* __collins_VP[]     = {"TO", "VBD", "VBN", "MD", "VBZ", "VB", "VBG", "VBP",
					    "AUX", "AUXG", "VP", "ADJP", "NN", "NNS", "NP"};
    static const char* __collins_WHADJP[] = {"CC", "WRB", "JJ", "ADJP"};
    static const char* __collins_WHADVP[] = {"CC", "WRB"};
    static const char* __collins_WHNP[]   = {"WDT", "WP", "WP$", "WHADJP", "WHPP", "WHNP"};
    static const char* __collins_WHPP[]   = {"IN", "TO", "FW"};
    static const char* __collins_X[]      = {};
    static const char* __collins_NP0[]    = {"NN", "NNP", "NNPS", "NNS", "NX", "POS", "JJR"};
    static const char* __collins_NP1[]    = {"NP"};
    static const char* __collins_NP2[]    = {"$", "ADJP", "PRN"};
    static const char* __collins_NP3[]    = {"CD"};
    static const char* __collins_NP4[]    = {"JJ", "JJS", "RB", "QP"};
    static const char* __collins_TYPO[]   = {};
    static const char* __collins_EDITED[] = {};
    static const char* __collins_XS[]     = {"IN"};
    
    Collins::Collins() : HeadFinder("collins")
    {
      categories["[ADJP]"]     = boost::assign::list_of<category_list_type>(left,  assign_category(__collins_ADJP));
      categories["[ADVP]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__collins_ADVP));
      categories["[CONJP]"]    = boost::assign::list_of<category_list_type>(right, assign_category(__collins_CONJP));
      categories["[FRAG]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__collins_FRAG));
      categories["[INTJ]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__collins_INTJ));
      categories["[LST]"]      = boost::assign::list_of<category_list_type>(right, assign_category(__collins_LST));
      categories["[NAC]"]      = boost::assign::list_of<category_list_type>(left, assign_category(__collins_NAC));
      categories["[NX]"]       = boost::assign::list_of<category_list_type>(left, assign_category(__collins_NX));
      categories["[PP]"]       = boost::assign::list_of<category_list_type>(right, assign_category(__collins_PP));
      categories["[PRN]"]      = boost::assign::list_of<category_list_type>(left, assign_category(__collins_PRN));
      categories["[PRT]"]      = boost::assign::list_of<category_list_type>(right, assign_category(__collins_PRT));
      categories["[QP]"]       = boost::assign::list_of<category_list_type>(left, assign_category(__collins_QP));
      categories["[RRC]"]      = boost::assign::list_of<category_list_type>(right, assign_category(__collins_RRC));
      categories["[S]"]        = boost::assign::list_of<category_list_type>(left, assign_category(__collins_S));
      categories["[SBAR]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__collins_SBAR));
      categories["[SBARQ]"]    = boost::assign::list_of<category_list_type>(left, assign_category(__collins_SBARQ));
      categories["[SINV]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__collins_SINV));
      categories["[SQ]"]       = boost::assign::list_of<category_list_type>(left, assign_category(__collins_SQ));
      categories["[UCP]"]      = boost::assign::list_of<category_list_type>(right, assign_category(__collins_UCP));
      categories["[VP]"]       = boost::assign::list_of<category_list_type>(left, assign_category(__collins_VP));
      categories["[WHADJP]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__collins_WHADJP));
      categories["[WHADVP]"]   = boost::assign::list_of<category_list_type>(right, assign_category(__collins_WHADVP));
      categories["[WHNP]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__collins_WHNP));
      categories["[WHPP]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__collins_WHPP));
      categories["[X]"]        = boost::assign::list_of<category_list_type>(right, assign_category(__collins_X));
      categories["[NP]"]       = boost::assign::list_of<category_list_type>(rightdis, assign_category(__collins_NP0))
	(left, assign_category(__collins_NP1))
	(rightdis, assign_category(__collins_NP2))
	(right, assign_category(__collins_NP3))
	(rightdis, assign_category(__collins_NP4));
      categories["[TYPO]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__collins_TYPO));
      categories["[EDITED]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__collins_EDITED));
      categories["[XS]"]       = boost::assign::list_of<category_list_type>(right, assign_category(__collins_XS));
    }
    
  };
};
