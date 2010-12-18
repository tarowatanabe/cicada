//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/assign/list_of.hpp>

#include "collins_modified.hpp"

namespace cicada
{
  namespace head
  {
    static const char* __modcollins_puncts[] = {"''", "``", "-LRB-", "-RRB-", ".", ":", ",", "PERIOD", "COLON", "COMMA"};
    
    static const char* __modcollins_ADJP0[]  = {"$"};
    static const char* __modcollins_ADJP1[]  = {"NNS", "NN", "JJ", "QP", "VBN", "VBG"};
    static const char* __modcollins_ADJP2[]  = {"ADJP"};
    static const char* __modcollins_ADJP3[]  = {"JJP", "JJR", "JJS", "DT", "RB", "RBR", "CD", "IN", "VBD"};
    static const char* __modcollins_ADJP4[]  = {"ADVP", "NP"};
    static const char* __modcollins_JJP[]    = {"NNS", "NN", "$", "QP", "JJ", "VBN", "VBG", "ADJP",
						"JJP", "JJR", "NP", "JJS", "DT", "FW", "RBR", "RBS",
						"SBAR", "RB"};
    static const char* __modcollins_ADVP0[]  = {"ADVP", "IN"};
    static const char* __modcollins_ADVP1[]  = {"RB", "RBR", "RBS", "JJ", "JJR", "JJS"};
    static const char* __modcollins_ADVP2[]  = {"RP", "DT", "NN", "CD", "NP", "VBN", "NNP", "CC", "FW", "NNS", "ADJP", "NML"};
    static const char* __modcollins_CONJP[]  = {"CC", "RB", "IN"};
    //static const char* __modcollins_FRAG[] = {} // use puncts...
    static const char* __modcollins_INTJ[]   = {};
    static const char* __modcollins_LST[]    = {"LS", ":", "COMMA"};
    static const char* __modcollins_NAC[]    = {"NN", "NNS", "NML", "NNP", "NNPS", "NP", "NAC", "EX",
						"$", "CD", "QP", "PRP", "VBG", "JJ", "JJS", "JJR",
						"ADJP", "JJP", "FW"};
    static const char* __modcollins_NX[]     = {"NP", "NX"};
    static const char* __modcollins_PP0[]    = {"IN", "TO", "VBG", "VBN", "RP", "FW", "JJ"};
    static const char* __modcollins_PP1[]    = {"PP"};
    static const char* __modcollins_PRN[]    = {"VP", "NP", "PP", "S", "SINV", "SBAR", "ADJP", "JJP",
						"ADVP", "INTJ", "WHNP", "NAC", "VBP", "JJ", "NN", "NNP"};
    static const char* __modcollins_PRT[]    = {"RP"};
    static const char* __modcollins_QP[]     = {"$", "IN", "NNS", "NN", "JJ", "CD", "PDT", "DT",
						"RB", "NCD", "QP", "JJR", "JJS"};
    static const char* __modcollins_RRC[]    = {"VP", "NP", "ADVP", "ADJP", "JJP", "PP"};
    static const char* __modcollins_S[]      = {"TO", "VP", "S", "FRAG", "SBAR", "ADJP", "JJP", "UCP", "NP"};
    static const char* __modcollins_SBAR[]   = {"WHNP", "WHPP", "WHADVP", "WHADJP", "IN", "DT", "S", "SQ",
						"SINV", "SBAR", "FRAG"};
    static const char* __modcollins_SBARQ[]  = {"SQ", "S", "SINV", "SBARQ", "FRAG", "SBAR"};
    static const char* __modcollins_SINV[]   = {"VBZ", "VBD", "VBP", "VB", "MD", "VP", "S", "SINV",
						"ADJP", "JJP", "NP"};
    static const char* __modcollins_SQ[]     = {"VBZ", "VBD", "VBP", "VB", "MD", "AUX", "AUXG", "VP",
						"SQ"};
    //static const char* __modcollins_UCP[]    = {}; // use puncts
    static const char* __modcollins_VP[]     = {"TO", "VBD", "VBN", "MD", "VBZ", "VB", "VBG", "VBP",
						"VP", "AUX", "AUXG", "ADJP", "JJP", "NN", "NNS", "JJ",
						"NP", "NNP"};
    static const char* __modcollins_WHADJP[] = {"WRB", "WHADVP", "RB", "JJ", "ADJP", "JJP", "JJR"};
    static const char* __modcollins_WHADVP[] = {"WRB", "WHADVP"};
    static const char* __modcollins_WHNP[]   = {"WDT", "WP", "WP$", "WHADJP", "WHPP", "WHNP"};
    static const char* __modcollins_WHPP[]   = {"IN", "TO", "FW"};
    static const char* __modcollins_X[]      = {"S", "VP", "ADJP", "JJP", "NP", "SBAR", "PP", "X"};
    static const char* __modcollins_NP0[]    = {"NN", "NNP", "NNPS", "NNS", "NX", "POS", "JJR"};
    static const char* __modcollins_NP1[]    = {"NP", "NML", "PRP"};
    static const char* __modcollins_NP2[]    = {"$", "ADJP", "JJP", "PRN", "FW"};
    static const char* __modcollins_NP3[]    = {"CD"};
    static const char* __modcollins_NP4[]    = {"JJ", "JJS", "RB", "QP", "DT", "WDT", "RBR", "ADVP"};
    static const char* __modcollins_NML0[]   = {"NN", "NNP", "NNPS", "NNS", "NX", "POS", "JJR"};
    static const char* __modcollins_NML1[]   = {"NP", "NML", "PRP"};
    static const char* __modcollins_NML2[]   = {"$", "ADJP", "JJP", "PRN"};
    static const char* __modcollins_NML3[]   = {"CD"};
    static const char* __modcollins_NML4[]   =  {"JJ", "JJS", "RB", "QP", "DT", "WDT", "RBR", "ADVP"};
    static const char* __modcollins_POSSP[]  = {"POS"};
    static const char* __modcollins_ROOT[]   = {"S", "SQ", "SINV", "SBAR", "FRAG"};
    static const char* __modcollins_TYPO[]   = {"NN", "NP", "NML", "NNP", "NNPS", "TO", "VBD", "VBN",
						"MD", "VBZ", "VB", "VBG", "VBP", "VP", "ADJP", "JJP",
						"FRAG"};
    static const char* __modcollins_ADV[]    = {"RB", "RBR", "RBS", "FW", "ADVP", "TO", "CD", "JJR",
						"JJ", "IN", "NP", "NML", "JJS", "NN"};
    static const char* __modcollins_EDITED[] = {};
    static const char* __modcollins_META[]   = {};
    static const char* __modcollins_XS[]     = {"IN"};
    
    CollinsModified::CollinsModified() : HeadFinder("collins-modified")
    {
      punctuations = assign_category(__modcollins_puncts);
      std::sort(punctuations.begin(), punctuations.end());
      
      categories["[ADJP]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_ADJP0))
	(rightdis,    assign_category(__modcollins_ADJP1))
	(left,        assign_category(__modcollins_ADJP2))
	(rightdis,    assign_category(__modcollins_ADJP3))
	(left,        assign_category(__modcollins_ADJP4))
	(rightexcept, assign_category(__modcollins_puncts));
      categories["[JJP]"]    = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_JJP));
      categories["[ADVP]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_ADVP0))
	(rightdis,    assign_category(__modcollins_ADVP1))
	(rightdis,    assign_category(__modcollins_ADVP2))
	(rightexcept, assign_category(__modcollins_puncts));
      categories["[CONJP]"]  = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_CONJP));
      categories["[FRAG]"]   = boost::assign::list_of<category_list_type>(rightexcept, assign_category(__modcollins_puncts)); // FRAG..
      categories["[INTJ]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_INTJ));
      categories["[LST]"]    = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_LST));
      categories["[NAC]"]    = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_NAC));
      categories["[NX]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_NX));
      categories["[PP]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_PP0))
	(right, assign_category(__modcollins_PP1));
      categories["[PRN]"]    = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_PRN))
	(leftexcept, assign_category(__modcollins_puncts));
      categories["[PRT]"]    = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_PRT));
      categories["[QP]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_QP));
      categories["[RRC]"]    = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_RRC));
      categories["[S]"]      = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_S));
      categories["[SBAR]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_SBAR));
      categories["[SBARQ]"]  = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_SBARQ));
      categories["[SINV]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_SINV));
      categories["[SQ]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_SQ));
      categories["[UCP]"]    = boost::assign::list_of<category_list_type>(rightexcept, assign_category(__modcollins_puncts)); // UCP
      categories["[VP]"]     = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_VP));
      categories["[WHADJP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_WHADJP));
      categories["[WHADVP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_WHADVP));
      categories["[WHNP]"]   = boost::assign::list_of<category_list_type>(left, assign_category(__modcollins_WHNP));
      categories["[WHPP]"]   = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_WHPP));
      categories["[X]"]      = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_X))
	(rightexcept, assign_category(__modcollins_puncts));
      categories["[NP]"]     = boost::assign::list_of<category_list_type>(rightdis, assign_category(__modcollins_NP0))
	(left,        assign_category(__modcollins_NP1))
	(rightdis,    assign_category(__modcollins_NP2))
	(right,       assign_category(__modcollins_NP3))
	(rightdis,    assign_category(__modcollins_NP4))
	(rightexcept, assign_category(__modcollins_puncts));
      categories["[NML]"]    = boost::assign::list_of<category_list_type>(rightdis, assign_category(__modcollins_NML0))
	(left,        assign_category(__modcollins_NML1))
	(rightdis,    assign_category(__modcollins_NML2))
	(right,       assign_category(__modcollins_NML3))
	(rightdis,    assign_category(__modcollins_NML4))
	(rightexcept, assign_category(__modcollins_puncts));
      categories["[POSSP]"]  = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_POSSP));
      categories["[ROOT]"]   = boost::assign::list_of<category_list_type>(left,  assign_category(__modcollins_ROOT));
      categories["[TYPO]"]   = boost::assign::list_of<category_list_type>(left,  assign_category(__modcollins_TYPO));
      categories["[ADV]"]    = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_ADV));
      categories["[EDITED]"] = boost::assign::list_of<category_list_type>(left,  assign_category(__modcollins_EDITED));
      categories["[META]"]   = boost::assign::list_of<category_list_type>(left,  assign_category(__modcollins_META));
      categories["[XS]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__modcollins_XS));
    }
    
  };
};
