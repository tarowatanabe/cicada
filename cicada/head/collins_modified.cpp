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
    //static const char* __modcollins_FRAG[] = {}
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
    static const char* __modcollins_UCP[]    = {};
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
      
      punctuations = assign_category(__collins_puncts);
      std::sort(punctuations.begin(), punctuations.end());
      
      categories["[ADJP]"]   = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_ADJP));
      categories["[ADVP]"]   = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_ADVP));
      categories["[CONJP]"]  = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_CONJP));
      categories["[FRAG]"]   = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_FRAG));
      categories["[INTJ]"]   = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_INTJ));
      categories["[LST]"]    = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_LST));
      categories["[NAC]"]    = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_NAC));
      categories["[NX]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_NX));
      categories["[PP]"]     = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_PP));
      categories["[PRN]"]    = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_PRN));
      categories["[PRT]"]    = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_PRT));
      categories["[QP]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_QP));
      categories["[RRC]"]    = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_RRC));
      categories["[S]"]      = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_S));
      categories["[SBAR]"]   = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SBAR));
      categories["[SBARQ]"]  = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SBARQ));
      categories["[SINV]"]   = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SINV));
      categories["[SQ]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SQ));
      categories["[UCP]"]    = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_UCP));
      categories["[VP]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_VP));
      categories["[WHADJP]"] = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_WHADJP));
      categories["[WHADVP]"] = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_WHADVP));
      categories["[WHNP]"]   = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_WHNP));
      categories["[WHPP]"]   = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_WHPP));
      categories["[X]"]      = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_X));
      categories["[NP]"]     = boost::assign::list_of<category_list_type>(rightdis, assign_category(__collins_NP0))
	(left,     assign_category(__collins_NP1))
	(rightdis, assign_category(__collins_NP2))
	(right,    assign_category(__collins_NP3))
	(rightdis, assign_category(__collins_NP4));
      categories["[TYPO]"]   = boost::assign::list_of<category_list_type>(left,  assign_category(__collins_TYPO));
      categories["[EDITED]"] = boost::assign::list_of<category_list_type>(left,  assign_category(__collins_EDITED));
      categories["[XS]"]     = boost::assign::list_of<category_list_type>(right, assign_category(__collins_XS));
    }
    
  };
};
