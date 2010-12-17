//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/assign/list_of.hpp>

#include "collins.hpp"

namespace cicada
{
  namespace head
  {
    static const char* __chinese_punct[] = {"PU"};
    
    static const char* __chinese_ROOT[] = {"IP"};
    static const char* __chinese_PAIR[] = {"IP"};
    static const char* __chinese_ADJP[] = {"JJ", "ADJP"};
    static const char* __chinese_ADVP[] = {"AD", "CS", "ADVP", "JJ"};
    static const char* __chinese_CLP[]  = {"M", "CLP"};
    static const char* __chinese_CP[]   = {"DEC", "WHNP", "WHPP"};
    static const char* __chinese_DNP[]  = {"DEG", "DEC"};
    static const char* __chinese_DP[]   = {"DT", "DP"};
    static const char* __chinese_DVP[]  = {"DEV", "DEC"};
    static const char* __chinese_FRAG[] = {"VV", "NN"};
    static const char* __chinese_INTJ[] = {"INTJ", "IJ", "SP"};
    static const char* __chinese_IP[]   = {"VP", "IP"};
    static const char* __chinese_LCP[]  = {"LC", "LCP"};
    static const char* __chinese_LST[]  = {"CD", "PU"};
    static const char* __chinese_NP[]   = {"NN", "NR", "NT", "NP", "PN", "CP"};
    static const char* __chinese_PP[]   = {"P", "PP"};
    static const char* __chinese_PRN0[] = {"NP", "VP", "IP", "QP", "PP", "ADJP", "CLP", "LCP"};
    static const char* __chinese_PRN1[] = {"NN", "NR", "NT", "FW"};
    static const char* __chinese_QP[]   = {"QP", "CLP", "CD", "OD", "NP", "NT", "M"};
    static const char* __chinese_UCP[]  = {};
    static const char* __chinese_VP[]   = {"VP", "VCD", "VPT", "VV", "VCP", "VA", "VC", "VE", "IP", "VSB", "VCP", "VRD", "VNV"};
    static const char* __chinese_VCD[]  = {"VCD", "VV", "VA", "VC", "VE"};
    static const char* __chinese_VCP[]  = {"VCD", "VV", "VA", "VC", "VE"};
    static const char* __chinese_VRD[]  = {"VCD", "VRD", "VV", "VA", "VC", "VE"};
    static const char* __chinese_VSB[]  = {"VCD", "VSB", "VV", "VA", "VC", "VE"};
    static const char* __chinese_VNV[]  = {"VV", "VA", "VC", "VE"};
    static const char* __chinese_VPT[]  = {"VV", "VA", "VC", "VE"};
    static const char* __chinese_CD[]   = {"CD"};
    static const char* __chinese_NN[]   = {"NN"};
    static const char* __chinese_NR[]   = {"NR"};
    
    static const char* __chinese_VV[]   = {};
    static const char* __chinese_VA[]   = {};
    static const char* __chinese_VC[]   = {};
    static const char* __chinese_VE[]   = {};
    static const char* __chinese_FLR[]  = {};
    
    Chinese::Chinese() : HeadFinder("chinese")
    {
      categories["[ADJP]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_ADJP));
      categories["[ADVP]"]     = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_ADVP));
      categories["[CONJP]"]    = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_CONJP));
      categories["[FRAG]"]     = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_FRAG));
      categories["[INTJ]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_INTJ));
      categories["[LST]"]      = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_LST));
      categories["[NAC]"]      = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_NAC));
      categories["[NX]"]       = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_NX));
      categories["[PP]"]       = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_PP));
      categories["[PRN]"]      = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_PRN));
      categories["[PRT]"]      = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_PRT));
      categories["[QP]"]       = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_QP));
      categories["[RRC]"]      = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_RRC));
      categories["[S]"]        = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_S));
      categories["[SBAR]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SBAR));
      categories["[SBARQ]"]    = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SBARQ));
      categories["[SINV]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SINV));
      categories["[SQ]"]       = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_SQ));
      categories["[UCP]"]      = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_UCP));
      categories["[VP]"]       = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_VP));
      categories["[WHADJP]"]   = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_WHADJP));
      categories["[WHADVP]"]   = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_WHADVP));
      categories["[WHNP]"]     = boost::assign::list_of<category_list_type>(left,     assign_category(__collins_WHNP));
      categories["[WHPP]"]     = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_WHPP));
      categories["[X]"]        = boost::assign::list_of<category_list_type>(right,    assign_category(__collins_X));
      categories["[NP]"]       = boost::assign::list_of<category_list_type>(rightdis, assign_category(__collins_NP0))
	(left, assign_category(__collins_NP1))
	(rightdis, assign_category(__collins_NP2))
	(right, assign_category(__collins_NP3))
	(rightdis, assign_category(__collins_NP4));
      categories["[TYPO]"]     = boost::assign::list_of<category_list_type>(left,  assign_category(__collins_TYPO));
      categories["[EDITED]"]   = boost::assign::list_of<category_list_type>(left,  assign_category(__collins_EDITED));
      categories["[XS]"]       = boost::assign::list_of<category_list_type>(right, assign_category(__collins_XS));
    }
    
  };
};
