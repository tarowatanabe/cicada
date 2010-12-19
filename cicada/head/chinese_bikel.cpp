//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/assign/list_of.hpp>

#include "chinese_bikel.hpp"

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
    
    ChineseBikel::ChineseBikel() : HeadFinder("chinese-bikel")
    {
      categories["[ROOT]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_ROOT));
      categories["[PAIR]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_PAIR));
      categories["[ADJP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_ADJP));
      categories["[ADVP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_ADVP));
      categories["[CLP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_CLP));
      categories["[CP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_CP))
	(rightexcept, assign_category(__chinese_punct));
      categories["[DNP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_DNP))
	(rightexcept, assign_category(__chinese_punct));
      categories["[DP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_DP));
      categories["[DVP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_DVP));
      categories["[FRAG]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_FRAG))
	(rightexcept, assign_category(__chinese_punct));
      categories["[INTJ]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_INTJ));
      categories["[IP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_IP))
	(rightexcept, assign_category(__chinese_punct));
      categories["[LCP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_LCP));
      categories["[LST]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_LST));
      categories["[NP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_NP));
      categories["[PP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_PP));
      categories["[PRN]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_PRN0))
	(rightdis, assign_category(__chinese_PRN1));
      categories["[QP]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_QP));
      categories["[UCP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_UCP));
      categories["[VP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VP))
	(leftexcept, assign_category(__chinese_punct));
      categories["[VCD]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VCD));
      categories["[VCP]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VCP));
      categories["[VRD]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VRD));
      categories["[VSB]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_VSB));
      categories["[VNV]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VNV));
      categories["[VPT]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VPT));
      categories["[CD]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_CD));
      categories["[NN]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_NN));
      categories["[NR]"] = boost::assign::list_of<category_list_type>(right, assign_category(__chinese_NR));
      categories["[VV]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VV));
      categories["[VA]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VA));
      categories["[VC]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VC));
      categories["[VE]"] = boost::assign::list_of<category_list_type>(left, assign_category(__chinese_VE));
      categories["[FLR]"] = boost::assign::list_of<category_list_type>(rightexcept, assign_category(__chinese_punct));
    }
    
  };
};
