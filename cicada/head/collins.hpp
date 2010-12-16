// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__HEAD__COLLINS__HPP__
#define __CICADA__HEAD__COLLINS__HPP__ 1

#include <boost/assign/list_of.hpp>
#include <boost/assign/list_inserter.hpp>

namespace cicada
{
  namespace head
  {
    
    class Collins : public cicada::HeadFinder
    {
    public:
      Collins() : HeadFinder("collins")
      {
	namespace ba = boost::assign;
	
	categories["[ADJP]"]  = ba::list_of<category_list_type>(left, ba::list_of("[NNS]")("[QP]")("[NN]")("[$]")
								("[ADVP]")("[JJ]")("[VBN]")("[VBG]")
								("[ADJP]")("[JJR]")("[NP]")("[JJS]")
								("[DT]")("[FW]")("[RBR]")("[RBS]")
								("[SBAR]")("[RB]"));
	
	categories["[ADVP]"]  = ba::list_of<category_list_type>(right, ba::list_of("[RB]")("[RBR]")("[RBS]")("[FW]")
								("[ADVP]")("[TO]")("[CD]")("[JJR]")
								("[JJ]")("[IN]")("[NP]")("[JJS]")
								("[NN]"));
	categories["[CONJP]"] = ba::list_of<category_list_type>(right, ba::list_of("[CC]")("[RB]")("[IN]"));
	categories["[FRAG]"]  = ba::list_of<category_list_type>(right, category_set_type());
	categories["[INTJ]"]  = ba::list_of<category_list_type>(left, category_set_type());
	categories["[LST]"]   = ba::list_of<category_list_type>(right, ba::list_of("[LS]")("[:]")("[COLON]"));
	categories["[NAC]"]   = ba::list_of<category_list_type>(left, ba::list_of("[NN]")("[NNS]")("[NNP]")("[NNPS]")
								("[NP]")("[NAC]")("[EX]")("[$]")
								("[CD]")("[QP]")("[PRP]")("[VBG]")
								("[JJ]")("[JJS]")("[JJR]")("[ADJP]")
								("[FW]"));
	categories["[NX]"]    = ba::list_of<category_list_type>(left, category_set_type());
	categories["[PP]"]    = ba::list_of<category_list_type>(right, ba::list_of("[IN]")("[TO]")("[VBG]")("[VBN]")
								("[RP]")("[FW]"));
	categories["[PRN]"]   = ba::list_of<category_list_type>(left, category_set_type());
	categories["[PRT]"]   = ba::list_of<category_list_type>(right, ba::list_of("[RP]"));
	categories["[QP]"]    = ba::list_of<category_list_type>(left, ba::list_of("[$]")("[IN]")("[NNS]")("[NN]")
								("[JJ]")("[RB]")("[DT]")("[CD]")
								("[NCD]")("[QP]")("[JJR]")("[JJS]"));
	categories["[PPC]"]   = ba::list_of<category_list_type>(right, ba::list_of("[VP]")("[NP]")("[ADVP]")("[ADJP]")
								("[PP]"));
	categories["[S]"]     = ba::list_of<category_list_type>(left, ba::list_of("[TO]")("[IN]")("[VP]")("[S]")
								("[SBAR]")("[ADJP]")("[UCP]")("[NP]"));
	categories["[SBAR]"]  = ba::list_of<category_list_type>(left, ba::list_of("[WHNP]")("[WHPP]")("[WHADVP]")("[WHADJP]")
								("[IN]")("[DT]")("[S]")("[SQ]")
								("[SINV]")("[SBAR]")("[FRAG]"));
	categories["[SBARQ]"] = ba::list_of<category_list_type>(left, ba::list_of("[SQ]")("[S]")("[SINV]")("[SBARQ]")
								("[FRAG]"));
	categories["[SINV]"]  = ba::list_of<category_list_type>(left, ba::list_of("[VBZ]")("[VBD]")("[VBP]")("[VB]")
								("[MD]")("[VP]")("[S]")("[SINV]")
								("[ADVP]")("[NP]"));
	categories["[SQ]"]    = ba::list_of<category_list_type>(left, ba::list_of("[VBZ]")("[VBD]")("[VBP]")("[VB]")
								("[MD]")("[VP]")("[SQ]"));
	categories["[UCP]"]   = ba::list_of<category_list_type>(right, category_set_type());
      }
      
      size_type find_marked_head(const rule_type& rule, const symbol_type& parent) const { return size_type(-1); }
      size_type find_head(const rule_type& rule, const symbol_type& parent) const
      {
	
	
      }
      
    private:
      category_info_type categories;
    };
  };
};

#endif
