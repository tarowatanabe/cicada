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
	
	categories["[ADJP]"] = ba::list_of<category_list_type>(left, ba::list_of("[NNS]")("[QP]")("[NN]")("[$]")("[ADVP]")("[JJ]")("[VBN]")("[VBG]")("[ADJP]")("[JJR]")("[NP]")("[JJS]")("[DT]")("[FW]")("[RBR]")("[RBS]")("[SBAR]")("[RB]"));
	categories["[ADVP]"] = ba::list_of<category_list_type>(left, ba::list_of("[RB]")("[RBR]")("[RBS]")("[FW]")("[ADVP]")("[TO]")("[CD]")("[JJR]")("[JJ]")("[IN]")("[NP]")("[JJS]")("[NN]"));
	
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
