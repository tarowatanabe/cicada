//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <boost/assign/list_of.hpp>

#include "collins.hpp"

namespace cicada
{
  namespace head
  {
    static String[] leftExceptPunct = {"leftexcept", "PU"};
    static String[] rightExceptPunct = {"rightexcept", "PU"};
    
    static const char* __chinese_punct[] = {"PU"};
    
    static const char* __chinese_ROOT[] = {"IP"};
    static const char* __chinese_PAIR[] = {"IP"};
    static const char* __chinese_ADJP[] = {"JJ", "ADJP"};
    static const char* __chinese_ADVP[] = {"AD", "CS", "ADVP", "JJ"};
    static const char* __chinese_CLP[] = {"M", "CLP"};
    static const char* __chinese_CP[] = {"DEC", "WHNP", "WHPP"};
    static const char* __chinese_DNP[] = {"DEG", "DEC"};
    
    
    nonTerminalInfo.put("DP", new String[][]{{left, "DT", "DP"}}); // there's one instance of DP adjunction
    nonTerminalInfo.put("DVP", new String[][]{{right, "DEV", "DEC"}}); // DVP always has DEV under it
    nonTerminalInfo.put("FRAG", new String[][]{{right, "VV", "NN"}, rightExceptPunct}); //FRAG seems only to be used for bits at the beginnings of articles: "Xinwenshe<DATE>" and "(wan)"
    nonTerminalInfo.put("INTJ", new String[][]{{right, "INTJ", "IJ", "SP"}});
    nonTerminalInfo.put("IP", new String[][]{{left, "VP", "IP"}, rightExceptPunct});  // CDM July 2010 following email from Pi-Chuan changed preference to VP over IP: IP can be -SBJ, -OBJ, or -ADV, and shouldn't be head
    nonTerminalInfo.put("LCP", new String[][]{{right, "LC", "LCP"}}); // there's a bit of LCP adjunction
    nonTerminalInfo.put("LST", new String[][]{{right, "CD", "PU"}}); // covers all examples
    nonTerminalInfo.put("NP", new String[][]{{right, "NN", "NR", "NT", "NP", "PN", "CP"}}); // Basic heads are NN/NR/NT/NP; PN is pronoun.  Some NPs are nominalized relative clauses without overt nominal material; these are NP->CP unary rewrites.  Finally, note that this doesn't give any special treatment of coordination.
    nonTerminalInfo.put("PP", new String[][]{{left, "P", "PP"}}); // in the manual there's an example of VV heading PP but I couldn't find such an example with tgrep2
    // cdm 2006: PRN changed to not choose punctuation.  Helped parsing (if not significantly)
    // nonTerminalInfo.put("PRN", new String[][]{{left, "PU"}}); //presumably left/right doesn't matter
    nonTerminalInfo.put("PRN", new String[][]{{left, "NP", "VP", "IP", "QP", "PP", "ADJP", "CLP", "LCP"}, {rightdis, "NN", "NR", "NT", "FW"}});
    // cdm 2006: QP: add OD -- occurs some; occasionally NP, NT, M; parsing performance no-op
    nonTerminalInfo.put("QP", new String[][]{{right, "QP", "CLP", "CD", "OD", "NP", "NT", "M"}}); // there's some QP adjunction
    // add OD?
    nonTerminalInfo.put("UCP", new String[][]{{left, }}); //an alternative would be "PU","CC"
    nonTerminalInfo.put("VP", new String[][]{{left, "VP", "VCD", "VPT", "VV", "VCP", "VA", "VC", "VE", "IP", "VSB", "VCP", "VRD", "VNV"}, leftExceptPunct}); //note that ba and long bei introduce IP-OBJ small clauses; short bei introduces VP
    // add BA, LB, as needed

    // verb compounds
    nonTerminalInfo.put("VCD", new String[][]{{left, "VCD", "VV", "VA", "VC", "VE"}}); // could easily be right instead
    nonTerminalInfo.put("VCP", new String[][]{{left, "VCD", "VV", "VA", "VC", "VE"}}); // not much info from documentation
    nonTerminalInfo.put("VRD", new String[][]{{left, "VCD", "VRD", "VV", "VA", "VC", "VE"}}); // definitely left
    nonTerminalInfo.put("VSB", new String[][]{{right, "VCD", "VSB", "VV", "VA", "VC", "VE"}}); // definitely right, though some examples look questionably classified (na2lai2 zhi1fu4)
    nonTerminalInfo.put("VNV", new String[][]{{left, "VV", "VA", "VC", "VE"}}); // left/right doesn't matter
    nonTerminalInfo.put("VPT", new String[][]{{left, "VV", "VA", "VC", "VE"}}); // activity verb is to the left

    // some POS tags apparently sit where phrases are supposed to be
    nonTerminalInfo.put("CD", new String[][]{{right, "CD"}});
    nonTerminalInfo.put("NN", new String[][]{{right, "NN"}});
    nonTerminalInfo.put("NR", new String[][]{{right, "NR"}});

    // I'm adding these POS tags to do primitive morphology for character-level
    // parsing.  It shouldn't affect anything else because heads of preterminals are not
    // generally queried - GMA
    nonTerminalInfo.put("VV", new String[][]{{left}});
    nonTerminalInfo.put("VA", new String[][]{{left}});
    nonTerminalInfo.put("VC", new String[][]{{left}});
    nonTerminalInfo.put("VE", new String[][]{{left}});

    // new for ctb6.
    nonTerminalInfo.put("FLR", new String[][]{rightExceptPunct});

    
    
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
