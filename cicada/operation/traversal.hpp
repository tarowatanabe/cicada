// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__TRAVERSAL__HPP__
#define __CICADA__OPERATION__TRAVERSAL__HPP__ 1

#include <vector>

#include <cicada/alignment.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/rule.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring.hpp>

#include <utils/sgi_hash_set.hpp>

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace operation
  {
    struct kbest_traversal_edges
    {
      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Rule       rule_type;
      
      typedef hypergraph_type::feature_set_type feature_set_type;
      typedef std::vector<hypergraph_type::id_type, std::allocator<hypergraph_type::id_type> > edge_set_type;

      typedef boost::tuple<edge_set_type, feature_set_type> value_type;
  
      template <typename Edge, typename Iterator>
      void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
      {
	boost::get<0>(yield).clear();

	boost::get<0>(yield).push_back(edge.id);
	boost::get<1>(yield) = edge.features;
    
	// collect edge and features
	for (/**/; first != last; ++ first) {
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*first).begin(), boost::get<0>(*first).end());
	  boost::get<1>(yield) += boost::get<1>(*first);
	}
      }
    };

    struct kbest_alignment_traversal
    {
      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Rule       rule_type;
      typedef cicada::Alignment  alignment_type;
      typedef cicada::Vocab      vocab_type;
      
      typedef hypergraph_type::feature_set_type   feature_set_type;
      typedef hypergraph_type::attribute_set_type attribute_set_type;

      typedef attribute_set_type::attribute_type attribute_type;
  
      typedef boost::tuple<alignment_type, feature_set_type> value_type;

      kbest_alignment_traversal() :
	attr_source_position("source-position"),
	attr_target_position("target-position") {}

      attribute_type attr_source_position;
      attribute_type attr_target_position;

      struct __point : public boost::static_visitor<int>
      {
	int operator()(const attribute_set_type::int_type& x) const { return x; }
	int operator()(const attribute_set_type::float_type& x) const { return -1; }
	int operator()(const attribute_set_type::string_type& x) const { return -1; }
      };

      template <typename Edge, typename Iterator>
      void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
      {
	boost::get<0>(yield).clear();
	boost::get<1>(yield) = edge.features;
	
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_source_position);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);

	if (siter != edge.attributes.end() && titer != edge.attributes.end()) {
	  const int source_pos = boost::apply_visitor(__point(), siter->second);
	  const int target_pos = boost::apply_visitor(__point(), titer->second);
	  
	  if (source_pos >= 0 && target_pos >= 0)
	    boost::get<0>(yield).push_back(std::make_pair(source_pos, target_pos));
	}
	
	// collect features...
	for (/**/; first != last; ++ first) {
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*first).begin(), boost::get<0>(*first).end());
	  boost::get<1>(yield) += boost::get<1>(*first);
	}
      }
    };


    struct kbest_traversal
    {
      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Rule       rule_type;
      typedef cicada::Sentence   sentence_type;
      typedef cicada::Vocab      vocab_type;
      
      typedef hypergraph_type::feature_set_type   feature_set_type;
      typedef hypergraph_type::attribute_set_type attribute_set_type;
  
      typedef boost::tuple<sentence_type, feature_set_type> value_type;

      kbest_traversal() {}
      kbest_traversal(const std::string& __insertion_prefix) : insertion_prefix(__insertion_prefix) {}

      struct __inserted : public boost::static_visitor<bool>
      {
	bool operator()(const attribute_set_type::int_type& x) const { return x; }
	bool operator()(const attribute_set_type::float_type& x) const { return false; }
	bool operator()(const attribute_set_type::string_type& x) const { return false; }
      };

      template <typename Edge, typename Iterator>
      void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
      {
	boost::get<0>(yield).clear();
	boost::get<1>(yield) = edge.features;
	
	bool is_insertion = false;
	if (! insertion_prefix.empty()) {
	  attribute_set_type::const_iterator aiter = edge.attributes.find("insertion");
	  if (aiter != edge.attributes.end())
	    is_insertion = boost::apply_visitor(__inserted(), aiter->second);
	}

	if (! is_insertion) {
	  int non_terminal_pos = 0;
	  rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	    if (titer->is_non_terminal()) {
	      int pos = titer->non_terminal_index() - 1;
	      if (pos < 0)
		pos = non_terminal_pos;
	      ++ non_terminal_pos;
	      
	      boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + pos)).begin(), boost::get<0>(*(first + pos)).end());
	    } else if (*titer != vocab_type::EPSILON)
	      boost::get<0>(yield).push_back(*titer);
	} else {
	  int non_terminal_pos = 0;
	  rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	    if (titer->is_non_terminal()) {
	      int pos = titer->non_terminal_index() - 1;
	      if (pos < 0)
		pos = non_terminal_pos;
	      ++ non_terminal_pos;
	      
	      boost::get<0>(yield).insert(boost::get<0>(yield).end(), boost::get<0>(*(first + pos)).begin(), boost::get<0>(*(first + pos)).end());
	    } else if (*titer != vocab_type::EPSILON)
	      boost::get<0>(yield).push_back(insertion_prefix + static_cast<const std::string&>(*titer));
	}
    
	// collect features...
	for (/**/; first != last; ++ first)
	  boost::get<1>(yield) += boost::get<1>(*first);
      }

      std::string insertion_prefix;
    };

    struct kbest_alignment_filter
    {
      template <typename Node, typename Yield>
      bool operator()(const Node& node, const Yield& yield) const
      {
	return false;
      }
    };
    

    struct kbest_filter
    {
      template <typename Node, typename Yield>
      bool operator()(const Node& node, const Yield& yield) const
      {
	return false;
      }
    };

    struct kbest_alignment_filter_unique
    {
      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Alignment  alignment_type;
      
#ifdef HAVE_TR1_UNORDERED_SET
      typedef std::tr1::unordered_set<alignment_type, boost::hash<alignment_type>, std::equal_to<alignment_type>, std::allocator<alignment_type> > unique_type;
#else
      typedef sgi::hash_set<alignment_type, boost::hash<alignment_type>, std::equal_to<alignment_type>, std::allocator<alignment_type> > unique_type;
#endif
      typedef std::vector<unique_type, std::allocator<unique_type> > unique_set_type;
 

      kbest_alignment_filter_unique(const hypergraph_type& graph) : uniques(graph.nodes.size()) {}
  
      template <typename Node, typename Yield>
      bool operator()(const Node& node, const Yield& yield) const
      {
	unique_set_type& aligns = const_cast<unique_set_type&>(uniques);
	unique_type::iterator iter = aligns[node.id].find(boost::get<0>(yield));
	if (iter == aligns[node.id].end()) {
	  aligns[node.id].insert(boost::get<0>(yield));
	  return false;
	} else
	  return true;
      }

      unique_set_type uniques;
    };

    struct kbest_filter_unique
    {
      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Sentence sentence_type;
      
#ifdef HAVE_TR1_UNORDERED_SET
      typedef std::tr1::unordered_set<sentence_type, boost::hash<sentence_type>, std::equal_to<sentence_type>, std::allocator<sentence_type> > unique_type;
#else
      typedef sgi::hash_set<sentence_type, boost::hash<sentence_type>, std::equal_to<sentence_type>, std::allocator<sentence_type> > unique_type;
#endif
      typedef std::vector<unique_type, std::allocator<unique_type> > unique_set_type;
 

      kbest_filter_unique(const hypergraph_type& graph) : uniques(graph.nodes.size()) {}
  
      template <typename Node, typename Yield>
      bool operator()(const Node& node, const Yield& yield) const
      {
	unique_set_type& sents = const_cast<unique_set_type&>(uniques);
	unique_type::iterator iter = sents[node.id].find(boost::get<0>(yield));
	if (iter == sents[node.id].end()) {
	  sents[node.id].insert(boost::get<0>(yield));
	  return false;
	} else
	  return true;
      }

      unique_set_type uniques;
    };

  };
};

#endif
