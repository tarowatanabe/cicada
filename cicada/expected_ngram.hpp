// -*- mode: c++ -*-

#ifndef __CICADA__EXPECTED_NGRAM__HPP__
#define __CICADA__EXPECTED_NGRAM__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/inside_outside.hpp>

namespace cicada
{
  
  template <typename Function, typename Counts>
  class ExpectedNGram
  {
  public:
    typedef HyperGraph hypergraph_type;
    
    typedef hypergraph_type::node_type node_type;
    typedef hypergraph_type::edge_type edge_type;
    typedef hypergraph_type::rule_type rule_type;
    
    typedef rule_type::vocab_type      vocab_type;
    typedef rule_type::symbol_type     symbol_type;
    typedef rule_type::symbol_set_type symbol_set_type;
    
    
    typedef Function function_type;
    typedef typename function_type::value_type weight_type;

    typedef symbol_type word_type;

    typedef std::vector<word_type, std::allocator<word_type> > buffer_type;
    
    typedef symbol_set_type context_type;
    
    typedef std::pair<context_type, context_type> context_pair_type;
    typedef std::vector<context_pair_type, std::allocator<context_pair_type> > context_pair_set_type;

    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;

    struct ExtractSource
    {
      template <typename Edge>
      const context_type& operator()(const Edge& edge) const
      {
	return edge.rule->source;
      }
    };

    struct ExtractTarget
    {
      template <typename Edge>
      const context_type& operator()(const Edge& edge) const
      {
	return edge.rule->target;
      }
    };
    
  public:
    ExpectedNGram(Function __function,
		  const int __order,
		  const bool __bos_eos=false)
      : function(__function),
	order(__order),
	bos_eos(__bos_eos) {}

    
    template <typename Iterator>
    void collect_counts(Iterator first, Iterator iter, Iterator last, const weight_type& weight, Counts& counts)
    {
      const int context_size = order - 1;
      
      first = std::max(iter - context_size, first);
	
      for (/**/; first != iter; ++ first)
	for (Iterator iter2 = iter; iter2 != std::min(first + order, last); ++ iter2)
	  counts[typename Counts::key_type(first, iter2 + 1)] += weight;
    }


    template <typename Iterator>
    void collect_counts(Iterator first, Iterator last, const weight_type& weight, Counts& counts)
    {
      for (/**/; first != last; ++ first)
	for (Iterator iter = first; iter != std::min(first + order, last); ++ iter)
	  counts[typename Counts::key_type(first, iter + 1)] += weight;
    }
    
    template <typename Extract>
    void operator()(const hypergraph_type& graph, Counts& counts, Extract extract)
    {
      // we assume that the graph is already disambiguated in terms of ngram contexts..
      // if not, you should generate graph with additional feature-function, "variational" or "mbr"

      const int context_size = order - 1;

      buffer_type buffer;
      
      context_pair_set_type contexts(graph.nodes.size());
      weight_set_type       inside(graph.nodes.size());
      weight_set_type       outside(graph.nodes.size());
      
      cicada::inside(graph, inside, function);
      cicada::outside(graph, inside, outside, function);

      const weight_type weight_goal = inside[graph.goal];
      
      hypergraph_type::node_set_type::const_iterator niter_end = graph.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = graph.nodes.begin(); niter != niter_end; ++ niter) {
	const node_type& node = *niter;
	const bool is_goal = (node.id == graph.goal);
	
	node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
	for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	  const edge_type& edge = graph.edges[*eiter];
	  
	  weight_type weight = function(edge) * inside[node.id] / weight_goal;
	  edge_type::node_set_type::const_iterator niter_end = edge.tails.end();
	  for (edge_type::node_set_type::const_iterator niter = edge.tails.begin(); niter != niter_end; ++ niter)
	    weight *= outside[*niter];
	  
	  buffer.clear();
	  const context_type& context = extract(edge);
	  
	  if (edge.tails.empty()) {
	    context_type::const_iterator citer_end = context.end();
	    for (context_type::const_iterator citer = context.begin(); citer != citer_end; ++ citer)
	      if (*citer != vocab_type::EPSILON)
		buffer.push_back(*citer);
	    
	    collect_counts(buffer.begin(), buffer.end(), weight, counts);
	    
	    if (contexts[node.id].first.empty()) {
	      if (buffer.size() <= context_size)
		contexts[node.id] = std::make_pair(context_type(buffer.begin(), buffer.end()),
						   context_type());
	      else
		contexts[node.id] = std::make_pair(context_type(buffer.begin(), buffer.begin() + context_size),
						   context_type(buffer.end() - context_size, buffer.end()));
	    }
	    
	    if (is_goal && bos_eos) {
	      buffer.insert(buffer.begin(), vocab_type::BOS);
	      buffer.insert(buffer.end(), vocab_type::EOS);
	      
	      collect_counts(buffer.begin(), buffer.begin() + 1, buffer.end(), weight, counts);
	      collect_counts(buffer.begin(), buffer.end() - 1, buffer.end(), weight, counts);
	    }
	    
	  } else {
	    buffer.reserve(context.size() + edge.tails.size() * order * 2);
	    
	    int star_first = -1;
	    int star_last  = -1;

	    buffer_type::iterator biter_first = buffer.begin();
	    buffer_type::iterator biter       = buffer.begin();
	    
	    int non_terminal_pos = 0;
	    context_type::const_iterator citer_end = context.end();
	    for (context_type::const_iterator citer = context.begin(); citer != citer_end; ++ citer) {
	      if (citer->is_non_terminal()) {
		
		// collect ngram counts
		collect_counts(biter_first, biter, buffer.end(), weight, counts);
		biter = buffer.end();
		
		int antecedent_index = citer->non_terminal_index() - 1;
		if (antecedent_index < 0)
		  antecedent_index = non_terminal_pos;
		++ non_terminal_pos;

		const context_pair_type& context_pair = contexts[edge.tails[antecedent_index]];
		
		buffer.insert(buffer.end(), context_pair.first.begin(), context_pair.first.end());
		
		// collect ngram counts
		collect_counts(biter_first, biter, buffer.end(), weight, counts);
		biter = buffer.end();
		
		if (! context_pair.second.empty()) {
		  biter_first = buffer.end();
		  
		  star_last = buffer.size();
		  if (star_first < 0)
		    star_first = buffer.size();
		  
		  buffer.insert(buffer.end(), context_pair.second.begin(), context_pair.second.end());
		}
		
	      } else if (*citer != vocab_type::EPSILON)
		buffer.push_back(*citer);
	    }
	    
	    if (biter != buffer.end()) {
	      collect_counts(biter_first, biter, buffer.end(), weight, counts);
	      biter = buffer.end();
	    }
	    
	    if (contexts[node.id].first.empty()) {
	      
	      if (star_first >= 0) {
		const int prefix_size = utils::bithack::min(star_first, context_size);
		const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
		
		contexts[node.id]= std::make_pair(context_type(buffer.begin(), buffer.begin() + prefix_size),
						  context_type(buffer.end() - suffix_size, buffer.end()));
		
	      } else {
		if (buffer.size() <= context_size)
		  contexts[node.id] = std::make_pair(context_type(buffer.begin(), buffer.end()),
						     context_type());
		else
		  contexts[node.id] = std::make_pair(context_type(buffer.begin(), buffer.begin() + context_size),
						     context_type(buffer.end() - context_size, buffer.end()));
	      }
	    }
	    
	    if (is_goal && bos_eos) {
	      if (star_first >= 0) {
		const int prefix_size = utils::bithack::min(star_first, context_size);
		const int suffix_size = utils::bithack::min(int(buffer.size() - star_last), context_size);
		
		buffer.insert(buffer.begin(), vocab_type::BOS);
		buffer.insert(buffer.end(), vocab_type::EOS);

		collect_counts(buffer.begin(), buffer.begin() + 1, buffer.begin() + 1 + prefix_size, weight, counts);
		collect_counts(buffer.end() - suffix_size - 1, buffer.end() - 1, buffer.end(), weight, counts);
	      } else {
		buffer.insert(buffer.begin(), vocab_type::BOS);
		buffer.insert(buffer.end(), vocab_type::EOS);

		collect_counts(buffer.begin(), buffer.begin() + 1, buffer.end(), weight, counts);
		collect_counts(buffer.begin(), buffer.end() - 1, buffer.end(), weight, counts);
	      }
	    }
	  }
	}
      }
      
    }
    
    function_type function;
    int  order;
    bool bos_eos;
  };
  
  
  template <typename Function, typename Counts>
  void expected_ngram(const HyperGraph& graph, Function function, Counts& counts, const int order, const bool bos_eos=false, const bool yield_source=false)
  {
    ExpectedNGram<Function, Counts> __expected(function, order, bos_eos);
    
    if (yield_source)
      __expected(graph, counts, typename ExpectedNGram<Function, Counts>::ExtractSource());
    else
      __expected(graph, counts, typename ExpectedNGram<Function, Counts>::ExtractTarget());
  }
  
};

#endif
