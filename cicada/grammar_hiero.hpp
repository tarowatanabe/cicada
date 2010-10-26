// -*- mode: c++ -*-

#ifndef __CICADA__GRAMMAR_HIERO__HPP__
#define __CICADA__GRAMMAR_HIERO__HPP__ 1

// very siple mutable grammar class..

#include <string>
#include <vector>

#include <cicada/grammar_mutable.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  
  class GrammarGlue : public GrammarMutable
  {
  public:
    GrammarGlue(const symbol_type& goal, const symbol_type& non_terminal, const bool __straight, const bool __inverted)
      : straight(__straight), inverted(__inverted)
    {
      rule_ptr_type rule_unary(new rule_type(goal,
					     rule_type::symbol_set_type(1, non_terminal.non_terminal(1)),
					     rule_type::symbol_set_type(1, non_terminal.non_terminal(1)),
					     1));
      insert(rule_unary);

      if (straight) {
	std::vector<symbol_type, std::allocator<symbol_type> > phrase(2);
	phrase.front() = goal.non_terminal(1);
	phrase.back()  = non_terminal.non_terminal(2);
	
	rule_ptr_type rule(new rule_type(goal,
					 rule_type::symbol_set_type(phrase.begin(), phrase.end()),
					 rule_type::symbol_set_type(phrase.begin(), phrase.end()),
					 2));
	rule->features["glue-straight-penalty"] = -1;
	
	insert(rule);
      }
      
      if (inverted) {
	std::vector<symbol_type, std::allocator<symbol_type> > phrase1(2);
	std::vector<symbol_type, std::allocator<symbol_type> > phrase2(2);
	
	phrase1.front() = goal.non_terminal(1);
	phrase1.back()  = non_terminal.non_terminal(2);
	
	phrase2.front() = non_terminal.non_terminal(2);
	phrase2.back()  = goal.non_terminal(1);
	
	rule_ptr_type rule(new rule_type(goal,
					 rule_type::symbol_set_type(phrase1.begin(), phrase1.end()),
					 rule_type::symbol_set_type(phrase2.begin(), phrase2.end()),
					 2));
	rule->features["glue-inverted-penalty"] = -1;
	
	insert(rule);
      }
    }
    
    bool valid_span(int first, int last, int distance) const
    {
      // how to check inverted only...?
      if (straight && ! inverted)
	return first == 0;
      else
	return true;
    }
    
  private:
    bool straight;
    bool inverted;
  };

  class GrammarInsertion : public GrammarMutable
  {
  public:
    typedef Lattice    lattice_type;
    typedef HyperGraph hypergraph_type;

  private:
    typedef std::vector<bool, std::allocator<bool> > pos_set_type;
    typedef std::vector<pos_set_type, std::allocator<pos_set_type> > pos_pair_set_type;
    
  public:
    GrammarInsertion(const hypergraph_type& graph, const symbol_type& non_terminal)
    {
      typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > symbol_set_type;
      
      symbol_set_type symbols;
      symbols.set_empty_key(symbol_type());
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) 
	if (eiter->rule) {
	  const rule_type& rule = *(eiter->rule);
	  
	  rule_type::symbol_set_type::const_iterator siter_end = rule.source.end();
	  for (rule_type::symbol_set_type::const_iterator siter = rule.source.begin(); siter != siter_end; ++ siter) 
	    if (*siter != vocab_type::EPSILON && siter->is_terminal() && symbols.find(*siter) == symbols.end()) {
	      
	      rule_ptr_type rule(new rule_type(non_terminal,
					       rule_type::symbol_set_type(1, *siter),
					       rule_type::symbol_set_type(1, *siter)));
	      rule->features["insertion-penalty"] = - 1.0;
	      
	      insert(rule);
	      
	      symbols.insert(*siter);
	    }
	}
    }

    GrammarInsertion(const lattice_type& lattice, const symbol_type& non_terminal)
      : positions(lattice.size(), pos_set_type(lattice.size() + 1, false))
    {
      typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > symbol_set_type;
      
      symbol_set_type symbols;
      symbols.set_empty_key(symbol_type());

      for (size_t first = 0; first != lattice.size(); ++ first) {
	const lattice_type::arc_set_type& arcs = lattice[first];
	
	if (arcs.empty())
	  positions[first].clear();
	
	lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter) {
	  const size_t last = first + aiter->distance;
	  
	  positions[first][last] = true;

	  if (aiter->label != vocab_type::EPSILON && symbols.find(aiter->label) == symbols.end()) {
	  
	    rule_ptr_type rule(new rule_type(non_terminal,
					     rule_type::symbol_set_type(1, aiter->label),
					     rule_type::symbol_set_type(1, aiter->label)));
	    rule->features["insertion-penalty"] = - 1.0;
	    
	    insert(rule);
	    
	    symbols.insert(aiter->label);
	  }
	}
      }
    }

    bool valid_span(int first, int last, int distance) const
    {
      return positions.empty() || (! positions[first].empty() && (first == last || positions[first][last]));
    }

  private:    
    pos_pair_set_type positions;
  };
  
  
  class GrammarDeletion : public GrammarMutable
  {
  public:
    typedef Lattice    lattice_type;
    typedef HyperGraph hypergraph_type;

  private:
    typedef std::vector<bool, std::allocator<bool> > pos_set_type;
    typedef std::vector<pos_set_type, std::allocator<pos_set_type> > pos_pair_set_type;

  public:
    GrammarDeletion(const hypergraph_type& graph, const symbol_type& non_terminal)
    {
      typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > symbol_set_type;
      
      symbol_set_type symbols;
      symbols.set_empty_key(symbol_type());
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = graph.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = graph.edges.begin(); eiter != eiter_end; ++ eiter) 
	if (eiter->rule) {
	  const rule_type& rule = *(eiter->rule);
	  
	  rule_type::symbol_set_type::const_iterator siter_end = rule.source.end();
	  for (rule_type::symbol_set_type::const_iterator siter = rule.source.begin(); siter != siter_end; ++ siter) 
	    if (*siter != vocab_type::EPSILON && siter->is_terminal() && symbols.find(*siter) == symbols.end()) {
	      
	      rule_ptr_type rule(new rule_type(non_terminal,
					       rule_type::symbol_set_type(1, *siter),
					       rule_type::symbol_set_type(1, vocab_type::EPSILON)));
	      rule->features["deletion-penalty"] = - 1.0;
	      
	      insert(rule);
	      
	      symbols.insert(*siter);
	    }
	}
    }

    GrammarDeletion(const lattice_type& lattice, const symbol_type& non_terminal)
      : positions(lattice.size(), pos_set_type(lattice.size() + 1, false))
    {
      typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > symbol_set_type;
      
      symbol_set_type symbols;
      symbols.set_empty_key(symbol_type());
      
      for (size_t first = 0; first != lattice.size(); ++ first) {
	const lattice_type::arc_set_type& arcs = lattice[first];

	if (arcs.empty())
	  positions[first].clear();
		
	lattice_type::arc_set_type::const_iterator aiter_end = arcs.end();
	for (lattice_type::arc_set_type::const_iterator aiter = arcs.begin(); aiter != aiter_end; ++ aiter) {
	  const size_t last = first + aiter->distance;
	  
	  positions[first][last] = true;

	  if (aiter->label != vocab_type::EPSILON && symbols.find(aiter->label) == symbols.end()) {
	  
	    rule_ptr_type rule(new rule_type(non_terminal,
					     rule_type::symbol_set_type(1, aiter->label),
					     rule_type::symbol_set_type(1, vocab_type::EPSILON)));
	    rule->features["deletion-penalty"] = - 1.0;
	    
	    insert(rule);
	    
	    symbols.insert(aiter->label);
	  }
	}
      }
    }
    
    bool valid_span(int first, int last, int distance) const
    {
      return positions.empty() || (! positions[first].empty() && (first == last || positions[first][last]));
    }

  private:    
    pos_pair_set_type positions;
  };
};

#endif
