
#include <map>

#include "cicada/feature_bleu.hpp"
#include "cicada/parameter.hpp"

#include "utils/hashmurmur.hpp"
#include "utils/compact_trie.hpp"
#include "utils/indexed_set.hpp"

namespace cicada
{
  namespace feature
  {

    class BleuImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      // this implementation specific...
      typedef uint32_t id_type;
      typedef uint32_t count_type;

      typedef symbol_type word_type;
      
      typedef std::allocator<std::pair<const word_type, count_type> >  ngram_allocator_type;
      typedef utils::compact_trie<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>, ngram_allocator_type> ngram_set_type;
      
      struct Node
      {
	Node() : word(), parent(id_type(-1)), order(0) {}
	
	word_type word;
	id_type   parent;
	int       order;
      };
      typedef Node node_type;
      typedef std::vector<node_type, std::allocator<node_type> > node_set_type;
      typedef std::vector<int, std::allocator<int> > size_set_type;
    
      typedef utils::simple_vector<count_type, std::allocator<count_type> > count_set_type;
      
      struct count_set_hash : public utils::hashmurmur<size_t>
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	size_t operator()(const count_set_type& x) const
	{
	  return hasher_type::operator()(x.begin(), x.end(), 0);
	}
      };
      
      typedef utils::indexed_set<count_set_type, count_set_hash, std::equal_to<count_set_type>, std::allocator<count_set_type> > states_count_set_type;


    public:
      BleuImpl(const int __order,
	       const bool __exact)
	: order(__order), exact(__exact) {}
      
      void clear()
      {
	// for ngram handling
	ngrams.clear();
	nodes.clear();
	sizes.clear();
	
	// for count set representation
	states_counts.clear();
	
	source_size = 0;
      }
      
      void insert(const int __source_size, const sentence_type& sentence)
      {
	typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
	source_size = __source_size;
	
	counts_type counts;
	sentence_type::const_iterator siter_end = sentence.end();
	for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  ngram_set_type::id_type id = ngrams.root();
	  
	  int n = 1;
	  for (sentence_type::const_iterator iter = siter; iter != std::min(siter + order, siter_end); ++ iter, ++ n) {
	    const ngram_set_type::id_type id_next = ngrams.insert(id, *iter);
	    
	    ++ counts[id_next];
	    
	    if (id_next >= nodes.size())
	      nodes.resize(id_next + 1, node_type());
	    
	    nodes[id_next].word = *iter;
	    nodes[id_next].parent = id;
	    nodes[id_next].order  = n;
	    
	    id = id_next;
	  }
	}
	
	// collect clipped ngram counts
	counts_type::const_iterator citer_end = counts.end();
	for (counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
	  ngrams[citer->first] = utils::bithack::max(ngrams[citer->first], citer->second);
	
	// keep sizes...
	sizes.push_back(sentence.size());
	std::sort(sizes.begin(), sizes.end());
      }
      
    private:
      ngram_set_type ngrams;
      node_set_type  nodes;
      size_set_type  sizes;
      
      states_count_set_type states_counts;
      
      int source_size;

      int order;
      bool exact;
    };
    
    inline
    bool true_false(const std::string& token)
    {
      if (strcasecmp(token.c_str(), "true") == 0)
	return true;
      if (strcasecmp(token.c_str(), "yes") == 0)
	return true;
      if (atoi(token.c_str()) > 0)
	return true;
      return false;
    }
    
    Bleu::Bleu(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "bleu")
	throw std::runtime_error("this is not Bleu feature: " + parameter);
      
      const int order = (param.find("order") != param.end() ? boost::lexical_cast<int>(param.find("order")->second) : 4);
      const bool exact = (param.find("exact") != param.end() ? true_false(param.find("exact")->second) : false);
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, exact));
      
      // two-side context + length (hypothesis/reference) + counts-id (hypothesis/reference)
      base_type::__state_size = sizeof(symbol_type) * order * 2 + sizeof(int) * 2 + sizeof(impl_type::id_type) * 2;
      base_type::__feature_name = "bleu";
      
      pimpl = bleu_impl.release();
    }
    
    Bleu::~Bleu() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    void Bleu::operator()(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates) const
    {
      
    }
    
    void Bleu::operator()(state_ptr_type& state,
			  feature_set_type& features) const
    {
      
    }
    
    void Bleu::clear()
    {
      
    }
    
    void Bleu::insert(const int source_size, const sentence_type& sentence)
    {
      
      
    }

  };
};
