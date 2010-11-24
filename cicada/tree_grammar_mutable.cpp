
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "tree_grammar_mutable.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/compact_trie_dense_set.hpp"
#include "utils/compress_stream.hpp"

#include "utils/bithack.hpp"
#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/tempfile.hpp"
#include "utils/group_aligned_code.hpp"
#include "utils/byte_aligned_code.hpp"
#include "utils/simple_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/arc_list.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/filesystem.hpp>

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  struct TreeGrammarMutableImpl : public utils::hashmurmur<uint64_t>
  {
    friend class TreeGrammarMutable;

    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol  symbol_type;
    typedef cicada::Symbol  word_type;
    typedef cicada::Feature feature_type;
    typedef cicada::Vocab   vocab_type;
    
    
    typedef TreeTransducer::rule_type          rule_type;
    typedef TreeTransducer::rule_ptr_type      rule_ptr_type;
    typedef TreeTransducer::rule_pair_type     rule_pair_type;
    typedef TreeTransducer::rule_pair_set_type rule_pair_set_type;
    
    typedef TreeTransducer::feature_set_type feature_set_type;
    
    typedef utils::compact_trie_dense_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					  std::allocator<symbol_type > > edge_trie_type;
    typedef edge_trie_type::id_type edge_id_type;
    
    typedef utils::compact_trie_dense<edge_id_type, rule_pair_set_type, boost::hash<edge_id_type>, std::equal_to<edge_id_type>,
				      std::allocator<std::pair<const edge_id_type, rule_pair_set_type> > > trie_type;
    typedef trie_type::id_type id_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
    
    TreeGrammarMutableImpl(const std::string& parameter)
      : trie(edge_id_type(-1)), edges(symbol_type()) { read(parameter); }
    
    edge_id_type edge(const symbol_type* first, const symbol_type* last) const
    {
      return edges.find(first, last);
    }
    
    id_type root() const { return trie.root(); }
    id_type next(id_type node, const edge_id_type& edge) const { return trie.find(node, edge); }
    bool has_next(id_type node) const { return ! trie.empty(node); }
    const rule_pair_set_type& rules(id_type node) { return trie[node]; }
    
    void insert(const std::string& pattern);
    void insert(const rule_pair_type& rule);
    
  public:
    void clear() { trie.clear();  edges.clear(); }
    void read(const std::string& parameter);
    
  public:
    trie_type      trie;
    edge_trie_type edges;
  };
  
  //
  // we will parse only scores part and rely on tree-rule to parse partial string...
  //
  
  
  void TreeGrammarMutableImpl::read(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
	  
    const parameter_type param(parameter);
    
    const path_type path = param.name();
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no grammar file") + param.name());
    
    
    
  }
  

  void TreeGrammarMutableImpl::insert(const std::string& pattern)
  {
    
    
  }
  
  inline
  int tree_depth(const TreeRule& tree, const int depth)
  {
    // pre-order traversal
    int max_depth = depth;
    for (TreeRule::const_iterator aiter = tree.begin(); aiter != tree.end(); ++ aiter)
      max_depth = utils::bithack::max(max_depth, tree_depth(*aiter, depth + 1));
    
    return max_depth;
  }

  inline
  void tree_add_epsilon(TreeRule& tree, const int max_depth, const int depth)
  {
    // pre-order traversal
    if (tree.empty()) {
      TreeRule* curr = &tree;
      
      for (int i = depth; i != max_depth; ++ i) {
	curr->antecedents = TreeRule::antecedent_set_type(1, TreeRule(Vocab::EPSILON));
	curr = &curr->front();
      }
    } else {
      for (TreeRule::iterator aiter = tree.begin(); aiter != tree.end(); ++ aiter)
	tree_add_epsilon(*aiter, max_depth, depth + 1);
    }
  }

  template <typename Path>
  inline
  void tree_to_hyperpath(const TreeRule& tree, Path& path, const int depth)
  {
    path[depth].push_back(tree.label);

    if (! tree.empty()) {
      for (TreeRule::const_iterator aiter = tree.begin(); aiter != tree.end(); ++ aiter)
	tree_to_hyperpath(*aiter, path, depth + 1);
      path[depth + 1].push_back(Vocab::NONE);
    }
  }
  
  template <typename Path>
  inline
  void tree_to_hyperpath(const TreeRule& tree, Path& path)
  {
    // compute max-depth by pre-order
    const int max_depth = tree_depth(tree, 0);
    
    // add epsilon annotation by pre-order
    TreeRule tree_epsilon(tree);
    tree_add_epsilon(tree_epsilon, max_depth, 0);
    
    // convert into hyperpath by pre-order
    path.clear();
    path.resize(max_depth + 1);
    
    tree_to_hyperpath(tree_epsilon, path, 0);
    path[0].push_back(Vocab::NONE);
  }
  
  void TreeGrammarMutableImpl::insert(const rule_pair_type& rule_pair)
  {
    
    
    
  }

  TreeGrammarMutable::TreeGrammarMutable(const std::string& parameter)
    : pimpl(new impl_type(parameter)) {}
  
  TreeGrammarMutable::~TreeGrammarMutable() { std::auto_ptr<impl_type> tmp(pimpl); }
  
  TreeGrammarMutable::TreeGrammarMutable(const TreeGrammarMutable& x)
    : pimpl(new impl_type(*x.pimpl)) {}

  TreeGrammarMutable& TreeGrammarMutable::operator=(const TreeGrammarMutable& x)
  {
    *pimpl = *x.pimpl;
    return *this;
  }
  
  TreeGrammarMutable::transducer_ptr_type TreeGrammarMutable::clone() const
  {
    return transducer_ptr_type(new TreeGrammarMutable(*this));
  }
  
  TreeGrammarMutable::edge_id_type TreeGrammarMutable::edge(const symbol_type& symbol) const
  {
    return edge(&symbol, (&symbol) + 1);
  }
  
  TreeGrammarMutable::edge_id_type TreeGrammarMutable::edge(const symbol_set_type& symbols) const
  {
    return edge(&(*symbols.begin()), &(*symbols.end()));
  }
  
  TreeGrammarMutable::edge_id_type TreeGrammarMutable::edge(const symbol_type* first, const symbol_type* last) const
  {
    impl_type::edge_id_type id = pimpl->edge(first, last);
    
    return (id != pimpl->edges.npos() ? id : edge_id_type(-1));
  }


  TreeGrammarMutable::id_type TreeGrammarMutable::root() const
  {
    return pimpl->root();
  }
  
  TreeGrammarMutable::id_type TreeGrammarMutable::next(const id_type& node, const edge_id_type& edge) const
  {
    return pimpl->next(node, edge);
  }
  
  bool TreeGrammarMutable::has_next(const id_type& node) const
  {
    return pimpl->has_next(node);
  }
  
  const TreeGrammarMutable::rule_pair_set_type& TreeGrammarMutable::rules(const id_type& node) const
  {
    static const rule_pair_set_type __empty;
    return (node == pimpl->root() ? __empty : pimpl->rules(node));
  }
  
   
  void TreeGrammarMutable::read(const std::string& parameter)
  {
    pimpl->read(parameter);
  }

  void TreeGrammarMutable::clear()
  {
    pimpl->clear();
  }
  
  void TreeGrammarMutable::insert(const std::string& pattern)
  {
    pimpl->insert(pattern);
  }

  void TreeGrammarMutable::insert(const rule_pair_type& rule_pair)
  {
    pimpl->insert(rule_pair);
  }
  
};
