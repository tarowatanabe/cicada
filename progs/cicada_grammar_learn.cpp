//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// learn latent annotation grammar from treebank
//
// First, we will implement via EM-algorithm
// 1. read all the parse tree
// 2. left-binarization (or, all/right binarization...?)
// 3. initialize table by maximum-likelihood estimates...
// 5. Iterate
//    6.  Split
//    7.  EM-iterations
//    8.  Merge
//    9.  EM-iterations
//

#include <stdexcept>
#include <vector>
#include <deque>

#include <cicada/hypergraph.hpp>
#include <cicada/inside_ouside.hpp>
#include <cicada/semiring.hpp>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>

#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/compress_stream.hpp>
#include <utils/resource.hpp>

typedef cicada::HyperGraph hypergraph_type;
typedef hypergraph_type::rule_type     rule_type;
typedef hypergraph_type::rule_ptr_type rule_ptr_type;

typedef rule_type::symbol_type     symbol_type;
typedef rule_type::symbol_set_type symbol_set_type;

typedef hypergraph_type::feature_set_type   feature_set_type;
typedef hypergraph_type::attribute_set_type attribute_set_type;

typedef feature_set_type::feature_type     feature_type;
typedef attribute_set_type::attribute_type attribute_type;

typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

// use of google dense_map for holding const rule_type*, not rule_ptr_type!

template <typename Tp>
struct ptr_hash
{
  size_t operator()(const Tp* x) const
  {
    return boost::hash_value(*x);
  }
  
  size_t operator()(const boost::shared_ptr<Tp>* x) const
  {
    return boost::hash_value(*x);
  }

};

template <typename Tp>
struct ptr_equal
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x == y || *x == *y;
  }
  
  bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
  {
    return x == y || *x == *y;
  }
};

typedef double count_type;
typedef double prob_type;
typedef double logprob_type;

typedef cicada::semiring::Logprob<double> weight_type;

class Grammar : public google::dense_hash_map<rule_ptr_type, count_type, boost::hash<rule_ptr_type>, std::equal_to<rule_ptr_type> >
{
public:
  typedef google::dense_hash_map<rule_ptr_type, count_type, boost::hash<rule_ptr_type>, std::equal_to<rule_ptr_type> > count_set_type;
  
public:
  Grammar() : count_set_type() { count_set_type::set_empty_key(rule_ptr_type()); }
};

typedef Grammar grammar_type;

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<symbol_type, grammar_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, grammar_type> > > count_set_type;
#else
  typedef sgi::hash_map<symbol_type, grammar_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			std::allocator<std::pair<const symbol_type, grammar_type> > > count_set_type;
#endif


path_set_type input_files;

int max_iteration = 6;         // max split-merge iterations
int max_iteration_split = 20;  // max EM-iterations for split
int max_iteration_merge = 20;  // max EM-iterations for merge

bool binarize_left = false;
bool binarize_right = false;
bool binarize_all = false;

// naive variational bayes for smoothing... otherwise, dirichlet prior
bool variational_bayes_mode = false;

double prior          = 0.01;

double merge_ratio = 0.5;
double cutoff_threshold = 1e-10;
int    cutoff_unk_terminal = 2;

int threads = 1;

int debug = 0;

int main(int argc, char** argv)
{
  grammar_type grammar;
  
}

symbol_type annotate_symbol(const symbol_type& symbol, const int bitpos, const bool bit=true)
{
  if (symbol.is_non_terminal()) {
    namespace xpressive = boost::xpressive;
    
    typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
    typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
    
    static const symbol_type root("[ROOT]");
    
    if (symbol == root) return root;
    
    static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (xpressive::s2= -+xpressive::_d);
    
    const utils::piece piece = sumbol.non_terminal_strip();
    const int mask = 1 << bitpos;
    
    pmatch what;
    if (xpressive::regex_math(piece), what, re) {
      const int value = (utils::lexical_cast<int>(what[2]) & (~mask)) | (-bit & mask);
      return '[' + what[1] + '@' + utils::lexical_cast<std::string>(value) + ']';
    } else
      return '[' + piece + '@' + utils::lexical_cast<std::string>(-bit & mask) + ']';
  } else
    return symbol;
}

struct attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
{
  attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
  attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
  attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
};

void grammar_merge(hypergraph_set_type& treebanks, gramamr_type& grammar, const int bits)
{
  // compute "loss" incurred by splitting...

  typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
  
  weight_set_type inside;
  weight_set_type outside;

  static const attribute_type attr_node("node");
  
  hypergraph_set_type::const_iterator titer_end = treebanks.end();
  for (hypergraph_set_type::const_iterator titer = treebanks.begin(); titer != titer_end; ++ titer) {
    const hypergraph_type& treebank = *titer;
    
    inside.clear();
    outside.clear();
    inside.resize(treebank.nodes.size());
    outside.resize(treebank.nodes.size());
    
    cicada::inside(treebank, inside, weight_function(grammar));
    cicada::outside(treebank, inside, outside, weight_function(grammar));

    hypergraph_type::node_set_type::const_iterator niter_end = treebank.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = treebank.nodes.begin(); niter != niter_end; ++ niter) {
      const hypergraph_type::node_type& node = *niter;
      const hypergraph_type::edge_type& edge = treebank.edges[node.edges.front()];
      
      const symbol_type lhs = edge.rule->lhs;
      
      attribute_set_type::const_iterator aiter = edge.attributes.find(attr_node);
      if (aiter == edge.attributes.end())
	throw std::runtime_error("no node attribute?");
      const int node_id_prev = boost::apply_visitor(attribute_integer(), iter->second);
      if (node_id_prev < 0)
	throw std::runtime_error("invalid node attribute?");
      
      const weight_type weight = inside[node.id] * outside[node.id];
      
    }
  }
  
  // given the loss, ....
}

void grammar_split(hypergraph_set_type& treebanks, gramamr_type& grammar, const int bits)
{
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<rule_ptr_type, ptr_hash<rule_type>, ptr_equal<rule_type>, std::allocator<const rule_ptr_type> > rule_set_type;
#else
  typedef sgi::hash_set<rule_ptr_type, ptr_hash<rule_type>, ptr_equal<rule_type>, std::allocator<rule_ptr_type> > rule_set_type;
#endif

  typedef std::vector<int, std::allocator<int> > index_set_type;
  typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, hypergraph_type::id_type> > > node_set_type;
#else
  typedef sgi::hash_map<symbol_type, hypergraph_type::id_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			std::allocator<std::pair<const symbol_type, hypergraph_type::id_type> > > node_set_type;
#endif
  typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;
  
  // split treebanks...
  // we will control by "bits"
    
  rule_set_type rules;
  node_map_type node_map;
  hypergraph_type treebank_new;

  index_set_type  index;
  index_set_type  index_end;
  symbol_set_type symbols;
  symbol_set_type symbols_new;

  static const symbol_type root("[ROOT]");
  static const attribute_type attr_node("node");
  
  // construct treebanks...
  hypergraph_set_type::iterator titer_end = treebanks.end();
  for (hypergraph_set_type::iterator titer = treebanks.begin(); titer != titer_end; ++ titer) {
    hypergraph_type& treebank = *titer;
    treebank_new.clear();
    
    //
    // we will create node for original node-id + new-symbol
    //
    
    node_map.clear();
    node_map.resize(treebank.nodes.size());
    
    hypergraph_type::node_set_type::const_iterator niter_end = treebank.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = treebank.nodes.begin(); niter != niter_end; ++ niter) {
      const hypergraph_type::node_type& node = *niter;
      
      hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = treebank.edges[*eiter];
	const rule_ptr_type& rule = edge.rule;
	
	symbols.clear();
	symbols.push_back(rule->lhs);
	symbols.insert(symbols.end(), rule->rhs.begin(), rule->rhs.end());
	
	symbols_new.clear();
	symbols_new.insert(symbols_new.end(), symbols.begin(), symbols.end());
	
	j.clear();
	j.resize(rule->rhs.size() + 1, 0);
	j_end.resize(rule->rhs.size() + 1);
	
	hypergraph_type::edge_type::node_set_type tails(edge.tails.size());
	
	for (size_t i = 0; i != symbols.size(); ++ i)
	  j_end[i] = utils::bithack::branch(symbols[i] == root, 1, utils::bithack::branch(symbols[i].is_non_terminal(), 2, 0));
	
	for (;;) {
	  // construct rule
	  for (size_t i = 0; i != symbols.size(); ++ i)
	    if (j_end[i])
	      symbols_new[i] = annotate_symbol(symbols[i], bits, j[i]);
	  
	  const rule_ptr_type rule = *rules.insert(rule_type::create(rule_type(symbols_new.front(), symbols_new.begin() + 1, symbols_new.end()))).first;
	  
	  // construct edge
	  std::pair<node_set_type::iterator, bool> head = node_map[edge.head].insert(std::make_pair(symbols_new.front(), 0));
	  if (head.second) {
	    head.first->second = treebank_new.add_node().id;
	    
	    // handling goal... assuming penntreebank style...
	    if (node.id == treebank.goal)
	      treebank_new.goal = head.first->second;
	  }
	  
	  size_t pos = 0;
	  for (size_t i = 1; i != symbols_new.size(); ++ i)
	    if (j_end[i]) {
	      node_set_type::const_iterator iter = node_map[edge.tails[pos]].find(symbols_new[i]);
	      if (iter == node_map[edge.tails[pos]].end())
		throw std::runtime_error("invalid node...?");
	      
	      tails[pos] = iter->second;
	      ++ pos;
	    }
	  
	  hypergraph_type::edge_type& edge_new = treebank_new.add_edge(tails.begin(), tails.end());
	  edge_new.rule = rule;
	  edge_new.attributes[attr_node] = attribute_set_type::int_type(edge.head);
	  
	  treebank_new.connect_edge(edge_new.id, head.first->second);
	  
	  size_t index = 0;
	  for (/**/; index != j.size(); ++ index) 
	    if (j_end[index]) {
	      ++ j[index];
	      if (j[index] < j_end[index]) break;
	      j[index] = 0;
	    }
	  
	  if (index == j.size()) break;
	}
      }   
    }

    if (treebank_new.is_valid())
      treebank_new.topologically_sort();
    
    treebank.swap(treebank_new);
    treebank_new.clear();
  }
  
  // split grammar...
  count_set_type counts;
  
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter) {
    const rule_ptr_type& rule = giter->first;
    
    symbols.clear();
    symbols.push_back(rule->lhs);
    symbols.insert(symbols.end(), rule->rhs.begin(), rule->rhs.end());

    symbols_new.clear();
    symbols_new.insert(symbols_new.end(), symbols.begin(), symbols.end());
    
    j.clear();
    j.resize(rule->rhs.size() + 1, 0);
    j_end.resize(rule->rhs.size() + 1);
    
    for (size_t i = 0; i != symbols.size(); ++ i)
      j_end[i] = utils::bithack::branch(symbols[i] == root, 1, utils::bithack::branch(symbols[i].is_non_terminal(), 2, 0));
    
    for (;;) {
      for (size_t i = 0; i != symbols.size(); ++ i)
	if (j_end[i])
	  symbols_new[i] = annotate_symbol(symbols[i], bits, j[i]);
      
      const rule_ptr_type rule = *rules.insert(rule_type::create(rule_type(symbols_new.front(), symbols_new.begin() + 1, symbols_new.end()))).first;
      
      counts[rule->lhs][rule] += utils::mathop::exp(giter->second);
      
      size_t index = 0;
      for (/**/; index != j.size(); ++ index) 
	if (j_end[index]) {
	  ++ j[index];
	  if (j[index] < j_end[index]) break;
	  j[index] = 0;
	}
      
      if (index == j.size()) break;
    }
  }
  
  // maximization
  if (variational_bayes_mode)
    maximize_gramamr(counts, grammar, MaximizeBayes());
  else
    maximize_grammar(counts, grammar, Maximize());
}

struct zero_function
{
  template <typename Edge>
  weight_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<weight_type>::exp(0.0);
  }
  
  const gramamr_type& grammar;
};

struct weight_function
{
  wegiht_function(const grammar_type& __grammar) : grammar(__grammar) {}
  
  template <typename Edge>
  weight_type operator()(const Edge& edge) const
  {
    gramamr_type::const_iterator giter = grammar.find(edge.rule);
    if (giter == grammar.end())
      throw std::runtime_error("invalid rule");

    return cicada::semiring::traits<weight_type>::exp(giter->second);
  }
  
  const gramamr_type& grammar;
};

struct accumulator_type
{
  accumulator_type(const hypergraph_type& __treebank,
		   count_set_type& __counts)
    : treebank(__treebank), counts(__counts) {}
  
  count_type& operator[](const hypergraph_type::id_type& x)
  {
    const rule_ptr_type& rule = treebank.edges[x].rule;
    
    return counts[rule->lhs][rule];
  }
  
  const hypergraph_type& treebank;
  count_set_type& counts;
};

template <typename Function>
void grammar_learn(const hypergraph_set_type& treebanks, gramamr_type& grammar, Function function)
{
  // TODO: support threading!
  
  typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
  
  weight_set_type inside;
  weight_set_type outside;

  count_set_type counts;

  weight_type logprob(cicada::semiring::traits<weight_type>::one());
  
  // expectation
  hypergraph_set_type::const_iterator titer_end = treebanks.end();
  for (hypergraph_set_type::const_iterator titer = treebanks.begin(); titer != titer_end; ++ titer) {
    const hypergraph_type& treebank = *titer;
    
    inside.clear();
    outside.clear();
    inside.resize(treebank.nodes.size());
    outside.resize(treebank.nodes.size());
    
    accumulator_type accumulator(treebank, counts);
    
    cicada::inside_outside(treebank, inside, outside, accumulator, function, function);
    
    if (debug >= 2)
      std::cerr << "inside: " << cicada::semiring::log(inside.back()) << std::endl;
    
    logprob *= inside.back();
  }
  
  if (debug)
    std::cerr << "log-likelihood: " << cicada::semiring::log(logprob) << std::endl;
  
  // maximization
  if (variational_bayes_mode)
    maximize_gramamr(counts, grammar, MaximizeBayes());
  else
    maximize_grammar(counts, grammar, Maximize());
}


struct Maximize
{
  void operator()(const grammar_type& counts, grammar_type& grammar) const
  {
    double sum = 0.0;
    grammar_type::const_iterator citer_end = counts.end();
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      sum += citer->second + prior;
    
    const double logsum = utils::mathop::log(sum);
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      grammar[citer->first] = utils::mathop::log(citer->second + prior) - logsum;
  }
};

struct MaximizeBayes
{
  void operator()(const grammar_type& counts, grammar_type& grammar) const
  {
    double sum = 0.0;
    grammar_type::const_iterator citer_end = counts.end();
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      sum += citer->second + prior;
    
    const double logsum = utils::mathop::digamma(sum);
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      grammar[citer->first] = utils::mathop::digamma(citer->second + prior) - logsum;
  }  
};

template <typename Maximizer>
void maximize_grammar(const count_set_type& counts, grammar_type& grammar, Maximizer maximizer)
{
  grammar.clear();
  
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
    maximize(citer->second, grammar);
}

void read_treebank(const path_set_type& files, hypergraph_set_type& treebanks)
{
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<rule_ptr_type, ptr_hash<rule_type>, ptr_equal<rule_type>, std::allocator<const rule_ptr_type> > rule_set_type;
#else
  typedef sgi::hash_set<rule_ptr_type, ptr_hash<rule_type>, ptr_equal<rule_type>, std::allocator<rule_ptr_type> > rule_set_type;
#endif
  
  hypergrap_type treebank;
  rule_set_type rules;
  
  path_set_type::const_iterator fiter_end = files.end();
  for (path_set_type::const_iterator fiter = files.begin(); fiter != fiter_end; ++ fiter) {
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (is >> treebank) {
      if (binarize_left)
	cicada::binarize_left(treebank, 0);
      else if (binarize_right)
	cicada::binarize_right(treebank, 0);
      else if (binarize_all)
	cicada::binarize_all(treebank);
      
      // use of unique rules...
      hypergraph_type::edge_set_type::iterator eiter_end = treebank.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = treebank.edges.begin(); eiter != eiter_end; ++ eiter)
	eiter->rule = *(rules.insert(eiter->rule).first);
      
      treebanks.push_back(hypergraph_type());
      treebanks.back().swap(treebank);
    }
  }
}
