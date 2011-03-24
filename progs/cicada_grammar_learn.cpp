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
// TODO: Do we differentiate terminal-rules and rules?
// 

#include <stdexcept>
#include <vector>
#include <deque>

#include <cicada/hypergraph.hpp>
#include <cicada/inside_outside.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
#include <cicada/binarize.hpp>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/xpressive/xpressive.hpp>

#include <utils/bithack.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/compress_stream.hpp>
#include <utils/resource.hpp>
#include <utils/mathop.hpp>
#include <utils/lexical_cast.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/array_power2.hpp>

#include <google/dense_hash_map>
#include <google/dense_hash_set>

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
struct ptr_hash : public boost::hash<Tp>
{
  typedef boost::hash<Tp> hasher_type;

  size_t operator()(const Tp* x) const
  {
    return (x ? hasher_type::operator()(*x) : size_t(0));
  }
  
  size_t operator()(const boost::shared_ptr<Tp>& x) const
  {
    return (x ? hasher_type::operator()(*x) : size_t(0));
  }

};

template <typename Tp>
struct ptr_equal
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x == y || (x && y && *x == *y);
  }
  
  bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
  {
    return x == y || (x && y && *x == *y);
  }
};

typedef double count_type;
typedef double prob_type;
typedef double logprob_type;

typedef cicada::semiring::Logprob<double> weight_type;

class Grammar : public google::dense_hash_map<rule_ptr_type, count_type, ptr_hash<rule_type>, ptr_equal<rule_type> >
{
public:
  typedef google::dense_hash_map<rule_ptr_type, count_type, ptr_hash<rule_type>, ptr_equal<rule_type> > count_set_type;
  
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
path_type     output_grammar_file = "-";
path_type     output_unknown_file = "-";

int max_iteration = 6;         // max split-merge iterations
int max_iteration_split = 20;  // max EM-iterations for split
int max_iteration_merge = 20;  // max EM-iterations for merge

bool binarize_left = false;
bool binarize_right = false;
bool binarize_all = false;

// naive variational bayes for smoothing... otherwise, dirichlet prior
bool variational_bayes_mode = false;

double prior          = 0.01;
double prior_terminal = 0.01;

bool keep_terminal = false;
double merge_ratio = 0.5;

double cutoff_threshold = 1e-20;
int    cutoff_unknown = 10;

int threads = 1;

int debug = 0;

template <typename Generator>
void grammar_merge(hypergraph_set_type& treebanks, grammar_type& grammar, const int bits, Generator& generator);

template <typename Generator>
void grammar_split(hypergraph_set_type& treebanks, grammar_type& grammar, const int bits, Generator& generator);

template <typename Function>
double grammar_learn(const hypergraph_set_type& treebanks, grammar_type& grammar, Function function);

template <typename Maximizer>
void grammar_maximize(const count_set_type& counts, grammar_type& grammar, Maximizer maximizer);

void write_grammar(const path_type& file, const grammar_type& grammar);
void read_treebank(const path_set_type& files, hypergraph_set_type& treebanks);

void options(int argc, char** argv);

struct zero_function
{
  typedef weight_type value_type;
  
  template <typename Edge>
  weight_type operator()(const Edge& edge) const
  {
    return cicada::semiring::traits<weight_type>::exp(0.0);
  }
};

struct weight_function
{
  typedef weight_type value_type;

  weight_function(const grammar_type& __grammar) : grammar(__grammar) {}
  
  template <typename Edge>
  weight_type operator()(const Edge& edge) const
  {
    grammar_type::const_iterator giter = grammar.find(edge.rule);
    if (giter == grammar.end())
      throw std::runtime_error("invalid rule");

    return cicada::semiring::traits<weight_type>::exp(giter->second);
  }
  
  const grammar_type& grammar;
};

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (merge_ratio <= 0.0 || 1.0 <= merge_ratio)
      throw std::runtime_error("invalid merge ratio");
      
    if (int(binarize_left) + binarize_right + binarize_all > 1)
      throw std::runtime_error("specify either binarize-{left,right,all}");

    if (int(binarize_left) + binarize_right + binarize_all == 0)
      binarize_left = true;

    if (input_files.empty())
      input_files.push_back("-");
    
    threads = utils::bithack::max(threads, 1);
    
    hypergraph_set_type treebanks;
    read_treebank(input_files, treebanks);
    
    grammar_type grammar;
    grammar_learn(treebanks, grammar, zero_function());
    
    boost::mt19937 generator;
    generator.seed(time(0) * getpid());

    if (debug)
      std::cerr << "grammar size: " << grammar.size() << std::endl;
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      
      if (debug)
	std::cerr << "iteration: " << (iter + 1) << std::endl;
      
      // split...
      {
	const utils::resource split_start;
	grammar_split(treebanks, grammar, iter, generator);
	const utils::resource split_end;
	
	if (debug)
	  std::cerr << "split: " << "grammar size: " << grammar.size() << std::endl;

	if (debug)
	  std::cerr << "cpu time: " << (split_end.cpu_time() - split_start.cpu_time())
		    << " user time: " << (split_end.user_time() - split_start.user_time())
		    << std::endl;
	
	double logprob = 0.0;
	for (int i = 0; i < max_iteration_split; ++ i) {
	  if (debug)
	    std::cerr << "split iteration: " << (i + 1) << std::endl;
	  
	  const utils::resource learn_start;
	  const double logprob_curr = grammar_learn(treebanks, grammar, weight_function(grammar));
	  const utils::resource learn_end;

	  if (debug)
	    std::cerr << "cpu time: " << (learn_end.cpu_time() - learn_start.cpu_time())
		      << " user time: " << (learn_end.user_time() - learn_start.user_time())
		      << std::endl;
	  
	  if (i && logprob_curr < logprob) break;
	  
	  logprob = logprob_curr;
	}
      }
      
      // merge..
      {
	const utils::resource merge_start;
	grammar_merge(treebanks, grammar, iter, generator);
	const utils::resource merge_end;
	
	if (debug)
	  std::cerr << "merge: " << "grammar size: " << grammar.size() << std::endl;

	if (debug)
	  std::cerr << "cpu time: " << (merge_end.cpu_time() - merge_start.cpu_time())
		    << " user time: " << (merge_end.user_time() - merge_start.user_time())
		    << std::endl;
	
	double logprob = 0.0;
	for (int i = 0; i < max_iteration_merge; ++ i) {
	  if (debug)
	    std::cerr << "merge iteration: " << (i + 1) << std::endl;
	  
	  const utils::resource learn_start;
	  const double logprob_curr = grammar_learn(treebanks, grammar, weight_function(grammar));
	  const utils::resource learn_end;
	  
	  if (debug)
	    std::cerr << "cpu time: " << (learn_end.cpu_time() - learn_start.cpu_time())
		      << " user time: " << (learn_end.user_time() - learn_start.user_time())
		      << std::endl;
	  
	  if (i && logprob_curr < logprob) break;
	  
	  logprob = logprob_curr;
	}
      }
    }
    
    // output grammar...
    write_grammar(output_grammar_file, grammar);
    
    // output unknown handler...
    
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


symbol_type annotate_symbol(const symbol_type& symbol, const int bitpos, const bool bit)
{
  if (symbol.is_non_terminal()) {
    namespace xpressive = boost::xpressive;
    
    typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
    typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
    
    static const symbol_type root("[ROOT]");
    
    if (symbol == root) return root;
    
    static pregex re = (xpressive::s1= +(~xpressive::_s)) >> '@' >> (xpressive::s2= -+xpressive::_d);
    
    const utils::piece piece = symbol.non_terminal_strip();
    const int mask = 1 << bitpos;
    
    pmatch what;
    if (xpressive::regex_match(piece, what, re)) {
      const int value = (utils::lexical_cast<int>(what[2]) & (~mask)) | (-bit & mask);
      return '[' + what[1] + '@' + utils::lexical_cast<std::string>(value) + ']';
    } else
      return '[' + piece + '@' + utils::lexical_cast<std::string>(-bit & mask) + ']';
  } else
    return symbol;
}

struct Annotator : public utils::hashmurmur<size_t>
{
  typedef utils::hashmurmur<size_t> hasher_type;
  
  Annotator(const int __bits) : bits(__bits) {}
  
  struct Cache
  {
    symbol_type symbol;
    symbol_type annotated;
    bool bit;
    
    Cache() : symbol(), annotated(), bit(false) {}
  };
  typedef Cache cache_type;
  typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;
  
  const symbol_type& annotate(const symbol_type& symbol, const bool bit)
  {
    const size_t cache_pos = hasher_type::operator()(symbol.id(), bit) & (caches.size() - 1);
    cache_type& cache = caches[cache_pos];
    if (cache.symbol != symbol || cache.bit != bit) {
      cache.symbol = symbol;
      cache.bit = bit;
      cache.annotated = annotate_symbol(symbol, bits, bit);
    }
    return cache.annotated;
  }
  
  cache_set_type caches;
  const int bits;
};

struct attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
{
  attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
  attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
  attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
};

template <typename Tp>
struct greater_ptr_second
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

template <typename Tp>
struct less_ptr_second
{
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second < y->second;
  }
};

struct filter_pruned
{
  typedef std::vector<bool, std::allocator<bool> > removed_type;

  const removed_type& removed;
  
  filter_pruned(const removed_type& __removed) : removed(__removed) {}
  
  template <typename Edge>
  bool operator()(const Edge& edge) const
  {
    return removed[edge.id];
  }
};


struct Maximize
{
  void operator()(const grammar_type& counts, grammar_type& grammar) const
  {
    bool is_terminal = false;
    double sum = 0.0;
    grammar_type::const_iterator citer_end = counts.end();
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      sum += citer->second;
      is_terminal |= citer->first->rhs.front().is_terminal();
    }

    const double prior_lhs = (is_terminal ? prior_terminal : prior);
    const double logsum = utils::mathop::log(sum + prior_lhs * counts.size());
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      grammar[citer->first] = utils::mathop::log(citer->second + prior_lhs) - logsum;
  }
};

struct MaximizeBayes
{
  void operator()(const grammar_type& counts, grammar_type& grammar) const
  {
    bool is_terminal = false;
    double sum = 0.0;
    grammar_type::const_iterator citer_end = counts.end();
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      sum += citer->second;
      is_terminal |= citer->first->rhs.front().is_terminal();
    }
    
    const double prior_lhs = (is_terminal ? prior_terminal : prior);
    const double logsum = utils::mathop::digamma(sum + prior_lhs * counts.size());
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      grammar[citer->first] = utils::mathop::digamma(citer->second + prior_lhs) - logsum;
  }  
};

template <typename Scale>
struct TaskMergeScale
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  typedef Scale scale_set_type;
  
  TaskMergeScale(const hypergraph_set_type& __treebanks,
		 const grammar_type& __grammar,
		 queue_type& __queue)
    : treebanks(__treebanks),
      grammar(__grammar),
      queue(__queue),
      scale()
  {
    scale.set_empty_key(symbol_type());
  }
  
  void operator()()
  {
    // we will statistics that are requried for P(A_1 | A) and P(A_2 | A)
    
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    
    weight_set_type inside;
    weight_set_type outside;
    
    int id = 0;
    for (;;) {
      queue.pop(id);
      if (id < 0) break;
      
      const hypergraph_type& treebank = treebanks[id];
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
      
      cicada::inside(treebank, inside, weight_function(grammar));
      cicada::outside(treebank, inside, outside, weight_function(grammar));
      
      const weight_type weight_total = inside.back();

      hypergraph_type::node_set_type::const_iterator niter_end = treebank.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = treebank.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	const hypergraph_type::edge_type& edge = treebank.edges[node.edges.front()];
	const symbol_type lhs = edge.rule->lhs;
	
	scale[lhs] += inside[node.id] * outside[node.id] / weight_total;
      }
    }
  }
  
  const hypergraph_set_type& treebanks;
  const grammar_type& grammar;
  queue_type& queue;

  scale_set_type scale;
};

template <typename Loss, typename Scale>
struct TaskMergeLoss : public Annotator
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;

  typedef Loss loss_set_type;
  typedef Scale scale_set_type;

  TaskMergeLoss(const hypergraph_set_type& __treebanks,
		const grammar_type& __grammar,
		const scale_set_type& __scale,
		const int& __bits,
		queue_type& __queue)
    : Annotator(__bits),
      treebanks(__treebanks),
      grammar(__grammar),
      scale(__scale),
      queue(__queue),
      loss()
  {
    loss.set_empty_key(symbol_type());
  }
  
  void operator()()
  {
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    typedef std::pair<symbol_type, hypergraph_type::id_type> symbol_id_type;
    typedef std::vector<symbol_id_type, std::allocator<symbol_id_type> > symbol_id_set_type;
    typedef std::vector<symbol_id_set_type, std::allocator<symbol_id_set_type> > symbol_id_map_type;
    
    weight_set_type inside;
    weight_set_type outside;
    symbol_id_map_type symbols;
    
    const attribute_type attr_node("node");
    
    int id = 0;
    for (;;) {
      queue.pop(id);
      if (id < 0) break;
      
      const hypergraph_type& treebank = treebanks[id];
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
    
      cicada::inside(treebank, inside, weight_function(grammar));
      cicada::outside(treebank, inside, outside, weight_function(grammar));
    
      symbols.clear();
      hypergraph_type::node_set_type::const_iterator niter_end = treebank.nodes.end();
      for (hypergraph_type::node_set_type::const_iterator niter = treebank.nodes.begin(); niter != niter_end; ++ niter) {
	const hypergraph_type::node_type& node = *niter;
	const hypergraph_type::edge_type& edge = treebank.edges[node.edges.front()];
      
	const symbol_type lhs = edge.rule->lhs;
      
	attribute_set_type::const_iterator aiter = edge.attributes.find(attr_node);
	if (aiter == edge.attributes.end())
	  throw std::runtime_error("no node attribute?");
	const int node_id_prev = boost::apply_visitor(attribute_integer(), aiter->second);
	if (node_id_prev < 0)
	  throw std::runtime_error("invalid node attribute?");
      
	if (node_id_prev >= static_cast<int>(symbols.size()))
	  symbols.resize(node_id_prev + 1);
      
	symbols[node_id_prev].push_back(symbol_id_type(lhs, node.id));
      }
    
      // now, collect loss...
      const weight_type weight_total = inside.back();
      
      symbol_id_map_type::const_iterator siter_end = symbols.end();
      for (symbol_id_map_type::const_iterator siter = symbols.begin(); siter != siter_end; ++ siter) {
	if (siter->size() == 2) {
	  weight_type prob_split;
	  weight_type inside_merge;
	  weight_type outside_merge;
	  count_type  scale_norm = 0.0;
	  
	  // is it correct?
	  symbol_id_set_type::const_iterator iter_end = siter->end();
	  for (symbol_id_set_type::const_iterator iter = siter->begin(); iter != iter_end; ++ iter) {
	    
	    typename scale_set_type::const_iterator witer = scale.find(iter->first);
	    if (witer == scale.end())
	      throw std::runtime_error("no scale? " + static_cast<const std::string&>(iter->first));

	    prob_split += inside[iter->second] * outside[iter->second];
	    inside_merge += inside[iter->second] * weight_type(witer->second);
	    outside_merge += outside[iter->second];
	    
	    scale_norm += witer->second;
	  }
	  
	  const weight_type loss_node = (inside_merge * outside_merge / weight_type(scale_norm)) / prob_split;
	  
	  std::pair<typename loss_set_type::iterator, bool> result = loss.insert(std::make_pair(annotate(siter->front().first, false), loss_node));
	  if (! result.second)
	    result.first->second *= loss_node;
	} else if (siter->size() > 2)
	  throw std::runtime_error("more than two splitting?");
      }
    }
  }
  
  const hypergraph_set_type& treebanks;
  const grammar_type& grammar;
  const scale_set_type& scale;
  
  queue_type& queue;
  
  loss_set_type loss;
};

template <typename Merged>
struct TaskMergeTreebank
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  TaskMergeTreebank(hypergraph_set_type& __treebanks,
		    const Merged& __merged,
		    queue_type& __queue)
    : treebanks(__treebanks),
      merged(__merged),
      queue(__queue) {}
  
  void operator()()
  {
    hypergraph_type treebank_new;
    filter_pruned::removed_type removed;

    int id = 0;
    for (;;) {
      queue.pop(id);
      if (id < 0) break;
      
      hypergraph_type& treebank = treebanks[id];
      
      removed.clear();
      removed.resize(treebank.edges.size(), false);
      
      hypergraph_type::edge_set_type::iterator eiter_end = treebank.edges.end();
      for (hypergraph_type::edge_set_type::iterator eiter = treebank.edges.begin(); eiter != eiter_end; ++ eiter) {
	hypergraph_type::edge_type& edge = *eiter;
	
	const symbol_type lhs = edge.rule->lhs;
	if (merged.find(lhs) != merged.end())
	  removed[edge.id] = true;
	else {
	  symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	    if (riter->is_non_terminal() && merged.find(*riter) != merged.end())
	      removed[edge.id] = true;
	}
      }
      
      cicada::topologically_sort(treebank, treebank_new, filter_pruned(removed));
      
      treebank.swap(treebank_new);
    }
  }
  
  hypergraph_set_type& treebanks;
  const Merged& merged;
  queue_type& queue;
};

template <typename Merged>
struct TaskMergeGrammar : public Annotator
{
  typedef utils::lockfree_list_queue<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > queue_type;
  
  TaskMergeGrammar(const int __bits,
		   const Merged& __merged,
		   queue_type& __queue)
    : Annotator(__bits),
      merged(__merged),
      queue(__queue) {}
  
  void operator()()
  {
    const grammar_type::value_type* ptr = 0;
    for (;;) {
      queue.pop(ptr);
      if (! ptr) break;
      
      const rule_ptr_type& rule = ptr->first;
      
      symbol_type lhs = rule->lhs;
      if (merged.find(lhs) != merged.end())
	lhs = annotate(lhs, false);
      
      symbol_set_type symbols(rule->rhs);
      
      symbol_set_type::iterator siter_end = symbols.end();
      for (symbol_set_type::iterator siter = symbols.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal() && merged.find(*siter) != merged.end())
	  *siter = annotate(*siter, false);
      
      counts[lhs][rule_type::create(rule_type(lhs, symbols))] += utils::mathop::exp(ptr->second);
    }
  }
  
  const Merged& merged;
  queue_type& queue;
  count_set_type counts;
};


template <typename Generator>
void grammar_merge(hypergraph_set_type& treebanks, grammar_type& grammar, const int bits, Generator& generator)
{
  typedef google::dense_hash_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > merged_set_type;
  typedef google::dense_hash_map<symbol_type, count_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > scale_set_type;
  typedef google::dense_hash_map<symbol_type, weight_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > loss_set_type;

  typedef TaskMergeScale<scale_set_type>               task_scale_type;
  typedef TaskMergeLoss<loss_set_type, scale_set_type> task_loss_type;
  typedef TaskMergeTreebank<merged_set_type>           task_treebank_type;
  typedef TaskMergeGrammar<merged_set_type>            task_grammar_type;
  
  typedef std::vector<task_scale_type, std::allocator<task_scale_type> >       task_scale_set_type;
  typedef std::vector<task_loss_type, std::allocator<task_loss_type> >         task_loss_set_type;
  typedef std::vector<task_treebank_type, std::allocator<task_treebank_type> > task_treebank_set_type;
  typedef std::vector<task_grammar_type, std::allocator<task_grammar_type> >   task_grammar_set_type;
  
  typedef typename task_scale_type::queue_type    queue_scale_type;
  typedef typename task_loss_type::queue_type     queue_loss_type;
  typedef typename task_treebank_type::queue_type queue_treebank_type;
  typedef typename task_grammar_type::queue_type  queue_grammar_type;
  
  typedef std::vector<const loss_set_type::value_type*, std::allocator<const loss_set_type::value_type*> > sorted_type;
  

  // MapReduce to compute scaling
  queue_scale_type queue_scale;
  task_scale_set_type tasks_scale(threads, task_scale_type(treebanks, grammar, queue_scale));
  
  boost::thread_group workers_scale;
  for (int i = 0; i != threads; ++ i)
    workers_scale.add_thread(new boost::thread(boost::ref(tasks_scale[i])));
  
  for (size_t i = 0; i != treebanks.size(); ++ i)
    queue_scale.push(i);
  
  for (int i = 0; i != threads; ++ i)
    queue_scale.push(-1);
  
  workers_scale.join_all();
  
  scale_set_type scale;
  scale.set_empty_key(symbol_type());
  
  for (int i = 0; i != threads; ++ i) {
    if (scale.empty())
      scale.swap(tasks_scale[i].scale);
    else {
      scale_set_type::const_iterator siter_end = tasks_scale[i].scale.end();
      for (scale_set_type::const_iterator siter = tasks_scale[i].scale.begin(); siter != siter_end; ++ siter)
	scale[siter->first] += siter->second;
    }
  }
  
  // MapReduce to compute loss
  queue_loss_type queue_loss;
  task_loss_set_type tasks_loss(threads, task_loss_type(treebanks, grammar, scale, bits, queue_loss));
  
  boost::thread_group workers_loss;
  for (int i = 0; i != threads; ++ i)
    workers_loss.add_thread(new boost::thread(boost::ref(tasks_loss[i])));
  
  for (size_t i = 0; i != treebanks.size(); ++ i)
    queue_loss.push(i);
  
  for (int i = 0; i != threads; ++ i)
    queue_loss.push(-1);
  
  workers_loss.join_all();
  
  loss_set_type      loss;
  loss.set_empty_key(symbol_type());

  for (int i = 0; i != threads; ++ i) {
    if (loss.empty())
      loss.swap(tasks_loss[i].loss);
    else {
      loss_set_type::const_iterator liter_end = tasks_loss[i].loss.end();
      for (loss_set_type::const_iterator liter = tasks_loss[i].loss.begin(); liter != liter_end; ++ liter) {
	std::pair<loss_set_type::iterator, bool> result = loss.insert(*liter);
	if (! result.second)
	  result.first->second *= liter->second;
      }
    }
  }
  
  // sort wrt gain of merging == loss of splitting...
  sorted_type sorted;
  sorted.reserve(loss.size());
  
  loss_set_type::const_iterator liter_end = loss.end();
  for (loss_set_type::const_iterator liter = loss.begin(); liter != liter_end; ++ liter)
    sorted.push_back(&(*liter));
  
  const size_t sorted_size = utils::bithack::min(utils::bithack::max(1, int(merge_ratio * sorted.size())), int(sorted.size() - 1));
  std::nth_element(sorted.begin(), sorted.begin() + sorted_size, sorted.end(), greater_ptr_second<loss_set_type::value_type>());
  
  const weight_type threshold = sorted[sorted_size]->second;
  
  merged_set_type merged;
  merged.set_empty_key(symbol_type());
  
  sorted_type::const_iterator siter_end = sorted.end();
  for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end && (*siter)->second >= threshold; ++ siter) {
    //std::cerr << "merge: " << (*siter)->first << " gain: " << (*siter)->second << std::endl;
    merged.insert(annotate_symbol((*siter)->first, bits, true));
  }
  
  // MapReduce to merge treeebanks
  queue_treebank_type queue_treebank;

  boost::thread_group workers_treebank;
  for (int i = 0; i != threads; ++ i)
    workers_treebank.add_thread(new boost::thread(task_treebank_type(treebanks, merged, queue_treebank)));
  
  for (size_t i = 0; i != treebanks.size(); ++ i)
    queue_treebank.push(i);
  
  for (int i = 0; i != threads; ++ i)
    queue_treebank.push(-1);

  workers_treebank.join_all();
  
  // MapReduce to merge grammar
  queue_grammar_type queue_grammar;
  task_grammar_set_type tasks_grammar(threads, task_grammar_type(bits, merged, queue_grammar));
  
  boost::thread_group workers_grammar;
  for (int i = 0; i != threads; ++ i)
    workers_grammar.add_thread(new boost::thread(boost::ref(tasks_grammar[i])));
  
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    queue_grammar.push(&(*giter));
  
  
  for (int i = 0; i != threads; ++ i)
    queue_grammar.push(0);
  
  workers_grammar.join_all();
  
  count_set_type counts;
  
  for (int i = 0; i != threads; ++ i) {
    if (counts.empty())
      counts.swap(tasks_grammar[i].counts);
    else {
      count_set_type::const_iterator citer_end = tasks_grammar[i].counts.end();
      for (count_set_type::const_iterator citer = tasks_grammar[i].counts.begin(); citer != citer_end; ++ citer) {
	grammar_type& grammar = counts[citer->first];
	
	grammar_type::const_iterator giter_end = citer->second.end();
	for (grammar_type::const_iterator giter = citer->second.begin(); giter != giter_end; ++ giter)
	  grammar[giter->first] += giter->second;
      }
    }
  }
  
  if (debug)
    std::cerr << "# of symbols: " << counts.size() << std::endl;
  
  // maximization
  if (variational_bayes_mode)
    grammar_maximize(counts, grammar, MaximizeBayes());
  else
    grammar_maximize(counts, grammar, Maximize());
}

struct TaskSplitTreebank : public Annotator
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  typedef utils::hashmurmur<size_t> hasher_type;
  
  TaskSplitTreebank(hypergraph_set_type& __treebanks,
		    const int __bits,
		    queue_type& __queue)
    : Annotator(__bits),
      treebanks(__treebanks),
      queue(__queue)
  {}
  
  void operator()()
  {
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
    
    node_map_type   node_map;
    hypergraph_type treebank_new;
    
    index_set_type  j;
    index_set_type  j_end;
    symbol_set_type symbols;
    symbol_set_type symbols_new;
    
    const symbol_type root("[ROOT]");
    const attribute_type attr_node("node");
    
    int id = 0;
    for (;;) {
      queue.pop(id);
      if (id < 0) break;
      
      hypergraph_type& treebank = treebanks[id];
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
		symbols_new[i] = annotate(symbols[i], j[i]);
	    
	    const rule_ptr_type rule = rule_type::create(rule_type(symbols_new.front(), symbols_new.begin() + 1, symbols_new.end()));
	  
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
    
  }
  
  hypergraph_set_type& treebanks;
  queue_type& queue;
};

template <typename Generator>
struct TaskSplitGrammar : public Annotator
{
  typedef utils::lockfree_list_queue<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > queue_type;

  TaskSplitGrammar(Generator __generator,
		   const int& __bits,
		   queue_type& __queue)
    : Annotator(__bits),
      generator(__generator),
      queue(__queue) {}

  void operator()()
  {
    typedef std::vector<int, std::allocator<int> > index_set_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    
    index_set_type  j;
    index_set_type  j_end;
    symbol_set_type symbols;
    symbol_set_type symbols_new;
    
    const symbol_type root("[ROOT]");
    
    const grammar_type::value_type* ptr = 0;
    
    for (;;) {
      queue.pop(ptr);
      if (! ptr) break;
      
      const rule_ptr_type& rule = ptr->first;
      
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
	    symbols_new[i] = annotate(symbols[i], j[i]);
      
	const rule_ptr_type rule = rule_type::create(rule_type(symbols_new.front(), symbols_new.begin() + 1, symbols_new.end()));
      
	// we will add 1% of randomness...
	counts[rule->lhs][rule] += utils::mathop::exp(ptr->second) * boost::uniform_real<double>(0.99, 1.01)(generator);
	
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
  
  Generator generator;
  queue_type& queue;
  count_set_type counts;
};

template <typename Generator>
void grammar_split(hypergraph_set_type& treebanks, grammar_type& grammar, const int bits, Generator& generator)
{
  typedef TaskSplitTreebank           task_treebank_type;
  typedef TaskSplitGrammar<Generator> task_grammar_type;

  typedef std::vector<task_grammar_type, std::allocator<task_grammar_type> > task_grammar_set_type;
  
  typedef typename task_treebank_type::queue_type queue_treebank_type;
  typedef typename task_grammar_type::queue_type  queue_grammar_type;
  
  queue_treebank_type queue_treebank;
  
  boost::thread_group workers_treebank;
  for (int i = 0; i != threads; ++ i)
    workers_treebank.add_thread(new boost::thread(task_treebank_type(treebanks, bits, queue_treebank)));

  for (size_t i = 0; i != treebanks.size(); ++ i)
    queue_treebank.push(i);
  
  for (int i = 0; i != threads; ++ i)
    queue_treebank.push(-1);

  workers_treebank.join_all();

  queue_grammar_type  queue_grammar;
  task_grammar_set_type tasks_grammar(threads, task_grammar_type(generator, bits, queue_grammar));
  
  boost::thread_group workers_grammar;
  for (int i = 0; i != threads; ++ i)
    workers_grammar.add_thread(new boost::thread(boost::ref(tasks_grammar[i])));
  
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    queue_grammar.push(&(*giter));
  
  for (int i = 0; i != threads; ++ i)
    queue_grammar.push(0);
  
  workers_grammar.join_all();

  // split grammar...
  count_set_type counts;
  
  for (int i = 0; i != threads; ++ i) {
    if (counts.empty())
      counts.swap(tasks_grammar[i].counts);
    else {
      count_set_type::const_iterator citer_end = tasks_grammar[i].counts.end();
      for (count_set_type::const_iterator citer = tasks_grammar[i].counts.begin(); citer != citer_end; ++ citer) {
	grammar_type& grammar = counts[citer->first];
	
	grammar_type::const_iterator giter_end = citer->second.end();
	for (grammar_type::const_iterator giter = citer->second.begin(); giter != giter_end; ++ giter)
	  grammar[giter->first] += giter->second;
      }
    }
  }

  if (debug)
    std::cerr << "# of symbols: " << counts.size() << std::endl;
  
  // maximization
  if (variational_bayes_mode)
    grammar_maximize(counts, grammar, MaximizeBayes());
  else
    grammar_maximize(counts, grammar, Maximize());
}


struct accumulator_type
{
  typedef weight_type value_type;

  accumulator_type(const weight_type& __weight_total,
		   const hypergraph_type& __treebank,
		   count_set_type& __counts)
    : weight_total(__weight_total), treebank(__treebank), counts(__counts) {}

  struct Count
  {
    Count(count_type& __count, const weight_type& __weight) : count(__count), weight(__weight) {}
    
    Count& operator+=(const weight_type& value)
    {
      count += value / weight;
      return *this;
    }
    
    count_type& count;
    const weight_type& weight;
  };
  
  Count operator[](const hypergraph_type::id_type& x)
  {
    const rule_ptr_type& rule = treebank.edges[x].rule;
    
    return Count(counts[rule->lhs][rule], weight_total);
  }
  
  const weight_type& weight_total;
  const hypergraph_type& treebank;
  count_set_type& counts;
};

template <typename Function>
struct TaskLearn
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  TaskLearn(const hypergraph_set_type& __treebanks, queue_type& __queue, Function __function)
    : treebanks(__treebanks),
      queue(__queue),
      function(__function),
      logprob(cicada::semiring::traits<weight_type>::one()),
      counts() {}
  
  void operator()()
  {
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    
    weight_set_type inside;
    weight_set_type outside;
    
    int id = 0;
    
    for (;;) {
      queue.pop(id);
      if (id < 0) break;
      
      const hypergraph_type& treebank = treebanks[id];
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
      
      accumulator_type accumulator(inside.back(), treebank, counts);
      
      cicada::inside_outside(treebank, inside, outside, accumulator, function, function);
      
      if (debug >= 2)
	std::cerr << "inside: " << cicada::semiring::log(inside.back()) << std::endl;
      
      logprob *= inside.back();
    }
  }
  
  const hypergraph_set_type& treebanks;
  queue_type&    queue;
  Function       function;
  
  weight_type    logprob;
  count_set_type counts;
};

template <typename Function>
double grammar_learn(const hypergraph_set_type& treebanks, grammar_type& grammar, Function function)
{
  typedef TaskLearn<Function> task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  typedef typename task_type::queue_type queue_type;

  queue_type queue;
  task_set_type tasks(threads, task_type(treebanks, queue, function));

  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));

  for (size_t i = 0; i != treebanks.size(); ++ i)
    queue.push(i);
  
  for (int i = 0; i != threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
  
  weight_type logprob(cicada::semiring::traits<weight_type>::one());
  count_set_type counts;
  for (int i = 0; i != threads; ++ i) {
    logprob *= tasks[i].logprob;
    
    if (counts.empty())
      counts.swap(tasks[i].counts);
    else {
      count_set_type::const_iterator liter_end = tasks[i].counts.end();
      for (count_set_type::const_iterator liter = tasks[i].counts.begin(); liter != liter_end; ++ liter) {
	grammar_type& grammar = counts[liter->first];

	grammar_type::const_iterator giter_end = liter->second.end();
	for (grammar_type::const_iterator giter = liter->second.begin(); giter != giter_end; ++ giter)
	  grammar[giter->first] += giter->second;
      }
    }
  }
  
  if (debug)
    std::cerr << "log-likelihood: " << cicada::semiring::log(logprob)
	      << " # of symbols: " << counts.size()
	      << std::endl;
  
  // maximization
  if (variational_bayes_mode)
    grammar_maximize(counts, grammar, MaximizeBayes());
  else
    grammar_maximize(counts, grammar, Maximize());
  
  return cicada::semiring::log(logprob);
}

template <typename Maximizer>
struct TaskMaximize
{
  typedef utils::lockfree_list_queue<const count_set_type::value_type*, std::allocator<const count_set_type::value_type*> > queue_type;
  
  TaskMaximize(queue_type& __queue,
	       Maximizer __maximizer)
    : queue(__queue),
      maximizer(__maximizer) {}
  
  void operator()()
  {
    const count_set_type::value_type* ptr = 0;
    
    for (;;) {
      queue.pop(ptr);
      if (! ptr) break;
      
      maximizer(ptr->second, grammar);
    }
  }
  
  queue_type& queue;
  grammar_type grammar;
  Maximizer    maximizer;
};


template <typename Maximizer>
void grammar_maximize(const count_set_type& counts, grammar_type& grammar, Maximizer maximizer)
{
  typedef TaskMaximize<Maximizer> task_type;
  typedef typename task_type::queue_type queue_type;
  
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  
  queue_type queue;
  task_set_type tasks(threads, task_type(queue, maximizer));
  
  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
    queue.push(&(*citer));
  
  for (int i = 0; i != threads; ++ i)
    queue.push(0);
  
  workers.join_all();
  
  grammar.clear();
  for (int i = 0; i != threads; ++ i)
    grammar.insert(tasks[i].grammar.begin(), tasks[i].grammar.end());
}

void write_grammar(const path_type& file, const grammar_type& grammar)
{
  typedef std::vector<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > sorted_type;

  if (grammar.empty()) return;
  
  if (0.0 < cutoff_threshold && cutoff_threshold < 1.0) {
    count_set_type counts;
    sorted_type sorted;
    
    grammar_type::const_iterator giter_end = grammar.end();
    for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
      counts[giter->first->lhs][giter->first] = giter->second;

    utils::compress_ostream os(file, 1024 * 1024);
    
    count_set_type::const_iterator citer_end = counts.end();
    for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      const grammar_type& grammar = citer->second;
      
      sorted.clear();
      
      grammar_type::const_iterator giter_end = grammar.end();
      for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
	sorted.push_back(&(*giter));
      
      std::sort(sorted.begin(), sorted.end(), greater_ptr_second<grammar_type::value_type>());
      
      const double logprob_max = sorted.front()->second;
      const double logprob_threshold = logprob_max + utils::mathop::log(cutoff_threshold);
      
      sorted_type::const_iterator siter_end = sorted.end();
      for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end && (*siter)->second >= logprob_threshold; ++ siter)
	os << *((*siter)->first) << " ||| ||| " << (*siter)->second << '\n';
    }
  } else {
    // we will dump in sorted order...
    sorted_type sorted;
    sorted.reserve(grammar.size());
    
    grammar_type::const_iterator giter_end = grammar.end();
    for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
      sorted.push_back(&(*giter));
    
    std::sort(sorted.begin(), sorted.end(), greater_ptr_second<grammar_type::value_type>());
    
    utils::compress_ostream os(file, 1024 * 1024);
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << *((*siter)->first) << " ||| ||| " << (*siter)->second << '\n';
  }
}

void read_treebank(const path_set_type& files, hypergraph_set_type& treebanks)
{
  hypergraph_type treebank;
  
  path_set_type::const_iterator fiter_end = files.end();
  for (path_set_type::const_iterator fiter = files.begin(); fiter != fiter_end; ++ fiter) {
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (is >> treebank) {
      if (! treebank.is_valid()) continue;

      if (binarize_left)
	cicada::binarize_left(treebank, 0);
      else if (binarize_right)
	cicada::binarize_right(treebank, 0);
      else if (binarize_all)
	cicada::binarize_all(treebank);
      
      treebanks.push_back(hypergraph_type());
      treebanks.back().swap(treebank);
    }
  }
  
  if (debug)
    std::cerr << "# of treebank: " << treebanks.size() << std::endl;

}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  
  desc.add_options()
    ("input",  po::value<path_set_type>(&input_files), "input treebank")
    ("output-grammar", po::value<path_type>(&output_grammar_file), "output grammar")
    ("output-unknown", po::value<path_type>(&output_unknown_file), "output unknown rules")
    
    ("max-iteration",       po::value<int>(&max_iteration)->default_value(max_iteration),             "maximum split/merge iterations")
    ("max-iteration-split", po::value<int>(&max_iteration_split)->default_value(max_iteration_split), "maximum EM iterations after split")
    ("max-iteration-merge", po::value<int>(&max_iteration_merge)->default_value(max_iteration_merge), "maximum EM iterations after merge")
    
    ("binarize-left",  po::bool_switch(&binarize_left),  "left binarization")
    ("binarize-right", po::bool_switch(&binarize_right), "right binarization")
    ("binarize-all",   po::bool_switch(&binarize_all),   "all binarization")
    
    ("variational-bayes", po::bool_switch(&variational_bayes_mode), "variational Bayes estimates")
    
    ("prior",           po::value<double>(&prior)->default_value(prior),                   "Dirichlet prior")
    ("prior-terminal",  po::value<double>(&prior_terminal)->default_value(prior_terminal), "Dirichlet prior for terminal rule")

    
    ("keep-terminal", po::bool_switch(&keep_terminal),                             "keep terminal rule (and do not split/merge)")
    ("merge-ratio",   po::value<double>(&merge_ratio)->default_value(merge_ratio), "merging ratio")
    
    ("cutoff-threshold", po::value<double>(&cutoff_threshold)->default_value(cutoff_threshold), "dump with beam-threshold (<= 0.0 implies no beam)")
    ("cutoff-unknown",   po::value<int>(&cutoff_unknown)->default_value(cutoff_unknown),        "cut-off threshold for unknown word (<=1 implies no cutoff)")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::variables_map variables;
  
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
