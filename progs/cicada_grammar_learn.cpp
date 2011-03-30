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
//
// TODO: sample unknown word...
// 

#include <stdexcept>
#include <vector>
#include <deque>

#include <cicada/vocab.hpp>
#include <cicada/hypergraph.hpp>
#include <cicada/inside_outside.hpp>
#include <cicada/semiring.hpp>
#include <cicada/sort.hpp>
#include <cicada/binarize.hpp>
#include <cicada/signature.hpp>
#include <cicada/tokenizer.hpp>
#include <cicada/sentence.hpp>

#include <boost/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/bind.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <boost/math/special_functions/expm1.hpp>

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

typedef cicada::Signature signature_type;
typedef cicada::Tokenizer tokenizer_type;

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

typedef symbol_set_type ngram_type;
class NGramCounts : public google::dense_hash_map<ngram_type, count_type, boost::hash<ngram_type>, std::equal_to<ngram_type> >
{
public:
  typedef google::dense_hash_map<ngram_type, count_type, boost::hash<ngram_type>, std::equal_to<ngram_type> > count_set_type;

  NGramCounts() : count_set_type() { count_set_type::set_empty_key(ngram_type()); }
};
typedef NGramCounts ngram_count_set_type;

class WordCounts : public google::dense_hash_map<symbol_type, count_type, boost::hash<symbol_type>, std::equal_to<symbol_type> >
{
public:
  typedef google::dense_hash_map<symbol_type, count_type, boost::hash<symbol_type>, std::equal_to<symbol_type> > count_set_type;

  WordCounts() : count_set_type() { count_set_type::set_empty_key(symbol_type()); }
};
typedef WordCounts word_count_set_type;

path_set_type input_files;
path_type     output_grammar_file = "-";
path_type     output_lexicon_file;
path_type     output_character_file;

int max_iteration = 6;         // max split-merge iterations
int max_iteration_split = 20;  // max EM-iterations for split
int max_iteration_merge = 20;  // max EM-iterations for merge

bool binarize_left = false;
bool binarize_right = false;
bool binarize_all = false;


double prior_rule      = 0.1;
double prior_lexicon   = 0.01;
double prior_signature = 0.01;
double prior_character = 0.01;

double merge_ratio = 0.5;
double unknown_ratio = 0.5;
double unknown_threshold = 20;

std::string signature = "";
bool signature_list = false;

double cutoff_rule = 1e-30;
double cutoff_lexicon = 1e-40;
double cutoff_character = 0;

int threads = 1;

int debug = 0;

template <typename Generator, typename Maximizer>
void grammar_merge(hypergraph_set_type& treebanks,
		   grammar_type& grammar,
		   const int bits,
		   Generator& generator,
		   Maximizer maximizer);

template <typename Generator, typename Maximizer>
void grammar_split(hypergraph_set_type& treebanks,
		   grammar_type& grammar,
		   const int bits,
		   Generator& generator,
		   Maximizer maximizer);

template <typename Function, typename Maximizer>
double grammar_learn(const hypergraph_set_type& treebanks,
		     grammar_type& grammar,
		     Function function,
		     Maximizer maximier);

template <typename Function>
void lexicon_learn(const hypergraph_set_type& treebanks,
		   grammar_type& lexicon,
		   Function function);

template <typename Function>
void characters_learn(const hypergraph_set_type& treebanks,
		      ngram_count_set_type& model,
		      ngram_count_set_type& backoff,
		      Function function);

template <typename Maximizer>
void grammar_maximize(const count_set_type& counts,
		      grammar_type& grammar,
		      Maximizer maximizer);

void write_characters(const path_type& file,
		      const ngram_count_set_type& model,
		      const ngram_count_set_type& backoff,
		      const double cutoff);
void write_grammar(const path_type& file,
		   const grammar_type& grammar);

void read_treebank(const path_set_type& files,
		   hypergraph_set_type& treebanks);

void grammar_prune(grammar_type& grammar, const double cutoff);
void lexicon_prune(grammar_type& grammar, const double cutoff);

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

struct Maximize
{
  void operator()(const grammar_type& counts, grammar_type& grammar) const
  {
    // simle maximizer...
    double sum = 0.0;
    grammar_type::const_iterator citer_end = counts.end();
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      sum += citer->second;
    
    const double logsum = cicada::semiring::log(weight_type(sum));
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
      grammar[citer->first] = cicada::semiring::log(weight_type(citer->second)) - logsum;
  }
};

struct MaximizeBayes : public utils::hashmurmur<size_t>
{
  typedef utils::hashmurmur<size_t> hasher_type;
  
  MaximizeBayes(const grammar_type& __base) : base(__base) {}
  
  typedef std::vector<prob_type, std::allocator<prob_type> > prob_set_type;
  
  struct Cache
  {
    symbol_type symbol;
    symbol_type coarse;
    
    Cache() : symbol(), coarse() {}
  };
  typedef Cache cache_type;
  typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;

  const symbol_type& coarse(const symbol_type& symbol) const
  {
    if (symbol.is_terminal()) return symbol;

    const size_t cache_pos = hasher_type::operator()(symbol.id()) & (caches.size() - 1);
    cache_type& cache = const_cast<cache_type&>(caches[cache_pos]);
    if (cache.symbol != symbol) {
      cache.symbol = symbol;
      
      const std::string str = symbol.non_terminal_strip();
      std::string::size_type pos = str.rfind('@');
      if (pos != std::string::npos)
	cache.coarse = '[' + str.substr(0, pos) + ']';
      else
	cache.coarse = symbol;
    }
    return cache.coarse;
  }
  
  
  prob_set_type  __probs;
  cache_set_type caches;
  const grammar_type& base;
  
  void operator()(const grammar_type& counts, grammar_type& grammar) const
  {
    using namespace boost::math::policies;
    typedef policy<domain_error<errno_on_error>,
		   pole_error<errno_on_error>,
		   overflow_error<errno_on_error>,
		   rounding_error<errno_on_error>,
		   evaluation_error<errno_on_error> > policy_type;
    
    if (counts.empty()) return;
    
    const bool is_terminal = counts.begin()->first->rhs.front().is_terminal();
    const double prior = (is_terminal ? prior_lexicon : prior_rule);

    prob_set_type& probs = const_cast<prob_set_type&>(__probs);
    
    probs.resize(counts.size());
    
    double sum = 0.0;
    prob_set_type::iterator piter = probs.begin();
    grammar_type::const_iterator citer_end = counts.end();
    for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ piter) {
      const symbol_type lhs = coarse(citer->first->lhs);
      symbol_set_type rhs(citer->first->rhs);
      symbol_set_type::iterator riter_end = rhs.end();
      for (symbol_set_type::iterator riter = rhs.begin(); riter != riter_end; ++ riter)
	*riter = coarse(*riter);
      
      const rule_ptr_type rule_coarse(rule_type::create(rule_type(lhs, rhs)));
      
      grammar_type::const_iterator biter = base.find(rule_coarse);
      if (biter == base.end())
	throw std::runtime_error("no base?");
      
      *piter = utils::mathop::exp(biter->second);
      
      sum += citer->second + prior * (*piter);
    }
    
    for (;;) {
      weight_type logprob_sum;
      const double logsum = utils::mathop::digamma(sum);
      
      prob_set_type::iterator piter = probs.begin();
      for (grammar_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ piter) {
	const double logprob = utils::mathop::digamma(citer->second + prior * *piter) - logsum;
	
	grammar[citer->first] = logprob;
	logprob_sum += cicada::semiring::traits<weight_type>::exp(logprob);
      }
      
      const double discount = - boost::math::expm1(cicada::semiring::log(logprob_sum), policy_type());
      if (discount > 0.0) break;
      
      ++ sum;
    }
  }
};

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (signature_list) {
      std::cout << signature_type::lists();
      return 0;
    }
    
    if (merge_ratio <= 0.0 || 1.0 <= merge_ratio)
      throw std::runtime_error("invalid merge ratio");
      
    if (int(binarize_left) + binarize_right + binarize_all > 1)
      throw std::runtime_error("specify either binarize-{left,right,all}");

    const signature_type* sig = (! signature.empty() ? &signature_type::create(signature) : 0);
    
    if (! output_character_file.empty()) {
      if (! sig)
	throw std::runtime_error("character estimation requires signature");
      
      if (output_lexicon_file.empty())
	throw std::runtime_error("we will dump character file, but no lexicon file");
    }
    
    if (int(binarize_left) + binarize_right + binarize_all == 0)
      binarize_left = true;

    if (input_files.empty())
      input_files.push_back("-");
    
    threads = utils::bithack::max(threads, 1);
    
    hypergraph_set_type treebanks;
    read_treebank(input_files, treebanks);

    grammar_type grammar;
    grammar_learn(treebanks, grammar, zero_function(), Maximize());
    
    boost::mt19937 generator;
    generator.seed(time(0) * getpid());
    
    if (debug)
      std::cerr << "grammar size: " << grammar.size() << std::endl;
    
    grammar_type base(grammar);
    
    for (int iter = 0; iter < max_iteration; ++ iter) {
      
      if (debug)
	std::cerr << "iteration: " << (iter + 1) << std::endl;
      
      // split...
      {
	const utils::resource split_start;
	grammar_split(treebanks, grammar, iter, generator, MaximizeBayes(base));
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
	  const double logprob_curr = grammar_learn(treebanks, grammar, weight_function(grammar), MaximizeBayes(base));
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
	grammar_merge(treebanks, grammar, iter, generator, MaximizeBayes(base));
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
	  const double logprob_curr = grammar_learn(treebanks, grammar, weight_function(grammar), MaximizeBayes(base));
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

    if (! output_lexicon_file.empty()) {
      // first, split into two
      grammar_type rules;
      grammar_type lexicon;
      
      grammar_type::const_iterator giter_end = grammar.end();
      for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter) {
	if (giter->first->rhs.size() == 1 && giter->first->rhs.front().is_terminal())
	  lexicon.insert(*giter);
	else
	  rules.insert(*giter);
      }
      
      if (! output_character_file.empty()) {
	ngram_count_set_type model;
	ngram_count_set_type backoff;
	
	characters_learn(treebanks, model, backoff, weight_function(grammar));
	
	write_characters(output_character_file, model, backoff, cutoff_character);
      } 
      
      if (sig)
	lexicon_learn(treebanks, lexicon, weight_function(grammar));

      if (0.0 < cutoff_rule && cutoff_rule < 1.0)
	grammar_prune(rules, cutoff_rule);
      if (0.0 < cutoff_lexicon && cutoff_lexicon < 1.0)
	lexicon_prune(lexicon, cutoff_lexicon);
      
      write_grammar(output_grammar_file, rules);
      write_grammar(output_lexicon_file, lexicon);
    } else {
      if (0.0 < cutoff_rule && cutoff_rule < 1.0)
	grammar_prune(grammar, cutoff_rule);
      
      write_grammar(output_grammar_file, grammar);
    }
  }
  catch (std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

bool is_fixed_non_terminal(const symbol_type& symbol)
{ 
  static const symbol_type root("[ROOT]");
  
  return symbol.is_non_terminal() && symbol == root;
};

symbol_type annotate_symbol(const symbol_type& symbol, const int bitpos, const bool bit)
{
  if (symbol.is_non_terminal()) {
    if (is_fixed_non_terminal(symbol)) return symbol;

    namespace xpressive = boost::xpressive;
    
    typedef xpressive::basic_regex<utils::piece::const_iterator> pregex;
    typedef xpressive::match_results<utils::piece::const_iterator> pmatch;
    
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



template <typename Scale>
struct TaskMergeScale
{
  typedef utils::lockfree_list_queue<const hypergraph_type*, std::allocator<const hypergraph_type*> > queue_type;
  typedef Scale scale_set_type;
  
  TaskMergeScale(const grammar_type& __grammar,
		 queue_type& __queue)
    : grammar(__grammar),
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
    
    const hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      const hypergraph_type& treebank = *__treebank;
      
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
  
  const grammar_type& grammar;
  queue_type& queue;

  scale_set_type scale;
};

template <typename Loss, typename Scale>
struct TaskMergeLoss : public Annotator
{
  typedef utils::lockfree_list_queue<const hypergraph_type*, std::allocator<const hypergraph_type*> > queue_type;

  typedef Loss loss_set_type;
  typedef Scale scale_set_type;

  TaskMergeLoss(const grammar_type& __grammar,
		const scale_set_type& __scale,
		const int& __bits,
		queue_type& __queue)
    : Annotator(__bits),
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
    
    const hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      const hypergraph_type& treebank = *__treebank;
      
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
	  
	  std::pair<typename loss_set_type::iterator, bool> result = loss.insert(std::make_pair(annotate(siter->front().first, true), loss_node));
	  if (! result.second)
	    result.first->second *= loss_node;
	} else if (siter->size() > 2)
	  throw std::runtime_error("more than two splitting?");
      }
    }
  }
  
  const grammar_type& grammar;
  const scale_set_type& scale;
  
  queue_type& queue;
  
  loss_set_type loss;
};

template <typename Merged>
struct TaskMergeTreebank
{
  typedef utils::lockfree_list_queue<hypergraph_type*, std::allocator<hypergraph_type*> > queue_type;
  
  TaskMergeTreebank(const Merged& __merged,
		    queue_type& __queue)
    : merged(__merged),
      queue(__queue) {}
  
  void operator()()
  {
    hypergraph_type treebank_new;
    filter_pruned::removed_type removed;
    
    hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      hypergraph_type& treebank = *__treebank;
      
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
      treebank_new.clear();
      hypergraph_type(treebank).swap(treebank);
    }
  }
  
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
      
      bool annotated = false;
      symbol_type lhs = rule->lhs;
      if (merged.find(lhs) != merged.end()) {
	lhs = annotate(lhs, false);
	annotated = true;
      }
      
      symbol_set_type symbols(rule->rhs);
      
      symbol_set_type::iterator siter_end = symbols.end();
      for (symbol_set_type::iterator siter = symbols.begin(); siter != siter_end; ++ siter)
	if (siter->is_non_terminal() && merged.find(*siter) != merged.end()) {
	  *siter = annotate(*siter, false);
	  annotated = true;
	}
      
      if (annotated)
	counts[lhs][rule_type::create(rule_type(lhs, symbols))] += utils::mathop::exp(ptr->second);
      else
	counts[rule->lhs][rule] += utils::mathop::exp(ptr->second);
    }
  }
  
  const Merged& merged;
  queue_type& queue;
  count_set_type counts;
};


template <typename Generator, typename Maximizer>
void grammar_merge(hypergraph_set_type& treebanks,
		   grammar_type& grammar,
		   const int bits,
		   Generator& generator,
		   Maximizer maximizer)
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
  task_scale_set_type tasks_scale(threads, task_scale_type(grammar, queue_scale));
  
  boost::thread_group workers_scale;
  for (int i = 0; i != threads; ++ i)
    workers_scale.add_thread(new boost::thread(boost::ref(tasks_scale[i])));
  
  hypergraph_set_type::iterator titer_end = treebanks.end();
  for (hypergraph_set_type::iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue_scale.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue_scale.push(0);
  
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
  task_loss_set_type tasks_loss(threads, task_loss_type(grammar, scale, bits, queue_loss));
  
  boost::thread_group workers_loss;
  for (int i = 0; i != threads; ++ i)
    workers_loss.add_thread(new boost::thread(boost::ref(tasks_loss[i])));
  
  for (hypergraph_set_type::iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue_loss.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue_loss.push(0);
  
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
  
  const size_t sorted_size = utils::bithack::min(utils::bithack::max(size_t(1), size_t(merge_ratio * sorted.size())), size_t(sorted.size() - 1));
  std::nth_element(sorted.begin(), sorted.begin() + sorted_size, sorted.end(), greater_ptr_second<loss_set_type::value_type>());
  
  const weight_type threshold = sorted[sorted_size]->second;
  
  merged_set_type merged;
  merged.set_empty_key(symbol_type());
  
  if (debug >= 2)
    std::cerr << "threshold: " << threshold << std::endl;
  
  // do we stricktly erase...?
  sorted_type::const_iterator siter_end = sorted.end();
  for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end && (*siter)->second >= threshold; ++ siter) {
    if (debug >= 2)
      std::cerr << "merge: " << (*siter)->first << " gain: " << (*siter)->second << std::endl;
    merged.insert((*siter)->first);
  }
  
  if (debug)
    std::cerr << "merged: " << merged.size() << " split: " << (sorted.size() - merged.size()) << std::endl;
  
  // MapReduce to merge treeebanks
  queue_treebank_type queue_treebank;

  boost::thread_group workers_treebank;
  for (int i = 0; i != threads; ++ i)
    workers_treebank.add_thread(new boost::thread(task_treebank_type(merged, queue_treebank)));
  
  for (hypergraph_set_type::iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue_treebank.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue_treebank.push(0);

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
  
  grammar_maximize(counts, grammar, maximizer);
}

struct TaskSplitTreebank : public Annotator
{
  typedef utils::lockfree_list_queue<hypergraph_type*, std::allocator<hypergraph_type*> > queue_type;
  
  typedef utils::hashmurmur<size_t> hasher_type;
  
  TaskSplitTreebank(const int __bits,
		    queue_type& __queue,
		    const grammar_type& __grammar)
    : Annotator(__bits),
      queue(__queue),
      grammar(__grammar)
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
    
    const attribute_type attr_node("node");
    
    hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      hypergraph_type& treebank = *__treebank;
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
	    j_end[i] = utils::bithack::branch(symbols[i].is_non_terminal(), utils::bithack::branch(is_fixed_non_terminal(symbols[i]), 1, 2), 0);
	
	  for (;;) {
	    // construct rule
	    for (size_t i = 0; i != symbols.size(); ++ i)
	      if (j_end[i])
		symbols_new[i] = annotate(symbols[i], j[i]);

	    rule_ptr_type rule = rule_type::create(rule_type(symbols_new.front(), symbols_new.begin() + 1, symbols_new.end()));
	    grammar_type::const_iterator giter = grammar.find(rule);
	    if (giter != grammar.end())
	      rule = giter->first;
	  
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
      hypergraph_type(treebank).swap(treebank);
    }
  }
  
  queue_type& queue;
  const grammar_type& grammar;
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
	j_end[i] = utils::bithack::branch(symbols[i].is_non_terminal(), utils::bithack::branch(is_fixed_non_terminal(symbols[i]), 1, 2), 0);
    
      for (;;) {
	for (size_t i = 0; i != symbols.size(); ++ i)
	  if (j_end[i])
	    symbols_new[i] = annotate(symbols[i], j[i]);
      
	const rule_ptr_type rule = rule_type::create(rule_type(symbols_new.front(), symbols_new.begin() + 1, symbols_new.end()));
      
	// we will add 1% of randomness...
	counts[rule->lhs][rule] = utils::mathop::exp(ptr->second) * boost::uniform_real<double>(0.99, 1.01)(generator);
	
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

template <typename Generator, typename Maximizer>
void grammar_split(hypergraph_set_type& treebanks,
		   grammar_type& grammar,
		   const int bits,
		   Generator& generator,
		   Maximizer maximizer)
{
  typedef TaskSplitTreebank           task_treebank_type;
  typedef TaskSplitGrammar<Generator> task_grammar_type;

  typedef std::vector<task_grammar_type, std::allocator<task_grammar_type> > task_grammar_set_type;
  
  typedef typename task_treebank_type::queue_type queue_treebank_type;
  typedef typename task_grammar_type::queue_type  queue_grammar_type;

  
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

  tasks_grammar.clear();
  
  if (debug)
    std::cerr << "# of symbols: " << counts.size() << std::endl;
  
  // maximization
  grammar_maximize(counts, grammar, maximizer);

  counts.clear();
  
  queue_treebank_type queue_treebank;
  
  boost::thread_group workers_treebank;
  for (int i = 0; i != threads; ++ i)
    workers_treebank.add_thread(new boost::thread(task_treebank_type(bits, queue_treebank, grammar)));
  
  hypergraph_set_type::iterator titer_end = treebanks.end();
  for (hypergraph_set_type::iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue_treebank.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue_treebank.push(0);

  workers_treebank.join_all();
}


template <typename Function>
struct TaskLearn
{
  typedef utils::lockfree_list_queue<const hypergraph_type*, std::allocator<const hypergraph_type*> > queue_type;
  
  TaskLearn(queue_type& __queue, Function __function)
    : queue(__queue),
      function(__function),
      logprob(cicada::semiring::traits<weight_type>::one()),
      counts() {}
  
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


  void operator()()
  {
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    
    weight_set_type inside;
    weight_set_type outside;
    
    const hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      const hypergraph_type& treebank = *__treebank;
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
      
      accumulator_type accumulator(inside.back(), treebank, counts);
      
      cicada::inside_outside(treebank, inside, outside, accumulator, function, function);
      
      if (debug >= 3)
	std::cerr << "inside: " << cicada::semiring::log(inside.back()) << std::endl;
      
      logprob *= inside.back();
    }
  }
  
  queue_type&    queue;
  Function       function;
  
  weight_type    logprob;
  count_set_type counts;
};

template <typename Function, typename Maximizer>
double grammar_learn(const hypergraph_set_type& treebanks,
		     grammar_type& grammar,
		     Function function,
		     Maximizer maximizer)
{
  typedef TaskLearn<Function> task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  typedef typename task_type::queue_type queue_type;

  queue_type queue;
  task_set_type tasks(threads, task_type(queue, function));

  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));

  hypergraph_set_type::const_iterator titer_end = treebanks.end();
  for (hypergraph_set_type::const_iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue.push(0);
  
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
  grammar_maximize(counts, grammar, maximizer);
  
  return cicada::semiring::log(logprob);
}

template <typename Function>
struct TaskLexiconFrequency
{
  typedef utils::lockfree_list_queue<const hypergraph_type*, std::allocator<const hypergraph_type*> > queue_type;

  TaskLexiconFrequency(queue_type& __queue, Function __function)
    : queue(__queue),
      function(__function),
      counts() { }

  void operator()()
  {
        typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    
    weight_set_type inside;
    weight_set_type outside;
    weight_set_type scores;
    
    ngram_type trigram(3);
    ngram_type bigram(2);
    
    const signature_type& __signature = signature_type::create(signature);
    
    const hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      const hypergraph_type& treebank = *__treebank;
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
      
      scores.clear();
      scores.resize(treebank.edges.size());
      
      cicada::inside_outside(treebank, inside, outside, scores, function, function);
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = treebank.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = treebank.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	if (! edge.tails.empty()) continue;
	
	counts[edge.rule->rhs.front()] += scores[edge.id] / inside.back();
      }
    }
  }
  
  queue_type&    queue;
  Function       function;
  
  word_count_set_type counts;
};

template <typename Function>
struct TaskLexiconCount
{
  typedef utils::lockfree_list_queue<const hypergraph_type*, std::allocator<const hypergraph_type*> > queue_type;
  
  // we will collect, tag-signature-word, signature-word, word
  
  TaskLexiconCount(const word_count_set_type& __word_counts, queue_type& __queue, Function __function)
    : word_counts(__word_counts),
      queue(__queue),
      function(__function),
      counts(),
      counts_sig() { }
  
  void operator()()
  {
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    
    weight_set_type inside;
    weight_set_type outside;
    weight_set_type scores;
    
    ngram_type trigram(3);
    ngram_type bigram(2);
    
    const signature_type& __signature = signature_type::create(signature);
    
    const hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      const hypergraph_type& treebank = *__treebank;
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
      
      scores.clear();
      scores.resize(treebank.edges.size());
      
      cicada::inside_outside(treebank, inside, outside, scores, function, function);
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = treebank.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = treebank.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	if (! edge.tails.empty()) continue;
	
	const rule_type& rule = *edge.rule;
	
	// assume penntreebank style...
	bigram[0] = rule.lhs;
	bigram[1] = __signature(rule.rhs.front());

	trigram[0] = bigram[1];
	trigram[1] = bigram[0];
	trigram[2] = rule.rhs.front();
	
	const count_type count = scores[edge.id] / inside.back();
	
	counts[trigram] += count;
	counts[ngram_type(trigram.begin() + 1, trigram.end())] += count;
	counts[ngram_type(trigram.begin() + 2, trigram.end())] += count;
	
	counts_sig[bigram] += count;
	counts_sig[ngram_type(bigram.begin() + 1, bigram.end())] += count;

	word_count_set_type::const_iterator witer = word_counts.find(rule.rhs.front());
	if (witer == word_counts.end())
	  throw std::runtime_error("invalid word???");
	
	if (witer->second <= unknown_threshold)
	  counts_unknown[bigram] += count;
      }
    }
  }
  
  const word_count_set_type& word_counts;
  queue_type&    queue;
  Function       function;
  
  ngram_count_set_type counts; // counts of tag-signature-word
  ngram_count_set_type counts_sig; // counts of tag-signature
  ngram_count_set_type counts_unknown; // count of unknown pair of tag-signature
};

struct LexiconEstimate
{
  typedef std::vector<logprob_type, std::allocator<logprob_type> > logprob_set_type;
  
  typedef std::vector<const ngram_count_set_type::value_type*, std::allocator<const ngram_count_set_type::value_type*> > ngram_set_type;
  typedef google::dense_hash_map<ngram_type, ngram_set_type, boost::hash<ngram_type>, std::equal_to<ngram_type> > ngram_count_map_type;
  
  LexiconEstimate(const double& __prior, const int __order) : prior(__prior), order(__order) {}
  
  double operator()(const ngram_count_set_type& counts, ngram_count_set_type& model, ngram_count_set_type& backoff)
  {
    using namespace boost::math::policies;
    typedef policy<domain_error<errno_on_error>,
		   pole_error<errno_on_error>,
		   overflow_error<errno_on_error>,
		   rounding_error<errno_on_error>,
		   evaluation_error<errno_on_error> > policy_type;
  
    double total = 0.0;
    size_t vocab_size = 0;
    
    ngram_set_type ngrams_local;
    {
      ngram_count_set_type::const_iterator niter_end = counts.end();
      for (ngram_count_set_type::const_iterator niter = counts.begin(); niter != niter_end; ++ niter)
	if (niter->first.size() == 1) {
	  total += niter->second;
	  ++ vocab_size;
	  ngrams_local.push_back(&(*niter));
	}
    }
    
    logprob_set_type logprobs_local(ngrams_local.size());
    
    // we will loop, increment total until we have enough mass discounted...
    double discount = 0.0;
    for (;;) {
      discount = 0.0;
      
      double logprob_sum = boost::numeric::bounds<double>::lowest();
      const double lognorm = utils::mathop::digamma(prior_lexicon * vocab_size + total);
      
      logprob_set_type::iterator liter = logprobs_local.begin();
      ngram_set_type::const_iterator niter_end = ngrams_local.end();
      for (ngram_set_type::const_iterator niter = ngrams_local.begin(); niter != niter_end; ++ niter, ++ liter) {
	const double logprob = utils::mathop::digamma(prior_lexicon + (*niter)->second) - lognorm;
	logprob_sum = utils::mathop::logsum(logprob_sum, logprob);
	*liter = logprob;
      }
      
      discount = - boost::math::expm1(logprob_sum, policy_type());
      
      if (discount > 0.0) break;
      ++ total;
    }
    
    {
      // copy logprob into actual storage..
      logprob_set_type::iterator liter = logprobs_local.begin();
      ngram_set_type::const_iterator niter_end = ngrams_local.end();
      for (ngram_set_type::const_iterator niter = ngrams_local.begin(); niter != niter_end; ++ niter, ++ liter)
	model[(*niter)->first] = *liter;
    }
    
    const double logprob_unk = utils::mathop::log(discount);

    ngram_count_map_type ngrams;
    ngrams.set_empty_key(ngram_type());
        
    for (int n = 2; n <= order; ++ n) {
      ngrams.clear();
      
      {
	ngram_count_set_type::const_iterator niter_end = counts.end();
	for (ngram_count_set_type::const_iterator niter = counts.begin(); niter != niter_end; ++ niter)
	  if (static_cast<int>(niter->first.size()) == n) 
	    ngrams[ngram_type(niter->first.begin(), niter->first.begin() + n - 1)].push_back(&(*niter));
      }
      
      ngram_count_map_type::const_iterator citer_end = ngrams.end();
      for (ngram_count_map_type::const_iterator citer = ngrams.begin(); citer != citer_end; ++ citer) {
	const ngram_set_type& ngrams_local = citer->second;
	logprobs_local.resize(ngrams_local.size());
	
	double total = 0.0;
	double logsum_lower = boost::numeric::bounds<double>::lowest();
	ngram_set_type::const_iterator niter_end = ngrams_local.end();
	for (ngram_set_type::const_iterator niter = ngrams_local.begin(); niter != niter_end; ++ niter) {
	  total += (*niter)->second;
	  
	  ngram_count_set_type::const_iterator liter = model.find(ngram_type((*niter)->first.end() - n + 1, (*niter)->first.end()));
	  if (liter == model.end())
	    throw std::runtime_error("invalid lower order count: " + utils::lexical_cast<std::string>((*niter)->first));
	  
	  logsum_lower = utils::mathop::logsum(logsum_lower, liter->second);
	}
	
	const double discount_lower = - boost::math::expm1(logsum_lower, policy_type());
	double discount = 0.0;
	
	for (;;) {
	  discount = 0.0;
	
	  double logprob_sum = boost::numeric::bounds<double>::lowest();
	  const double lognorm = utils::mathop::digamma(prior_lexicon * ngrams_local.size() + total);
	  
	  logprob_set_type::iterator liter = logprobs_local.begin();
	  for (ngram_set_type::const_iterator niter = ngrams_local.begin(); niter != niter_end; ++ niter, ++ liter) {
	    const double logprob = utils::mathop::digamma(prior_lexicon + (*niter)->second) - lognorm;
	    logprob_sum = utils::mathop::logsum(logprob_sum, logprob);
	    *liter = logprob;
	  }
	  
	  discount = - boost::math::expm1(logprob_sum, policy_type());
	  
	  if (discount > 0.0) break;
	  ++ total;
	}
	
	// copy logprob into actual storage..
	logprob_set_type::const_iterator liter = logprobs_local.begin();
	for (ngram_set_type::const_iterator niter = ngrams_local.begin(); niter != niter_end; ++ niter, ++ liter)
	  model[(*niter)->first] = *liter;
	
	backoff[citer->first] = utils::mathop::log(discount) -  utils::mathop::log(discount_lower);
      }
    }

    return logprob_unk;
  }
  
  const double prior;
  const int order;
};

template <typename Function>
void lexicon_learn(const hypergraph_set_type& treebanks,
		   grammar_type& lexicon,
		   Function function)
{
  // we will learn a trigram of tag-signature-word, but dump tag-signature only...
  // 
  // a trick is: since we are learning OOV probability, tag-signature-word trigram
  // will always backoff (with tag-signature backoff penalty) to signature-word which will
  // always backoff (with signature backoff penalty) to uniform distribution...
  //
  // Thus, we simply preserve tag-signature bigram probability with the penalties 
  // consisting of tag-sinature backoff + signarue backoff + uniform distribution
  
  typedef TaskLexiconFrequency<Function> task_frequency_type;
  typedef TaskLexiconCount<Function>     task_count_type;
  
  typedef std::vector<task_frequency_type, std::allocator<task_frequency_type> > task_frequency_set_type;
  typedef std::vector<task_count_type, std::allocator<task_count_type> >         task_count_set_type;
  
  typedef typename task_frequency_type::queue_type queue_frequency_type;
  typedef typename task_count_type::queue_type     queue_count_type;

  queue_frequency_type queue_frequency;
  task_frequency_set_type tasks_frequency(threads, task_frequency_type(queue_frequency, function));
  
  boost::thread_group workers_frequency;
  for (int i = 0; i != threads; ++ i)
    workers_frequency.add_thread(new boost::thread(boost::ref(tasks_frequency[i])));
  
  hypergraph_set_type::const_iterator titer_end = treebanks.end();
  for (hypergraph_set_type::const_iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue_frequency.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue_frequency.push(0);
  
  workers_frequency.join_all();
  
  word_count_set_type word_counts;
  for (int i = 0; i != threads; ++ i) {
    if (word_counts.empty())
      word_counts.swap(tasks_frequency[i].counts);
    else {
      word_count_set_type::const_iterator witer_end = tasks_frequency[i].counts.end();
      for (word_count_set_type::const_iterator witer = tasks_frequency[i].counts.begin(); witer != witer_end; ++ witer)
	word_counts[witer->first] += witer->second;
    }
  }

  queue_count_type queue_count;
  task_count_set_type tasks_count(threads, task_count_type(word_counts, queue_count, function));
  
  boost::thread_group workers_count;
  for (int i = 0; i != threads; ++ i)
    workers_count.add_thread(new boost::thread(boost::ref(tasks_count[i])));
  
  for (hypergraph_set_type::const_iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue_count.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue_count.push(0);
  
  workers_count.join_all();
  
  ngram_count_set_type counts;
  ngram_count_set_type counts_sig;
  ngram_count_set_type counts_unknown;
  for (int i = 0; i != threads; ++ i) {
    if (counts.empty())
      counts.swap(tasks_count[i].counts);
    else {
      ngram_count_set_type::const_iterator niter_end = tasks_count[i].counts.end();
      for (ngram_count_set_type::const_iterator niter = tasks_count[i].counts.begin(); niter != niter_end; ++ niter)
	counts[niter->first] += niter->second;
    }
    
    if (counts_sig.empty())
      counts_sig.swap(tasks_count[i].counts_sig);
    else {
      ngram_count_set_type::const_iterator niter_end = tasks_count[i].counts_sig.end();
      for (ngram_count_set_type::const_iterator niter = tasks_count[i].counts_sig.begin(); niter != niter_end; ++ niter)
	counts_sig[niter->first] += niter->second;
    }
    
    if (counts_unknown.empty())
      counts_unknown.swap(tasks_count[i].counts_unknown);
    else {
      ngram_count_set_type::const_iterator niter_end = tasks_count[i].counts_unknown.end();
      for (ngram_count_set_type::const_iterator niter = tasks_count[i].counts_unknown.begin(); niter != niter_end; ++ niter)
	counts_unknown[niter->first] += niter->second;
    }
  }

  ngram_count_set_type model;
  ngram_count_set_type backoff;

  ngram_count_set_type model_sig;
  ngram_count_set_type backoff_sig;
  
  // estimate for tag-sig-word
  const double logprob_unk = LexiconEstimate(prior_lexicon, 3)(counts, model, backoff);
  
  //std::cerr << "logprob-unk: " << logprob_unk << std::endl;

  // estimate for tag-sig
  const double logprob_unk_sig = LexiconEstimate(prior_signature, 2)(counts_sig, model_sig, backoff_sig);

  //std::cerr << "logprob-unk-sig: " << logprob_unk_sig << std::endl;
  
  // finished computation...
  // actually, we do not need full-trigram!
  // we need: tag-signature-<UNK>
  // which will be computed by bakoff(tag-signature) + backoff(signature) + unk which is stored in backoffs!
  
  // also, we need tag-signature probability
  // and probability mass left out from KNOWN lexicon...
  
  //
  // From lexicon, compute discounted mass in fully observed counts, that is (probably) available in lexicon...
  //

  ngram_count_set_type model_tag;
  grammar_type::const_iterator liter_end = lexicon.end();
  for (grammar_type::const_iterator liter = lexicon.begin(); liter != liter_end; ++ liter) {
    const symbol_type& lhs = liter->first->lhs;
    
    std::pair<ngram_count_set_type::iterator, bool> result = model_tag.insert(std::make_pair(ngram_type(1, lhs), liter->second));
    if (! result.second)
      result.first->second = utils::mathop::logsum(result.first->second, liter->second);
  }

  ngram_type bigram(2);
  
  ngram_count_set_type::const_iterator biter_end = backoff.end();
  for (ngram_count_set_type::const_iterator biter = backoff.begin(); biter != biter_end; ++ biter) 
    if (biter->first.size() == 2) {
      using namespace boost::math::policies;
      typedef policy<domain_error<errno_on_error>,
	pole_error<errno_on_error>,
	overflow_error<errno_on_error>,
	rounding_error<errno_on_error>,
	evaluation_error<errno_on_error> > policy_type;
      
      // swap backoff context..!
      bigram[0] = biter->first[1];
      bigram[1] = biter->first[0];
      
      // check if this is really unknown rule...
      if (counts_unknown.find(bigram) == counts_unknown.end()) continue;
      
      ngram_count_set_type::const_iterator tag_iter = model_tag.find(ngram_type(1, bigram.front()));
      if (tag_iter == model_tag.end())
	throw std::runtime_error("invalid tag model!?");
      
      ngram_count_set_type::const_iterator sig_iter = model_sig.find(bigram);
      if (sig_iter == model_sig.end())
	throw std::runtime_error("invalid signature model!?");
      
      ngram_count_set_type::const_iterator siter = backoff.find(ngram_type(1, bigram.front()));
      if (siter == backoff.end())
	throw std::runtime_error("invalid backoffs!");
      
      const logprob_type score_backoff = utils::mathop::log(- boost::math::expm1(tag_iter->second, policy_type()));
      const logprob_type score_bigram  = sig_iter->second;
      const logprob_type score_trigram = biter->second + siter->second + logprob_unk;
      
      //std::cerr << "backoff: " << score_backoff << " bigram: " << score_bigram << " trigram: " << score_trigram << std::endl;
      //std::cerr << "biter: " << biter->second << " siter: " << siter->second << std::endl;

      const logprob_type score = score_trigram + score_bigram + score_backoff;
      
      lexicon.insert(std::make_pair(rule_type::create(rule_type(bigram.front(), rule_type::symbol_set_type(1, bigram.back()))), score));
    }
}

template <typename Function>
struct TaskCharacterCount
{
  typedef utils::lockfree_list_queue<const hypergraph_type*, std::allocator<const hypergraph_type*> > queue_type;
  
  // we will collect, tag-signature-word, signature-word, word
  
  TaskCharacterCount(queue_type& __queue, Function __function)
    : queue(__queue),
      function(__function),
      counts() { }
  
  void operator()()
  {
    typedef std::vector<weight_type, std::allocator<weight_type> > weight_set_type;
    typedef cicada::Sentence phrase_type;
    
    weight_set_type inside;
    weight_set_type outside;
    weight_set_type scores;
    
    ngram_type ngram(3);

    phrase_type phrase(1);
    phrase_type tokenized;
    
    const signature_type& __signature = signature_type::create(signature);
    const tokenizer_type& __tokenizer = tokenizer_type::create("character");
    
    const hypergraph_type* __treebank = 0;
    for (;;) {
      queue.pop(__treebank);
      if (! __treebank) break;
      
      const hypergraph_type& treebank = *__treebank;
      
      inside.clear();
      outside.clear();
      inside.resize(treebank.nodes.size());
      outside.resize(treebank.nodes.size());
      
      scores.clear();
      scores.resize(treebank.edges.size());
      
      cicada::inside_outside(treebank, inside, outside, scores, function, function);
      
      hypergraph_type::edge_set_type::const_iterator eiter_end = treebank.edges.end();
      for (hypergraph_type::edge_set_type::const_iterator eiter = treebank.edges.begin(); eiter != eiter_end; ++ eiter) {
	const hypergraph_type::edge_type& edge = *eiter;
	if (! edge.tails.empty()) continue;
	
	const rule_type& rule = *edge.rule;

	const count_type count = scores[edge.id] / inside.back();
	
	// assume penntreebank style...
	ngram[0] = __signature(rule.rhs.front());
	ngram[1] = rule.lhs;
	phrase.front() = rule.rhs.front();
	
	__tokenizer(phrase, tokenized);
	phrase_type::const_iterator titer_end = tokenized.end();
	for (phrase_type::const_iterator titer = tokenized.begin(); titer != titer_end; ++ titer) {
	  ngram[2] = *titer;
	  
	  counts[ngram] += count;
	  counts[ngram_type(ngram.begin() + 1, ngram.end())] += count;
	  counts[ngram_type(ngram.begin() + 2, ngram.end())] += count;
	}
      }
    }
  }
  
  queue_type&    queue;
  Function       function;
  
  ngram_count_set_type counts; // counts of tag-signature-word
};

template <typename Function>
void characters_learn(const hypergraph_set_type& treebanks,
		      ngram_count_set_type& model,
		      ngram_count_set_type& backoff,
		      Function function)
{
  typedef cicada::Vocab vocab_type;

  typedef TaskCharacterCount<Function> task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  typedef typename task_type::queue_type queue_type;

  queue_type queue;
  task_set_type tasks(threads, task_type(queue, function));

  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  hypergraph_set_type::const_iterator titer_end = treebanks.end();
  for (hypergraph_set_type::const_iterator titer = treebanks.begin(); titer != titer_end; ++ titer)
    queue.push(&(*titer));
  
  for (int i = 0; i != threads; ++ i)
    queue.push(0);
  
  workers.join_all();
  
  ngram_count_set_type counts;
  for (int i = 0; i != threads; ++ i) {
    if (counts.empty())
      counts.swap(tasks[i].counts);
    else {
      ngram_count_set_type::const_iterator niter_end = tasks[i].counts.end();
      for (ngram_count_set_type::const_iterator niter = tasks[i].counts.begin(); niter != niter_end; ++ niter)
	counts[niter->first] += niter->second;
    }
  }
  
  // estimate for tag-sig-word
  const double logprob_unk = LexiconEstimate(prior_character, 3)(counts, model, backoff);
  
  model[ngram_type(1, vocab_type::UNK)] = logprob_unk;
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
void grammar_maximize(const count_set_type& counts,
		      grammar_type& grammar,
		      Maximizer maximizer)
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
  for (int i = 0; i != threads; ++ i) {
    if (grammar.empty())
      grammar.swap(tasks[i].grammar);
    else
      grammar.insert(tasks[i].grammar.begin(), tasks[i].grammar.end());
  }
}

inline
bool is_sgml_tag(const symbol_type& symbol)
{
  const size_t size = symbol.size();
  return size != 0 && symbol[0] == '<' && symbol[size - 1] == '>';
}

void write_characters(const path_type& file,
		      const ngram_count_set_type& model,
		      const ngram_count_set_type& backoff,
		      const double cutoff)
{
  typedef std::vector<const ngram_count_set_type::value_type*, std::allocator<const ngram_count_set_type::value_type*> > sorted_type;
  typedef google::dense_hash_map<ngram_type, sorted_type, boost::hash<ngram_type>, std::equal_to<ngram_type> > sorted_map_type;
  
  utils::compress_ostream os(file, 1024 * 1024);
  os.precision(10);
  
  sorted_map_type sorted;
  sorted.set_empty_key(ngram_type(1, symbol_type()));
  
  ngram_count_set_type::const_iterator biter_end = backoff.end();
  for (ngram_count_set_type::const_iterator biter = backoff.begin(); biter != biter_end; ++ biter)
    sorted[ngram_type(biter->first.begin(), biter->first.end() - 1)].push_back(&(*biter));
  
  for (int order = 1; order <= 2; ++ order) {
    sorted_map_type::iterator biter_end = sorted.end();
    for (sorted_map_type::iterator biter = sorted.begin(); biter != biter_end; ++ biter)
      if (static_cast<int>(biter->first.size()) == order - 1) {
	std::sort(biter->second.begin(), biter->second.end(), greater_ptr_second<ngram_count_set_type::value_type>());
	
	sorted_type::const_iterator siter_end = biter->second.end();
	for (sorted_type::const_iterator siter = biter->second.begin(); siter != siter_end; ++ siter)
	  os << "backoff: " << (*siter)->first << ' ' << (*siter)->second << '\n';
      }
  }
  
  sorted.clear();
  ngram_count_set_type::const_iterator miter_end = model.end();
  for (ngram_count_set_type::const_iterator miter = model.begin(); miter != miter_end; ++ miter)
    sorted[ngram_type(miter->first.begin(), miter->first.end() - 1)].push_back(&(*miter));
  
  for (int order = 1; order <= 3; ++ order) {
    sorted_map_type::iterator miter_end = sorted.end();
    for (sorted_map_type::iterator miter = sorted.begin(); miter != miter_end; ++ miter)
      if (static_cast<int>(miter->first.size()) == order - 1) {
	std::sort(miter->second.begin(), miter->second.end(), greater_ptr_second<ngram_count_set_type::value_type>());
	
	sorted_type::const_iterator siter_end = miter->second.end();
	for (sorted_type::const_iterator siter = miter->second.begin(); siter != siter_end; ++ siter)
	  os << "model: " << (*siter)->first << ' ' << (*siter)->second << '\n';
      }
  }
}

void grammar_prune(grammar_type& grammar, const double cutoff)
{
  typedef std::vector<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > sorted_type;
  
  typedef std::pair<rule_ptr_type, logprob_type> rule_logprob_type;
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<symbol_type, rule_logprob_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, rule_logprob_type> > > reachable_set_type;
#else
  typedef sgi::hash_map<symbol_type, rule_logprob_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
			std::allocator<std::pair<const symbol_type, rule_logprob_type> > > reachable_set_type;
#endif
  // we will first compute "reachable" rules...
  // reachable label -> rule mapping
  
  count_set_type counts;
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    counts[giter->first->lhs].insert(*giter);
  
  reachable_set_type reachables;
  reachable_set_type reachables_next;
  
  const symbol_type goal("[ROOT]");
  
  reachables[goal];
  
  for (;;) {
    bool equilibrate = true;
    
    reachables_next.clear();
    reachables_next = reachables;
    
    reachable_set_type::const_iterator riter_end = reachables.end();
    for (reachable_set_type::const_iterator riter = reachables.begin(); riter != riter_end; ++ riter) {
      count_set_type::const_iterator citer = counts.find(riter->first);
      if (citer == counts.end()) continue; // ignore lexical rule, preterminals
      
      const grammar_type& grammar = citer->second;
      
      grammar_type::const_iterator giter_end = grammar.end();
      for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter) {
      
	symbol_set_type::const_iterator siter_end = giter->first->rhs.end();
	for (symbol_set_type::const_iterator siter = giter->first->rhs.begin(); siter != siter_end; ++ siter)
	  if (siter->is_non_terminal()) {
	    // we will keep the best rule....
	    std::pair<reachable_set_type::iterator, bool> result = reachables_next.insert(std::make_pair(*siter, *giter));
	    if (result.second)
	      equilibrate = false;
	    else if (giter->second > result.first->second.second)
	      result.first->second = *giter;
	  }
      }
    }
    
    reachables.swap(reachables_next);
    
    if (equilibrate) break;
  }
  
  // reachables are set of rules we "must" preserve, and should not be pruned away...
  
  grammar.clear();
  reachable_set_type::const_iterator riter_end = reachables.end();
  for (reachable_set_type::const_iterator riter = reachables.begin(); riter != riter_end; ++ riter)
    if (riter->second.first != rule_ptr_type())
      grammar.insert(riter->second);
  
  const double logcutoff = utils::mathop::log(cutoff);
  sorted_type sorted;
  
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
    const grammar_type& grammar_local = citer->second;
    
    sorted.clear();
    grammar_type::const_iterator giter_end = grammar_local.end();
    for (grammar_type::const_iterator giter = grammar_local.begin(); giter != giter_end; ++ giter)
      sorted.push_back(&(*giter));
    
    std::sort(sorted.begin(), sorted.end(), greater_ptr_second<grammar_type::value_type>());
    
    const double logprob_max = sorted.front()->second;
    const double logprob_threshold = logprob_max + logcutoff;
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end && (*siter)->second >= logprob_threshold; ++ siter)
      grammar.insert(*(*siter));
  }
}

void lexicon_prune(grammar_type& grammar, const double cutoff)
{
  typedef std::vector<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > sorted_type;

  count_set_type counts;
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    counts[giter->first->lhs].insert(*giter);
  
  grammar.clear();
  
  const double logcutoff = utils::mathop::log(cutoff);
  sorted_type sorted;
  
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
    const grammar_type& grammar_local = citer->second;
    
    sorted.clear();
    grammar_type::const_iterator giter_end = grammar_local.end();
    for (grammar_type::const_iterator giter = grammar_local.begin(); giter != giter_end; ++ giter)
      sorted.push_back(&(*giter));
    
    std::sort(sorted.begin(), sorted.end(), greater_ptr_second<grammar_type::value_type>());
    
    const double logprob_max = sorted.front()->second;
    const double logprob_threshold = logprob_max + logcutoff;
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end && (*siter)->second >= logprob_threshold; ++ siter)
      grammar.insert(*(*siter));
  }
}

void write_grammar(const path_type& file,
		   const grammar_type& grammar)
{
  typedef std::vector<const grammar_type::value_type*, std::allocator<const grammar_type::value_type*> > sorted_type;

  if (grammar.empty()) return;
  
  count_set_type counts;
  sorted_type sorted;
    
  grammar_type::const_iterator giter_end = grammar.end();
  for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
    counts[giter->first->lhs][giter->first] = giter->second;

  utils::compress_ostream os(file, 1024 * 1024);
  os.precision(10);
    
  count_set_type::const_iterator citer_end = counts.end();
  for (count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
    const grammar_type& grammar = citer->second;
      
    sorted.clear();
      
    grammar_type::const_iterator giter_end = grammar.end();
    for (grammar_type::const_iterator giter = grammar.begin(); giter != giter_end; ++ giter)
      sorted.push_back(&(*giter));
      
    std::sort(sorted.begin(), sorted.end(), greater_ptr_second<grammar_type::value_type>());
    
    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      os << *((*siter)->first) << " ||| ||| " << (*siter)->second << '\n';
  }
}


void read_treebank(const path_set_type& files,
		   hypergraph_set_type& treebanks)
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
    ("output-grammar",   po::value<path_type>(&output_grammar_file),   "output grammar")
    ("output-lexicon",   po::value<path_type>(&output_lexicon_file),   "output lexical rules")
    ("output-character", po::value<path_type>(&output_character_file), "output character model")
    
    ("max-iteration",       po::value<int>(&max_iteration)->default_value(max_iteration),             "maximum split/merge iterations")
    ("max-iteration-split", po::value<int>(&max_iteration_split)->default_value(max_iteration_split), "maximum EM iterations after split")
    ("max-iteration-merge", po::value<int>(&max_iteration_merge)->default_value(max_iteration_merge), "maximum EM iterations after merge")
    
    ("binarize-left",  po::bool_switch(&binarize_left),  "left binarization")
    ("binarize-right", po::bool_switch(&binarize_right), "right binarization")
    ("binarize-all",   po::bool_switch(&binarize_all),   "all binarization")
    
    ("prior-rule",      po::value<double>(&prior_rule)->default_value(prior_rule),           "Dirichlet prior for rules")
    ("prior-lexicon",   po::value<double>(&prior_lexicon)->default_value(prior_lexicon),     "Dirichlet prior for lexical rule")
    ("prior-signature", po::value<double>(&prior_signature)->default_value(prior_signature), "Dirichlet prior for signature")
    ("prior-character", po::value<double>(&prior_character)->default_value(prior_character), "Dirichlet prior for character")

    ("cutoff-rule",      po::value<double>(&cutoff_rule)->default_value(cutoff_rule),           "cutoff for rules")
    ("cutoff-lexicon",   po::value<double>(&cutoff_lexicon)->default_value(cutoff_lexicon),     "cutoff for lexical rule")
    ("cutoff-character", po::value<double>(&cutoff_character)->default_value(cutoff_character), "cutoff for character")
    
    ("merge-ratio",   po::value<double>(&merge_ratio)->default_value(merge_ratio),     "merging ratio")
    
    ("unknown-ratio",     po::value<double>(&unknown_ratio)->default_value(unknown_ratio),         "unknown word ratio")
    ("unknown-threshold", po::value<double>(&unknown_threshold)->default_value(unknown_threshold), "unknown word threshold")
    
    ("signature",      po::value<std::string>(&signature), "signature for unknown word")
    ("signature-list", po::bool_switch(&signature_list),   "list of signatures")
    
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
