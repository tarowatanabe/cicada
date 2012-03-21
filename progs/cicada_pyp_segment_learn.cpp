//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-segmenter!
//
// @InProceedings{mochihashi-yamada-ueda:2009:ACLIJCNLP,
//  author    = {Mochihashi, Daichi  and  Yamada, Takeshi  and  Ueda, Naonori},
//  title     = {Bayesian Unsupervised Word Segmentation with Nested Pitman-Yor Language Modeling},
//  booktitle = {Proceedings of the Joint Conference of the 47th Annual Meeting of the ACL and the 4th International Joint Conference on Natural Language Processing of the AFNLP},
//  month     = {August},
//  year      = {2009},
//  address   = {Suntec, Singapore},
//  publisher = {Association for Computational Linguistics},
//  pages     = {100--108},
//  url       = {http://www.aclweb.org/anthology/P/P09/P09-1012}
// }

//
// we will use a tabular-based DP, thus, the order must be fixed!
// (or, do we augment extra contexts using vectors????)
//
// we will provide two ngram LM structure, 
// one for character-based one (we will use utils::piece) and ngram-based one also using utils::piece as a storage...
// For utils::piece details, seed the implementation in cicada_translit_learn_pyp.cpp
// 


#include <map>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/semiring/logprob.hpp>

#include "utils/chunk_vector.hpp"
#include "utils/utf8.hpp"
#inlcude "utils/array_power2.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/chinese_restaurant_process.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

typedef cicada::Vocab     vocab_type;

// PYP Word model.... this is basically the same as the PYPLM with additional length model...
struct PYP
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef std::string sentence_type;
  
  typedef utils::piece piece_type;
  typedef utils::piece segment_type;
  typedef utils::piece word_type;
};

struct PYPWord
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::sentence_type sentence_type;
  typedef PYP::piece_type    piece_type;
  typedef PYP::segment_type  segment_type;
  typedef PYP::word_type     word_type;
  
  typedef uint32_t  id_type;
  
  struct Node
  {
    typedef utils::chinese_restaurant_process<word_type, boost::hash<word_type>, std::equal_to<word_type>,
					      std::allocator<word_type > > table_type;
  
    Node() : table(), parent(id_type(-1)), order(0)  {}
    
    table_type table;
    id_type parent;
    int     order;
  };
  typedef Node node_type;
  
  typedef utils::compact_trie_dense<word_type, node_type, boost::hash<word_type>, std::equal_to<word_type>,
				    std::allocator<std::pair<const word_type, node_type> > > trie_type;
  
  
  typedef std::vector<double, std::allocator<double> > parameter_set_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > node_set_type;
  typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;

  typedef std::vector<segment_type, std::allocator<segment_type> > buffer_type;
  
  static segment_type BOS()
  {
    static std::string __bos = "\n";
    return __bos;
  }

  static segment_type EOS()
  {
    static std::string __eos = "\n";
    return __eos;
  }
  
  struct length_base_type
  {
    typedef utils::array_power2<double, 64, std::allocator<double> > cache_type;

    length_base_type(const double& __lambda,
		     const double& __strength_shape,
		     const double& __strength_rate)
      : lambda(__lambda),
	strength_shape(__strength_shape),
	strength_rate(__strength_rate)
    {
      initialize(__lambda);
    }
    
    double logprob(const size_type size) const
    {
      if (size < cache.size())
	return cache[size];
      else
	return utis::mathop::log_poisson(size, lambda);
    }

    template <typename Sampler>
    void sample_parameters(const double& word_size, const double& table_size, Sample& sampler)
    {
      initialize(sampler.gamma(strength_shape + word_size, strength_rate + table_size));
    }
    
    void initialize(const double& __lambda)
    {
      cache.clear();
      lambda = __lambda;
      
      cache[0] = - std::numeric_limits<double>::infinity();
      for (size_type size = 1; size != cache.size(); ++ size)
	cache[size] = utis::mathop::log_poisson(size, lambda);
    }

    cache_type cache;
    
    double lambda;
    
    double strength_shape;
    double strength_rate;
  };
  
  PYPWord(const int order,
	  const double __p0,
	  const double __discount,
	  const double __strength,
	  const double __discount_alpha,
	  const double __discount_beta,
	  const double __strength_shape,
	  const double __strength_rate,
	  const double __lambda,
	  const double __lambda_strength,
	  const double __lambda_rate)
    : trie(word_type()),
      nodes(order),
      discount(order, __discount),
      strength(order, __strength),
      discount_alpha(__discount_alpha),
      discount_beta(__discount_beta),
      strength_shape(__strength_shape),
      strength_rate(__strength_rate),
      p0(__p0),
      counts0(0),
      base(__lambda, __lambda_strength, __lambda_rate)
  {
    // unitialize root table...
    root.parent = id_type(-1);
    root.order = 0;
    root.table = node_type::table_type(discount[0], strength[0]);
  }
  
  template <typename Iterator, typename Sampler>
  bool increment(const segment_type& segment, Iterator first, Iterator last, Sampler& sampler, const double temperature=1.0)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    int order = 1;
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
      const id_type node_prev = node;
      
      node = trie.insert(node, *iter);
      
      node_type& trie_node = trie[node];
      
      if (! trie_node.order) {
	trie_node.parent = node_prev;
	trie_node.order = order;
	trie_node.table = node_type::table_type(discount[order], strength[order]);
	
	nodes[order].push_back(node);
      }
    }
    
    if (node == trie.root()) {
      if (root.table.increment(segment, p0, sampler, temperature))
	++ counts0;
      else
	return false;
    } else {
      const double backoff = prob(segment, trie[node].parent);
      
      // we will also increment lower-order when new table is created!
      if (trie[node].table.increment(segment, backoff, sampler, temperature))
	increment(segment, trie[node].parent, sampler, temperature);
      else
	return false;
    }
    
    return true;
  }

  template <typename Iterator, typename Sampler>
  bool decrement(const segment_type& segment, Iterator first, Iterator last, Sampler& sampler)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter)
      node = trie.find(node, *iter);
    
    if (node == trie.root()) {
      if (root.table.decrement(segment, sampler))
	-- counts0;
      else
	return false;
    } else {
      if (trie[node].table.decrement(segment, sampler))
	decrement(segment, trie[node].parent, sampler);
      else
	return false;
    }
    
    return true;
  }
  

  double prob(const segment_type& segment, const id_type& node) const
  {
    if (node == trie.root())
      return root.table.prob(segment, p0);
    else {
      const double p = prob(segment, trie[node].parent);
      
      if (trie[node].table.empty())
	return p;
      else
	return trie[node].table.prob(segment, p);
    }
  }
  
  template <typename Iterator>
  double prob(const segment_type& segment, Iterator first, Iterator last)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    // this may happen when accessing 0-gram!
    if (! (first <= last))
      return p0;
    
    double p = root.table.prob(word, p0);
      
    // we will traverse from the back!
    reverse_iterator begin(last);
    reverse_iterator end(first);
      
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter) {
      node = trie.find(node, *iter);
	
      if (node == trie_type::npos() || trie[node].table.empty())
	return p;
      else
	p = trie[node].table.prob(word, p);
    }
      
    return p;
  }
  
  template <typename Sampler>
  void increment(const segment_type& segment, Sample& sampler, const double temperature=1.0)
  {
    const size_type conext_size = discount.size() - 1;
    
    buffer.clear();
    buffer.push_back(BOS());
    
    segment_type::const_iterator siter_end = segment.end();
    for (segment_type::const_iterator siter = segment.begin(); siter != siter_end; /**/) {
      const size_type char_size = utils::utf8_size(*siter);
      const segment_type seg(siter, siter + char_size);
      
      increment(seg, buffer.begin(), std::min(buffer.begin(), buffer.end() - context_size), sampler, temperature);
      
      buffer.push_back(seg);
      siter += char_size;
    }
    
    increment(EOS(), buffer.begin(), std::min(buffer.begin(), buffer.end() - context_size), sampler, temperature);
  }
  
  template <typename Sampler>
  void decrement(const segment_type& segment, Sample& sampler)
  {
    const size_type conext_size = discount.size() - 1;

    buffer.clear();
    buffer.push_back(BOS());
    
    segment_type::const_iterator siter_end = segment.end();
    for (segment_type::const_iterator siter = segment.begin(); siter != siter_end; /**/) {
      const size_type char_size = utils::utf8_size(*siter);
      const segment_type seg(siter, siter + char_size);
      
      decrement(seg, buffer.begin(), std::min(buffer.begin(), buffer.end() - context_size), sampler);
      
      buffer.push_back(seg);
      siter += char_size;
    }
    
    decrement(EOS(), buffer.begin(), std::min(buffer.begin(), buffer.end() - context_size), sampler);
  }
  
  double prob(const segment_type& segment)
  {
    const size_type conext_size = discount.size() - 1;
    
    buffer.clear();
    buffer.push_back(BOS());
    
    double loprob = 0.0;
    
    segment_type::const_iterator siter_end = segment.end();
    for (segment_type::const_iterator siter = segment.begin(); siter != siter_end; /**/) {
      const size_type char_size = utils::utf8_size(*siter);
      const segment_type seg(siter, siter + char_size);
      
      logprob += std::log(prob(seg, buffer.begin(), std::min(buffer.begin(), buffer.end() - context_size)));
      
      buffer.push_back(seg);
      siter += char_size;
    }
    
    logprob += std::log(prob(EOS(), buffer.begin(), std::min(buffer.begin(), buffer.end() - context_size)));
    
    // exclude BOS/EOS
    logprob += base.logprob(buffer.size() - 2);
    
    return std::exp(logprob);
  }

  // TODO: add log_likelihood
  //     : experiment with sampled length parameters... (Do we really need this...? Is it simply a normalization constant...?)
  //

  double log_likelihood() const
  {
    //double logprob = std::log(p0) * counts0;
    
    double logprob = 0.0;
    for (size_type order = 0; order != discount.size(); ++ order)
      logprob += log_likelihood(order, discount[order], strength[order]);
    
    return logprob;
  }
  
  
  double log_likelihood(const int order, const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    double logprob = (utils::mathop::log_beta_density(discount, discount_alpha, discount_beta)
		      + utils::mathop::log_gamma_density(strength + discount, strength_shape, strength_rate));
    
    if (order == 0)
      return logprob + (! root.table.empty() ? root.table.log_likelihood(discount, strength) : 0.0);
    else {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (! trie[*niter].table.empty())
	  logprob += trie[*niter].table.log_likelihood(discount, strength);
      
      return logprob;
    }
  }

  struct DiscountSampler
  {
    DiscountSampler(const PYPWord& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPWord& pyplm;
    int order;
    
    double operator()(const double& proposed_discount) const
    {
      return pyplm.log_likelihood(order, proposed_discount, pyplm.strength[order]);
    }
  };
  
  struct StrengthSampler
  {
    StrengthSampler(const PYPWord& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPWord& pyplm;
    int order;
    
    double operator()(const double& proposed_strength) const
    {
      return pyplm.log_likelihood(order, pyplm.discount[order], proposed_strength);
    }
  };
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    for (size_type order = 0; order != discount.size(); ++ order) {
      for (int iter = 0; iter != num_loop; ++ iter) {
	strength[order] = sample_strength(order, sampler, discount[order], strength[order]);
	
	discount[order] = sample_discount(order, sampler, discount[order], strength[order]);
      }
      
      strength[order] = sample_strength(order, sampler, discount[order], strength[order]);
    }
  }

  template <typename Sampler>
  double sample_strength(const int order, Sampler& sampler, const double& discount, const double& strength) const
  {
    double x = 0.0;
    double y = 0.0;

    if (order == 0) {
      x += root.table.sample_log_x(sampler, discount, strength);
      y += root.table.sample_y(sampler, discount, strength);
    } else {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (! trie[*niter].table.empty()) {
	  x += trie[*niter].table.sample_log_x(sampler, discount, strength);
	  y += trie[*niter].table.sample_y(sampler, discount, strength);
	}
    }
    
    return sampler.gamma(strength_shape + y, strength_rate - x);
  }
  
  template <typename Sampler>
  double sample_discount(const int order, Sampler& sampler, const double& discount, const double& strength) const
  {
    double y = 0.0;
    double z = 0.0;
    
    if (order == 0) {
      y += root.table.sample_y_inv(sampler, discount, strength);
      z += root.table.sample_z_inv(sampler, discount, strength);
    } else {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (! trie[*niter].table.empty()) {
	  y += trie[*niter].table.sample_y_inv(sampler, discount, strength);
	  z += trie[*niter].table.sample_z_inv(sampler, discount, strength);
	}
    }
    
    return sampler.beta(discount_alpha + y, discount_beta + z);
  }
  
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    for (size_type order = 0; order != discount.size(); ++ order) {
      DiscountSampler discount_sampler(*this, order);
      StrengthSampler strength_sampler(*this, order);

      for (int iter = 0; iter != num_loop; ++ iter) {
	strength[order] = utils::slice_sampler(strength_sampler,
					       strength[order],
					       sampler,
					       - discount[order] + std::numeric_limits<double>::min(),
					       std::numeric_limits<double>::infinity(),
					       0.0,
					       num_iterations,
					       100 * num_iterations);
	
	discount[order] = utils::slice_sampler(discount_sampler,
					       discount[order],
					       sampler,
					       (strength[order] < 0.0 ? - strength[order] : 0.0) + std::numeric_limits<double>::min(),
					       1.0,
					       0.0,
					       num_iterations,
					       100 * num_iterations);
      }
      
      strength[order] = utils::slice_sampler(strength_sampler,
					     strength[order],
					     sampler,
					     - discount[order] + std::numeric_limits<double>::min(),
					     std::numeric_limits<double>::infinity(),
					     0.0,
					     num_iterations,
					     100 * num_iterations);
      
      if (order == 0) {
	root.table.discount() = discount[order];
	root.table.strength() = strength[order];

	root.table.verify_parameters();
      } else {
	node_set_type::const_iterator niter_end = nodes[order].end();
	for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
	  trie[*niter].table.discount() = discount[order];
	  trie[*niter].table.strength() = strength[order];
	  
	  trie[*niter].table.verify_parameters();
	}
      }
    }
  }
  
public: 
  trie_type trie;
  node_type root;
  node_map_type nodes;
  
  parameter_set_type discount;
  parameter_set_type strength;
  
  double discount_alpha;
  double discount_beta;
  double strength_shape;
  double strength_rate;
  
  double    p0;
  size_type counts0;
  
  length_base_type base;
  
  buffer_type buffer;
};

struct PYPLM
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::sentence_type sentence_type;
  typedef PYP::piece_type    piece_type;
  typedef PYP::segment_type  segment_type;

  typedef segment_type word_type;

  typedef uint32_t  id_type;
  
  typedef boost::filesystem::path path_type;
  
  struct Node
  {
    typedef utils::chinese_restaurant_process<word_type, boost::hash<word_type>, std::equal_to<word_type>,
					      std::allocator<word_type > > table_type;
  
    Node() : table(), parent(id_type(-1)), order(0)  {}
  
    table_type table;
    id_type parent;
    int     order;
  };
  typedef Node node_type;
  
  typedef utils::compact_trie_dense<word_type, node_type, boost::hash<word_type>, std::equal_to<word_type>,
				    std::allocator<std::pair<const word_type, node_type> > > trie_type;
  
  typedef std::vector<double, std::allocator<double> > parameter_set_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > node_set_type;
  typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;

  
  PYPLM(PYPWord& __base,
	const int order,
	const double __discount,
	const double __strength,
	const double __discount_alpha,
	const double __discount_beta,
	const double __strength_shape,
	const double __strength_rate)
    : base(__base),
      trie(word_type()),
      nodes(order),
      discount(order, __discount),
      strength(order, __strength),
      discount_alpha(__discount_alpha),
      discount_beta(__discount_beta),
      strength_shape(__strength_shape),
      strength_rate(__strength_rate),
      counts0(0)
  {
    // unitialize root table...
    root.parent = id_type(-1);
    root.order = 0;
    root.table = node_type::table_type(discount[0], strength[0]);
  }

  template <typename Iterator>
  id_type insert(Iterator first, Iterator last)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    int order = 1;
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
      const id_type node_prev = node;
      
      node = trie.insert(node, *iter);
      
      node_type& trie_node = trie[node];
      
      if (! trie_node.order) {
	trie_node.parent = node_prev;
	trie_node.order = order;
	trie_node.table = node_type::table_type(discount[order], strength[order]);
	
	nodes[order].push_back(node);
      }
    }
    
    return node;
  }
  
  
  template <typename Sampler>
  bool increment(const word_type& word, const id_type& node, Sampler& sampler, const double temperature=1.0)
  {
    if (node == trie.root()) {
      if (root.table.increment(word, p0, sampler, temperature))
	++ counts0;
      else
	return false;
    } else {
      const double backoff = prob(word, trie[node].parent);
      
      // we will also increment lower-order when new table is created!
      if (trie[node].table.increment(word, backoff, sampler, temperature))
	increment(word, trie[node].parent, sampler, temperature);
      else
	return false;
    }
      
    return true;
  }
  
  template <typename Sampler>
  bool decrement(const word_type& word, const id_type& node, Sampler& sampler)
  {
    if (node == trie.root()) {
      if (root.table.decrement(word, sampler))
	-- counts0;
      else
	return false;
    } else {
      if (trie[node].table.decrement(word, sampler))
	decrement(word, trie[node].parent, sampler);
      else
	return false;
    }
      
    return true;
  }
  
  double prob(const word_type& word, const id_type& node) const
  {
    if (node == trie.root())
      return root.table.prob(word, p0);
    else {
      const double p = prob(word, trie[node].parent);
      
      if (trie[node].table.empty())
	return p;
      else
	return trie[node].table.prob(word, p);
    }
  }

  template <typename Iterator>
  double prob(const word_type& word, Iterator first, Iterator last)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    // this may happen when accessing 0-gram!
    if (! (first <= last))
      return p0;
    
    double p = root.table.prob(word, p0);
    
    // we will traverse from the back!
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter) {
      node = trie.find(node, *iter);
      
      if (node == trie_type::npos() || trie[node].table.empty())
	return p;
      else
	p = trie[node].table.prob(word, p);
    }
    
    return p;
  }
  
  double log_likelihood() const
  {
    double logprob = 0.0;
    for (size_type order = 0; order != discount.size(); ++ order)
      logprob += log_likelihood(order, discount[order], strength[order]);
    
    return logprob;
  }
  
  double log_likelihood(const int order, const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    double logprob = (utils::mathop::log_beta_density(discount, discount_alpha, discount_beta)
		      + utils::mathop::log_gamma_density(strength + discount, strength_shape, strength_rate));
    
    if (order == 0)
      return logprob + (! root.table.empty() ? root.table.log_likelihood(discount, strength) : 0.0);
    else {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (! trie[*niter].table.empty())
	  logprob += trie[*niter].table.log_likelihood(discount, strength);
      
      return logprob;
    }
  }
  
  struct DiscountSampler
  {
    DiscountSampler(const PYPLM& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPLM& pyplm;
    int order;
    
    double operator()(const double& proposed_discount) const
    {
      return pyplm.log_likelihood(order, proposed_discount, pyplm.strength[order]);
    }
  };
  
  struct StrengthSampler
  {
    StrengthSampler(const PYPLM& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPLM& pyplm;
    int order;
    
    double operator()(const double& proposed_strength) const
    {
      return pyplm.log_likelihood(order, pyplm.discount[order], proposed_strength);
    }
  };
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    double word_size = 0.0;
    double table_size = 0.0;
    for (size_type order = 0; order != nodes.size(); ++ order) {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (trie.empty(*niter)) {
	  // we will collect information from nodes without children...
	  
	  typename table_type::const_iterator titer_end = trie[*niter].table.end();
	  for (typename table_type::const_iterator titer = trie[*niter].table.begin(); titer != titer_end; ++ titer)
	    word_size += titer->first.size() * titer->second.size_table();
	  
	  table_size += trie[*niter].table.size_table();
	}
    }
    
    base.base.sample_parameters(word_size, table_size, sampler);
    
    base.sample_parameters(sampler, num_loop, num_iterations);
    
    for (size_type order = 0; order != discount.size(); ++ order) {
      
      for (int iter = 0; iter != num_loop; ++ iter) {
	strength[order] = sample_strength(order, sampler, discount[order], strength[order]);
	
	discount[order] = sample_discount(order, sampler, discount[order], strength[order]);
      }
      
      strength[order] = sample_strength(order, sampler, discount[order], strength[order]);
    }
  }

  template <typename Sampler>
  double sample_strength(const int order, Sampler& sampler, const double& discount, const double& strength) const
  {
    double x = 0.0;
    double y = 0.0;

    if (order == 0) {
      x += root.table.sample_log_x(sampler, discount, strength);
      y += root.table.sample_y(sampler, discount, strength);
    } else {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (! trie[*niter].table.empty()) {
	  x += trie[*niter].table.sample_log_x(sampler, discount, strength);
	  y += trie[*niter].table.sample_y(sampler, discount, strength);
	}
    }
    
    return sampler.gamma(strength_shape + y, strength_rate - x);
  }
  
  template <typename Sampler>
  double sample_discount(const int order, Sampler& sampler, const double& discount, const double& strength) const
  {
    double y = 0.0;
    double z = 0.0;
    
    if (order == 0) {
      y += root.table.sample_y_inv(sampler, discount, strength);
      z += root.table.sample_z_inv(sampler, discount, strength);
    } else {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (! trie[*niter].table.empty()) {
	  y += trie[*niter].table.sample_y_inv(sampler, discount, strength);
	  z += trie[*niter].table.sample_z_inv(sampler, discount, strength);
	}
    }
    
    return sampler.beta(discount_alpha + y, discount_beta + z);
  }
  
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    double word_size = 0.0;
    double table_size = 0.0;
    for (size_type order = 0; order != nodes.size(); ++ order) {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (trie.empty(*niter)) {
	  // we will collect information from nodes without children...
	  
	  typename table_type::const_iterator titer_end = trie[*niter].table.end();
	  for (typename table_type::const_iterator titer = trie[*niter].table.begin(); titer != titer_end; ++ titer)
	    word_size += titer->first.size() * titer->second.size_table();
	  
	  table_size += trie[*niter].table.size_table();
	}
    }
    
    base.base.sample_parameters(word_size, table_size, sampler);

    base.slice_sample_parameters(sampler, num_loop, num_iterations);

    for (size_type order = 0; order != discount.size(); ++ order) {
      DiscountSampler discount_sampler(*this, order);
      StrengthSampler strength_sampler(*this, order);

      for (int iter = 0; iter != num_loop; ++ iter) {
	strength[order] = utils::slice_sampler(strength_sampler,
					       strength[order],
					       sampler,
					       - discount[order] + std::numeric_limits<double>::min(),
					       std::numeric_limits<double>::infinity(),
					       0.0,
					       num_iterations,
					       100 * num_iterations);
	
	discount[order] = utils::slice_sampler(discount_sampler,
					       discount[order],
					       sampler,
					       (strength[order] < 0.0 ? - strength[order] : 0.0) + std::numeric_limits<double>::min(),
					       1.0,
					       0.0,
					       num_iterations,
					       100 * num_iterations);
      }
      
      strength[order] = utils::slice_sampler(strength_sampler,
					     strength[order],
					     sampler,
					     - discount[order] + std::numeric_limits<double>::min(),
					     std::numeric_limits<double>::infinity(),
					     0.0,
					     num_iterations,
					     100 * num_iterations);
      
      if (order == 0) {
	root.table.discount() = discount[order];
	root.table.strength() = strength[order];

	root.table.verify_parameters();
      } else {
	node_set_type::const_iterator niter_end = nodes[order].end();
	for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
	  trie[*niter].table.discount() = discount[order];
	  trie[*niter].table.strength() = strength[order];
	  
	  trie[*niter].table.verify_parameters();
	}
      }
    }
  }

  void write(const path_type& path)
  {
    typedef utils::repository repository_type;
    typedef std::vector<id_type, std::allocator<id_type> > node_set_type;

    typedef uint64_t count_type;

    typedef boost::fusion::tuple<id_type, count_type, count_type> data_type;

    typedef std::map<word_type, data_type, std::less<word_type>, std::allocator<std::pair<const word_type, data_type> >  > word_set_type;
    typedef utils::succinct_vector<std::allocator<int32_t> > position_set_type;
    typedef std::vector<count_type, std::allocator<count_type> > offset_set_type;
      
    repository_type rep(path, repository_type::write);
      
    rep["order"] = boost::lexical_cast<std::string>(discount.size());

    rep["counts0"] = boost::lexical_cast<std::string>(counts0);
    
    rep["discount-alpha"] = boost::lexical_cast<std::string>(discount_alpha);
    rep["discount-beta"]  = boost::lexical_cast<std::string>(discount_beta);
    rep["strength-shape"] = boost::lexical_cast<std::string>(strength_shape);
    rep["strength-rate"]  = boost::lexical_cast<std::string>(strength_rate);
          
    for (size_type order = 0; order != discount.size(); ++ order) {
      rep["discount" + boost::lexical_cast<std::string>(order)] = boost::lexical_cast<std::string>(discount[order]);
      rep["strength" + boost::lexical_cast<std::string>(order)] = boost::lexical_cast<std::string>(strength[order]);
    }
    
    // we will compute on-memory for faster indexing... (and assuming small data!)
    
    boost::iostreams::filtering_ostream os_index;
    boost::iostreams::filtering_ostream os_count;
    boost::iostreams::filtering_ostream os_total;
    
    os_index.push(utils::packed_sink<word_type::id_type, std::allocator<word_type::id_type> >(rep.path("index")));
    os_count.push(utils::packed_sink<count_type, std::allocator<count_type> >(rep.path("count")));
    os_total.push(utils::packed_sink<count_type, std::allocator<count_type> >(rep.path("total")));
    
    // dump total counts!!!!!!!
      
    position_set_type positions;
    offset_set_type   offsets(1, 0);
    count_type        offset = 0;
    
    node_set_type nodes;
    node_set_type nodes_next;
    
    word_set_type words;
    
    // unigram!
    {
      const count_type count_customer = root.table.size_customer();
      const count_type count_table    = root.table.size_table();
      
      os_total.write((char*) &count_customer, sizeof(count_type));
      os_total.write((char*) &count_table, sizeof(count_type));

      node_type::table_type::const_iterator titer_end = root.table.end();
      for (node_type::table_type::const_iterator titer = root.table.begin(); titer != titer_end; ++ titer)
	words.insert(std::make_pair(titer->first, data_type(trie_type::npos(),
							    titer->second.size_customer(),
							    titer->second.size_table())));
      
      trie_type::const_iterator iter_end = trie.end();
      for (trie_type::const_iterator iter = trie.begin(); iter != iter_end; ++ iter)
	if (! trie[iter->second].table.empty()) {
	  std::pair<word_set_type::iterator, bool> result = words.insert(std::make_pair(iter->first, data_type(iter->second, 0, 0)));
	  if (! result.second)
	    boost::fusion::get<0>(result.first->second) = iter->second;
	}
      
      word_set_type::const_iterator witer_end = words.end();
      for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) {
	nodes.push_back(boost::fusion::get<0>(witer->second));
	
	word_type::id_type word_id = witer->first.id();
	os_index.write((char*) &word_id, sizeof(word_type::id_type));
	os_count.write((char*) &boost::fusion::get<1>(witer->second), sizeof(count_type));
	os_count.write((char*) &boost::fusion::get<2>(witer->second), sizeof(count_type));
      }
      
      offset += words.size();
      offsets.push_back(offset);
    }
    
    
    for (size_type order = 1; order < discount.size(); ++ order) {
      nodes_next.clear();
      
      node_set_type::const_iterator niter_end = nodes.end();
      for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	if (*niter == trie_type::npos()) {
	  const count_type count_customer = 0;
	  const count_type count_table    = 0;
	  
	  os_total.write((char*) &count_customer, sizeof(count_type));
	  os_total.write((char*) &count_table, sizeof(count_type));
	  
	  positions.set(positions.size(), false); // we will set bit!
	} else {
	  const node_type& node = trie[*niter];
	  trie_type::const_iterator iter_begin = trie.begin(*niter);
	  trie_type::const_iterator iter_end   = trie.end(*niter);
	  
	  const count_type count_customer = node.table.size_customer();
	  const count_type count_table    = node.table.size_table();
	  
	  os_total.write((char*) &count_customer, sizeof(count_type));
	  os_total.write((char*) &count_table, sizeof(count_type));

	  words.clear();
	  
	  // we have an entry in the model
	  node_type::table_type::const_iterator titer_end = node.table.end();
	  for (node_type::table_type::const_iterator titer = node.table.begin(); titer != titer_end; ++ titer)
	    words.insert(std::make_pair(titer->first, data_type(trie_type::npos(),
								titer->second.size_customer(),
								titer->second.size_table())));
	  
	  // or, we have an entry with the longer contexts..
	  for (trie_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter)
	    if (! trie[iter->second].table.empty()) {
	      std::pair<word_set_type::iterator, bool> result = words.insert(std::make_pair(iter->first, data_type(iter->second, 0, 0)));
	      if (! result.second)
		boost::fusion::get<0>(result.first->second) = iter->second;
	    }
	  
	  word_set_type::const_iterator witer_end = words.end();
	  for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) {
	    nodes_next.push_back(boost::fusion::get<0>(witer->second));
	    
	    word_type::id_type word_id = witer->first.id();
	    os_index.write((char*) &word_id, sizeof(word_type::id_type));
	    os_count.write((char*) &boost::fusion::get<1>(witer->second), sizeof(count_type));
	    os_count.write((char*) &boost::fusion::get<2>(witer->second), sizeof(count_type));
	    
	    positions.set(positions.size(), true);
	  }
	  positions.set(positions.size(), false);
	  
	  offset += words.size();
	}
      }
      
      offsets.push_back(offset);
      
      nodes.swap(nodes_next);
    }
    
    // dump position
    positions.write(rep.path("position"));
    
    // dump offsets
    for (size_type order = 1; order != offsets.size(); ++ order)
      rep[boost::lexical_cast<std::string>(order) + "-gram-offset"] = boost::lexical_cast<std::string>(offsets[order]);
    
    // vocabulary...
    word_type::write(rep.path("vocab"));
  }

public: 
  PYPWord& base;

  trie_type trie;
  node_type root;
  node_map_type nodes;
  
  parameter_set_type discount;
  parameter_set_type strength;
  
  double discount_alpha;
  double discount_beta;
  double strength_shape;
  double strength_rate;
  
  size_type counts0;
};

struct PYPGraph
{
  typedef PYP::size_type       size_type;
  typedef PYP::difference_type difference_type;

  typedef PYP::sentence_type sentence_type;
  typedef PYP::piece_type    piece_type;
  typedef PYP::word_type     word_type;


  typedef uint32_t id_type;
  
  typedef std::vectotr<word_type, std::allocator<word_type> > derivation_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;

  struct segment_type
  {
    word_type       word;
    difference_type first;
    difference_type last;
    
    segment_type() : word(), first(0), last(0) {}
    segment_type(const word_type& __word, const difference_type& __first, const difference_type& __last)
      : word(__word), first(__first), last(__last) {}

    friend
    bool operator==(const segment_type& x, const segment_type& y)
    {
      return x.first == y.first && x.last == y.last;
    }
    
    friend
    size_t hash_value(segment_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(first, last);
    }
  };
  
  typedef utils::simple_vector<segment_type, std::allocator<segment_type> > segment_set_type;
  
  struct segment_set_hash_type : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    size_t operator()(const segment_set_type& x) const
    {
      size_t seed = 0;
      segment_set_type::const_iterator iter_end = x.end();
      for (segment_set_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	seed = hasher_type::operator()(first, hasher_type::operator()(last, seed));
      return seed;
    }
  };
  
  typedef std::equal_to<segment_set_type> segment_set_equal_type;
  
  typedef utils::unordered_map<segment_set_type, id_type, segment_set_hash_type, segment_set_equal_type,
			       std::allocator<std::pair<const segment_set_type, id_type> > >::type segment_set_unique_type;
  
  typedef utils::chunk_vector<segment_set_type, 4096/sizeof(segment_set_type), std::allocator<segment_set_type> > segment_set_map_type;
  
  typedef std::vector<logprob_type, std::allocator<logprob_type> > alpha_type;
  typedef std::vector<size_type, std::allocator<size_type> > position_set_type;
  
  typedef std::vector<logprob_type, std::allocator<logprob_type> > logprob_set_type;
  typedef std::vector<prob_type, std::allocator<prob_type> >       prob_set_type;
  
  struct node_type
  {
    typedef std::vector<id_type, std::allocator<id_type> > edge_set_type;
    
    edge_set_type edges;
    id_type       id;
    
    node_type() : edges(), id() {}
  };
  
  struct edge_type
  { 
    segment_type  segment;
    logprob_type  prob;
    
    id_type       id;
    id_type       head;
    id_type       tail;
    
    edge_type() : segment(), prob(), id(), head(), tail()  {}
    edge_type(const segment_type& __segment,
	      const logprob_type& __prob) 
      : segment(__segment), prob(__prob), id(), head(), tail() {}
  };
  
  typedef utils::chunk_vector<node_type, 4096 /sizeof(node_type), std::allocator<node_type> > node_set_type;
  typedef utils::chunk_vector<edge_type, 4096 /sizeof(edge_type), std::allocator<edge_type> > edge_set_type;

  node_type& add_node()
  {
    const id_type node_id = nodes.size();
    
    nodes.push_back(node_type());
    nodes.back().id = node_id;
    
    return nodes.back();
  }

  edge_type& add_edge()
  {
    const id_type edge_id = edges.size();
    
    edges.push_back(edge_type());
    edges.back().id = edge_id;
    
    return edges.back();
  }
  
  void connect_edge(const id_type edge, const id_type head)
  {
    edges[edge].head = head;
    nodes[head].edges.push_back(edge);
  };
  
  void initialize(const sentence_type& sentence)
  {
    positions.clear();
    positions.push_back(0);
    
    sentence_type::const_iterator siter_begin = sentence.begin();
    sentence_type::const_iterator siter_end   = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; /**/) {
      siter += utils::utf8_size(*siter);
      
      positions.push_back(siter - siter_begin);
    }
    
    nodes.clear();
    edges.clear();
    
    segments.clear();
    alpha.clear();
    states.clear();
    
    node_type& start = nodes.add_node();
    
    segments.push_back(segment_set_type(1, segment_type(vocab_type::BOS, -1, 0)));
    alpha.push_back(cicada::semiring::traits<logprob_type>::one());
  }

  logprob_type forward(const sentence_type& sentence, PYPLM& lm)
  {
    initialize(sentence);
    
    const size_type size = positions.size();
    
    
    
    
    
  }
  
  node_set_type nodes;
  edge_set_type edges;
  
  segment_set_map_type segments;
  position_set_type    positions;
  alpha_type           alpha;
  
  segment_set_unique_type states;
};


typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef utils::sampler<boost::mt19937> sampler_type;

path_set_type train_files;
path_set_type test_files;
path_type     output_file;

int order = 3;
int spell_order = 4;
int spell_length = 10;
int samples = 30;
int baby_steps = 0;
int anneal_steps = 0;
int resample_rate = 1;
int resample_iterations = 2;
bool slice_sampling = false;

double discount = 0.9;
double strength = 1;

double discount_prior_alpha = 1.0;
double discount_prior_beta  = 1.0;
double strength_prior_shape = 1.0;
double strength_prior_rate  = 1.0;

double spell_discount = 0.9;
double spell_strength = 1.0;

double spell_discount_prior_alpha = 1.0;
double spell_discount_prior_beta  = 1.0;
double spell_strength_prior_shape = 1.0;
double spell_strength_prior_rate  = 1.0;

double lambda = 4;
double lambda_alpha = 0.2;
double lambda_beta = 0.1;

int threads = 1;
int debug = 0;

size_t vocabulary_size(const path_set_type& files);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    if (order <= 0)
      throw std::runtime_error("order must be positive");
    
    if (samples < 0)
      throw std::runtime_error("# of samples must be positive");
            
    if (resample_rate <= 0)
      throw std::runtime_error("resample rate must be >= 1");

    if (train_files.empty())
      throw std::runtime_error("no training data?");
    
    if (! slice_sampling && strength < 0.0)
      throw std::runtime_error("negative strength w/o slice sampling is not supported!");

    sampler_type sampler;
    const size_t num_vocab = vocabulary_size(train_files);
    
    PYPLM lm(order,
	     1.0 / num_vocab,
	     discount,
	     strength,
	     discount_prior_alpha,
	     discount_prior_beta,
	     strength_prior_shape,
	     strength_prior_rate,
	     infinite);
    
    // we will precompute <word, node> pair...
    typedef boost::fusion::tuple<word_type, PYPLM::id_type, int> data_type;
    typedef std::vector<data_type, std::allocator<data_type> > data_set_type;
    typedef std::vector<bool, std::allocator<bool> > non_oov_type;

    typedef utils::unordered_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_unique_type;
    typedef utils::chunk_vector<word_unique_type, 4096 / sizeof(word_unique_type), std::allocator<word_unique_type> > word_unique_set_type;
    typedef std::vector<size_t, std::allocator<size_t> > index_type;

    data_set_type        training;
    word_unique_set_type uniques;
    non_oov_type non_oov;
    
    if (vocab_type::BOS.id() >= non_oov.size())
      non_oov.resize(vocab_type::BOS.id() + 1, false);
    if (vocab_type::EOS.id() >= non_oov.size())
      non_oov.resize(vocab_type::EOS.id() + 1, false);
    
    non_oov[vocab_type::BOS.id()] = true;
    non_oov[vocab_type::EOS.id()] = true;
    
    for (path_set_type::const_iterator fiter = train_files.begin(); fiter != train_files.end(); ++ fiter) {
      utils::compress_istream is(*fiter, 1024 * 1024);
      
      sentence_type sentence;
      sentence_type ngram(1, vocab_type::BOS);
      
      while (is >> sentence) {
	ngram.resize(1);
	
	sentence_type::const_iterator siter_end = sentence.end();
	for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  training.push_back(data_type(*siter, lm.insert(std::max(ngram.begin(), ngram.end() - order + 1), ngram.end()), 0));
	  
	  if (boost::fusion::get<1>(training.back()) >= uniques.size())
	    uniques.resize(boost::fusion::get<1>(training.back()) + 1);
	  uniques[boost::fusion::get<1>(training.back())].insert(boost::fusion::get<0>(training.back()));
	  
	  ngram.push_back(*siter);
	  
	  if (siter->id() >= non_oov.size())
	    non_oov.resize(siter->id() + 1, false);
	  non_oov[siter->id()] = true;
	}
	
	training.push_back(data_type(vocab_type::EOS, lm.insert(std::max(ngram.begin(), ngram.end() - order + 1), ngram.end()), 0));
	
	if (boost::fusion::get<1>(training.back()) >= uniques.size())
	  uniques.resize(boost::fusion::get<1>(training.back()) + 1);
	uniques[boost::fusion::get<1>(training.back())].insert(boost::fusion::get<0>(training.back()));
      }
    }

    if (training.empty())
      throw std::runtime_error("no training data?");
    
    data_set_type(training).swap(training);

    
    // assign rank...
    {
      data_set_type::iterator titer_end = training.end();
      for (data_set_type::iterator titer = training.begin(); titer != titer_end; ++ titer)
	boost::fusion::get<2>(*titer) = uniques[boost::fusion::get<1>(*titer)].size();
      
      uniques.clear();
    }
    
    // sort by rank
    std::sort(training.begin(), training.end(), less_data<data_type>());
    
    // compute index
    index_type index(1, 0);
    
    data_set_type::const_iterator titer_begin = training.begin();
    data_set_type::const_iterator titer_end   = training.end();
    for (data_set_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
      if (boost::fusion::get<2>(*titer) != boost::fusion::get<2>(training[index.back()]))
	index.push_back(titer - titer_begin);
    
    index.push_back(training.size());
    
    if (debug >= 2)
      std::cerr << "# of baby step levels: " << (index.size() - 1) << std::endl;
    
    // sample parameters, first...
    if (slice_sampling)
      lm.slice_sample_parameters(sampler, resample_iterations);
    else
      lm.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2)
      for (int n = 0; n != order; ++ n)
	std::cerr << "order=" << n << " discount=" << lm.discount[n] << " strength=" << lm.strength[n] << std::endl;
    
    {
      data_set_type::iterator titer_end = training.end();
      for (data_set_type::iterator titer = training.begin(); titer != titer_end; ++ titer)
	boost::fusion::get<2>(*titer) = false;
    }
    
    size_t baby_index = 0;
    size_t baby_iter = utils::bithack::branch(baby_steps > 0, size_t(0), index.back());
    const size_t baby_last = index.back();
    const size_t baby_size = (baby_steps > 0 ? (index.back() + (baby_steps - 1)) / baby_steps : size_t(0));
    
    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);
      
    data_set_type training_samples;
    
    if (baby_iter == baby_last)
      training_samples = training;
    else
      training_samples.reserve(training.size());
      
    bool sampling = false;
    int sample_iter = 0;
    
    // then, learn!
    for (size_t iter = 0; sample_iter != samples; ++ iter, sample_iter += sampling) {

      bool baby_finished = true;
      if (baby_iter != baby_last) {
	baby_finished = false;
	
	const size_t baby_next = utils::bithack::min(baby_iter + baby_size, baby_last);
	
	while (baby_iter < baby_next) {
	  std::copy(training.begin() + index[baby_index], training.begin() + index[baby_index + 1], std::back_inserter(training_samples));
	  
	  baby_iter = index[baby_index + 1];
	  ++ baby_index;
	}
	
	if (debug >= 2)
	  std::cerr << "baby: " << training_samples.size() << std::endl;
      }

      double temperature = 1.0;
      bool anneal_finished = true;
      if (anneal_iter != anneal_last) {
	anneal_finished = false;
	temperature = double(anneal_last - anneal_iter) + 1;
	
	++ anneal_iter;

	if (debug >= 2)
	  std::cerr << "temperature: " << temperature << std::endl;
      }
      
      sampling = baby_finished && anneal_finished;
      
      if (debug) {
	if (sampling)
	  std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
	else
	  std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
      }
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(training_samples.begin(), training_samples.end(), gen);

      if (infinite)  {
	data_set_type::iterator titer_end = training_samples.end();
	for (data_set_type::iterator titer = training_samples.begin(); titer != titer_end; ++ titer) {
	  if (boost::fusion::get<2>(*titer))
	    lm.decrement(boost::fusion::get<0>(*titer), boost::fusion::get<1>(*titer), boost::fusion::get<2>(*titer), sampler);
	  
	  lm.increment(boost::fusion::get<0>(*titer), boost::fusion::get<1>(*titer), boost::fusion::get<2>(*titer), sampler, temperature);
	}
	
	if (debug >= 2) {
	  std::cerr << "penetration count" << std::endl;
	  size_t total = 0;
	  for (size_t n = 0; n != lm.orders.size(); ++ n) {
	    std::cerr << "order=" << n << " a=" << lm.orders[n].first << " b=" << lm.orders[n].second << std::endl;
	    total += lm.orders[n].first;
	  }
	  std::cerr << "total=" << total << std::endl;
	}

      } else {
	data_set_type::iterator titer_end = training_samples.end();
	for (data_set_type::iterator titer = training_samples.begin(); titer != titer_end; ++ titer) {
	  if (boost::fusion::get<2>(*titer))
	    lm.decrement(boost::fusion::get<0>(*titer), boost::fusion::get<1>(*titer), sampler);
	  else
	    boost::fusion::get<2>(*titer) = true;
	  
	  lm.increment(boost::fusion::get<0>(*titer), boost::fusion::get<1>(*titer), sampler, temperature);
	}
      }
      
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  lm.slice_sample_parameters(sampler, resample_iterations);
	else
	  lm.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2)
	  for (int n = 0; n != order; ++ n)
	    std::cerr << "order=" << n << " discount=" << lm.discount[n] << " strength=" << lm.strength[n] << std::endl;
      }
	
      if (debug)
	std::cerr << "log-likelihood: " << lm.log_likelihood() << std::endl;
    }
    
    // clear training data
    training.clear();
    data_set_type(training).swap(training);
    
    // we will dump LM... now, define a format!
    if (! output_file.empty())
      lm.write(output_file);
    
    // testing!
    if (! test_files.empty()) {
      sentence_type sentence;
      sentence_type ngram(1, vocab_type::BOS);

      double logprob_total = 0.0;
      size_t num_word = 0;
      size_t num_oov = 0;
      size_t num_sentence = 0;
      
      for (path_set_type::const_iterator fiter = test_files.begin(); fiter != test_files.end(); ++ fiter) {
	utils::compress_istream is(*fiter, 1024 * 1024);
	
	while (is >> sentence) {
	  ngram.resize(1);
	  
	  sentence_type::const_iterator siter_end = sentence.end();
	  for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	    const bool is_oov = ! (siter->id() < non_oov.size() && non_oov[siter->id()]);
	    const double prob = lm.prob(*siter, std::max(ngram.begin(), ngram.end() - order + 1), ngram.end());
	    
	    if (! is_oov)
	      logprob_total += std::log(prob);
	    
	    num_oov += is_oov;
	    
	    ngram.push_back(*siter);
	  }

	  const double prob = lm.prob(vocab_type::EOS, std::max(ngram.begin(), ngram.end() - order + 1), ngram.end());
	  logprob_total += std::log(prob);
	  
	  num_word += sentence.size();
	  ++ num_sentence;
	}
      }
      
      std::cerr << "# of sentences: " << num_sentence
		<< " # of words: " << num_word
		<< " # of OOV: " << num_oov
		<< " order: " << order
		<< std::endl;
      
      std::cerr << "logprob = " << logprob_total << std::endl;
      std::cerr << "ppl     = " << utils::mathop::exp(- logprob_total / (num_word - num_oov + num_sentence)) << std::endl;
      std::cerr << "ppl1    = " << utils::mathop::exp(- logprob_total / (num_word - num_oov)) << std::endl;
    }
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct string_hash : public utils::hashmurmur<size_t>
{
  typedef utils::hashmurmur<size_t> hasher_type;
  
  size_t operator()(const std::string& x) const
  {
    return hasher_type::operator()(x.begin(), x.end(), 0);
  }
};

template <typename Tp>
struct greater_second
{
  bool operator()(const Tp& x, const Tp& y) const
  {
    return x.second > y.second;
  }
  
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

size_t vocabulary_size(const path_set_type& files)
{
  typedef uint64_t count_type;
  typedef utils::unordered_map<std::string, count_type, string_hash, std::equal_to<std::string>, std::allocator<std::pair<const std::string, count_type> > >::type vocab_type;

  typedef std::vector<const vocab_type::value_type*, std::allocator<const vocab_type::value_type*> > sorted_type;
  
  vocab_type vocab;
  std::string word;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    if (! boost::filesystem::exists(*fiter))
      throw std::runtime_error("no file? " + fiter->string());
    
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (is >> word) 
      ++ vocab[word];
  }

  sorted_type sorted;
  sorted.reserve(vocab.size());
  
  vocab_type::const_iterator viter_end = vocab.end();
  for (vocab_type::const_iterator viter = vocab.begin(); viter != viter_end; ++ viter)
    sorted.push_back(&(*viter));
  
  std::sort(sorted.begin(), sorted.end(), greater_second<vocab_type::value_type>());

  sorted_type::const_iterator siter_end = sorted.end();
  for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter) {
    word_type __tmptmp((*siter)->first);
  }
  
  return vocab.size();
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("train", po::value<path_set_type>(&train_files)->multitoken(), "train file(s)")
    ("test",  po::value<path_set_type>(&test_files)->multitoken(),  "test file(s)")
    ("output", po::value<path_type>(&output_file), "output file")
    
    ("order", po::value<int>(&order)->default_value(order), "max ngram order")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    ("infinite",            po::bool_switch(&infinite),                                               "infinite n-gram language model")
    
    ("discount",       po::value<double>(&discount)->default_value(discount),                         "discount ~ Beta(alpha,beta)")
    ("discount-alpha", po::value<double>(&discount_prior_alpha)->default_value(discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("discount-beta",  po::value<double>(&discount_prior_beta)->default_value(discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("strength",       po::value<double>(&strength)->default_value(strength),                         "strength ~ Gamma(shape,rate)")
    ("strength-shape", po::value<double>(&strength_prior_shape)->default_value(strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("strength-rate",  po::value<double>(&strength_prior_rate)->default_value(strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}

