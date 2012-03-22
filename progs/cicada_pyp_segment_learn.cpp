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
#include "utils/array_power2.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/chinese_restaurant_process.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/dense_hash_set.hpp"
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

#include <unicode/uchar.h>
#include <unicode/uscript.h>

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
    struct poisson_type
    {
      typedef utils::array_power2<double, 32, std::allocator<double> > cache_type;

      poisson_type() { initialize(1.0); }
      poisson_type(const double& __lambda) { initialize(__lambda); }
      
      double logprob(const size_type size) const
      {
	if (size < cache.size())
	  return cache[size];
	else
	  return utils::mathop::log_poisson(size, lambda);
      }
      
      void initialize(const double& __lambda)
      {
	cache.clear();
	lambda = __lambda;
	
	cache[0] = - std::numeric_limits<double>::infinity();
	for (size_type size = 1; size != cache.size(); ++ size)
	  cache[size] = utils::mathop::log_poisson(size, lambda);
      }
      
      cache_type cache;
      double lambda;
    };
    
    typedef utils::unordered_map<id_type, poisson_type, utils::hashmurmur<size_t>, std::equal_to<id_type>,
				 std::allocator<std::pair<const id_type, poisson_type> > >::type poisson_set_type;
    
    template <typename Iterator>
    length_base_type(Iterator first, Iterator last,
		     const double& __lambda,
		     const double& __strength_shape,
		     const double& __strength_rate)
      : strength_shape(__strength_shape),
	strength_rate(__strength_rate)
    {
      for (/**/; first != last; ++ first)
	poissons[code_class(*first)] = poisson_type(__lambda);
      
      poisson = poisson_type(__lambda);
    }
    
    double logprob(const segment_type& segment, const size_type size) const
    {
      const id_type id = code_class(segment);
      
      poisson_set_type::const_iterator piter = poissons.find(id);
      if (piter != poissons.end())
	return piter->second.logprob(size);
      else
	return poisson.logprob(size);
    }
    
    template <typename Iterator, typename Sampler>
    void sample_parameters(Iterator first, Iterator last, Sampler& sampler)
    {
      typedef std::pair<double, double> count_type;
      typedef utils::unordered_map<id_type, count_type, utils::hashmurmur<size_t>, std::equal_to<id_type>,
				   std::allocator<std::pair<const id_type, count_type> > >::type count_set_type;

    
      count_type     total;
      count_set_type counts;
      for (/**/; first != last; ++ first) {
	count_type& count = counts[code_class(first->first)];

	count.first  += first->first.size() * first->second;
	count.second += first->second;
	total.first  += first->first.size() * first->second;
	total.second += first->second;
      }
      
      typename count_set_type::const_iterator citer_end = counts.end();
      for (typename count_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
	poissons[citer->first].initialize(sampler.gamma(citer->second.first + strength_shape,
							citer->second.second + strength_rate));
      
      poisson.initialize(sampler.gamma(total.first + strength_shape,
				       total.second + strength_rate));
    }
    
    id_type code_class(const segment_type& segment) const
    {
      const UChar32 code = utils::utf8_code(segment.begin());
      
      UErrorCode status = U_ZERO_ERROR;
      const id_type script = uscript_getScript(code, &status);
      const id_type category = u_getIntPropertyValue(code, UCHAR_GENERAL_CATEGORY);
      
      return ((category & 0xffff) << 16) | (script & 0xffff);
    }
    
    friend
    std::ostream& operator<<(std::ostream& os, const length_base_type& base)
    {
      poisson_set_type::const_iterator piter_end = base.poissons.end();
      for (poisson_set_type::const_iterator piter = base.poissons.begin(); piter != piter_end; ++ piter)
        os << "script: " << ((piter->first >> 16) & 0xffff) << " cat: " << (piter->first & 0xffff) << " lambda: " << piter->second.lambda << std::endl;
      
      os << "lambda: " << base.poisson.lambda << std::endl;
      
      return os;
    }
    

    poisson_set_type poissons;
    poisson_type     poisson;
    
    double strength_shape;
    double strength_rate;
  };
  
  template <typename Iterator>
  PYPWord(Iterator first, Iterator last,
	  const int order,
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
      base(first, last, __lambda, __lambda_strength, __lambda_rate)
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
    
    return increment(segment, node, sampler, temperature);
  }

  template <typename Sampler>
  bool increment(const segment_type& segment, const id_type& node, Sampler& sampler, const double temperature=1.0)
  {
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
    
    return decrement(segment, node, sampler);
  }
  
  template <typename Sampler>
  bool decrement(const segment_type& segment, const id_type& node, Sampler& sampler)
  {
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
    
    double p = root.table.prob(segment, p0);
      
    // we will traverse from the back!
    reverse_iterator begin(last);
    reverse_iterator end(first);
      
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter) {
      node = trie.find(node, *iter);
	
      if (node == trie_type::npos() || trie[node].table.empty())
	return p;
      else
	p = trie[node].table.prob(segment, p);
    }
      
    return p;
  }
  
  template <typename Sampler>
  void increment(const segment_type& segment, Sampler& sampler, const double temperature=1.0)
  {
    const size_type context_size = discount.size() - 1;
    
    buffer.clear();
    buffer.push_back(BOS());
    
    segment_type::const_iterator siter_end = segment.end();
    for (segment_type::const_iterator siter = segment.begin(); siter != siter_end; /**/) {
      const size_type char_size = utils::utf8_size(*siter);
      const segment_type seg(siter, siter + char_size);
      
      increment(seg, std::max(buffer.begin(), buffer.end() - context_size), buffer.end(), sampler, temperature);
      
      buffer.push_back(seg);
      siter += char_size;
    }
    
    increment(EOS(), std::max(buffer.begin(), buffer.end() - context_size), buffer.end(), sampler, temperature);
  }
  
  template <typename Sampler>
  void decrement(const segment_type& segment, Sampler& sampler)
  {
    const size_type context_size = discount.size() - 1;

    buffer.clear();
    buffer.push_back(BOS());
    
    segment_type::const_iterator siter_end = segment.end();
    for (segment_type::const_iterator siter = segment.begin(); siter != siter_end; /**/) {
      const size_type char_size = utils::utf8_size(*siter);
      const segment_type seg(siter, siter + char_size);
      
      decrement(seg, std::max(buffer.begin(), buffer.end() - context_size), buffer.end(), sampler);
      
      buffer.push_back(seg);
      siter += char_size;
    }
    
    decrement(EOS(), std::max(buffer.begin(), buffer.end() - context_size), buffer.end(), sampler);
  }
  
  double prob(const segment_type& segment)
  {
    const size_type context_size = discount.size() - 1;
    
    buffer.clear();
    buffer.push_back(BOS());
    
    double logprob = 0.0;
    
    segment_type::const_iterator siter_end = segment.end();
    for (segment_type::const_iterator siter = segment.begin(); siter != siter_end; /**/) {
      const size_type char_size = utils::utf8_size(*siter);
      const segment_type seg(siter, siter + char_size);
      
      logprob += std::log(prob(seg, std::max(buffer.begin(), buffer.end() - context_size), buffer.end()));
      
      buffer.push_back(seg);
      siter += char_size;
    }
    
    logprob += std::log(prob(EOS(), std::max(buffer.begin(), buffer.end() - context_size), buffer.end()));
    
    // exclude BOS...
    logprob += base.logprob(segment, buffer.size() - 1);
    
    return std::exp(logprob);
  }

  // TODO: add log_likelihood
  //     : experiment with sampled length parameters... (Do we really need this...? Is it simply a normalization constant...?)
  //

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
  
  template <typename Iterator, typename Sampler>
  bool increment(const word_type& word, Iterator first, Iterator last, Sampler& sampler, const double temperature=1.0)
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

    // first, we will insert word...
    base.increment(word, sampler, temperature);
    
    // then, ngram...
    return increment(word, node, sampler, temperature);
  }

  template <typename Iterator, typename Sampler>
  bool decrement(const word_type& word, Iterator first, Iterator last, Sampler& sampler)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter)
      node = trie.find(node, *iter);
    
    // first, we will decrement word...
    base.decrement(word, sampler);
    
    // then, ngram...
    return decrement(word, node, sampler);
  }
  
  template <typename Sampler>
  bool increment(const word_type& word, const id_type& node, Sampler& sampler, const double temperature=1.0)
  {
    if (node == trie.root()) {
      if (root.table.increment(word, base.prob(word), sampler, temperature))
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
      return root.table.prob(word, base.prob(word));
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
    const double p0 = base.prob(word);
    
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
    sample_length(sampler);
    
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
    sample_length(sampler);
    
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
  
  template <typename Sampler>
  void sample_length(Sampler& sampler)
  {
    typedef size_type count_type;
    typedef utils::dense_hash_map<segment_type, count_type, boost::hash<segment_type>, std::equal_to<segment_type>,
				  std::allocator<std::pair<const segment_type, count_type> > >::type count_set_type;

    count_set_type counts;
    counts.set_empty_key(segment_type());
    
    for (size_type order = 0; order != nodes.size(); ++ order) {
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
	if (trie.empty(*niter)) {
	  // we will collect information from nodes without children...
	  
	  typename node_type::table_type::const_iterator titer_end = trie[*niter].table.end();
	  for (typename node_type::table_type::const_iterator titer = trie[*niter].table.begin(); titer != titer_end; ++ titer)
	    counts[titer->first] += titer->second.size_table();
	}
    }
    
    base.base.sample_parameters(counts.begin(), counts.end(), sampler);
  }
  

  void write(const path_type& path)
  {
    
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
  typedef PYP::segment_type  segment_type;
  typedef PYP::word_type     word_type;

  typedef uint32_t id_type;
  
  typedef std::vector<word_type, std::allocator<word_type> > derivation_type;
  
  typedef cicada::semiring::Logprob<double> logprob_type;
  typedef double prob_type;
  
  
  typedef utils::simple_vector<segment_type, std::allocator<segment_type> > segment_set_type;
  
  struct segment_set_hash_type : public utils::hashmurmur<size_t>
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    size_t operator()(const segment_set_type& x) const
    {
      size_t seed = 0;
      segment_set_type::const_iterator iter_end = x.end();
      for (segment_set_type::const_iterator iter = x.begin(); iter != iter_end; ++ iter)
	seed = hasher_type::operator()(iter->begin(), iter->end(), seed);
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

  typedef std::vector<id_type, std::allocator<id_type> > frontier_type;
  typedef std::vector<frontier_type, std::allocator<frontier_type> > frontier_set_type;
  typedef std::vector<bool, std::allocator<bool> > boundary_set_type;

  node_type& add_node()
  {
    const id_type node_id = nodes.size();
    
    nodes.push_back(node_type());
    nodes.back().id = node_id;
    
    return nodes.back();
  }

  edge_type& add_edge(const id_type head, const id_type tail)
  {
    const id_type edge_id = edges.size();
    
    edges.push_back(edge_type());
    
    edge_type& edge = edges.back();
    
    edge.id = edge_id;
    edge.head = head;
    edge.tail = tail;
    
    nodes[head].edges.push_back(edge_id);
    
    return edge;
  }
  
  void connect_edge(const id_type edge, const id_type head)
  {
    edges[edge].head = head;
    nodes[head].edges.push_back(edge);
  };
  
  void initialize(const sentence_type& sentence)
  {
    boundaries.clear();
    positions.clear();
    positions.push_back(0);
    
    sentence_type::const_iterator siter_begin = sentence.begin();
    sentence_type::const_iterator siter_end   = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; /**/) {
      boundaries.push_back(u_isUWhiteSpace(utils::utf8_code(siter)));
      
      siter += utils::utf8_size(*siter);
      
      positions.push_back(siter - siter_begin);
    }
    
    nodes.clear();
    edges.clear();
    
    frontiers.clear();
    segments.clear();
    alpha.clear();
    states.clear();
    
    frontiers.resize(positions.size());
    
    node_type& start = add_node();
    segments.push_back(segment_set_type());
    alpha.push_back(cicada::semiring::traits<logprob_type>::one());
    
    node_type& node = add_node();
    frontiers.front().push_back(node.id);
    segments.push_back(segment_set_type(1, static_cast<const std::string&>(vocab_type::BOS)));
    alpha.push_back(cicada::semiring::traits<logprob_type>::one());
    
    edge_type& edge = add_edge(node.id, start.id);
    edge.segment = static_cast<const std::string&>(vocab_type::BOS);
    edge.prob = cicada::semiring::traits<logprob_type>::one();
  }

  logprob_type forward(const sentence_type& sentence, PYPLM& lm, const difference_type length_max, const bool ignore_boundary)
  {
    initialize(sentence);
    
    const difference_type context_size = lm.discount.size() - 1;
    const difference_type size = positions.size();
    
    for (difference_type last = 1; last != size; ++ last) {
      
      if (! ignore_boundary && boundaries[last - 1]) {
	frontiers[last].swap(frontiers[last - 1]);
	
	continue;
      }

      states.clear();
      
      frontier_type& curr = frontiers[last];
      
      for (difference_type first = last - 1; first >= last - utils::bithack::min(last, length_max); -- first) {
	const frontier_type& prevs = frontiers[first];

	if (prevs.empty()) break;
	
	const segment_type seg(sentence.begin() + positions[first], sentence.begin() + positions[last]);
	
	frontier_type::const_iterator piter_end = prevs.end();
	for (frontier_type::const_iterator piter = prevs.begin(); piter != piter_end; ++ piter) {
	  const node_type& prev = nodes[*piter];
	  
	  const segment_set_type& segs_prev = segments[*piter];
	  segment_set_type segs(segs_prev);
	  segs.push_back(seg);
	  
	  const segment_set_type context(std::max(segs.begin(), segs.end() - context_size), segs.end());
	  
	  std::pair<segment_set_unique_type::iterator, bool> result = states.insert(std::make_pair(context, 0));
	  if (result.second) {
	    result.first->second = add_node().id;
	    
	    segments.push_back(context);
	    alpha.push_back(cicada::semiring::traits<logprob_type>::zero());
	    curr.push_back(result.first->second);
	  }
	  
	  edge_type& edge = add_edge(result.first->second, prev.id);
	  
	  edge.segment = seg;
	  edge.prob = lm.prob(seg, segs_prev.begin(), segs_prev.end());
	  
	  alpha[result.first->second] += edge.prob * alpha[prev.id];
	}
      }
      
      states.clear();
    }
    
    // final EOS...
    const segment_type seg(static_cast<const std::string&>(vocab_type::EOS));
    const node_type& goal = add_node();
    alpha.push_back(cicada::semiring::traits<logprob_type>::zero());
    
    const frontier_type& prevs = frontiers[size - 1];
    frontier_type::const_iterator piter_end = prevs.end();
    for (frontier_type::const_iterator piter = prevs.begin(); piter != piter_end; ++ piter) {
      const node_type& prev = nodes[*piter];
      const segment_set_type& segs_prev = segments[*piter];
      
      edge_type& edge = add_edge(goal.id, prev.id);
      
      edge.segment = seg;
      edge.prob = lm.prob(seg, segs_prev.begin(), segs_prev.end());
      
      alpha[goal.id] += edge.prob * alpha[prev.id];
    }
    
    return alpha[goal.id];
  }
  
  // backward sampling
  template <typename Sampler>
  logprob_type backward(Sampler& sampler, derivation_type& derivation)
  {
    // we will traverse from EOS, the last node...
    
    derivation.clear();
    logprob_type prob_derivation = cicada::semiring::traits<logprob_type>::one();
    
    id_type node_id = nodes.size() - 1;
    
    while (! nodes[node_id].edges.empty()) {
      const node_type::edge_set_type& prevs = nodes[node_id].edges;
      
      logprob_type logsum;
      logprobs.clear();
      node_type::edge_set_type::const_iterator piter_end = prevs.end();
      for (node_type::edge_set_type::const_iterator piter = prevs.begin(); piter != piter_end; ++ piter) {
	const edge_type& edge = edges[*piter];
	
	logprobs.push_back(alpha[edge.tail] * edge.prob);
	logsum += logprobs.back();
      }
      
      probs.clear();
      logprob_set_type::const_iterator liter_end = logprobs.end();
      for (logprob_set_type::const_iterator liter = logprobs.begin(); liter != liter_end; ++ liter)
	probs.push_back(*liter / logsum);
      
      const size_type pos = sampler.select(probs.begin(), probs.end()) - probs.begin();
      const edge_type& edge = edges[prevs[pos]];
      
      derivation.push_back(edge.segment);
      prob_derivation *= edge.prob;
      
      node_id = edge.tail;
    }
    
    std::reverse(derivation.begin(), derivation.end());

    
    return prob_derivation;
  }
  
  node_set_type nodes;
  edge_set_type edges;
  
  frontier_set_type    frontiers;
  segment_set_map_type segments;
  position_set_type    positions;
  boundary_set_type    boundaries;
  alpha_type           alpha;
  
  segment_set_unique_type states;

  logprob_set_type logprobs;
  prob_set_type    probs;
};


typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef utils::sampler<boost::mt19937> sampler_type;

typedef PYP::size_type       size_type;
typedef PYP::difference_type difference_type;

typedef PYP::sentence_type sentence_type;
typedef PYP::piece_type    piece_type;
typedef PYP::segment_type  segment_type;
typedef PYP::word_type     word_type;

typedef PYPGraph::derivation_type derivation_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > data_set_type;
typedef std::vector<derivation_type, std::allocator<derivation_type> > derivation_set_type;
typedef std::vector<size_type, std::allocator<size_type> > position_set_type;
typedef std::vector<segment_type, std::allocator<segment_type> > vocabulary_type;

path_set_type train_files;
path_set_type test_files;
path_type     output_file;

int order = 3;
int spell_order = 4;
int spell_length = 20;
bool ignore_boundary = false;
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

double spell_lambda = 8;
double spell_lambda_shape = 0.2;
double spell_lambda_rate = 0.1;

int threads = 1;
int debug = 0;

void options(int argc, char** argv);
void read_data(const path_set_type& paths, data_set_type& data);
void vocabulary(const data_set_type& data, vocabulary_type& vocab);

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

    
    data_set_type training;
    read_data(train_files, training);
    
    if (training.empty())
      throw std::runtime_error("no training data?");
    
    vocabulary_type vocab;
    vocabulary(training, vocab);
    
    if (debug)
      std::cerr << "char set size: " << vocab.size() << std::endl;
    
    derivation_set_type derivations(training.size());
    position_set_type positions(training.size());
    
    for (size_t i = 0; i != training.size(); ++ i)
      positions[i] = i;
    
    sampler_type sampler;

    PYPWord spell(vocab.begin(), vocab.end(),
		  spell_order,
		  1.0 / vocab.size(),
		  spell_discount,
		  spell_strength,
		  spell_discount_prior_alpha,
		  spell_discount_prior_beta,
		  spell_strength_prior_shape,
		  spell_strength_prior_rate,
		  spell_lambda,
		  spell_lambda_shape,
		  spell_lambda_rate);
    
    
    PYPLM lm(spell,
	     order,
	     discount,
	     strength,
	     discount_prior_alpha,
	     discount_prior_beta,
	     strength_prior_shape,
	     strength_prior_rate);
    
    PYPGraph graph;
    
    if (slice_sampling)
      lm.slice_sample_parameters(sampler, resample_iterations);
    else
      lm.sample_parameters(sampler, resample_iterations);
    
    if (debug >= 2) {
      for (size_t n = 0; n != lm.discount.size(); ++ n)
	std::cerr << "word order=" << n << " discount=" << lm.discount[n] << " strength=" << lm.strength[n] << std::endl;
      
      for (size_t n = 0; n != lm.base.discount.size(); ++ n)
	std::cerr << "spell order=" << n << " discount=" << lm.base.discount[n] << " strength=" << lm.base.strength[n] << std::endl;
      
      std::cerr << lm.base.base;
    }

    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);
    
    bool sampling = false;
    int sample_iter = 0;
    
    // then, learn!
    for (size_t iter = 0; sample_iter != samples; ++ iter, sample_iter += sampling) {
      
      double temperature = 1.0;
      bool anneal_finished = true;
      if (anneal_iter != anneal_last) {
	anneal_finished = false;
	temperature = double(anneal_last - anneal_iter) + 1;
	
	++ anneal_iter;
	
	if (debug >= 2)
	  std::cerr << "temperature: " << temperature << std::endl;
      }
      
      sampling = anneal_finished;
      
      if (debug) {
	if (sampling)
	  std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
	else
	  std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
      }
      
      
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(positions.begin(), positions.end(), gen);
      
      for (size_t i = 0; i != positions.size(); ++ i) {
	const size_t pos = positions[i];
	
	if (training[pos].empty()) continue;

	const size_t context_size = order - 1;

	if (debug >= 3)
	  std::cerr << "training=" << pos << std::endl;
	
	if (! derivations[pos].empty()) {
	  derivation_type::const_iterator diter_begin = derivations[pos].begin();
	  derivation_type::const_iterator diter_end   = derivations[pos].end();
	  
	  for (derivation_type::const_iterator diter = diter_begin + 1; diter != diter_end; ++ diter)
	    lm.decrement(*diter, std::max(diter_begin, diter - context_size), diter, sampler);
	}
	
	const PYPGraph::logprob_type logsum = graph.forward(training[pos], lm, spell_length, ignore_boundary);
	
	const PYPGraph::logprob_type logderivation = graph.backward(sampler, derivations[pos]);
	
	if (debug >= 3) {
	  std::cerr << "sum=" << logsum << " derivation=" << logderivation << std::endl;
	  
	  derivation_type::const_iterator diter_end = derivations[pos].end();
	  for (derivation_type::const_iterator diter = derivations[pos].begin(); diter != diter_end; ++ diter)
	    std::cerr << "word=\"" << *diter << "\"" << std::endl;
	}
	
	derivation_type::const_iterator diter_begin = derivations[pos].begin();
	derivation_type::const_iterator diter_end   = derivations[pos].end();
	
	for (derivation_type::const_iterator diter = diter_begin + 1; diter != diter_end; ++ diter)
	  lm.increment(*diter, std::max(diter_begin, diter - context_size), diter, sampler);
      }

      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  lm.slice_sample_parameters(sampler, resample_iterations);
	else
	  lm.sample_parameters(sampler, resample_iterations);
	
	if (debug >= 2) {
	  for (size_t n = 0; n != lm.discount.size(); ++ n)
	    std::cerr << "word order=" << n << " discount=" << lm.discount[n] << " strength=" << lm.strength[n] << std::endl;
	  
	  for (size_t n = 0; n != lm.base.discount.size(); ++ n)
	    std::cerr << "spell order=" << n << " discount=" << lm.base.discount[n] << " strength=" << lm.base.strength[n] << std::endl;
	  
	  std::cerr << lm.base.base;
	}
      }
	
      if (debug)
	std::cerr << "log-likelihood: " << lm.log_likelihood() << std::endl;
    }
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void read_data(const path_set_type& paths, data_set_type& data)
{
  data.clear();
  
  for (path_set_type::const_iterator fiter = paths.begin(); fiter != paths.end(); ++ fiter) { 
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    std::string line;
    while (std::getline(is, line))
      data.push_back(line);
  }
}

void vocabulary(const data_set_type& data, vocabulary_type& vocab)
{
  typedef utils::piece piece_type;
  typedef utils::dense_hash_set<piece_type, boost::hash<piece_type>, std::equal_to<piece_type>, std::allocator<piece_type> >::type vocab_type;
  
  vocab_type voc;
  voc.set_empty_key(piece_type());
  
  data_set_type::const_iterator diter_end = data.end();
  for (data_set_type::const_iterator diter = data.begin(); diter != diter_end; ++ diter) {
    const sentence_type& sentence = *diter;
    
    if (sentence.empty()) continue;
    
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; /**/) {
      const size_t char_width = utils::utf8_size(*siter);
      
      voc.insert(piece_type(siter, siter + char_width));
      
      siter += char_width;
    }
  }

  vocab.clear();
  vocab.insert(vocab.end(), voc.begin(), voc.end());
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
    
    ("order",        po::value<int>(&order)->default_value(order),               "max ngram order")
    ("spell-order",  po::value<int>(&spell_order)->default_value(spell_order),   "max spell ngram order")
    ("spell-length", po::value<int>(&spell_length)->default_value(spell_length), "max spell length")
    
    ("ignore-boundary", po::bool_switch(&ignore_boundary), "ignore word boundaries")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    
    ("discount",       po::value<double>(&discount)->default_value(discount),                         "discount ~ Beta(alpha,beta)")
    ("discount-alpha", po::value<double>(&discount_prior_alpha)->default_value(discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("discount-beta",  po::value<double>(&discount_prior_beta)->default_value(discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("strength",       po::value<double>(&strength)->default_value(strength),                         "strength ~ Gamma(shape,rate)")
    ("strength-shape", po::value<double>(&strength_prior_shape)->default_value(strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("strength-rate",  po::value<double>(&strength_prior_rate)->default_value(strength_prior_rate),   "strength ~ Gamma(shape,rate)")

    ("spell-discount",       po::value<double>(&spell_discount)->default_value(spell_discount),                         "discount ~ Beta(alpha,beta)")
    ("spell-discount-alpha", po::value<double>(&spell_discount_prior_alpha)->default_value(spell_discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("spell-discount-beta",  po::value<double>(&spell_discount_prior_beta)->default_value(spell_discount_prior_beta),   "discount ~ Beta(alpha,beta)")

    ("spell-strength",       po::value<double>(&spell_strength)->default_value(spell_strength),                         "strength ~ Gamma(shape,rate)")
    ("spell-strength-shape", po::value<double>(&spell_strength_prior_shape)->default_value(spell_strength_prior_shape), "strength ~ Gamma(shape,rate)")
    ("spell-strength-rate",  po::value<double>(&spell_strength_prior_rate)->default_value(spell_strength_prior_rate),   "strength ~ Gamma(shape,rate)")
    
    ("spell-lambda",       po::value<double>(&spell_lambda)->default_value(spell_lambda),             "lambda for spell")
    ("spell-lambda-shape", po::value<double>(&spell_lambda_shape)->default_value(spell_lambda_shape), "lambda ~ Gamma(shape,rate)")
    ("spell-lambda-rate",  po::value<double>(&spell_lambda_rate)->default_value(spell_lambda_rate),   "lambda ~ Gamma(shape,rate)")
    
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

