//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-LM!

//
// @InProceedings{teh:2006:COLACL,
//   author    = {Teh, Yee Whye},
//   title     = {A Hierarchical Bayesian Language Model Based On Pitman-Yor Processes},
//   booktitle = {Proceedings of the 21st International Conference on Computational Linguistics and 44th Annual Meeting of the Association for Computational Linguistics},
//   month     = {July},
//   year      = {2006},
//   address   = {Sydney, Australia},
//   publisher = {Association for Computational Linguistics},
//   pages     = {985--992},
//   url       = {http://www.aclweb.org/anthology/P06-1124},
//   doi       = {10.3115/1220175.1220299}
// }
//

// parallel training by:
// @article{DBLP:journals/taslp/HuangR10,
//   author    = {Songfang Huang and
//                Steve Renals},
//   title     = {Hierarchical Bayesian Language Models for Conversational
//                Speech Recognition},
//   journal   = {IEEE Transactions on Audio, Speech {\&} Language Processing},
//   volume    = {18},
//   number    = {8},
//   year      = {2010},
//   pages     = {1941-1954},
//   ee        = {http://dx.doi.org/10.1109/TASL.2010.2040782},
//   bibsource = {DBLP, http://dblp.uni-trier.de}
// }

#include <map>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/lexical_cast.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/pyp_parameter.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/restaurant.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/trie_compact.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"
#include "utils/rwticket.hpp"
#include "utils/atomicop.hpp"
#include "utils/piece.hpp"
#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

typedef cicada::Symbol    word_type;
typedef cicada::Sentence  sentence_type;
typedef cicada::Vocab     vocab_type;

struct PYPLM
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef uint32_t  id_type;
  
  typedef boost::filesystem::path path_type;
  
  struct Node
  {
    typedef utils::restaurant<word_type, boost::hash<word_type>, std::equal_to<word_type>,
			      std::allocator<word_type > > table_type;
  
    Node() : table(), parent(id_type(-1)), order(0) {}
    
    table_type table;
    id_type parent;
    int     order;
  };
  typedef Node node_type;
  
  typedef utils::trie_compact<word_type, node_type,
			      utils::unassigned<word_type>, 
			      boost::hash<word_type>, std::equal_to<word_type>,
			      std::allocator<std::pair<const word_type, node_type> > > trie_type;
  
  typedef std::vector<node_type*, std::allocator<node_type*> > node_set_type;
  typedef std::vector<node_set_type, std::allocator<node_set_type> > node_map_type;
  
  typedef utils::pyp_parameter parameter_type;
  typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;

  struct shard_type
  {
    shard_type() : trie()  {}
    
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
	}
      }
      
      return node;
    }
    
    template <typename Sampler>
    bool increment(const word_type& word, const id_type& node, Sampler& sampler, const double temperature=1.0)
    {
      if (node == trie.root()) {
	root.table.increment(word, p0, sampler, temperature);
	
	return true;
      } else {
	const double backoff = prob(word, trie[node].parent);
	
	// we will also increment lower-order when new table is created!
	if (trie[node].table.increment(word, backoff, sampler, temperature))
	  return increment(word, trie[node].parent, sampler, temperature);
	else
	  return false;
      }
    }
    
    template <typename Sampler>
    bool decrement(const word_type& word, const id_type& node, Sampler& sampler)
    {
      if (node == trie.root()) {
	root.table.decrement(word, sampler);
	
	return true;
      } else {
	if (trie[node].table.decrement(word, sampler))
	  return decrement(word, trie[node].parent, sampler);
	else
	  return false;
      }
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
    double prob(const word_type& word, Iterator first, Iterator last) const
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
    
    trie_type trie;
    
    node_type root;
    double p0;
  };
  
  typedef std::vector<shard_type, std::allocator<shard_type> > shard_set_type;
  
  PYPLM(const size_type num_shard,
	const int order,
	const parameter_type& param)
    : shards(num_shard),
      root(),
      p0(1.0),
      counts0(0),
      nodes(order),
      parameters(order, param)
  { }
  
  size_type shard(const word_type& word) const
  {
    return word.id() % shards.size();
  }
  
  template <typename Iterator>
  double prob(const word_type& word, Iterator first, Iterator last) const 
  {
    if (! (first <= last)) // this may happen when accessing 0-gram!
      return p0;
    else if (first == last) // no context
      return root.table.prob(word, p0); 
    else // we will shard wrt the last word in the context!
      return shards[shard(*(last - 1))].prob(word, first, last);
  }
  
  template <typename Sampler>
  void increment(const word_type& word, Sampler& sampler, const double temperature=1.0)
  {
    root.table.increment(word, p0, sampler, temperature);
  }
  
  template <typename Sampler>
  void decrement(const word_type& word, Sampler& sampler)
  {
    root.table.decrement(word, sampler);
  }
  
  void initialize(const double __p0)
  {
    // assign p0
    p0 = __p0;
    for (size_type i = 0; i != shards.size(); ++ i)
      shards[i].p0 = __p0;
    
    // assign nodes
    nodes.clear();
    nodes.reserve(parameters.size());
    nodes.resize(parameters.size());
    
    nodes[0].push_back(&root);
    for (size_type i = 0; i != shards.size(); ++ i)
      for (PYPLM::id_type id = 0; id != shards[i].trie.size(); ++ id) {
	PYPLM::node_type& node = shards[i].trie[id];
	
	nodes[node.order].push_back(&node);
      }
  }
  
  double log_likelihood() const
  {
    double logprob = std::log(p0) * counts0;
    
    for (size_type order = 0; order != parameters.size(); ++ order)
      logprob += log_likelihood(order, parameters[order].discount, parameters[order].strength);
    
    return logprob;
  }
  
  double log_likelihood(const int order, const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    double logprob = parameters[order].log_likelihood(discount, strength);
    
    node_set_type::const_iterator niter_end = nodes[order].end();
    for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
      if (! (*niter)->table.empty())
	logprob += (*niter)->table.log_likelihood(discount, strength);
    
    return logprob;
  }
  
  struct DiscountSampler
  {
    DiscountSampler(const PYPLM& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPLM& pyplm;
    int order;
    
    double operator()(const double& proposed_discount) const
    {
      return pyplm.log_likelihood(order, proposed_discount, pyplm.parameters[order].strength);
    }
  };
  
  struct StrengthSampler
  {
    StrengthSampler(const PYPLM& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPLM& pyplm;
    int order;
    
    double operator()(const double& proposed_strength) const
    {
      return pyplm.log_likelihood(order, pyplm.parameters[order].discount, proposed_strength);
    }
  };
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    counts0 = root.table.size_table();
    
    for (size_type i = 0; i != shards.size(); ++ i)
      shards[i].root.table = root.table;
    
    for (size_type order = 0; order != parameters.size(); ++ order) {
      
      for (int iter = 0; iter != num_loop; ++ iter) {
	parameters[order].strength = sample_strength(order, sampler, parameters[order].discount, parameters[order].strength);
	parameters[order].verify_parameters();
	
	parameters[order].discount = sample_discount(order, sampler, parameters[order].discount, parameters[order].strength);
	parameters[order].verify_parameters();
      }
      
      parameters[order].strength = sample_strength(order, sampler, parameters[order].discount, parameters[order].strength);
      parameters[order].verify_parameters();
      
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
	(*niter)->table.discount() = parameters[order].discount;
	(*niter)->table.strength() = parameters[order].strength;
      }
    }
  }

  template <typename Sampler>
  double sample_strength(const int order, Sampler& sampler, const double& discount, const double& strength) const
  {
    double x = 0.0;
    double y = 0.0;

    node_set_type::const_iterator niter_end = nodes[order].end();
    for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
      if (! (*niter)->table.empty()) {
	x += (*niter)->table.sample_log_x(sampler, discount, strength);
	y += (*niter)->table.sample_y(sampler, discount, strength);
      }
    
    return sampler.gamma(parameters[order].strength_shape + y, parameters[order].strength_rate - x);
  }
  
  template <typename Sampler>
  double sample_discount(const int order, Sampler& sampler, const double& discount, const double& strength) const
  {
    double y = 0.0;
    double z = 0.0;
    
    node_set_type::const_iterator niter_end = nodes[order].end();
    for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
      if (! (*niter)->table.empty()) {
	y += (*niter)->table.sample_y_inv(sampler, discount, strength);
	z += (*niter)->table.sample_z_inv(sampler, discount, strength);
      }
    
    return sampler.beta(parameters[order].discount_alpha + y, parameters[order].discount_beta + z);
  }
  
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    counts0 = root.table.size_table();
    
    for (size_type i = 0; i != shards.size(); ++ i)
      shards[i].root.table = root.table;
    
    for (size_type order = 0; order != parameters.size(); ++ order) {
      DiscountSampler discount_sampler(*this, order);
      StrengthSampler strength_sampler(*this, order);

      for (int iter = 0; iter != num_loop; ++ iter) {
	parameters[order].strength = utils::slice_sampler(strength_sampler,
							  parameters[order].strength,
							  sampler,
							  - parameters[order].discount + std::numeric_limits<double>::min(),
							  std::numeric_limits<double>::infinity(),
							  0.0,
							  num_iterations,
							  100 * num_iterations);
	parameters[order].verify_parameters();
	
	parameters[order].discount = utils::slice_sampler(discount_sampler,
							  parameters[order].discount,
							  sampler,
							  (parameters[order].strength < 0.0 ? - parameters[order].strength : 0.0) + std::numeric_limits<double>::min(),
							  1.0,
							  0.0,
							  num_iterations,
							  100 * num_iterations);
	parameters[order].verify_parameters();
      }
      
      parameters[order].strength = utils::slice_sampler(strength_sampler,
							parameters[order].strength,
							sampler,
							- parameters[order].discount + std::numeric_limits<double>::min(),
							std::numeric_limits<double>::infinity(),
							0.0,
							num_iterations,
							100 * num_iterations);
      parameters[order].verify_parameters();
      
      node_set_type::const_iterator niter_end = nodes[order].end();
      for (node_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
	(*niter)->table.discount() = parameters[order].discount;
	(*niter)->table.strength() = parameters[order].strength;
      }
    }
  }

  struct vocab_map_type
  {
    typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > mapping_type;

    vocab_map_type() : mapping(), inverse() {}

    word_type::id_type operator[](const word_type& word)
    {
      if (word.id() >= mapping.size())
	mapping.resize(word.id() + 1, word_type::id_type(-1));
      
      if (mapping[word.id()] == word_type::id_type(-1)) {
	mapping[word.id()] = inverse.size();
	inverse.push_back(word.id());
      }
      
      return mapping[word.id()];
    }

    mapping_type mapping;
    mapping_type inverse;
  };
  
  template <typename Tp>
  struct greater_pcustomer
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->second.size_customer() > y->second.size_customer();
    }
  };
  
  void write(const path_type& path)
  {
    typedef utils::repository repository_type;
    typedef std::pair<id_type, id_type> shard_node_type;
    typedef std::vector<shard_node_type, std::allocator<shard_node_type> > node_set_type;

    typedef uint64_t count_type;

    typedef boost::fusion::tuple<shard_node_type, count_type, count_type> data_type;

    typedef std::map<word_type::id_type, data_type, std::less<word_type::id_type>, std::allocator<std::pair<const word_type::id_type, data_type> >  > word_set_type;
    typedef utils::succinct_vector<std::allocator<int32_t> > position_set_type;
    typedef std::vector<count_type, std::allocator<count_type> > offset_set_type;
      
    repository_type rep(path, repository_type::write);
    
    rep["order"] = boost::lexical_cast<std::string>(parameters.size());
    
    rep["p0"]      = boost::lexical_cast<std::string>(p0);
    rep["counts0"] = boost::lexical_cast<std::string>(counts0);
    
    rep["discount-alpha"] = boost::lexical_cast<std::string>(parameters[0].discount_alpha);
    rep["discount-beta"]  = boost::lexical_cast<std::string>(parameters[0].discount_beta);
    rep["strength-shape"] = boost::lexical_cast<std::string>(parameters[0].strength_shape);
    rep["strength-rate"]  = boost::lexical_cast<std::string>(parameters[0].strength_rate);
    
    for (size_type order = 0; order != parameters.size(); ++ order) {
      rep["discount" + boost::lexical_cast<std::string>(order)] = boost::lexical_cast<std::string>(parameters[order].discount);
      rep["strength" + boost::lexical_cast<std::string>(order)] = boost::lexical_cast<std::string>(parameters[order].strength);
    }
    
    vocab_map_type vocab_map;
    
    vocab_map[vocab_type::BOS];
    vocab_map[vocab_type::EOS];
    
    {
      // compute vocabulary mapping...
      typedef std::vector<const node_type::table_type::value_type*, std::allocator<const node_type::table_type::value_type*> > sorted_type;
      
      sorted_type sorted;
      sorted.reserve(root.table.size());
      
      node_type::table_type::const_iterator titer_end = root.table.end();
      for (node_type::table_type::const_iterator titer = root.table.begin(); titer != titer_end; ++ titer)
	sorted.push_back(&(*titer));
      
      std::sort(sorted.begin(), sorted.end(), greater_pcustomer<node_type::table_type::value_type>());
      
      sorted_type::const_iterator siter_end = sorted.end();
      for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
	vocab_map[(*siter)->first];
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
	words.insert(std::make_pair(vocab_map[titer->first], data_type(std::make_pair(id_type(-1), trie_type::npos()),
								       titer->second.size_customer(),
								       titer->second.size_table())));
      
      for (size_type shard = 0; shard != shards.size(); ++ shard) {
	trie_type::const_iterator iter_end = shards[shard].trie.end();
	for (trie_type::const_iterator iter = shards[shard].trie.begin(); iter != iter_end; ++ iter)
	  if (! shards[shard].trie[iter->second].table.empty()) {
	    std::pair<word_set_type::iterator, bool> result = words.insert(std::make_pair(vocab_map[iter->first],
											  data_type(std::make_pair(shard, iter->second), 0, 0)));
	    if (! result.second)
	      boost::fusion::get<0>(result.first->second) = std::make_pair(shard, iter->second);
	  }
      }
      
      word_set_type::const_iterator witer_end = words.end();
      for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) {
	nodes.push_back(boost::fusion::get<0>(witer->second));
	
	os_index.write((char*) &witer->first, sizeof(word_type::id_type));
	os_count.write((char*) &boost::fusion::get<1>(witer->second), sizeof(count_type));
	os_count.write((char*) &boost::fusion::get<2>(witer->second), sizeof(count_type));
      }
      
      offset += words.size();
      offsets.push_back(offset);
    }
    
    
    for (size_type order = 1; order < parameters.size(); ++ order) {
      nodes_next.clear();
      
      node_set_type::const_iterator niter_end = nodes.end();
      for (node_set_type::const_iterator niter = nodes.begin(); niter != niter_end; ++ niter) {
	if (niter->second == trie_type::npos()) {
	  const count_type count_customer = 0;
	  const count_type count_table    = 0;
	  
	  os_total.write((char*) &count_customer, sizeof(count_type));
	  os_total.write((char*) &count_table, sizeof(count_type));
	  
	  positions.set(positions.size(), false); // we will set bit!
	} else {
	  const id_type shard = niter->first;
	  const trie_type& trie = shards[shard].trie;
	  const node_type& node = trie[niter->second];
	  
	  trie_type::const_iterator iter_begin = trie.begin(niter->second);
	  trie_type::const_iterator iter_end   = trie.end(niter->second);
	  
	  const count_type count_customer = node.table.size_customer();
	  const count_type count_table    = node.table.size_table();
	  
	  os_total.write((char*) &count_customer, sizeof(count_type));
	  os_total.write((char*) &count_table, sizeof(count_type));
	  
	  // we have an entry in the model
	  words.clear();
	  node_type::table_type::const_iterator titer_end = node.table.end();
	  for (node_type::table_type::const_iterator titer = node.table.begin(); titer != titer_end; ++ titer)
	    words.insert(std::make_pair(vocab_map[titer->first], data_type(std::make_pair(shard, trie_type::npos()),
									   titer->second.size_customer(),
									   titer->second.size_table())));
	  
	  // or, we have an entry with the longer contexts..
	  for (trie_type::const_iterator iter = iter_begin; iter != iter_end; ++ iter)
	    if (! trie[iter->second].table.empty()) {
	      std::pair<word_set_type::iterator, bool> result = words.insert(std::make_pair(vocab_map[iter->first],
											    data_type(std::make_pair(shard, iter->second), 0, 0)));
	      if (! result.second)
		boost::fusion::get<0>(result.first->second) = std::make_pair(shard, iter->second);
	    }
	  
	  word_set_type::const_iterator witer_end = words.end();
	  for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) {
	    nodes_next.push_back(boost::fusion::get<0>(witer->second));
	    
	    os_index.write((char*) &witer->first, sizeof(word_type::id_type));
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
    vocab_type vocab;
    vocab.open(rep.path("vocab"), vocab_map.inverse.size() >> 1);
    
    vocab_map_type::mapping_type::const_iterator iiter_end = vocab_map.inverse.end();
    for (vocab_map_type::mapping_type::const_iterator iiter = vocab_map.inverse.begin(); iiter != iiter_end; ++ iiter)
      vocab.insert(static_cast<const std::string&>(word_type(*iiter)));

    vocab.close();
  }

public: 
  shard_set_type shards;
  
  node_type  root;
  double     p0;
  size_type  counts0;
  
  node_map_type nodes;
  parameter_set_type parameters;
};

typedef PYPLM::size_type size_type;
typedef PYPLM::id_type   id_type;

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef utils::sampler<boost::mt19937> sampler_type;
typedef std::vector<sampler_type, std::allocator<sampler_type> > sampler_set_type;

// we will precompute <word, node> pair...
struct data_type
{
  data_type(const word_type& __word, const id_type& __node)
    : count(1), word(__word), node(__node), rank(0) {}
  data_type(const word_type& __word, const id_type& __node, const size_type& __count)
    : count(__count), word(__word), node(__node), rank(0) {}
  
  size_type      count;
  word_type      word;
  PYPLM::id_type node;
  int            rank;
};
typedef std::vector<data_type, std::allocator<data_type> > data_set_type;
typedef std::vector<data_set_type, std::allocator<data_set_type> > data_map_type;
typedef std::vector<bool, std::allocator<bool> > non_oov_type;

path_set_type train_files;
path_set_type train_count_files;
path_set_type test_files;
path_type     output_file;

int order = 4;
int samples = 1;
int burns = 10;
int baby_steps = 0;
int anneal_steps = 0;
int resample_rate = 1;
int resample_iterations = 1;
bool slice_sampling = false;

double discount_alpha = 1.0;
double discount_beta  = 1.0;
double strength_shape = 1.0;
double strength_rate  = 1.0;

int threads = 1;
int debug = 0;


void read_training(const path_set_type& corpus,
		   const path_set_type& counts,
		   data_map_type& training,
		   non_oov_type& non_oov,
		   PYPLM& model);

void learn(const data_map_type& training,
	   PYPLM& model,
	   sampler_set_type& samplers);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    
    if (order <= 1)
      throw std::runtime_error("order must be positive and greater than 2");
    
    if (samples < 0)
      throw std::runtime_error("# of samples must be positive");
            
    if (resample_rate <= 0)
      throw std::runtime_error("resample rate must be >= 1");

    if (train_files.empty() && train_count_files.empty())
      throw std::runtime_error("no training data?");
    
    sampler_type sampler;
    sampler_set_type samplers(threads, sampler);
    
    PYPLM model(threads,
		order,
		PYPLM::parameter_type(discount_alpha,
				      discount_beta,
				      strength_shape,
				      strength_rate));
    
    data_map_type training;
    non_oov_type  non_oov;
    
    read_training(train_files, train_count_files, training, non_oov, model);
    
    // - 1 for BOS
    model.initialize(1.0 / (std::count(non_oov.begin(), non_oov.end(), true) - 1));
    
    if (debug) {
      std::cerr << "p0=" << model.p0 << std::endl;
      for (size_t i = 0; i != model.shards.size(); ++ i)
	std::cerr << "shard=" << i << " node size: " << model.shards[i].trie.size() << std::endl;
    }
    
    learn(training, model, samplers);
    
    // clear training data
    training.clear();
    
    // we will dump LM... now, define a format!
    if (! output_file.empty())
      model.write(output_file);
    
    // testing!
    if (! test_files.empty()) {
      sentence_type sentence;
      sentence_type ngram(1, vocab_type::BOS);

      double logprob_total = 0.0;
      double logprob = 0.0;
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
	    const double prob = model.prob(*siter, std::max(ngram.begin(), ngram.end() - order + 1), ngram.end());
	    const double lp = std::log(prob);
	    
	    if (! is_oov)
	      logprob_total += lp;
	    logprob += lp;
	    
	    num_oov += is_oov;
	    
	    ngram.push_back(*siter);
	  }
	  
	  const double prob = model.prob(vocab_type::EOS, std::max(ngram.begin(), ngram.end() - order + 1), ngram.end());
	  const double lp = std::log(prob);
	  
	  logprob_total += lp;
	  logprob += lp;
	  
	  num_word += sentence.size();
	  ++ num_sentence;
	}
      }
      
      std::cerr << "# of sentences: " << num_sentence
		<< " # of words: " << num_word
		<< " # of OOV: " << num_oov
		<< " order: " << order
		<< std::endl;
      
      std::cerr << "logprob       = " << logprob_total << std::endl;
      std::cerr << "logprob(+oov) = " << logprob << std::endl;
      std::cerr << "ppl           = " << utils::mathop::exp(- logprob_total / (num_word - num_oov + num_sentence)) << std::endl;
      std::cerr << "ppl1          = " << utils::mathop::exp(- logprob_total / (num_word - num_oov)) << std::endl;
      std::cerr << "ppl(+oov)     = " << utils::mathop::exp(- logprob / (num_word + num_sentence)) << std::endl;
      std::cerr << "ppl1(+oov)    = " << utils::mathop::exp(- logprob / (num_word)) << std::endl;
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct LearnMapReduce
{
  typedef std::pair<word_type, int> word_count_type;
  
  typedef utils::lockfree_list_queue<word_count_type, std::allocator<word_count_type> > queue_type;
};

struct LearnMapper
{
  typedef LearnMapReduce::word_count_type word_count_type;
  typedef LearnMapReduce::queue_type      queue_type;
  
  LearnMapper(queue_type& __queue,
	      const data_set_type& __training,
	      PYPLM::shard_type& __model,
	      sampler_type& __sampler,
	      const bool __baby,
	      const double __temperature,
	      const bool __remove)
    : queue(__queue), training(__training), model(__model), sampler(__sampler), baby(__baby), temperature(__temperature), remove(__remove) {}

  template <typename Training>
  struct less_rank
  {
    less_rank(const Training& __training) : training(__training) {}
    
    bool operator()(const size_t& x, const size_t& y) const
    {
      return training[x].rank < training[y].rank;
    }
    
    const Training& training;
  };
  
  void operator()()
  {
    typedef std::vector<size_type, std::allocator<size_type> > position_set_type;
    typedef std::vector<int, std::allocator<int> > count_set_type;

    position_set_type positions(training.size());
    count_set_type counts;
    
    for (size_type i = 0; i != training.size(); ++ i)
      positions[i] = i;
    
    boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
    std::random_shuffle(positions.begin(), positions.end(), gen);

    if (baby)
      std::sort(positions.begin(), positions.end(), less_rank<data_set_type>(training));
    
    position_set_type::const_iterator piter_end = positions.end();
    for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter) {
      const word_type& word = training[*piter].word;
      
      if (word.id() >= counts.size())
	counts.resize(word.id() + 1, 0);
      
      int& counter = counts[word.id()];
      
      for (size_type i = 0; i != training[*piter].count; ++ i) {
	if (remove)
	  counter -= model.decrement(training[*piter].word, training[*piter].node, sampler);
	
	counter += model.increment(training[*piter].word, training[*piter].node, sampler, temperature);
      }
    }
    
    for (word_type::id_type id = 0; id != counts.size(); ++ id)
      if (counts[id])
	queue.push(std::make_pair(word_type(id), counts[id]));
    
    queue.push(std::make_pair(word_type(), 0));
  }
  
  queue_type&          queue;
  const data_set_type& training;
  PYPLM::shard_type&   model;
  sampler_type&        sampler;
  
  const bool           baby;
  const double         temperature;
  const bool           remove;
};

struct LearnReducer
{
  typedef LearnMapReduce::word_count_type word_count_type;
  typedef LearnMapReduce::queue_type      queue_type;
  
  LearnReducer(const size_type __size,
              queue_type& __queue,
              PYPLM& __model,
              sampler_type& __sampler,
              const double __temperature)
    : reduce_size(__size), queue(__queue), model(__model), sampler(__sampler), temperature(__temperature) {}
  
  void operator()()
  {
    size_type finished = 0;
    
    while (finished != reduce_size) {
      word_count_type word_count;
      queue.pop(word_count);
      
      if (word_count.first.empty()) {
       ++ finished;
       continue;
      }
      
      if (word_count.second < 0) {
       for (int i = 0; i != - word_count.second; ++ i)
         model.decrement(word_count.first, sampler);
      } else if (word_count.second > 0) {
       for (int i = 0; i != word_count.second; ++ i)
         model.increment(word_count.first, sampler, temperature);
      }
    }
  }
  
  const size_type      reduce_size;
  queue_type&          queue;
  PYPLM&               model;
  sampler_type&        sampler;
  
  const double         temperature;
};


void learn(const data_map_type& training,
	   PYPLM& model,
	   sampler_set_type& samplers)
{
  if (samplers.size() != training.size() || model.shards.size() != training.size())
    throw std::runtime_error("invalid shard");
  
  const size_t shards_size = model.shards.size();
  
  size_t anneal_iter = 0;
  const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);

  size_t baby_iter = 0;
  const size_t baby_last = utils::bithack::branch(baby_steps > 0, baby_steps, 0);
    
  size_t burn_iter = 0;
  const size_t burn_last = utils::bithack::branch(burns > 0, burns, 0);
  
  bool sampling = false;
  int sample_iter = 0;

  sampler_type sampler = samplers.front();
  
  // sample parameters, first...
  if (slice_sampling)
    model.slice_sample_parameters(sampler, resample_iterations);
  else
    model.sample_parameters(sampler, resample_iterations);
  
  if (debug >= 2)
    for (int n = 0; n != order; ++ n)
      std::cerr << "order=" << n << " discount=" << model.parameters[n].discount << " strength=" << model.parameters[n].strength << std::endl;
  
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
      
    bool baby_finished = true;
    if (baby_iter != baby_last) {
      ++ baby_iter;
      baby_finished = false;
    }
    
    bool burn_finished = true;
    if (burn_iter != burn_last) {
      ++ burn_iter;
      burn_finished = false;
    }
      
    sampling = anneal_finished && baby_finished && burn_finished;
    
    if (debug) {
      if (sampling)
	std::cerr << "sampling iteration: " << (iter + 1) << std::endl;
      else
	std::cerr << "burn-in iteration: " << (iter + 1) << std::endl;
    }

    LearnMapReduce::queue_type queue;
    
    boost::thread_group workers;
    for (size_type i = 0; i != shards_size; ++ i)
      workers.add_thread(new boost::thread(LearnMapper(queue,
						       training[i],
						       model.shards[i],
						       samplers[i],
						       ! baby_finished,
						       temperature,
						       iter)));
    
    // reduce!
    LearnReducer(shards_size, queue, model, sampler, temperature)();

    // join
    workers.join_all();
    
    if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
      if (slice_sampling)
	model.slice_sample_parameters(sampler, resample_iterations);
      else
	model.sample_parameters(sampler, resample_iterations);
      
      if (debug >= 2)
	for (int n = 0; n != order; ++ n)
	  std::cerr << "order=" << n << " discount=" << model.parameters[n].discount << " strength=" << model.parameters[n].strength << std::endl;
    }

    if (debug)
      std::cerr << "log-likelihood: " << model.log_likelihood() << std::endl;
  }
}

struct ReadMapReduce
{
  typedef utils::simple_vector<word_type, std::allocator<word_type> > ngram_type;

  typedef std::pair<ngram_type, size_type> ngram_count_type;
  
  typedef utils::lockfree_list_queue<ngram_count_type, std::allocator<ngram_count_type> > queue_type;
  typedef std::vector<queue_type, std::allocator<queue_type> > queue_set_type;
};

namespace std
{
  inline
  void swap(ReadMapReduce::ngram_count_type& x, ReadMapReduce::ngram_count_type& y)
  {
    x.first.swap(y.first);
    std::swap(x.second, y.second);
  }
};

struct ReadMapper
{
  typedef ReadMapReduce::ngram_type       ngram_type;
  typedef ReadMapReduce::ngram_count_type ngram_count_type;
  typedef ReadMapReduce::queue_type       queue_type;
  
  ReadMapper(queue_type& __queue,
	     data_set_type& __training,
	     PYPLM::shard_type& __model)
    : queue(__queue), training(__training), model(__model)  {}
  
  void operator()()
  {
    typedef utils::unordered_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_unique_type;
    typedef utils::chunk_vector<word_unique_type, 4096 / sizeof(word_unique_type), std::allocator<word_unique_type> > word_unique_set_type;
    
    word_unique_set_type uniques;
    ngram_count_type     ngram_count;
    
    for (;;) {
      queue.pop_swap(ngram_count);
      
      if (ngram_count.first.empty()) break;
      
      training.push_back(data_type(ngram_count.first.back(),
				   model.insert(ngram_count.first.begin(), ngram_count.first.end() - 1),
				   ngram_count.second));
      
      if (training.back().node >= uniques.size())
	uniques.resize(training.back().node + 1);
      uniques[training.back().node].insert(training.back().word);
    }
    
    data_set_type::iterator titer_end = training.end();
    for (data_set_type::iterator titer = training.begin(); titer != titer_end; ++ titer)
      titer->rank = uniques[titer->node].size();
  }
  
  queue_type&        queue;
  data_set_type&     training;
  PYPLM::shard_type& model;
};
  
inline
word_type escape_word(const utils::piece& __word)
{
  static const std::string& __BOS = static_cast<const std::string&>(vocab_type::BOS);
  static const std::string& __EOS = static_cast<const std::string&>(vocab_type::EOS);
  static const std::string& __UNK = static_cast<const std::string&>(vocab_type::UNK);
  
  const utils::ipiece word(__word);
  
  if (word == __BOS)
    return vocab_type::BOS;
  else if (word == __EOS)
    return vocab_type::EOS;
  else if (word == __UNK)
    return vocab_type::UNK;
  else
    return __word;
}

void read_training(const path_set_type& corpus_files,
		   const path_set_type& counts_files,
		   data_map_type& training,
		   non_oov_type& non_oov,
		   PYPLM& model)
{
  typedef ReadMapReduce::ngram_type       ngram_type;
  typedef ReadMapReduce::ngram_count_type ngram_count_type;

  const size_t shards_size = model.shards.size();

  training.clear();
  training.reserve(shards_size);
  training.resize(shards_size);
  
  ReadMapReduce::queue_set_type queues(shards_size);
  
  boost::thread_group workers;
  for (size_t i = 0; i != shards_size; ++ i)
    workers.add_thread(new boost::thread(ReadMapper(queues[i], training[i], model.shards[i])));

  if (vocab_type::BOS.id() >= non_oov.size())
    non_oov.resize(vocab_type::BOS.id() + 1, false);
  if (vocab_type::EOS.id() >= non_oov.size())
    non_oov.resize(vocab_type::EOS.id() + 1, false);

  ngram_count_type ngram_count;
  
  for (path_set_type::const_iterator fiter = corpus_files.begin(); fiter != corpus_files.end(); ++ fiter) {
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    sentence_type sentence;
    sentence_type ngram(1, vocab_type::BOS);
    
    while (is >> sentence) {
      ngram.resize(1);
      
      sentence_type::const_iterator siter_end = sentence.end();
      for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	const size_t shard = model.shard(ngram.back());
	ngram.push_back(*siter);
	
	ngram_count.first  = ngram_type(std::max(ngram.begin(), ngram.end() - order), ngram.end());
	ngram_count.second = 1;
	
	queues[shard].push_swap(ngram_count);
	
	if (siter->id() >= non_oov.size())
	  non_oov.resize(siter->id() + 1, false);
	non_oov[siter->id()] = true;
      }
      
      const size_t shard = model.shard(ngram.back());
      ngram.push_back(vocab_type::EOS);
      
      ngram_count.first  = ngram_type(std::max(ngram.begin(), ngram.end() - order), ngram.end());
      ngram_count.second = 1;
      
      queues[shard].push_swap(ngram_count);
    }
  }
  
  for (path_set_type::const_iterator fiter = counts_files.begin(); fiter != counts_files.end(); ++ fiter) {
    typedef std::vector<utils::piece, std::allocator<utils::piece> > tokens_type;
    typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
      
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    std::string line;
    tokens_type tokens;
    
    sentence_type ngram;
    
    while (std::getline(is, line)) {
      utils::piece line_piece(line);
      tokenizer_type tokenizer(line_piece);
	
      tokens.clear();
      tokens.insert(tokens.end(), tokenizer.begin(), tokenizer.end());
	
      // exclude unigram, exclude non-ordered, or not prefixed by BOS
      if (tokens.size() == 2) continue;
      if (static_cast<int>(tokens.size()) != order + 1 || escape_word(tokens.front()) != vocab_type::BOS) continue;
	
      ngram.clear();
	
      tokens_type::const_iterator titer_end = tokens.end() - 1;
      for (tokens_type::const_iterator titer = tokens.begin(); titer != titer_end; ++ titer)
	ngram.push_back(escape_word(*titer));
	
      ngram_count.first  = ngram_type(ngram.begin(), ngram.end());
      ngram_count.second = utils::lexical_cast<size_type>(tokens.back());
	
      queues[model.shard(ngram[ngram.size() - 2])].push_swap(ngram_count);
      
      if (ngram.back().id() >= non_oov.size())
	non_oov.resize(ngram.back().id() + 1, false);
      non_oov[ngram.back().id()] = true;
    }
  }
  
  for (size_t i = 0; i != shards_size; ++ i)
    queues[i].push(ngram_count_type());
  
  workers.join_all();
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("train",       po::value<path_set_type>(&train_files)->multitoken(), "train file(s)")
    ("train-count", po::value<path_set_type>(&train_count_files)->multitoken(), "train count file(s)")
    ("test",  po::value<path_set_type>(&test_files)->multitoken(),  "test file(s)")
    ("output", po::value<path_type>(&output_file), "output file")
    
    ("order", po::value<int>(&order)->default_value(order), "max ngram order")
    
    ("samples",             po::value<int>(&samples)->default_value(samples),                         "# of samples")
    ("burns",               po::value<int>(&burns)->default_value(burns),                             "# of burn-ins")
    ("baby-steps",          po::value<int>(&baby_steps)->default_value(baby_steps),                   "# of baby steps")
    ("anneal-steps",        po::value<int>(&anneal_steps)->default_value(anneal_steps),               "# of anneal steps")
    ("resample",            po::value<int>(&resample_rate)->default_value(resample_rate),             "hyperparameter resample rate")
    ("resample-iterations", po::value<int>(&resample_iterations)->default_value(resample_iterations), "hyperparameter resample iterations")
    
    ("slice",               po::bool_switch(&slice_sampling),                                         "slice sampling for hyperparameters")
    
    ("discount-alpha", po::value<double>(&discount_alpha)->default_value(discount_alpha), "discount ~ Beta(alpha,beta)")
    ("discount-beta",  po::value<double>(&discount_beta)->default_value(discount_beta),   "discount ~ Beta(alpha,beta)")

    ("strength-shape", po::value<double>(&strength_shape)->default_value(strength_shape), "strength ~ Gamma(shape,rate)")
    ("strength-rate",  po::value<double>(&strength_rate)->default_value(strength_rate),   "strength ~ Gamma(shape,rate)")
    
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

