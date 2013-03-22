//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// PYP-LM-mixture!

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
//
// and
//
// @inproceedings{WooTeh2009a,
//  Author =        {F. Wood and Y. W. Teh},
//  Title =         {A Hierarchical Nonparametric {B}ayesian Approach 
//                  to Statistical Language Model Domain Adaptation},
//  Booktitle =     {Proceedings of the International Conference on 
//                  Artificial Intelligence and Statistics},
//  Volume =        {12},
//  Year =          {2009}}
//

// TODO:
//
// fix so that we also include mapping from normalized word into surface form, defined in each latent LM
//

#include <map>
#include <iterator>
#include <numeric>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/stemmer.hpp>
#include <cicada/cluster.hpp>
#include <cicada/cluster_stemmer.hpp>

#include "utils/alloc_vector.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/pyp_parameter.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/restaurant.hpp"
#include "utils/restaurant_floor.hpp"
#include "utils/restaurant_vector.hpp"
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

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

typedef cicada::Symbol    word_type;
typedef cicada::Sentence  sentence_type;
typedef cicada::Vocab     vocab_type;
typedef cicada::Stemmer   stemmer_type;
typedef cicada::Cluster   cluster_type;
typedef cicada::ClusterStemmer normalizer_type;

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
    
    typedef utils::rwticket mutex_type;
  
    Node() : table(), parent(id_type(-1)), order(0), mutex()  {}
    Node(const Node& x) : table(x.table), parent(x.parent), order(x.order), mutex() {}
  
    table_type table;
    id_type parent;
    int     order;
  
    mutex_type mutex;
  };
  typedef Node node_type;
  
  typedef utils::trie_compact<word_type, node_type,
			      utils::unassigned<word_type>, 
			      boost::hash<word_type>, std::equal_to<word_type>,
			      std::allocator<std::pair<const word_type, node_type> > > trie_type;

  typedef std::vector<node_type*, std::allocator<node_type*> > node_ptr_set_type;
  typedef std::vector<node_ptr_set_type, std::allocator<node_ptr_set_type> > node_ptr_map_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > history_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  typedef std::vector<prob_set_type, std::allocator<double> > prob_map_type;

  typedef utils::pyp_parameter parameter_type;
  typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;
  
  PYPLM(const int order,
	const double __p0,
	const parameter_type& __parameter)
    : trie(),
      nodes(),
      parameters(order, __parameter),
      p0(__p0),
      counts0(0)
  { }

  template <typename Iterator>
  id_type insert(Iterator first, Iterator last, const normalizer_type& normalizer)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    int order = 1;
    id_type node = trie.root();
    for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
      const word_type normalized = normalizer(*iter);
      const id_type node_prev = node;
      
      node = trie.insert(node, normalized);
      
      node_type& trie_node = trie[node];
      
      if (! trie_node.order) {
	trie_node.parent = node_prev;
	trie_node.order  = order;
      }
    }
    
    return node;
  }
  
  template <typename Sampler>
  void increment(const word_type& word, const id_type& node, Sampler& sampler, const double temperature=1.0)
  {
    if (node == trie.root()) {
      node_type::mutex_type::scoped_writer_lock lock(root.mutex);
      
      counts0 += root.table.increment(word, p0, sampler, temperature);
    } else {
      const double backoff = prob(word, trie[node].parent);
      
      bool propagate = false;
      {
	node_type::mutex_type::scoped_writer_lock lock(trie[node].mutex);
	
	propagate = trie[node].table.increment(word, backoff, sampler, temperature);
      }
      
      if (propagate)
	increment(word, trie[node].parent, sampler, temperature);
    }
  }
  
  template <typename Sampler>
  void decrement(const word_type& word, const id_type& node, Sampler& sampler)
  {
    if (node == trie.root()) {
      node_type::mutex_type::scoped_writer_lock lock(root.mutex);
      
      counts0 -= root.table.decrement(word, sampler);
    } else {
      bool propagate = false;
      {
	node_type::mutex_type::scoped_writer_lock lock(trie[node].mutex);
	
	propagate = trie[node].table.decrement(word, sampler);
      }
      
      if (propagate)
	decrement(word, trie[node].parent, sampler);
    }
  }
  
  double prob(const word_type& word, const id_type& node) const
  {
    if (node == trie.root()) {
      node_type::mutex_type::scoped_reader_lock lock(const_cast<node_type&>(root).mutex);
      
      return root.table.prob(word, p0);
    } else {
      const double p = prob(word, trie[node].parent);
      
      node_type::mutex_type::scoped_reader_lock lock(const_cast<node_type&>(trie[node]).mutex);
      
      return (trie[node].table.empty() ? p : trie[node].table.prob(word, p));
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
    
    node_ptr_set_type::const_iterator niter_end = nodes[order].end();
    for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
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
    initialize();
    
    for (size_type order = 0; order != parameters.size(); ++ order) {
      
      for (int iter = 0; iter != num_loop; ++ iter) {
	parameters[order].strength = sample_strength(order, sampler, parameters[order].discount, parameters[order].strength);
	parameters[order].verify_parameters();
	
	parameters[order].discount = sample_discount(order, sampler, parameters[order].discount, parameters[order].strength);
	parameters[order].verify_parameters();
      }
      
      parameters[order].strength = sample_strength(order, sampler, parameters[order].discount, parameters[order].strength);
      parameters[order].verify_parameters();
      
      node_ptr_set_type::const_iterator niter_end = nodes[order].end();
      for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
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

    node_ptr_set_type::const_iterator niter_end = nodes[order].end();
    for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
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
    
    node_ptr_set_type::const_iterator niter_end = nodes[order].end();
    for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
      if (! (*niter)->table.empty()) {
	y += (*niter)->table.sample_y_inv(sampler, discount, strength);
	z += (*niter)->table.sample_z_inv(sampler, discount, strength);
      }
    
    return sampler.beta(parameters[order].discount_alpha + y, parameters[order].discount_beta + z);
  }
  
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    initialize();

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
      
      node_ptr_set_type::const_iterator niter_end = nodes[order].end();
      for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
	(*niter)->table.discount() = parameters[order].discount;
	(*niter)->table.strength() = parameters[order].strength;
      }
    }
  }

  void initialize()
  {
    if (nodes.size() == parameters.size()) return;
    
    nodes.clear();
    nodes.reserve(parameters.size());
    nodes.resize(parameters.size());
    
    nodes[0].push_back(&root);
    
    for (id_type id = 0; id != trie.size(); ++ id) {
      node_type& node = trie[id];
      
      nodes[node.order].push_back(&node);
    }
  }
  
  void write(const path_type& path, const normalizer_type& normalizer)
  {
    typedef utils::repository repository_type;
    typedef std::vector<id_type, std::allocator<id_type> > node_set_type;

    typedef uint64_t count_type;

    typedef boost::fusion::tuple<id_type, count_type, count_type> data_type;

    typedef std::map<word_type, data_type, std::less<word_type>, std::allocator<std::pair<const word_type, data_type> >  > word_set_type;
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
    
    // dump stemming algorithm or copy cluster file...
    if (normalizer.is_stemmer())
      rep["stemmer"] = normalizer.stemmer_algorithm();
    else {
      cluster_type cluster(normalizer.cluster_path());
      
      cluster.write(rep.path("cluster"));
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
    
    for (size_type order = 1; order < parameters.size(); ++ order) {
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
  }
  
public: 
  trie_type trie;
  node_type root;
  node_ptr_map_type nodes;
  
  parameter_set_type parameters;

  double    p0;
  size_type counts0;
};

struct PYPMixture
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  typedef uint32_t  id_type;
  
  typedef boost::filesystem::path path_type;

  typedef utils::pyp_parameter parameter_type;
  typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;

  struct Mixture
  {
    typedef utils::restaurant_vector<> table_type;
    typedef std::vector<double, std::allocator<double> > prob_set_type;
    typedef std::vector<size_type, std::allocator<size_type> > count_set_type;

    typedef utils::spinlock mutex_type;

    Mixture() : table(), probs(), counts(), counts0(0) {}
    Mixture(const parameter_type& parameter, const size_type size)
      : table(parameter),
	probs(size, 1.0 / size),
	counts(size, 0),
	counts0(0) {}

    Mixture(const Mixture& x)
      : table(x.table),
	probs(x.probs),
	counts(x.counts),
	counts0(x.counts0) {}
    
    template <typename Sampler>
    void increment(const size_type pos, Sampler& sampler, const double temperature)
    {
      mutex_type::scoped_try_lock lock(mutex);
      
      if (lock)
	counts0 += table.increment(pos, 1.0 / counts.size(), sampler, temperature);
      else
	utils::atomicop::fetch_and_add(counts[pos], size_type(1));
    }

    void clear()
    {
      // clear count table...
      table.clear();
      counts0 = 0;
      std::fill(counts.begin(), counts.end(), 0);
    }

    template <typename Sampler>
    void sample_parameters(Sampler& sampler, const double temperature, const int num_loop = 2, const int num_iterations = 8)
    {
      // transform counts into table-counts..
      const double p0 = 1.0 / counts.size();
      
      for (size_type value = 0; value != counts.size(); ++ value)
	for (size_type i = 0; i != counts[value]; ++ i)
	  counts0 += table.increment(value, p0, sampler, temperature);
      
      table.sample_parameters(sampler, num_loop, num_iterations);
      
      probs.resize(counts.size());
      for (size_type value = 0; value != counts.size(); ++ value)
	probs[value] = table.prob(value, p0);
    }

    template <typename Sampler>
    void slice_sample_parameters(Sampler& sampler, const double temperature, const int num_loop = 2, const int num_iterations = 8)
    {
      // transform counts into table-counts..
      const double p0 = 1.0 / counts.size();
      
      for (size_type value = 0; value != counts.size(); ++ value)
	for (size_type i = 0; i != counts[value]; ++ i)
	  counts0 += table.increment(value, p0, sampler, temperature);
      
      table.slice_sample_parameters(sampler, num_loop, num_iterations);
      
      probs.resize(counts.size());
      for (size_type value = 0; value != counts.size(); ++ value)
	probs[value] = table.prob(value, p0);
    }
    
    double log_likelihood() const
    {
      return table.log_likelihood() - std::log(double(counts.size())) * counts0;
    }
    
    table_type     table;
    prob_set_type  probs;
    count_set_type counts;
    size_type      counts0;

    mutex_type mutex;
  };
  
  typedef Mixture mixture_type;
  typedef std::vector<mixture_type, std::allocator<mixture_type> >  mixture_set_type;
  
  struct Node
  {
    typedef utils::restaurant_floor<16, word_type, boost::hash<word_type>, std::equal_to<word_type>,
				    std::allocator<word_type > > table_type;
  
    typedef utils::rwticket mutex_type;
  
    Node() : table(), parent(id_type(-1)), order(0), mutex()  {}
    Node(const Node& x) : table(x.table), parent(x.parent), order(x.order), mutex() {}
    
    table_type   table;
    id_type parent;
    int     order;
    
    mutex_type mutex;
  };
  
  typedef Node node_type;
  
  typedef utils::trie_compact<word_type, node_type,
			      utils::unassigned<word_type>, 
			      boost::hash<word_type>, std::equal_to<word_type>,
			      std::allocator<std::pair<const word_type, node_type> > > trie_type;

  typedef std::vector<node_type*, std::allocator<node_type*> > node_ptr_set_type;
  typedef std::vector<node_ptr_set_type, std::allocator<node_ptr_set_type> > node_ptr_map_type;

  typedef std::vector<id_type, std::allocator<id_type> > node_set_type;
  typedef std::vector<id_type, std::allocator<id_type> > history_type;
  typedef std::vector<double, std::allocator<double> > prob_set_type;
  typedef std::vector<prob_set_type, std::allocator<double> > prob_map_type;
    
  typedef PYPLM lm_type;
  typedef std::vector<lm_type, std::allocator<lm_type> > lm_set_type;
  
  typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;
  
  PYPMixture(lm_set_type& __models,
	     const int order,
	     const double __p0,
	     const parameter_type& parameter,
	     const parameter_type& parameter_mixture)
    : models(__models),
      mixtures(order, mixture_type(parameter_mixture, __models.size() + 1)),
      parameters(order, parameter),
      trie(),
      nodes(),
      p0(__p0),
      counts0(0)
  { }
  
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
	trie_node.parent  = node_prev;
	trie_node.order   = order;
      }
    }
    
    return node;
  }
  
  template <typename Iterator, typename Norms, typename Sampler>
  void increment(const word_type& word, Iterator first, Norms normalizer, Sampler& sampler, const double temperature=1.0)
  {
    if (*first == trie.root()) {
      prob_set_type probs(models.size() + 1, p0);
      
      for (size_type i = 0; i != models.size(); ++ i)
	probs[i + 1] = models[i].prob(word, first[i + 1]);
      
      std::pair<size_type, bool> result;
      
      {
	node_type::mutex_type::scoped_writer_lock lock(root.mutex);
	
	result = root.table.increment(word,
				      probs.begin(),
				      probs.end(),
				      mixtures[0].probs.begin(),
				      sampler,
				      temperature);
      }
      
      if (result.second) {
	if (result.first == 0)
	  utils::atomicop::fetch_and_add(counts0, size_type(1));
	else
	  models[result.first - 1].increment(word, first[result.first], sampler, temperature);
      }
      
      mixtures[0].increment(result.first, sampler, temperature);
      
    } else {
      node_type& node = trie[*first];
      
      node_set_type nodes(models.size() + 1, node.parent);
      prob_set_type probs(models.size() + 1, 0.0);
      
      for (size_type i = 0; i != models.size(); ++ i) {
	nodes[i + 1] = models[i].trie[first[i + 1]].parent;
	probs[i + 1] = models[i].prob(word, first[i + 1]);
      }
      
      probs.front() = prob(word, nodes.begin(), normalizer);

      std::pair<size_type, bool> result;
      
      {
	node_type::mutex_type::scoped_writer_lock lock(node.mutex);
	
	result = node.table.increment(word,
				      probs.begin(),
				      probs.end(),
				      mixtures[node.order].probs.begin(),
				      sampler,
				      temperature);
      }
      
      if (result.second) {
	if (result.first == 0)
	  increment(word, nodes.begin(), normalizer, sampler, temperature);
	else
	  models[result.first - 1].increment(word, first[result.first], sampler, temperature);
      }
      
      mixtures[node.order].increment(result.first, sampler, temperature);
    }
  }
  
  template <typename Iterator, typename Norms, typename Sampler>
  void decrement(const word_type& word, Iterator first, Norms normalizer, Sampler& sampler)
  {
    const bool is_root = (*first == trie.root());
    
    node_type& node = (is_root ? root : trie[*first]);

    std::pair<size_type, bool> result;

    {
      node_type::mutex_type::scoped_writer_lock lock(node.mutex);
      
      result = node.table.decrement(word, sampler);
    }
    
    if (result.second) {
      if (result.first == 0) {
	if (is_root)
	  utils::atomicop::fetch_and_add(counts0, size_type(-1));
	else {
	  // we need to decrement to parents!
	  node_set_type nodes(models.size() + 1, trie[*first].parent);
	  
	  for (size_type i = 0; i != models.size(); ++ i)
	    nodes[i + 1] = models[i].trie[first[i + 1]].parent;
	  
	  decrement(word, nodes.begin(), normalizer, sampler);
	}
      } else
	models[result.first - 1].decrement(word, first[result.first], sampler);
    }
  }
  
  template <typename Iterator, typename Norms>
  double prob(const word_type& word, Iterator first, Norms normalizer) const
  {
    if (*first == trie.root()) {
      prob_set_type probs(models.size() + 1, p0);
      
      for (size_type i = 0; i != models.size(); ++ i)
	probs[i + 1] = models[i].prob(word, first[i + 1]);
      
      node_type::mutex_type::scoped_reader_lock lock(const_cast<node_type&>(root).mutex);
      
      return root.table.prob(word, probs.begin(), probs.end(), mixtures[0].probs.begin());
    } else {
      const node_type& node = trie[*first];
      
      node_set_type nodes(models.size() + 1, node.parent);
      prob_set_type probs(models.size() + 1, 0.0);
      
      for (size_type i = 0; i != models.size(); ++ i) {
	nodes[i + 1] = models[i].trie[first[i + 1]].parent;
	probs[i + 1] = models[i].prob(word, first[i + 1]);
      }
      
      probs.front() = prob(word, nodes.begin(), normalizer);
      
      node_type::mutex_type::scoped_reader_lock lock(const_cast<node_type&>(node).mutex);
      
      return node.table.prob(word, probs.begin(), probs.end(), mixtures[node.order].probs.begin());
    }
  }
  
  template <typename Iterator, typename Normalizers>
  double prob(const word_type& word, Iterator first, Iterator last, Normalizers normalizer) const
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    // how do we compute this...?
    // we will traverse from the back!
    reverse_iterator begin(last);
    reverse_iterator end(first);
    
    node_set_type nodes(models.size() + 1, trie.root());
    prob_set_type probs(models.size() + 1, p0);
    
    for (size_type i = 0; i != models.size(); ++ i) {
      nodes[i + 1] = models[i].trie.root();
      probs[i + 1] = models[i].prob(word, nodes[i + 1]);
    }
    
    probs.front() = root.table.prob(word, probs.begin(), probs.end(), mixtures[0].probs.begin());
    
    int order = 1;
    for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
      
      // update by the latent models
      for (size_type i = 0; i != models.size(); ++ i) 
	if (order == 1 || nodes[i + 1] != models[i].trie.root()) {
	  nodes[i + 1] = models[i].trie.find(nodes[i + 1], normalizer[i](*iter));
	  
	  if (nodes[i + 1] != models[i].trie.root())
	    probs[i + 1] = models[i].trie[nodes[i + 1]].table.prob(word, probs[i + 1]);
	}
      
      // update by the model
      if (order == 1 || nodes.front() != trie.root())
	nodes.front() = trie.find(nodes.front(), *iter);
      
      if (nodes.front() != trie.root())
	probs.front() = trie[nodes.front()].table.prob(word, probs.begin(), probs.end(), mixtures[order].probs.begin());
      else
	probs.front() = std::inner_product(probs.begin(), probs.end(), mixtures[order].probs.begin(), 0.0);
    }
    
    return probs.front();
  }

  double log_likelihood() const
  {
    double logprob = std::log(p0) * counts0;

    for (size_type i = 0; i != models.size(); ++ i)
      logprob += models[i].log_likelihood();

    for (size_type order = 0; order != mixtures.size(); ++ order)
      logprob += mixtures[order].log_likelihood();
    
    for (size_type order = 0; order != parameters.size(); ++ order)
      logprob += log_likelihood(order, parameters[order].discount, parameters[order].strength);
    
    return logprob;
  }
  
  double log_likelihood(const int order, const double& discount, const double& strength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    double logprob = parameters[order].log_likelihood(discount, strength);
    
    node_ptr_set_type::const_iterator niter_end = nodes[order].end();
    for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
      if (! (*niter)->table.empty())
	logprob += (*niter)->table.log_likelihood(discount, strength);
    
    return logprob;
  }
  
  struct DiscountSampler
  {
    DiscountSampler(const PYPMixture& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPMixture& pyplm;
    int order;
    
    double operator()(const double& proposed_discount) const
    {
      return pyplm.log_likelihood(order, proposed_discount, pyplm.parameters[order].strength);
    }
  };
  
  struct StrengthSampler
  {
    StrengthSampler(const PYPMixture& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPMixture& pyplm;
    int order;
    
    double operator()(const double& proposed_strength) const
    {
      return pyplm.log_likelihood(order, pyplm.parameters[order].discount, proposed_strength);
    }
  };

  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const double temperature, const int num_loop = 2, const int num_iterations = 8)
  {
    initialize();
    
    for (size_type i = 0; i != models.size(); ++ i)
      models[i].slice_sample_parameters(sampler, num_loop, num_iterations);
    
    for (size_type order = 0; order != mixtures.size(); ++ order)
      mixtures[order].slice_sample_parameters(sampler, temperature, num_loop, num_iterations);
    
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
      
      node_ptr_set_type::const_iterator niter_end = nodes[order].end();
      for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
	(*niter)->table.discount() = parameters[order].discount;
	(*niter)->table.strength() = parameters[order].strength;
      }
    }
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const double temperature, const int num_loop = 2, const int num_iterations = 8)
  {
    initialize();
    
    for (size_type i = 0; i != models.size(); ++ i)
      models[i].sample_parameters(sampler, num_loop, num_iterations);
    
    for (size_type order = 0; order != mixtures.size(); ++ order)
      mixtures[order].sample_parameters(sampler, temperature, num_loop, num_iterations);
    
    for (size_type order = 0; order != parameters.size(); ++ order) {
      
      for (int iter = 0; iter != num_loop; ++ iter) {
	parameters[order].strength = sample_strength(order, sampler, parameters[order].discount, parameters[order].strength);
	parameters[order].verify_parameters();
	
	parameters[order].discount = sample_discount(order, sampler, parameters[order].discount, parameters[order].strength);
	parameters[order].verify_parameters();
      }
      
      parameters[order].strength = sample_strength(order, sampler, parameters[order].discount, parameters[order].strength);
      parameters[order].verify_parameters();
      
      node_ptr_set_type::const_iterator niter_end = nodes[order].end();
      for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter) {
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

    node_ptr_set_type::const_iterator niter_end = nodes[order].end();
    for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
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
    
    node_ptr_set_type::const_iterator niter_end = nodes[order].end();
    for (node_ptr_set_type::const_iterator niter = nodes[order].begin(); niter != niter_end; ++ niter)
      if (! (*niter)->table.empty()) {
	y += (*niter)->table.sample_y_inv(sampler, discount, strength);
	z += (*niter)->table.sample_z_inv(sampler, discount, strength);
      }
    
    return sampler.beta(parameters[order].discount_alpha + y, parameters[order].discount_beta + z);
  }

  void initialize()
  {
    if (nodes.size() == parameters.size()) return;
    
    nodes.clear();
    nodes.reserve(parameters.size());
    nodes.resize(parameters.size());
    
    nodes[0].push_back(&root);
    
    for (id_type id = 0; id != trie.size(); ++ id) {
      node_type& node = trie[id];
      
      nodes[node.order].push_back(&node);
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
    
    // How to dump mixtures...?
    
    
    
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
    
    for (size_type order = 1; order < parameters.size(); ++ order) {
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
  lm_set_type& models;
  mixture_set_type mixtures;
  
  parameter_set_type parameters;
  
  trie_type trie;
  node_type root;
  node_ptr_map_type nodes;
  
  double p0;
  size_type counts0;
};



typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef utils::sampler<boost::mt19937> sampler_type;

typedef PYPLM::size_type size_type;

// we will precompute <word, node> pair...
typedef utils::simple_vector<PYPLM::id_type, std::allocator<PYPLM::id_type> >  node_set_type;
typedef boost::fusion::tuple<word_type, node_set_type, int> data_type;
typedef std::vector<data_type, std::allocator<data_type> > data_set_type;

typedef std::vector<size_type, std::allocator<size_type> > position_set_type;
typedef std::vector<int, std::allocator<int> > derivation_set_type;

typedef std::string spec_type;
typedef std::vector<spec_type, std::allocator<spec_type> > spec_set_type;

struct Counter
{
  Counter() : counter(0) {}
  
  void increment()
  {
    utils::atomicop::fetch_and_add(counter, size_type(1));
  }
  
  void wait(size_type target)
  {
    for (;;) {
      utils::atomicop::memory_barrier();
      
      for (int i = 0; i < 64; ++ i) {
	const size_type curr = counter;
	
	if (curr == target)
	  return;
	else
	  boost::thread::yield();
      }
      
      struct timespec tm;
      tm.tv_sec = 0;
      tm.tv_nsec = 2000001;
      nanosleep(&tm, NULL);
    }
  }

  void clear() { counter = 0; }
  
  size_type counter;
};
typedef Counter counter_type;

struct Task
{
  typedef PYPLM::size_type       size_type;
  typedef PYPLM::difference_type difference_type;

  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type > > queue_type;

  
  Task(queue_type&   __mapper,
       counter_type& __reducer,
       const data_set_type& __training,
       const PYPMixture::normalizer_set_type& __normalizers,
       derivation_set_type& __derivations,
       PYPMixture& __model,
       const sampler_type& __sampler)
    : mapper(__mapper),
      reducer(__reducer),
      training(__training),
      normalizers_base(__normalizers),
      derivations(__derivations),
      model(__model),
      sampler(__sampler) {}
  
  void operator()()
  {
    try {
      size_type pos;
      
      PYPMixture::normalizer_set_type normalizers(normalizers_base);
    
      for (;;) {
	mapper.pop(pos);
      
	if (pos == size_type(-1)) break;
      
	if (derivations[pos])
	  model.decrement(boost::fusion::get<0>(training[pos]),
			  boost::fusion::get<1>(training[pos]).begin(),
			  normalizers.begin(),
			  sampler);
	else
	  derivations[pos] = true;
      
	model.increment(boost::fusion::get<0>(training[pos]),
			boost::fusion::get<1>(training[pos]).begin(),
			normalizers.begin(),
			sampler,
			temperature);
      
	//reducer.push(pos);
	reducer.increment();
      }
    }
    catch (const std::exception& err) {
      std::cerr << "error: " << err.what() << std::endl;
      throw err;
    }
  }
  
  queue_type&   mapper;
  counter_type& reducer;
  
  const data_set_type& training;
  const PYPMixture::normalizer_set_type& normalizers_base;
  derivation_set_type& derivations;
  
  PYPMixture& model;
  sampler_type sampler;
  
  double temperature;
};


template <typename Training>
struct less_rank
{
  less_rank(const Training& __training) : training(__training) {}

  bool operator()(const size_t& x, const size_t& y) const
  {
    return boost::fusion::get<2>(training[x]) < boost::fusion::get<2>(training[y]);
  }

  const Training& training;
};

path_set_type train_files;
path_set_type test_files;
path_type     output_file;

path_set_type cluster_files;
spec_set_type stemmer_specs;

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

double class_discount_alpha = 1.0;
double class_discount_beta  = 1.0;
double class_strength_shape = 1.0;
double class_strength_rate  = 1.0;

double lambda_discount_alpha = 1.0;
double lambda_discount_beta  = 1.0;
double lambda_strength_shape = 1.0;
double lambda_strength_rate  = 1.0;

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

    if (cluster_files.empty() && stemmer_specs.empty())
      throw std::runtime_error("no stemmers nor clusters");
    
    PYPMixture::normalizer_set_type normalizers(cluster_files.size() + stemmer_specs.size());
    
    for (size_t i = 0; i != cluster_files.size(); ++ i)
      normalizers[i] = normalizer_type(&cluster_type::create(cluster_files[i]));

    for (size_t i = 0; i != stemmer_specs.size(); ++ i)
      normalizers[i + cluster_files.size()] = normalizer_type(&stemmer_type::create(stemmer_specs[i]));

    if (normalizers.size() > 16)
      throw std::runtime_error("we support up to 16 mixtures!");
    
    sampler_type sampler;
    const size_t num_vocab = vocabulary_size(train_files);
    
    PYPMixture::lm_set_type latent(normalizers.size(), PYPLM(order,
							     1.0 / num_vocab,
							     PYPLM::parameter_type(discount_alpha,
										   discount_beta,
										   strength_shape,
										   strength_rate)));
    
    PYPMixture model(latent,
		     order,
		     1.0 / num_vocab,
		     PYPMixture::parameter_type(discount_alpha,
						discount_beta,
						strength_shape,
						strength_rate),
		     PYPMixture::parameter_type(lambda_discount_alpha,
						lambda_discount_beta,
						lambda_strength_shape,
						lambda_strength_rate));
        
    typedef std::vector<bool, std::allocator<bool> > non_oov_type;

    typedef utils::unordered_set<word_type, boost::hash<word_type>, std::equal_to<word_type>, std::allocator<word_type> >::type word_unique_type;
    typedef utils::chunk_vector<word_unique_type, 4096 / sizeof(word_unique_type), std::allocator<word_unique_type> > word_unique_set_type;
    typedef std::vector<size_t, std::allocator<size_t> > index_type;
    
    data_set_type        training;
    word_unique_set_type uniques;
    non_oov_type         non_oov;
    
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

      node_set_type nodes(normalizers.size() + 1);
      
      while (is >> sentence) {
	ngram.resize(1);
	
	sentence_type::const_iterator siter_end = sentence.end();
	for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  sentence_type::const_iterator niter_begin = std::max(ngram.begin(), ngram.end() - order + 1);
	  sentence_type::const_iterator niter_end   = ngram.end();
	  
	  nodes.front() = model.insert(niter_begin, niter_end);
	  for (size_t i = 0; i != normalizers.size(); ++ i)
	    nodes[i + 1] = latent[i].insert(niter_begin, niter_end, normalizers[i]);
	  
	  training.push_back(data_type(*siter, nodes, 0));
	  
#if 0
	  std::cerr << "ngram: " << ngram << std::endl;
	  std::cerr << "node: ";
	  std::copy(nodes.begin(), nodes.end(), std::ostream_iterator<int>(std::cerr, " "));
	  std::cerr << std::endl;
#endif
	  
	  if (boost::fusion::get<1>(training.back()).front() >= uniques.size())
	    uniques.resize(boost::fusion::get<1>(training.back()).front() + 1);
	  uniques[boost::fusion::get<1>(training.back()).front()].insert(boost::fusion::get<0>(training.back()));
	  
	  ngram.push_back(*siter);
	  
	  if (siter->id() >= non_oov.size())
	    non_oov.resize(siter->id() + 1, false);
	  non_oov[siter->id()] = true;
	}
	
	sentence_type::const_iterator niter_begin = std::max(ngram.begin(), ngram.end() - order + 1);
	sentence_type::const_iterator niter_end   = ngram.end();

	nodes.front() = model.insert(niter_begin, niter_end);
	for (size_t i = 0; i != normalizers.size(); ++ i)
	  nodes[i + 1] = latent[i].insert(niter_begin, niter_end, normalizers[i]);
	
	training.push_back(data_type(vocab_type::EOS, nodes, 0));
	
	if (boost::fusion::get<1>(training.back()).front() >= uniques.size())
	  uniques.resize(boost::fusion::get<1>(training.back()).front() + 1);
	uniques[boost::fusion::get<1>(training.back()).front()].insert(boost::fusion::get<0>(training.back()));
      }
    }
    
    if (training.empty())
      throw std::runtime_error("no training data?");
    
    data_set_type(training).swap(training);
    
    // assign rank...
    {
      data_set_type::iterator titer_end = training.end();
      for (data_set_type::iterator titer = training.begin(); titer != titer_end; ++ titer)
	boost::fusion::get<2>(*titer) = uniques[boost::fusion::get<1>(*titer).front()].size();
      
      uniques.clear();
    }
    
    // assign positons and orders
    position_set_type   positions(training.size());
    derivation_set_type derivations(training.size(), 0);
    for (size_type pos = 0; pos != positions.size(); ++ pos)
      positions[pos]= pos;
    
    // sample parameters, first...
    if (slice_sampling)
      model.slice_sample_parameters(sampler, 1.0, resample_iterations);
    else
      model.sample_parameters(sampler, 1.0, resample_iterations);
    
    if (debug >= 2) {
      for (int n = 0; n != order; ++ n)
	std::cerr << "order=" << n << " discount=" << model.parameters[n].discount << " strength=" << model.parameters[n].strength << std::endl;
      
      for (size_t i = 0; i != latent.size(); ++ i)
	for (int n = 0; n != order; ++ n)
	  std::cerr << "latent=" << i << " order=" << n << " discount=" << latent[i].parameters[n].discount << " strength=" << latent[i].parameters[n].strength << std::endl;
      
      for (int n = 0; n != order; ++ n) {
	std::cerr << "mixture order=" << n << " lambda=";
	std::copy(model.mixtures[n].probs.begin(), model.mixtures[n].probs.end(), std::ostream_iterator<double>(std::cerr, " "));
	std::cerr << std::endl;
      }
    }
    
    Task::queue_type queue_mapper;
    counter_type reducer;
    
    std::vector<Task, std::allocator<Task> > tasks(threads, Task(queue_mapper,
								 reducer,
								 training,
								 normalizers,
								 derivations,
								 model,
								 sampler));

    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    size_t anneal_iter = 0;
    const size_t anneal_last = utils::bithack::branch(anneal_steps > 0, anneal_steps, 0);
    
    size_t baby_iter = 0;
    const size_t baby_last = utils::bithack::branch(baby_steps > 0, baby_steps, 0);
    
    size_t burn_iter = 0;
    const size_t burn_last = utils::bithack::branch(burns > 0, burns, 0);
    
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

      // assign temperature...
      for (size_type i = 0; i != tasks.size(); ++ i)
	tasks[i].temperature = temperature;

      // clear mixtures counts...
      for (size_type n = 0; n != model.mixtures.size(); ++ n)
	model.mixtures[n].clear();
      
      // shuffle
      boost::random_number_generator<sampler_type::generator_type> gen(sampler.generator());
      std::random_shuffle(positions.begin(), positions.end(), gen);
      
      if (! baby_finished)
	std::sort(positions.begin(), positions.end(), less_rank<data_set_type>(training));

      reducer.clear();
      
      position_set_type::const_iterator piter_end = positions.end();
      for (position_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter)
	queue_mapper.push(*piter);

      reducer.wait(positions.size());
      
      if (static_cast<int>(iter) % resample_rate == resample_rate - 1) {
	if (slice_sampling)
	  model.slice_sample_parameters(sampler, temperature, resample_iterations);
	else
	  model.sample_parameters(sampler, temperature, resample_iterations);
	
	if (debug >= 2) {
	  for (int n = 0; n != order; ++ n)
	    std::cerr << "order=" << n << " discount=" << model.parameters[n].discount << " strength=" << model.parameters[n].strength << std::endl;

	  for (size_t i = 0; i != latent.size(); ++ i)
	    for (int n = 0; n != order; ++ n)
	      std::cerr << "latent=" << i << " order=" << n << " discount=" << latent[i].parameters[n].discount << " strength=" << latent[i].parameters[n].strength << std::endl;
	  
	  for (int n = 0; n != order; ++ n) {
	    std::cerr << "mixture order=" << n << " lambda=";
	    std::copy(model.mixtures[n].probs.begin(), model.mixtures[n].probs.end(), std::ostream_iterator<double>(std::cerr, " "));
	    std::cerr << std::endl;
	  }
	}
      }
	
      if (debug)
	std::cerr << "log-likelihood: " << model.log_likelihood() << std::endl;
    }
    
    for (int i = 0; i != threads; ++ i)
      queue_mapper.push(size_type(-1));
    
    workers.join_all();

    // clear training data
    training.clear();
    data_set_type(training).swap(training);

    positions.clear();
    position_set_type(positions).swap(positions);

    derivations.clear();
    derivation_set_type(derivations).swap(derivations);
    
    // we will dump LM... now, define a format!
    //if (! output_file.empty())
    //model.write(output_file);
    
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
	    const double prob = model.prob(*siter, std::max(ngram.begin(), ngram.end() - order + 1), ngram.end(), normalizers.begin());
	    const double lp = std::log(prob);
	    
	    if (! is_oov)
	      logprob_total += lp;
	    logprob += lp;
	    
	    num_oov += is_oov;
	    
	    ngram.push_back(*siter);
	  }

	  const double prob = model.prob(vocab_type::EOS, std::max(ngram.begin(), ngram.end() - order + 1), ngram.end(), normalizers.begin());
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

struct TaskVocab
{
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
  
  typedef uint64_t count_type;
  typedef utils::unordered_map<std::string, count_type, boost::hash<utils::piece>, std::equal_to<std::string>, std::allocator<std::pair<const std::string, count_type> > >::type vocab_type;
  
  typedef std::vector<const vocab_type::value_type*, std::allocator<const vocab_type::value_type*> > sorted_type;

  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type > > queue_type;
  
  TaskVocab(queue_type& __queue)
    : queue(__queue) {}

  void operator()()
  {
    path_type path;
    std::string word;
    
    for (;;) {
      queue.pop_swap(path);
      
      if (path.empty()) break;
      
      utils::compress_istream is(path, 1024 * 1024);
      
      while (is >> word) 
	++ vocab[word];
    }
  }
  
  queue_type& queue;
  vocab_type  vocab;
};

size_t vocabulary_size(const path_set_type& files)
{
  typedef TaskVocab task_type;
  
  task_type::queue_type queue;

  std::vector<task_type, std::allocator<task_type> > tasks(threads, task_type(queue));
  
  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    if (! boost::filesystem::exists(*fiter))
      throw std::runtime_error("no file? " + fiter->string());
    
    queue.push(*fiter);
  }

  for (int i = 0; i != threads; ++ i)
    queue.push(path_type());
  
  workers.join_all();
  
  task_type::vocab_type vocab;
  for (int i = 0; i != threads; ++ i) {
    if (vocab.empty())
      vocab.swap(tasks[i].vocab);
    else {
      task_type::vocab_type::const_iterator viter_end = tasks[i].vocab.end();
      for (task_type::vocab_type::const_iterator viter = tasks[i].vocab.begin(); viter != viter_end; ++ viter)
	vocab[viter->first] += viter->second;
    }
    
    tasks[i].vocab.clear();
  }
  
  task_type::sorted_type sorted;
  sorted.reserve(vocab.size());
  
  task_type::vocab_type::const_iterator viter_end = vocab.end();
  for (task_type::vocab_type::const_iterator viter = vocab.begin(); viter != viter_end; ++ viter)
    sorted.push_back(&(*viter));
  
  std::sort(sorted.begin(), sorted.end(), task_type::greater_second<task_type::vocab_type::value_type>());

  task_type::sorted_type::const_iterator siter_end = sorted.end();
  for (task_type::sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter) {
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

    ("cluster", po::value<path_set_type>(&cluster_files)->multitoken(), "cluster file(s)")
    ("stemmer", po::value<spec_set_type>(&stemmer_specs)->multitoken(), "stemmer specifications")
    
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
    
    ("class-discount-alpha", po::value<double>(&class_discount_alpha)->default_value(class_discount_alpha), "discount ~ Beta(alpha,beta)")
    ("class-discount-beta",  po::value<double>(&class_discount_beta)->default_value(class_discount_beta),   "discount ~ Beta(alpha,beta)")

    ("class-strength-shape", po::value<double>(&class_strength_shape)->default_value(class_strength_shape), "strength ~ Gamma(shape,rate)")
    ("class-strength-rate",  po::value<double>(&class_strength_rate)->default_value(class_strength_rate),   "strength ~ Gamma(shape,rate)")

    ("lambda-discount-alpha", po::value<double>(&lambda_discount_alpha)->default_value(lambda_discount_alpha), "discount ~ Beta(alpha,beta)")
    ("lambda-discount-beta",  po::value<double>(&lambda_discount_beta)->default_value(lambda_discount_beta),   "discount ~ Beta(alpha,beta)")

    ("lambda-strength-shape", po::value<double>(&lambda_strength_shape)->default_value(lambda_strength_shape), "strength ~ Gamma(shape,rate)")
    ("lambda-strength-rate",  po::value<double>(&lambda_strength_rate)->default_value(lambda_strength_rate),   "strength ~ Gamma(shape,rate)")

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

