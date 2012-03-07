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

#include <map>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/chinese_restaurant_process.hpp"
#include "utils/unordered_map.hpp"
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
    typedef utils::chinese_restaurant_process<word_type, boost::hash<word_type>, std::equal_to<word_type>,
					      std::allocator<word_type > > table_type;
  
    Node() : table(), order(0)  {}
  
    table_type table;
    int order;
  };
  typedef Node node_type;
  
  typedef utils::compact_trie_dense<word_type, node_type, boost::hash<word_type>, std::equal_to<word_type>,
				    std::allocator<std::pair<const word_type, node_type> > > trie_type;

  typedef std::vector<double, std::allocator<double> > parameter_set_type;

  PYPLM(const int order,
	const double __p0,
	const double __discount_alpha,
	const double __discount_beta,
	const double __strength_shape,
	const double __strength_rate)
    : trie(word_type()),
      discount(order, 0.9),
      strength(order, 1.0),
      discount_alpha(__discount_alpha),
      discount_beta(__discount_beta),
      strength_shape(__strength_shape),
      strength_rate(__strength_rate),
      p0(__p0),
      counts0(0)
  {
    // unitialize root table...
    root.order = 0;
    root.table = node_type::table_type(discount[0],
				       strength[0],
				       discount_alpha,
				       discount_beta,
				       strength_shape,
				       strength_rate);
  }

  template <typename Iterator, typename Sampler>
  void increment(const word_type& word, Iterator first, Iterator last, Sampler& sampler)
  {
    const double backoff = prob(word, first + 1, last);

    if (first == last) {
      if (root.table.increment(word, backoff, sampler))
	++ counts0;
    } else {
      typedef std::reverse_iterator<Iterator> reverse_iterator;
      
      reverse_iterator begin(last);
      reverse_iterator end(first);
      
      int order = 1;
      id_type node = trie.root();
      for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
	node = trie.insert(node, *iter);
	
	if (! trie[node].order) {
	  trie[node].order = order;
	  trie[node].table = node_type::table_type(discount[order],
						   strength[order],
						   discount_alpha,
						   discount_beta,
						   strength_shape,
						   strength_rate);
	}
      }
      
      // we will also increment lower-order!
      if (trie[node].table.increment(word, backoff, sampler))
	increment(word, first + 1, last, sampler);
    }
  }
  
  template <typename Iterator, typename Sampler>
  void decrement(const word_type& word, Iterator first, Iterator last, Sampler& sampler)
  {
    typedef std::reverse_iterator<Iterator> reverse_iterator;
    
    if (first == last) {
      if (root.table.decrement(word, sampler))
	-- counts0;
    } else {
      const id_type node = trie.find(reverse_iterator(last), reverse_iterator(first));
      
      if (node == trie_type::npos())
	throw std::runtime_error("word is not in the model...?");
      
      // we will also decrement lower-order!
      if (trie[node].table.decrement(word, sampler))
	decrement(word, first + 1, last, sampler);
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
      
      if (node == trie_type::npos())
	return p;
      else
	p = trie[node].table.prob(word, p);
    }
    
    return p;
  }

  double log_likelihood() const
  {
    double logprob = std::log(p0) * counts0;
    
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
      return logprob + root.table.log_likelihood(discount, strength);
    else {
      for (id_type i = 0; i != trie.size(); ++ i)
	if (trie[i].order == order)
	  logprob += trie[i].table.log_likelihood(discount, strength);
      return logprob;
    }
  }
  
  struct DiscountResampler
  {
    DiscountResampler(const PYPLM& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPLM& pyplm;
    int order;
    
    double operator()(const double& proposed_discount) const
    {
      return pyplm.log_likelihood(order, proposed_discount, pyplm.strength[order]);
    }
  };
  
  struct StrengthResampler
  {
    StrengthResampler(const PYPLM& __pyplm, const int __order) : pyplm(__pyplm), order(__order) {}
    
    const PYPLM& pyplm;
    int order;
    
    double operator()(const double& proposed_strength) const
    {
      return pyplm.log_likelihood(order, pyplm.discount[order], proposed_strength);
    }
  };
  
  template <typename Sampler>
  void resample_hyperparameters(Sampler& sampler,
				const size_type num_loop = 5,
				const size_type num_iterations = 10)
  {
    for (int order = discount.size() - 1; order >= 0; -- order) {
      DiscountResampler discount_resampler(*this, order);
      StrengthResampler strength_resampler(*this, order);
      
      for (size_type iter = 0; iter < num_loop; ++ iter) {
	strength[order] = utils::slice_sampler(strength_resampler,
					       strength[order],
					       sampler,
					       - discount[order] + std::numeric_limits<double>::min(),
					       std::numeric_limits<double>::infinity(),
					       0.0,
					       num_iterations,
					       100 * num_iterations);
	
	discount[order] = utils::slice_sampler(discount_resampler,
					       discount[order],
					       sampler,
					       (strength[order] < 0.0 ? - strength[order] : 0.0) + std::numeric_limits<double>::min(),
					       1.0,
					       0.0,
					       num_iterations,
					       100 * num_iterations);
      }
      
      strength[order] = utils::slice_sampler(strength_resampler,
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
      } else {
	for (id_type i = 0; i != trie.size(); ++ i)
	  if (trie[i].order == order) {
	    trie[i].table.discount() = discount[order];
	    trie[i].table.strength() = strength[order];
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
    rep["discount-alpha"] = boost::lexical_cast<std::string>(discount_alpha);
    rep["discount-beta"]  = boost::lexical_cast<std::string>(discount_beta);
    rep["strength-shape"] = boost::lexical_cast<std::string>(strength_shape);
    rep["strength-rate"]  = boost::lexical_cast<std::string>(strength_rate);
      
    rep["p0"]      = boost::lexical_cast<std::string>(p0);
    rep["counts0"] = boost::lexical_cast<std::string>(counts0);
    for (size_type order = 0; order != discount.size(); ++ order) {
      rep["discount" + boost::lexical_cast<std::string>(order)] = boost::lexical_cast<std::string>(discount[order]);
      rep["strength" + boost::lexical_cast<std::string>(order)] = boost::lexical_cast<std::string>(strength[order]);
    }

    // we will compute on-memory for faster indexing... (and assuming small data!)
    
    boost::iostreams::filtering_ostream os_index;
    boost::iostreams::filtering_ostream os_count;
    
    os_index.push(utils::packed_sink<word_type::id_type, std::allocator<word_type::id_type> >(rep.path("index")));
    os_count.push(utils::packed_sink<count_type, std::allocator<count_type> >(rep.path("count")));
      
    position_set_type positions;
    offset_set_type   offsets(1, 0);
    count_type        offset = 0;
    
    node_set_type nodes;
    node_set_type nodes_next;
    
    word_set_type words;
    
    // unigram!
    {
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
	if (*niter == trie_type::npos())
	  positions.set(positions.size(), false); // we will set bit!
	else {
	  const node_type& node = trie[*niter];
	  trie_type::const_iterator iter_begin = trie.begin(*niter);
	  trie_type::const_iterator iter_end   = trie.end(*niter);
	  
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

private: 
  trie_type trie;
  node_type root;
  
  parameter_set_type discount;
  parameter_set_type strength;
  
  double discount_alpha;
  double discount_beta;
  double strength_shape;
  double strength_rate;
  
  double    p0;
  size_type counts0;
};


typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef utils::sampler<boost::mt19937> sampler_type;

path_set_type train_files;
path_type     output_file;

int order = 3;
int samples = 300;
int resample_rate = 20;

double discount_prior_alpha = 1.0;
double discount_prior_beta  = 1.0;
double strength_prior_shape = 1.0;
double strength_prior_rate  = 1.0;

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
    
    if (samples <= 0)
      throw std::runtime_error("# of samples must be positive");

    if (train_files.empty())
      throw std::runtime_error("no training data?");

    sampler_type sampler;
    const size_t num_vocab = vocabulary_size(train_files);
    
    PYPLM lm(order,
	     1.0 / num_vocab,
	     discount_prior_alpha,
	     discount_prior_beta,
	     strength_prior_shape,
	     strength_prior_rate);
    
    for (int iter = 0; iter < samples; ++ iter) {
      if (debug)
	std::cerr << "iteration: " << iter << std::endl;
      
      sentence_type sentence;
      sentence_type ngram(order - 1, vocab_type::BOS);
      
      for (path_set_type::const_iterator fiter = train_files.begin(); fiter != train_files.end(); ++ fiter) {
	utils::compress_istream is(*fiter, 1024 * 1024);

	while (is >> sentence) {
	  ngram.resize(order - 1);
	  
	  sentence_type::const_iterator siter_end = sentence.end();
	  for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	    
	    if (iter)
	      lm.decrement(*siter, ngram.end() - order + 1, ngram.end(), sampler);
	    
	    lm.increment(*siter, ngram.end() - order + 1, ngram.end(), sampler);
	    
	    ngram.push_back(*siter);
	  }
	  
	  if (iter)
	    lm.decrement(vocab_type::EOS, ngram.end() - order + 1, ngram.end(), sampler);
	  
	  lm.increment(vocab_type::EOS, ngram.end() - order + 1, ngram.end(), sampler);
	}
      }
      
      if (iter % resample_rate == resample_rate - 1)
	lm.resample_hyperparameters(sampler);
      
      if (debug)
	std::cerr << "log-likelihood: " << lm.log_likelihood() << std::endl;
    }
    
    // we will dump LM... now, define a format!
    if (! output_file.empty())
      lm.write(output_file);
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
    ("output", po::value<path_type>(&output_file), "output file")
    
    ("order", po::value<int>(&order)->default_value(order), "max ngram order")
    
    ("samples",  po::value<int>(&samples)->default_value(samples),             "# of samples")
    ("resample", po::value<int>(&resample_rate)->default_value(resample_rate), "resampling rate")

    ("discount-alpha", po::value<double>(&discount_prior_alpha)->default_value(discount_prior_alpha), "discount ~ Beta(alpha,beta)")
    ("discount-beta",  po::value<double>(&discount_prior_beta)->default_value(discount_prior_beta),   "discount ~ Beta(alpha,beta)")
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

