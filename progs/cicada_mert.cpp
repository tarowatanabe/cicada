
//
// refset format:
// 0 |||  reference translatin for source sentence 0
// 0 |||  another reference
// 1 |||  reference translation for source sentence 1
//

#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <stdexcept>
#include <numeric>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"
#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"
#include "cicada/viterbi.hpp"

#include "cicada/eval.hpp"

#include "cicada/optimize/line_search.hpp"
#include "cicada/optimize/powell.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"

#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/thread.hpp>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Symbol   symbol_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef rule_type::feature_set_type    feature_set_type;
typedef cicada::WeightVector<double>   weight_set_type;
typedef feature_set_type::feature_type feature_type;

typedef std::vector<weight_set_type, std::allocator<weight_set_type> > weight_set_collection_type;

typedef std::vector<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

path_set_type tstset_files;
path_type     refset_file = "-";
path_type     output_file = "-";

path_type bound_lower_file;
path_type bound_upper_file;


path_set_type feature_weights_files;

std::string scorer_name = "bleu";
bool scorer_list = false;

int iteration = 10;
int samples_restarts   = 4;
int samples_directions = 10;

bool initial_average = false;

double tolerance = 1e-4;

bool regularize_l1 = false;
bool regularize_l2 = false;
double regularize_scale = 1.0; // inverse of C == 1.0 / C : where C is a constant of SVM^{light}

bool weight_normalize_l1 = false;
bool weight_normalize_l2 = false;

int threads = 4;

int debug = 0;

template <typename Iterator>
inline
void randomize(Iterator first, Iterator last, Iterator lower, Iterator upper)
{
  for (/**/; first != last; ++ first, ++ lower, ++ upper)
    *first = *lower + (double(random()) / RAND_MAX) * std::min(double(*upper - *lower), 1.0);
}

template <typename Iterator>
inline
void normalize_l2(Iterator first, Iterator last, const double radius)
{
  const double sum = std::inner_product(first, last, first, 0.0);
  
  if (sum != 0.0)
    std::transform(first, last, first, std::bind2nd(std::multiplies<double>(), radius / std::sqrt(sum)));
}

template <typename Iterator>
void normalize_l1(Iterator first, Iterator last, const double radius)
{
  double sum = 0.0;
  for (Iterator iter = first; iter != last; ++ iter)
    sum += std::fabs(*iter);
  
  if (sum != 0.0)
    std::transform(first, last, first, std::bind2nd(std::multiplies<double>(), radius / std::sqrt(sum)));
}

template <typename Iterator, typename BoundIterator>
inline
bool valid_bounds(Iterator first, Iterator last, BoundIterator lower, BoundIterator upper)
{
  for (/**/; first != last; ++ first, ++ lower, ++ upper)
    if (*first < *lower || *upper < *first)
      return false;
  return true;
}

void read_tstset(const path_set_type& files, hypergraph_set_type& graphs);
void read_refset(const path_type& file, scorer_document_type& scorers);

void options(int argc, char** argv);

struct EnvelopeComputer
{
  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;

  EnvelopeComputer(const scorer_document_type& __scorers,
		   const hypergraph_set_type&  __graphs)
    : scorers(__scorers),
      graphs(__graphs) {}

  void operator()(segment_document_type& segments, const weight_set_type& origin, const weight_set_type& direction) const;

  const scorer_document_type& scorers;
  const hypergraph_set_type&  graphs;
};

struct ViterbiComputer
{
  ViterbiComputer(const scorer_document_type& __scorers,
		  const hypergraph_set_type&  __graphs)
    : scorers(__scorers),
      graphs(__graphs) {}
  
  double operator()(const weight_set_type& weights) const;

  const scorer_document_type& scorers;
  const hypergraph_set_type& graphs;
};

template <typename Regularizer>
bool powell(const scorer_document_type& scorers,
	    const hypergraph_set_type& graphs,
	    const weight_set_type& bound_lower,
	    const weight_set_type& bound_upper,
	    Regularizer regularizer,
	    const double tolerance,
	    const int samples,
	    double& score,
	    weight_set_type& weights)
{
  cicada::optimize::Powell<EnvelopeComputer, ViterbiComputer, Regularizer> optimizer(EnvelopeComputer(scorers, graphs),
										     ViterbiComputer(scorers, graphs),
										     regularizer,
										     bound_lower,
										     bound_upper,
										     tolerance,
										     samples,
										     scorers.error_metric(),
										     debug);
  
  return optimizer(score, weights);
}

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);

    if (scorer_list) {
      std::cout << cicada::eval::Scorer::lists();
      return 0;
    }
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("you cannot use both of L1 and L2...");
    
    if (regularize_l1 || regularize_l2) {
      weight_normalize_l1 = false;
      weight_normalize_l2 = false;
      
      if (regularize_scale <= 0.0)
	throw std::runtime_error("the scaling for L1/L2 must be positive");
    }
    
    if (weight_normalize_l1 && weight_normalize_l2)
      throw std::runtime_error("you cannot use both of L1 and L2 for weight normalization...");

    threads = utils::bithack::max(threads, 1);
    
    // read reference set
    scorer_document_type scorers(scorer_name);
    
    read_refset(refset_file, scorers);

    if (debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;

    // read test set
    if (tstset_files.empty())
      tstset_files.push_back("-");

    if (debug)
      std::cerr << "reading hypergraphs" << std::endl;

    hypergraph_set_type graphs(scorers.size());
    
    read_tstset(tstset_files, graphs);

    // collect initial weights
    weight_set_collection_type weights;
    
    if (! feature_weights_files.empty()) {

      for (path_set_type::const_iterator fiter = feature_weights_files.begin(); fiter != feature_weights_files.end(); ++ fiter) {
	if (*fiter != "-" && ! boost::filesystem::exists(*fiter))
	  throw std::runtime_error("no file? " + fiter->file_string());
	
	utils::compress_istream is(*fiter);
	
	weights.push_back(weight_set_type());
	is >> weights.back();
      }
      
      if (initial_average && weights.size() > 1) {
	weight_set_type weight;
	
	weight_set_collection_type::const_iterator witer_end = weights.end();
	for (weight_set_collection_type::const_iterator witer = weights.begin(); witer != witer_end; ++ witer)
	  weight += *witer;
	
	weight *= (1.0 / weights.size());
	
	weights.push_back(weight);
      }
      
      std::set<weight_set_type, std::less<weight_set_type>, std::allocator<weight_set_type> > uniques;
      uniques.insert(weights.begin(), weights.end());
      weights.clear();
      weights.insert(weights.end(), uniques.begin(), uniques.end());
    } else {
      weights.push_back(weight_set_type());
      
      // all one weight...
      for (feature_type::id_type id = 0; id < feature_type::allocated(); ++ id)
	if (! feature_type(id).empty())
	  weights.back()[feature_type(id)] = 1.0;
    }
    
    // collect lower/upper bounds
    weight_set_type bound_lower;
    weight_set_type bound_upper;
    
    if (! bound_lower_file.empty())
      if (bound_lower_file == "-" || boost::filesystem::exists(bound_lower_file)) {
	utils::compress_istream is(bound_lower_file);
	is >> bound_lower;
      } else
	throw std::runtime_error("no lower-bound file?" + bound_lower_file.file_string());
    
    if (! bound_upper_file.empty())
      if (bound_upper_file == "-" || boost::filesystem::exists(bound_upper_file)) {
	utils::compress_istream is(bound_upper_file);
	is >> bound_upper;
      } else
	throw std::runtime_error("no upper-bound file?" + bound_upper_file.file_string());


    cicada::optimize::LineSearch::initialize_bound(bound_lower, bound_upper);
    
    double          optimum_objective = std::numeric_limits<double>::infinity();
    weight_set_type optimum_weights;
    
    if (debug)
      std::cerr << "start optimization" << std::endl;
    
    int sample = 0;
    for (weight_set_collection_type::const_iterator witer = weights.begin(); witer != weights.end(); ++ witer, ++ sample) {
      typedef cicada::optimize::LineSearch line_search_type;

      double          sample_objective = std::numeric_limits<double>::infinity();
      weight_set_type sample_weights = *witer;
      
      utils::resource opt_start;
      
      bool moved = false;
      if (regularize_l1)
	moved = powell(scorers,
		       graphs,
		       bound_lower,
		       bound_upper,
		       line_search_type::RegularizeL1(regularize_scale),
		       tolerance,
		       samples_directions,
		       sample_objective,
		       sample_weights);
      else if (regularize_l2)
	moved = powell(scorers,
		       graphs,
		       bound_lower,
		       bound_upper,
		       line_search_type::RegularizeL2(regularize_scale),
		       tolerance,
		       samples_directions,
		       sample_objective,
		       sample_weights);
      else
	moved = powell(scorers,
		       graphs,
		       bound_lower,
		       bound_upper,
		       line_search_type::RegularizeNone(regularize_scale),
		       tolerance,
		       samples_directions,
		       sample_objective,
		       sample_weights);
      
      utils::resource opt_end;
      
      if (debug)
	std::cerr << "cpu time: " << (opt_end.cpu_time() - opt_start.cpu_time()) << '\n'
		  << "user time: " << (opt_end.user_time() - opt_start.user_time()) << '\n';
      
      if (debug)
	std::cerr << "sample: " << (sample + 1) << " objective: " << sample_objective << std::endl
		  << sample_weights;
      
      if (moved && sample_objective < optimum_objective) {
	optimum_objective = sample_objective;
	optimum_weights = sample_weights;
      }
    }
    
    for (/**/; sample < samples_restarts + weights.size(); ++ sample) {
      typedef cicada::optimize::LineSearch line_search_type;
      
      double          sample_objective = std::numeric_limits<double>::infinity();
      weight_set_type sample_weights = weights.back();
      
      if (sample > 0) {
	// perform randomize...
	sample_weights = optimum_weights;
	
	while (1) {
	  randomize(sample_weights.begin(), sample_weights.end(), bound_lower.begin(), bound_upper.begin());
	  
	  if (weight_normalize_l1)
	    normalize_l1(sample_weights.begin(), sample_weights.end(), 1.0);
	  else
	    normalize_l2(sample_weights.begin(), sample_weights.end(), 1.0);
	  
	  if (valid_bounds(sample_weights.begin(), sample_weights.end(), bound_lower.begin(), bound_upper.begin()))
	    break;
	}
      }
      
      utils::resource opt_start;
      
      bool moved = false;
      if (regularize_l1)
	moved = powell(scorers,
		       graphs,
		       bound_lower,
		       bound_upper,
		       line_search_type::RegularizeL1(regularize_scale),
		       tolerance,
		       samples_directions,
		       sample_objective,
		       sample_weights);
      else if (regularize_l2)
	moved = powell(scorers,
		       graphs,
		       bound_lower,
		       bound_upper,
		       line_search_type::RegularizeL2(regularize_scale),
		       tolerance,
		       samples_directions,
		       sample_objective,
		       sample_weights);
      else
	moved = powell(scorers,
		       graphs,
		       bound_lower,
		       bound_upper,
		       line_search_type::RegularizeNone(regularize_scale),
		       tolerance,
		       samples_directions,
		       sample_objective,
		       sample_weights);
      
      utils::resource opt_end;
      
      if (debug)
	std::cerr << "cpu time: " << (opt_end.cpu_time() - opt_start.cpu_time()) << '\n'
		  << "user time: " << (opt_end.user_time() - opt_start.user_time()) << '\n';
      
      if (debug)
	std::cerr << "sample: " << (sample + 1) << " objective: " << sample_objective << std::endl
		  << sample_weights;
      
      if (moved && sample_objective < optimum_objective) {
	optimum_objective = sample_objective;
	optimum_weights = sample_weights;
      }
    }

    if (debug)
      std::cerr << "objective: " << optimum_objective << std::endl;
    
    if (weight_normalize_l1)
      normalize_l1(optimum_weights.begin(), optimum_weights.end(), std::sqrt(feature_type::allocated()));
    else if (weight_normalize_l2)
      normalize_l2(optimum_weights.begin(), optimum_weights.end(), std::sqrt(feature_type::allocated()));
    
    utils::compress_ostream os(output_file);
    os.precision(10);
    os << optimum_weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


struct EnvelopeTask
{
  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;

  typedef cicada::semiring::Envelope envelope_type;
  typedef std::vector<envelope_type, std::allocator<envelope_type> >  envelope_set_type;

  typedef utils::lockfree_list_queue<int, std::allocator<int> >  queue_type;

    
  EnvelopeTask(queue_type& __queue,
	       segment_document_type&      __segments,
	       const weight_set_type&      __origin,
	       const weight_set_type&      __direction,
	       const scorer_document_type& __scorers,
	       const hypergraph_set_type&  __graphs)
    : queue(__queue),
      segments(__segments),
      origin(__origin),
      direction(__direction),
      scorers(__scorers),
      graphs(__graphs) {}

  void operator()()
  {
    envelope_set_type envelopes;
    sentence_type     yield;

    int seg;
    
    while (1) {
      queue.pop(seg);
      if (seg < 0) break;

      
      
      envelopes.clear();
      envelopes.resize(graphs[seg].nodes.size());

      cicada::inside(graphs[seg], envelopes, cicada::semiring::EnvelopeFunction<weight_set_type>(origin, direction));
      
      const envelope_type& envelope = envelopes[graphs[seg].goal];
      const_cast<envelope_type&>(envelope).sort();
      
      envelope_type::const_iterator eiter_end = envelope.end();
      for (envelope_type::const_iterator eiter = envelope.begin(); eiter != eiter_end; ++ eiter) {
	const envelope_type::line_ptr_type& line = *eiter;
	
	line->yield(yield);
	
	scorer_type::score_ptr_type score = scorers[seg]->score(yield);
	
	if (debug >= 4)
	  std::cerr << "segment: " << seg << " x: " << line->x << std::endl;
	
	segments[seg].push_back(std::make_pair(line->x, score));
      }
    }
  }
  
  queue_type& queue;

  segment_document_type& segments;

  const weight_set_type& origin;
  const weight_set_type& direction;
  
  const scorer_document_type& scorers;
  const hypergraph_set_type&  graphs;
};

void EnvelopeComputer::operator()(segment_document_type& segments, const weight_set_type& origin, const weight_set_type& direction) const
{
  typedef EnvelopeTask task_type;
  typedef task_type::queue_type queue_type;

  segments.clear();
  segments.resize(graphs.size());
  
  queue_type queue(graphs.size());

  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, segments, origin, direction, scorers, graphs)));
  
  for (int seg = 0; seg < graphs.size(); ++ seg)
    if (graphs[seg].goal != hypergraph_type::invalid)
      queue.push(seg);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
}


struct ViterbiTask
{
  typedef cicada::semiring::Logprob<double> weight_type;
  
  struct viterbi_function
  {
    typedef rule_type::feature_set_type feature_set_type;

    typedef weight_type value_type;

    viterbi_function(const weight_set_type& __weights)
      : weights(__weights) {}

    const weight_set_type& weights;
  
    template <typename Edge>
    value_type operator()(const Edge& edge) const
    {
      return cicada::semiring::traits<value_type>::log(edge.features.dot(weights));
    }
  };

  struct viterbi_traversal
  {
    typedef sentence_type value_type;
  
    template <typename Edge, typename Iterator>
    void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
    {
      // extract target-yield, features
    
      yield.clear();
    
      rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
      for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
	if (titer->is_non_terminal()) {
	  const int pos = titer->non_terminal_index() - 1;
	  yield.insert(yield.end(), (*(*(first + pos))).begin(), (*(*(first + pos))).end());
	} else if (*titer != vocab_type::EPSILON)
	  yield.push_back(*titer);
    }
  };

  
  typedef std::pair<int, const hypergraph_type*> mapper_type;
  typedef std::pair<int, sentence_type>     reducer_type;

  typedef utils::lockfree_list_queue<mapper_type, std::allocator<mapper_type> >  queue_mapper_type;
  typedef utils::lockfree_list_queue<reducer_type, std::allocator<reducer_type> > queue_reducer_type;
  
  ViterbiTask(const weight_set_type& __weights,
	      queue_mapper_type& __queue_mapper,
	      queue_reducer_type& __queue_reducer)
    : weights(__weights),
      queue_mapper(__queue_mapper),
      queue_reducer(__queue_reducer) {}

  void operator()()
  {
    mapper_type  mapped;
    reducer_type reduced;
    
    weight_type  weight;
    
    while (1) {
      queue_mapper.pop(mapped);
      if (mapped.first < 0) break;
      
      if (! mapped.second)
	throw std::runtime_error("no graph?");
      
      reduced.first = mapped.first;
      
      viterbi(*mapped.second, reduced.second, weight, viterbi_traversal(), viterbi_function(weights));
      
      queue_reducer.push(reduced);
    }
  }
  
  const weight_set_type& weights;
  
  queue_mapper_type& queue_mapper;
  queue_reducer_type& queue_reducer;
};

double ViterbiComputer::operator()(const weight_set_type& weights) const
{
  typedef ViterbiTask task_type;
  
  typedef task_type::queue_mapper_type queue_mapper_type;
  typedef task_type::queue_reducer_type queue_reducer_type;
  
  typedef task_type::mapper_type  mapper_type;
  typedef task_type::reducer_type reducer_type;
  
  queue_mapper_type  queue_mapper(graphs.size());
  queue_reducer_type queue_reducer;

  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(weights, queue_mapper, queue_reducer)));

  int mapped_size = 0;
  for (int seg = 0; seg < graphs.size(); ++ seg)
    if (graphs[seg].goal != hypergraph_type::invalid) {
      queue_mapper.push(mapper_type(seg, &graphs[seg]));
      ++ mapped_size;
    }
  
  // termination...
  for (int i = 0; i < threads; ++ i)
    queue_mapper.push(mapper_type(-1, 0));

  scorer_type::score_ptr_type score;
  
  reducer_type reduced;
  for (int seg = 0; seg < mapped_size; ++ seg) {
    queue_reducer.pop(reduced);
    if (! score)
      score = scorers[reduced.first]->score(reduced.second);
    else
      *score += *scorers[reduced.first]->score(reduced.second);
  }
  
  workers.join_all();

  const double score_factor = (scorers.error_metric() ? 1.0 : - 1.0);
  
  return score->score().first * score_factor;
}

void read_tstset(const path_set_type& files, hypergraph_set_type& graphs)
{
  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (debug)
      std::cerr << "file: " << *titer << std::endl;
      
    if (boost::filesystem::is_directory(*titer)) {

      for (int i = 0; /**/; ++ i) {
	const path_type path = (*titer) / (boost::lexical_cast<std::string>(i) + ".gz");

	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	
	int id;
	std::string sep;
	hypergraph_type hypergraph;
      
	weight_set_type origin;
	weight_set_type direction;
      
	while (is >> id >> sep >> hypergraph) {
	
	  if (sep != "|||")
	    throw std::runtime_error("format error?");
	
	  if (id >= graphs.size())
	    throw std::runtime_error("tstset size exceeds refset size?" + boost::lexical_cast<std::string>(id));
	
	  graphs[id].unite(hypergraph);
	}
      }
    } else {
      
      utils::compress_istream is(*titer, 1024 * 1024);
      
      int id;
      std::string sep;
      hypergraph_type hypergraph;
      
      weight_set_type origin;
      weight_set_type direction;
      
      while (is >> id >> sep >> hypergraph) {
	
	if (sep != "|||")
	  throw std::runtime_error("format error?");
	
	if (id >= graphs.size())
	  throw std::runtime_error("tstset size exceeds refset size?" + boost::lexical_cast<std::string>(id));
	
	graphs[id].unite(hypergraph);
      }
    }
  }
  
  for (int id = 0; id < graphs.size(); ++ id)
    if (graphs[id].goal == hypergraph_type::invalid)
      std::cerr << "invalid graph at: " << id << std::endl;
}

void read_refset(const path_type& file, scorer_document_type& scorers)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;
    
  scorers.clear();
  
  utils::compress_istream is(file);

  std::string line;
  
  while (std::getline(is, line)) {
    tokenizer_type tokenizer(line);
    
    tokenizer_type::iterator iter = tokenizer.begin();
    if (iter == tokenizer.end()) continue;
    
    const int id = boost::lexical_cast<int>(*iter);
    ++ iter;
    
    if (iter == tokenizer.end()) continue;
    if (*iter != "|||") continue;
    ++ iter;
    
    if (id >= scorers.size())
      scorers.resize(id + 1);
    
    if (! scorers[id])
      scorers[id] = scorers.create();
    
    scorers[id]->insert(sentence_type(iter, tokenizer.end()));
  }
  
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(),          "test set file(s) (in hypergraph format)")
    ("refset",  po::value<path_type>(&refset_file)->default_value(refset_file), "reference set file")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("bound-lower", po::value<path_type>(&bound_lower_file),                    "lower bounds definition for feature weights")
    ("bound-upper", po::value<path_type>(&bound_upper_file),                    "upper bounds definition for feature weights")
    
    // feature weight files
    ("feature-weights",  po::value<path_set_type>(&feature_weights_files), "feature weights file(s)")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")
    
    ("iteration",          po::value<int>(&iteration),          "# of mert iteration")
    ("samples-restarts",   po::value<int>(&samples_restarts),   "# of random sampling for initial starting point")
    ("samples-directions", po::value<int>(&samples_directions), "# of ramdom sampling for directions")
    ("initial-average",    po::bool_switch(&initial_average),   "averaged initial parameters")
    
    ("tolerance", po::value<double>(&tolerance)->default_value(tolerance), "tolerance")
    
    ("regularize-l1",    po::bool_switch(&regularize_l1),      "regularization via L1")
    ("regularize-l2",    po::bool_switch(&regularize_l2),      "regularization via L2")
    ("regularize-scale", po::value<double>(&regularize_scale)->default_value(regularize_scale), "scaling for regularizer")
    
    ("normalize-l1",    po::bool_switch(&weight_normalize_l1), "weight normalization via L1 (not a regularizer...)")
    ("normalize-l2",    po::bool_switch(&weight_normalize_l2), "weight normalization via L2 (not a regularizer...)")


    ("threads", po::value<int>(&threads), "# of threads")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config", po::value<path_type>(), "configuration file")
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  if (variables.count("config")) {
    utils::compress_istream is(variables["config"].as<path_type>());
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
