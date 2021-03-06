//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

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
#include <algorithm>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"
#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"
#include "cicada/viterbi.hpp"
#include "cicada/sentence_vector.hpp"
#include "cicada/prune.hpp"

#include "cicada/operation/functional.hpp"
#include "cicada/operation/traversal.hpp"

#include "cicada/apply.hpp"
#include "cicada/model.hpp"

#include "cicada/feature/scorer.hpp"
#include "cicada/parameter.hpp"

#include "cicada/eval.hpp"

#include "cicada_output_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/random_seed.hpp"
#include "utils/filesystem.hpp"
#include "utils/getline.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>
#include <boost/thread.hpp>

#include "cicada_impl.hpp"
#include "cicada_text_impl.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Symbol   symbol_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef cicada::Model model_type;

typedef hypergraph_type::feature_set_type    feature_set_type;
typedef cicada::WeightVector<double>   weight_set_type;
typedef feature_set_type::feature_type feature_type;

typedef cicada::FeatureFunction feature_function_type;
typedef feature_function_type::feature_function_ptr_type feature_function_ptr_type;

typedef std::vector<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;


typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type> > feature_function_ptr_set_type;

typedef cicada::SentenceVector sentence_set_type;
typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type  score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;


path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";

bool forest_mode = false;
bool directory_mode = false;

std::string scorer_name = "bleu:order=4,exact=true";
int scorer_cube = 200;
double scorer_beam = 1e-5;
int max_iteration = 10;
int min_iteration = 5;
bool apply_exact = false;

int threads = 2;

int debug = 0;

void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 feature_function_ptr_set_type& features);
void read_refset(const path_set_type& file,
		 scorer_document_type& scorers,
		 sentence_document_type& sentences);
template <typename Generator>
double compute_oracles(const hypergraph_set_type& graphs,
		       const feature_function_ptr_set_type& features,
		       const scorer_document_type& scorers,
		       sentence_set_type& sentences,
		       hypergraph_set_type& forests,
		       Generator& generator);


void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    threads = utils::bithack::max(threads, 1);
    min_iteration = utils::bithack::min(min_iteration, max_iteration);
    
    // read reference set
    scorer_document_type   scorers(scorer_name);
    sentence_document_type sentences;
    
    read_refset(refset_files, scorers, sentences);
    
    if (debug)
      std::cerr << "# of references: " << sentences.size() << std::endl;

    // read test set
    if (tstset_files.empty())
      tstset_files.push_back("-");

    if (debug)
      std::cerr << "reading hypergraphs" << std::endl;

    hypergraph_set_type       graphs(sentences.size());
    feature_function_ptr_set_type features(sentences.size());
    
    read_tstset(tstset_files, graphs, sentences, features);
    
    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    sentence_set_type   oracles(sentences.size());
    hypergraph_set_type oracles_forest(sentences.size());
    const double objective = compute_oracles(graphs, features, scorers, oracles, oracles_forest, generator);
    
    if (debug)
      std::cerr << "oracle score: " << objective << std::endl;
    
    if (! output_file.empty()) {
      if (directory_mode) {
	
	prepare_directory(output_file);
	
	if (forest_mode) {
	  for (size_t id = 0; id != oracles_forest.size(); ++ id)
	    if (oracles_forest[id].is_valid()) {
	      utils::compress_ostream os(output_file / (utils::lexical_cast<std::string>(id) + ".gz"), 1024 * 1024);
	      os.precision(10);
	      
	      os << id << " ||| " << oracles_forest[id] << '\n';
	    }
	} else {
	  for (size_t id = 0; id != oracles.size(); ++ id)
	    if (! oracles[id].empty()) {
	      utils::compress_ostream os(output_file / (utils::lexical_cast<std::string>(id) + ".gz"), 1024 * 1024);
	      os.precision(10);
	      
	      os << id << " ||| " << oracles[id] << '\n';
	    }
	}
	
      } else {
	utils::compress_ostream os(output_file, 1024 * 1024);
	os.precision(10);
	
	if (forest_mode) {
	  for (size_t id = 0; id != oracles_forest.size(); ++ id)
	    if (oracles_forest[id].is_valid())
	      os << id << " ||| " << oracles_forest[id] << '\n';
	} else {
	  for (size_t id = 0; id != oracles.size(); ++ id)
	    if (! oracles[id].empty())
	      os << id << " ||| " << oracles[id] << '\n';
	}
      }
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct TaskSingle
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;

  typedef std::vector<double, std::allocator<double> > score_set_type;

  TaskSingle(queue_type&                          __queue,
	     const hypergraph_set_type&           __graphs,
	     const feature_function_ptr_set_type& __features,
	     const scorer_document_type&          __scorers,
	     score_set_type&                      __scores)
    : queue(__queue),
      graphs(__graphs),
      features(__features),
      scorers(__scorers),
      scores(__scores)
  {
  }
  
  void operator()()
  {
    // we will try maximize    
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);
    
    weight_set_type::feature_type feature_scorer;
    for (size_t i = 0; i != features.size(); ++ i)
      if (features[i]) {
	feature_scorer = features[i]->feature_name();
	break;
      }
    
    hypergraph_type graph_oracle;
    
    int id = 0;
    while (1) {
      queue.pop(id);
      if (id < 0) break;
      
      if (! graphs[id].is_valid()) continue;

      typedef cicada::semiring::Logprob<double> weight_type;
      
      model_type model;
      model.push_back(features[id]);
      
      if (apply_exact)
	cicada::apply_exact(model, graphs[id], graph_oracle);
      else
	cicada::apply_cube_prune(model, graphs[id], graph_oracle, cicada::operation::single_scaled_function<weight_type >(feature_scorer, score_factor), scorer_cube);
      
      // compute viterbi...
      weight_type weight;
      sentence_type sentence;
      cicada::viterbi(graph_oracle, sentence, weight, cicada::operation::sentence_traversal(), cicada::operation::single_scaled_function<weight_type >(feature_scorer, score_factor));
      
      // compute objective
      scores[id] = scorers[id]->score(sentence)->score() * score_factor;
    }
    
  }
  
  queue_type&                          queue;
  const hypergraph_set_type&           graphs;
  const feature_function_ptr_set_type& features;
  const scorer_document_type&          scorers;
  score_set_type&                      scores;
};

struct TaskOracle
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  TaskOracle(queue_type&                          __queue,
	     const hypergraph_set_type&           __graphs,
	     const feature_function_ptr_set_type& __features,
	     const scorer_document_type&          __scorers,
	     sentence_set_type&                   __sentences,
	     hypergraph_set_type&                 __forests,
	     score_ptr_set_type&                  __scores)
    : queue(__queue),
      graphs(__graphs),
      features(__features),
      scorers(__scorers),
      sentences(__sentences),
      forests(__forests),
      scores(__scores)
  {
    score_optimum.reset();
    
    score_ptr_set_type::const_iterator siter_end = scores.end();
    for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter) 
      if (*siter) {
	if (! score_optimum)
	  score_optimum = (*siter)->clone();
	else
	  *score_optimum += *(*siter);
      } 
  }
  
  
  void operator()()
  {
    // we will try maximize    
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);
    
    weight_set_type::feature_type feature_scorer;
    for (size_t i = 0; i != features.size(); ++ i)
      if (features[i]) {
	feature_scorer = features[i]->feature_name();
	break;
      }
    
    double objective_optimum = (score_optimum
				? score_optimum->score() * score_factor
				: - std::numeric_limits<double>::infinity());

    hypergraph_type graph_oracle;
    
    int id = 0;
    while (1) {
      queue.pop(id);
      if (id < 0) break;
      
      if (! graphs[id].is_valid()) continue;

      score_ptr_type score_curr;
      if (score_optimum)
	score_curr = score_optimum->clone();
      
      if (scores[id])
	*score_curr -= *scores[id];
      
      cicada::feature::Scorer* __scorer = dynamic_cast<cicada::feature::Scorer*>(features[id].get());
      
      if (__scorer)
	__scorer->assign(score_curr);

      typedef cicada::semiring::Logprob<double> weight_type;
      
      model_type model;
      model.push_back(features[id]);
      
      if (apply_exact)
	cicada::apply_exact(model, graphs[id], graph_oracle);
      else
	cicada::apply_cube_prune(model, graphs[id], graph_oracle, cicada::operation::single_scaled_function<weight_type >(feature_scorer, score_factor), scorer_cube);
      
      // compute viterbi...
      weight_type weight;
      sentence_type sentence;
      cicada::viterbi(graph_oracle, sentence, weight, cicada::operation::sentence_traversal(), cicada::operation::single_scaled_function<weight_type >(feature_scorer, score_factor));
      
      // compute pruned forest
      hypergraph_type forest;
      cicada::prune_beam(graph_oracle, forest, cicada::operation::single_scaled_function<cicada::semiring::Tropical<double> >(feature_scorer, score_factor), scorer_beam);
      
      // compute scores...
      score_ptr_type score_sample = scorers[id]->score(sentence);
      if (score_curr)
	*score_curr += *score_sample;
      else
	score_curr = score_sample;
      
      const double objective = score_curr->score() * score_factor;
      
      if (objective > objective_optimum || ! scores[id]) {
	score_optimum = score_curr;
	objective_optimum = objective;
	scores[id] = score_sample;
	
	sentences[id].swap(sentence);
	forests[id].swap(forest);

	// remove features...
	hypergraph_type::edge_set_type::iterator eiter_end = forests[id].edges.end();
	for (hypergraph_type::edge_set_type::iterator eiter = forests[id].edges.begin(); eiter != eiter_end; ++ eiter)
	  eiter->features.erase(feature_scorer);
      }
    }
    
  }
  
  score_ptr_type score_optimum;
  
  queue_type&                          queue;
  const hypergraph_set_type&           graphs;
  const feature_function_ptr_set_type& features;
  const scorer_document_type&          scorers;
  sentence_set_type&                   sentences;
  hypergraph_set_type&                 forests;
  score_ptr_set_type&                  scores;
};

template <typename Scores>
struct greater_score
{
  greater_score(const Scores& __scores) : scores(__scores) {}

  bool operator()(const int& x, const int& y) const
  {
    return scores[x] > scores[y];
  }
  
  const Scores& scores;
};

template <typename Generator>
double compute_oracles(const hypergraph_set_type& graphs,
		       const feature_function_ptr_set_type& features,
		       const scorer_document_type& scorers,
		       sentence_set_type& sentences,
		       hypergraph_set_type& forests,
		       Generator& generator)
{
  typedef TaskOracle            task_type;
  typedef task_type::queue_type queue_type;
  
  typedef boost::shared_ptr<task_type> task_ptr_type;
  typedef std::vector<task_ptr_type, std::allocator<task_ptr_type> > task_set_type;
  
  typedef std::vector<size_t, std::allocator<size_t> > id_set_type;

  score_ptr_set_type scores(graphs.size());

  id_set_type ids;
  for (size_t id = 0; id != graphs.size(); ++ id)
    if (graphs[id].is_valid())
      ids.push_back(id);

  score_ptr_type score_optimum;
  
  double objective_prev = - std::numeric_limits<double>::infinity();
  double objective_best = - std::numeric_limits<double>::infinity();
  
  sentence_set_type   sentences_best(sentences.size());
  hypergraph_set_type forests_best(forests.size());
  
  const bool error_metric = scorers.error_metric();
  const double score_factor = (error_metric ? - 1.0 : 1.0);

#if 0
  {
    typedef TaskSingle task_type;
    
    task_type::queue_type queue;

    task_type::score_set_type scores(graphs.size());
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(task_type(queue, graphs, features, scorers, scores)));
    
    id_set_type::const_iterator iiter_end = ids.end();
    for (id_set_type::const_iterator iiter = ids.begin(); iiter != iiter_end; ++ iiter)
      queue.push(*iiter);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    workers.join_all();
    
    // sort ids by the scores so that we can process form the less-errored hypotheses
    std::sort(ids.begin(), ids.end(), greater_score<task_type::score_set_type>(scores));
  }
#endif
  
  for (int iter = 0; iter < max_iteration; ++ iter) {
    if (debug)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    queue_type queue;
    
    task_set_type tasks(threads);
    for (int i = 0; i < threads; ++ i)
      tasks[i].reset(new task_type(queue, graphs, features, scorers, sentences, forests, scores));
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(*tasks[i])));
    
    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);
    
    id_set_type::const_iterator iiter_end = ids.end();
    for (id_set_type::const_iterator iiter = ids.begin(); iiter != iiter_end; ++ iiter)
      queue.push(*iiter);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
        
    workers.join_all();
    
    score_optimum.reset();
    score_ptr_set_type::const_iterator siter_end = scores.end();
    for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter) 
      if (*siter) {
	if (! score_optimum)
	  score_optimum = (*siter)->clone();
	else
	  *score_optimum += *(*siter);
      } 
    
    const double objective = score_optimum->score() * score_factor;
    if (debug)
      std::cerr << "oracle score: " << objective << std::endl;
    
    // if we found better objective, keep going!
    if (objective > objective_best) {
      objective_best = objective;
      sentences_best = sentences;
      forests_best   = forests;
    }
    
    // termination condition
    if (objective <= objective_prev && iter >= min_iteration)
      break;
    
    objective_prev = objective;
  }
  
  sentences.swap(sentences_best);
  forests.swap(forests_best);
  
  return objective_best;
}

typedef std::pair<path_type, std::string> path_line_type;

namespace std
{
  inline
  void swap(path_line_type& x, path_line_type& y)
  {
    x.first.swap(y.first);
    x.second.swap(y.second);
  }
};

struct ReadTstset
{
  typedef utils::lockfree_list_queue<path_line_type, std::allocator<path_line_type> > queue_type;
  
  ReadTstset(queue_type& __queue, const size_t size) : queue(__queue), graphs(size) {}

  void operator()()
  {
    size_t id;
    hypergraph_type hypergraph;
    
    path_line_type path_line;
    
    for (;;) {
      queue.pop_swap(path_line);
      
      if (path_line.first.empty() && path_line.second.empty()) break;
      
      if (! path_line.first.empty()) {
	utils::compress_istream is(path_line.first, 1024 * 1024);
	
	if (! utils::getline(is, path_line.second))
	  throw std::runtime_error("no line: " + path_line.first.string());
      }

      const std::string& line = path_line.second;
      
      if (line.empty()) continue;
      
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end  = line.end();
      
      if (! parse_id(id, iter, end))
	throw std::runtime_error("invalid id input: " + line);
      if (id >= graphs.size())
	throw std::runtime_error("tstset size exceeds refset size? " + line);
      
      if (! hypergraph.assign(iter, end))
	throw std::runtime_error("invalid graph format: " + line);
      if (iter != end)
	throw std::runtime_error("invalid id ||| graph format: " + line);

      if (graphs[id].is_valid())
	graphs[id].unite(hypergraph);
      else
	graphs[id].swap(hypergraph);
    }
  }
  
  queue_type& queue;
  hypergraph_set_type graphs;
};

void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 feature_function_ptr_set_type& features)
{
  typedef ReadTstset task_type;

  task_type::queue_type queue(64 * threads);
  
  std::vector<task_type, std::allocator<task_type> > tasks(threads, task_type(queue, graphs.size()));
  
  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (debug)
      std::cerr << "file: " << *titer << std::endl;
    
    if (boost::filesystem::is_directory(*titer)) {
      for (size_t i = 0; /**/; ++ i) {
	const path_type path = (*titer) / (utils::lexical_cast<std::string>(i) + ".gz");
	
	if (! boost::filesystem::exists(path)) break;

	queue.push(std::make_pair(path, std::string()));
      }
    } else {
      path_line_type path_line;
      
      utils::compress_istream is(*titer, 1024 * 1024);
      
      while (utils::getline(is, path_line.second))
	if (! path_line.second.empty())
	  queue.push_swap(path_line);
    }
  }
  
  // terminaltion
  for (int i = 0; i != threads; ++ i)
    queue.push(std::make_pair(path_type(), std::string()));
  
  workers.join_all();
  
  // merging...
  for (int i = 0; i != threads; ++ i) {
    for (size_t seg = 0; seg != graphs.size(); ++ seg)
      if (tasks[i].graphs[seg].is_valid()) {
	
	if (graphs[seg].is_valid())
	  graphs[seg].unite(tasks[i].graphs[seg]);
	else
	  graphs[seg].swap(tasks[i].graphs[seg]);
	
	tasks[i].graphs[seg].clear();
      }
    tasks[i].graphs.clear();
  }
  
  // finish!
  tasks.clear();
  
  if (debug)
    std::cerr << "assign scorer feature" << std::endl;
  
  typedef cicada::Parameter parameter_type;
    
  for (size_t id = 0; id != graphs.size(); ++ id) {
    if (graphs[id].goal == hypergraph_type::invalid)
      std::cerr << "invalid graph at: " << id << std::endl;
    else {
      features[id] = feature_function_type::create(scorer_name);
      
      cicada::feature::Scorer* __scorer = dynamic_cast<cicada::feature::Scorer*>(features[id].get());
      
      if (! __scorer)
	throw std::runtime_error("invalid scorer feature function...");

      static const cicada::Lattice       __lattice;
      static const cicada::SpanVector    __spans;
      static const cicada::NGramCountSet __ngram_counts;
      
      __scorer->assign(id, graphs[id], __lattice, __spans, sentences[id], __ngram_counts);
    }
  }
}

void read_refset(const path_set_type& files,
		 scorer_document_type& scorers,
		 sentence_document_type& sentences)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  if (files.empty())
    throw std::runtime_error("no reference files?");

  parser_type parser;
  id_sentence_type id_sentence;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->string());

    utils::compress_istream is(*fiter, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      id_sentence.second.clear();
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
	if (iter != iter_end)
	  throw std::runtime_error("refset parsing failed");
      
      const int& id = id_sentence.first;
      
      if (id >= static_cast<int>(scorers.size()))
	scorers.resize(id + 1);
      if (id >= static_cast<int>(sentences.size()))
	sentences.resize(id + 1);
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      sentences[id].push_back(id_sentence.second);
      scorers[id]->insert(sentences[id].back());
    }
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s) (in hypergraph format)")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("forest",    po::bool_switch(&forest_mode),    "output by forest")
    ("directory", po::bool_switch(&directory_mode), "output in directory")
        
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-cube", po::value<int>(&scorer_cube)->default_value(scorer_cube),         "cube pruning size")
    ("scorer-beam", po::value<double>(&scorer_beam)->default_value(scorer_beam),      "beam pruning size")
    
    ("max-iteration", po::value<int>(&max_iteration), "# of hill-climbing iteration")
    ("min-iteration", po::value<int>(&min_iteration), "# of hill-climbing iteration")
    
    ("apply-exact", po::bool_switch(&apply_exact), "exact application")
    
    ("threads", po::value<int>(&threads), "# of threads")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc_command << std::endl;
    exit(0);
  }
}
