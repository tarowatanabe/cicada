
//
// online learning using the margin between forest
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

#include "cicada/apply.hpp"
#include "cicada/model.hpp"

#include "cicada/feature/bleu.hpp"
#include "cicada/feature/bleu_linear.hpp"
#include "cicada/parameter.hpp"

#include "cicada/eval.hpp"

#include "cicada_impl.hpp"

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


typedef std::vector<path_type, std::allocator<path_type> > path_set_type;



typedef std::vector<feature_function_ptr_type, std::allocator<feature_function_ptr_type> > feature_function_ptr_set_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::scorer_ptr_type scorer_ptr_type;
typedef scorer_type::score_ptr_type  score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;


typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file = "-";
path_type output_file = "-";

bool input_id_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_directory_mode = false;

std::string symbol_goal         = vocab_type::S;
std::string symbol_non_terminal = vocab_type::X;

grammar_file_set_type grammar_mutable_files;
grammar_file_set_type grammar_static_files;

bool grammar_glue_straight = false;
bool grammar_glue_inverted = false;
bool grammar_insertion = false;
bool grammar_deletion = false;

feature_parameter_set_type feature_parameters;
bool feature_list = false;

op_set_type ops;
bool op_list = false;

std::string algorithm_name = "perceptron";
std::string scorer_name    = "bleu:order=4,exact=true";

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;
double loss_margin = 0.01;
double score_margin = 0.01;
double loss_scale = 100;

bool apply_exact = false;
int  cube_size = 200;

int threads = 4;

int debug = 0;

void optimize(weight_set_type& weights);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (input_lattice_mode && input_forest_mode)
      throw std::runtime_error("input can be sentence, lattice or forest");
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("you cannot use both of L1 and L2...");

    if (feature_list) {
      std::cout << cicada::FeatureFunction::lists();
      return 0;
    }

    if (op_list) {
      std::cout << OperationSet::lists();
      return 0;
    }

    srandom(time(0) * getpid());
    
    threads = utils::bithack::max(threads, 1);
    
    
    weight_set_type weights;
    
    optimize(weights);
    
    utils::compress_ostream os(output_file);
    os.precision(20);
    os << weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct Task
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;

  Task(queue_type& __queue)
    : queue(__queue), weights(), weights_accumulated(), updated(1) { initialize(); }
  
  queue_type&        queue;
  
  weight_set_type    weights;
  weight_set_type    weights_accumulated;
  size_t             updated;
  score_ptr_type     score;
  score_ptr_set_type scores;
  
  grammar_type grammar;
  model_type model;
  
  typedef cicada::semiring::Logprob<double> weight_type;
  
  struct weight_set_function
  {
    typedef cicada::semiring::Logprob<double> value_type;
    
    weight_set_function(const weight_set_type& __weights, const double& __scale)
      : weights(__weights), scale(__scale) {}
    
    const weight_set_type& weights;
    const double scale;
    
    template <typename Edge>
    value_type operator()(const Edge& x) const
    {
      return cicada::semiring::traits<value_type>::log(x.features.dot(weights) * scale);
    }
    
    value_type operator()(const feature_set_type& x) const
    {
      return cicada::semiring::traits<value_type>::log(x.dot(weights) * scale);
    }

  };

  struct target_traversal
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
	  yield.insert(yield.end(), (first + pos)->begin(), (first + pos)->end());
	} else if (*titer != vocab_type::EPSILON)
	  yield.push_back(*titer);
    }
  };

  struct diff_norm
  {
    double operator()(const double& x, const double& y) const{
      return (x - y) * (x - y);
    }
  };
  
  
  void operator()()
  {
    typedef boost::tuple<sentence_type, feature_set_type> yield_type;

    scorer_ptr_type           scorer(scorer_type::create(scorer_name));
    feature_function_ptr_type feature_function(feature_function_type::create(scorer_name));

    cicada::feature::Bleu* __bleu = dynamic_cast<cicada::feature::Bleu*>(feature_function.get());
    if (! __bleu)
      throw std::runtime_error("not supported feature");
    
    model_type model_bleu;
    model_bleu.push_back(feature_function);
    
    OperationSet operations(ops.begin(), ops.end(),
			    grammar,
			    model,
			    symbol_goal,
			    symbol_non_terminal,
			    grammar_insertion,
			    grammar_deletion,
			    true,
			    input_lattice_mode,
			    input_forest_mode,
			    true,
			    false,
			    debug);
    
    
    operations.assign(weights);

    hypergraph_type hypergraph_reward;
    hypergraph_type hypergraph_penalty;      
    
    std::string line;
    while (1) {
      queue.pop(line);
      if (line.empty()) break;
      
      operations(line);
      
      // operations.hypergraph contains result...
      const size_t& id = operations.id;
      const lattice_type& lattice = operations.lattice;
      const hypergraph_type& hypergraph = operations.hypergraph;
      const sentence_set_type& targets = operations.targets;

      std::cerr << "id: " << id << std::endl;

      // compute source-length
      int source_length = lattice.shortest_distance();
      if (hypergraph.is_valid()) {
	// we will enumerate forest structure... and collect min-size...
	std::vector<source_length_function::value_type, std::allocator<source_length_function::value_type> > lengths(hypergraph.nodes.size());
	
	cicada::inside(hypergraph, lengths, source_length_function());
	
	source_length = - log(lengths.back());
      }
      
      std::cerr << "source length: " << source_length << std::endl;

      // collect max-feature from hypergraph
      yield_type  yield_viterbi;
      weight_type weight_viterbi;
      cicada::viterbi(hypergraph, yield_viterbi, weight_viterbi, kbest_traversal(), weight_set_function(weights, 1.0));

      std::cerr << "viterbi: " << boost::get<0>(yield_viterbi) << std::endl;
            
      // update scores...
      if (id >= scores.size())
	scores.resize(id + 1);
      
      // remove "this" score
      if (score && scores[id])
	*score -= *scores[id];
      
      // create scorers...
      scorer->clear();
      __bleu->clear();
      sentence_set_type::const_iterator titer_end = targets.end();
      for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer) {
	scorer->insert(*titer);
	__bleu->insert(source_length, *titer);
      }
      __bleu->insert(score);
      
      // compute bleu-rewarded instance
      yield_type  yield_reward;
      weight_type weight_reward;
      
      weights[__bleu->feature_name()] =  loss_scale * source_length;
      cicada::apply_cube_prune(model_bleu, hypergraph, hypergraph_reward, weight_set_function(weights, 1.0), cube_size);
      cicada::viterbi(hypergraph_reward, yield_reward, weight_reward, kbest_traversal(), weight_set_function(weights, 1.0));

      std::cerr << "reward: " << boost::get<0>(yield_reward) << std::endl;
      
      // compute bleu-penalty hypergraph
      yield_type  yield_penalty;
      weight_type weight_penalty;
      
      weights[__bleu->feature_name()] = - loss_scale * source_length;
      cicada::apply_cube_prune(model_bleu, hypergraph, hypergraph_penalty, weight_set_function(weights, 1.0), cube_size);
      cicada::viterbi(hypergraph_penalty, yield_penalty, weight_penalty, kbest_traversal(), weight_set_function(weights, 1.0));

      std::cerr << "penalty: " << boost::get<0>(yield_penalty) << std::endl;
      
      weights[__bleu->feature_name()] = 0.0;
      
      score_ptr_type score_reward  = scorer->score(boost::get<0>(yield_reward));
      score_ptr_type score_penalty = scorer->score(boost::get<0>(yield_penalty));
      scores[id] = scorer->score(boost::get<0>(yield_viterbi));
      
      if (score) {
	*score_reward  += *score;
	*score_penalty += *score;
      }
      
      if (! score)
	score = scores[id]->clone();
      else
	*score += *scores[id];
      
      const std::pair<double, double> bleu_viterbi = score->score();
      const std::pair<double, double> bleu_reward  = score_reward->score();
      const std::pair<double, double> bleu_penalty = score_penalty->score();
      
      if (debug)
	std::cerr << "viterbi score: "  << bleu_viterbi.first << " penalty: " << bleu_viterbi.second << std::endl
		  << "oracle score: "   << bleu_reward.first  << " penalty: " << bleu_reward.second << std::endl
		  << "violated score: " << bleu_penalty.first << " penalty: " << bleu_penalty.second << std::endl;
      
      const double loss   = loss_scale * source_length * (bleu_reward.first - bleu_penalty.first);
      const double margin = boost::get<1>(yield_penalty).dot(weights) - boost::get<1>(yield_reward).dot(weights);
      const double norm   = boost::get<1>(yield_penalty).dot(boost::get<1>(yield_reward), diff_norm());
      
      const double alpha = std::max(0.0, std::min(1.0 / C, (loss - margin) / norm));
      
      if (loss - margin > 0.0) {
	std::cerr << "loss: " << loss << " margin: " << margin << " norm: " << norm << " alpha: " << alpha << std::endl;

	// update...
	feature_set_type::const_iterator riter_end = boost::get<1>(yield_reward).end();
	for (feature_set_type::const_iterator riter = boost::get<1>(yield_reward).begin(); riter != riter_end; ++ riter) {
	  weights[riter->first] += riter->second * alpha;
	  weights_accumulated[riter->first] += riter->second * alpha * updated;
	}
	
	feature_set_type::const_iterator piter_end = boost::get<1>(yield_penalty).end();
	for (feature_set_type::const_iterator piter = boost::get<1>(yield_penalty).begin(); piter != piter_end; ++ piter) {
	  weights[piter->first] -= piter->second * alpha;
	  weights_accumulated[piter->first] -= piter->second * alpha * updated;
	}

	++ updated;
      }
    }
  }
  
  void initialize()
  {
    // read grammars...
    for (grammar_file_set_type::const_iterator fiter = grammar_static_files.begin(); fiter != grammar_static_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarStatic(*fiter)));

    if (debug)
      std::cerr << "loaded mutable grammar: " << grammar.size() << std::endl;
    
    for (grammar_file_set_type::const_iterator fiter = grammar_mutable_files.begin(); fiter != grammar_mutable_files.end(); ++ fiter)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarMutable(*fiter)));

    if (debug)
      std::cerr << "loaded static grammar: " << grammar.size() << std::endl;
    
    if (grammar_glue_straight || grammar_glue_inverted)
      grammar.push_back(grammar_type::transducer_ptr_type(new cicada::GrammarGlue(symbol_goal,
										  symbol_non_terminal,
										  grammar_glue_straight,
										  grammar_glue_inverted)));

    if (debug)
      std::cerr << "grammar: " << grammar.size() << std::endl;
    
    // read features...
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));
    model.initialize();
  }  
};


void optimize(weight_set_type& weights)
{
  typedef Task                         task_type;
  typedef boost::shared_ptr<task_type> task_ptr_type;
  typedef std::vector<task_ptr_type, std::allocator<task_ptr_type> > task_ptr_set_type;

  typedef task_type::queue_type queue_type;
  
  typedef std::vector<std::string, std::allocator<std::string> > sample_set_type;
  
  // read all the training data...
  sample_set_type samples;
  if (input_directory_mode) {
    std::string line;
      
    boost::filesystem::directory_iterator iter_end;
    for (boost::filesystem::directory_iterator iter(input_file); iter != iter_end; ++ iter) {
      utils::compress_istream is(*iter, 1024 * 1024);
      
      if (std::getline(is, line) && ! line.empty())
	samples.push_back(line);
    }
  } else {
    utils::compress_istream is(input_file, 1024 * 1024);
    
    size_t id = 0;

    std::string line;
    while (std::getline(is, line))
      if (! line.empty()) {
	if (! input_id_mode)
	  samples.push_back(boost::lexical_cast<std::string>(id) + " ||| " + line);
	else
	  samples.push_back(line);
	++ id;
      }
  }

  if (debug)
    std::cerr << "# of samples: " << samples.size() << std::endl;
  
  queue_type queue(threads * 2);
  
  task_ptr_set_type tasks(threads);
  for (int i = 0; i < threads; ++ i)
    tasks[i].reset(new task_type(queue));

  weight_set_type weights_mixed;
  weight_set_type weights_accumulated;
  double norm_accumulated = 0;
  
  for (int iter = 0; iter < iteration; ++ iter) {
    
    if (debug)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(*tasks[i])));
    
    sample_set_type::const_iterator siter_end = samples.end();
    for (sample_set_type::const_iterator siter = samples.begin(); siter != siter_end; ++ siter)
      queue.push(*siter);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(std::string());
    
    
    workers.join_all();
    
    std::random_shuffle(samples.begin(), samples.end());
    
    // merge vector...
    weights.clear();
    weights_mixed.clear();
    double norm = 0.0;

    
    task_ptr_set_type::iterator titer_end = tasks.end();
    for (task_ptr_set_type::iterator titer = tasks.begin(); titer != titer_end; ++ titer) {
      weight_set_type& weights_local = (*titer)->weights;
      
      // mixing weights...
      weights_mixed += weights_local;
      
      weights_local *= (*titer)->updated;
      weights_local -= (*titer)->weights_accumulated;
      
      weights += weights_local;
      norm    += (*titer)->updated;
    }
    
    weights_accumulated += weights;
    norm_accumulated += norm;
    
    // mixing...
    weights_mixed *= (1.0 / tasks.size());
    
    for (task_ptr_set_type::iterator titer = tasks.begin(); titer != titer_end; ++ titer) {
      (*titer)->weights = weights_mixed;
      (*titer)->weights_accumulated.clear();
      (*titer)->updated = 1;
    }
  }
  
  weights = weights_accumulated;
  weights /= norm_accumulated;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    // options for input/output format
    ("input-id",         po::bool_switch(&input_id_mode),         "id-prefixed input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")
    
    // grammar
    ("goal",           po::value<std::string>(&symbol_goal)->default_value(symbol_goal),                 "goal symbol")
    ("non-terminal",   po::value<std::string>(&symbol_non_terminal)->default_value(symbol_non_terminal), "default non-terminal symbol")
    ("grammar",        po::value<grammar_file_set_type >(&grammar_mutable_files)->composing(),           "grammar file(s)")
    ("grammar-static", po::value<grammar_file_set_type >(&grammar_static_files)->composing(),            "static binary grammar file(s)")
    
    // special handling
    ("grammar-glue-straight", po::bool_switch(&grammar_glue_straight), "add straight hiero glue rule")
    ("grammar-glue-inverted", po::bool_switch(&grammar_glue_inverted), "add inverted hiero glue rule")
    ("grammar-insertion",     po::bool_switch(&grammar_insertion),     "source-to-target transfer grammar")
    ("grammar-deletion",      po::bool_switch(&grammar_deletion),      "source-to-<epsilon> transfer grammar")
    
    // models...
    ("feature-function",      po::value<feature_parameter_set_type >(&feature_parameters)->composing(), "feature function(s)")
    ("feature-function-list", po::bool_switch(&feature_list),                                           "list of available feature function(s)")

    //operatins...
    ("operation",      po::value<op_set_type>(&ops)->composing(), "operations")
    ("operation-list", po::bool_switch(&op_list),                 "list of available operation(s)")
    
    // learning related..
    ("algorithm",   po::value<std::string>(&algorithm_name)->default_value(algorithm_name),     "learning algorithm")
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    
    ("iteration",          po::value<int>(&iteration),          "# of mert iteration")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "regularization via L1")
    ("regularize-l2", po::bool_switch(&regularize_l2), "regularization via L2")
    ("C"            , po::value<double>(&C),           "regularization constant")
    
    ("loss-margin",   po::value<double>(&loss_margin)->default_value(loss_margin),   "loss margin for oracle forest")
    ("score-margin",  po::value<double>(&score_margin)->default_value(score_margin), "score margin for hypothesis forest")
    ("loss-scale",    po::value<double>(&loss_scale)->default_value(loss_scale),     "loss scaling")
    
    ("apply-exact", po::bool_switch(&apply_exact), "exact feature applicatin w/o pruning")
    ("cube-size",   po::value<int>(&cube_size),    "cube-pruning size")
    
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
