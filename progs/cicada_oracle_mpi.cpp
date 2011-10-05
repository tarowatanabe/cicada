//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// refset format:
// 0 |||  reference translatin for source sentence 0
// 0 |||  another reference
// 1 |||  reference translation for source sentence 1
//


#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

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

#include "cicada/feature/scorer.hpp"
#include "cicada/parameter.hpp"
#include "cicada/prune.hpp"
#include "cicada/eval.hpp"

#include "cicada/operation/functional.hpp"
#include "cicada/operation/traversal.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/resource.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random.hpp>

#include <boost/thread.hpp>

#include "lbfgs.h"

#include "utils/base64.hpp"
#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"
#include "utils/mpi_device_bcast.hpp"
#include "utils/mpi_stream.hpp"
#include "utils/mpi_stream_simple.hpp"

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

int max_iteration = 10;
int min_iteration = 5;
bool apply_exact = false;
int cube_size = 200;
double beam_size = 1e-5;

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

void bcast_sentences(sentence_set_type& sentences, hypergraph_set_type& forests);
void bcast_weights(const int rank, weight_set_type& weights);
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  try {
    options(argc, argv);
    
    // read test set
    if (tstset_files.empty())
      throw std::runtime_error("no test set?");

    min_iteration = utils::bithack::min(min_iteration, max_iteration);
        
    // read reference set
    scorer_document_type   scorers(scorer_name);
    sentence_document_type sentences;
    
    read_refset(refset_files, scorers, sentences);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "# of references: " << sentences.size() << std::endl;
    
    if (mpi_rank == 0 && debug)
      std::cerr << "reading hypergraphs" << std::endl;
    
    hypergraph_set_type       graphs(sentences.size());
    feature_function_ptr_set_type features(sentences.size());
    
    read_tstset(tstset_files, graphs, sentences, features);

    if (mpi_rank == 0 && debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    boost::mt19937 generator;
    generator.seed(time(0) * getpid());
    
    sentence_set_type oracles(sentences.size());
    hypergraph_set_type oracles_forest(sentences.size());
    const double objective = compute_oracles(graphs, features, scorers, oracles, oracles_forest, generator);

    if (mpi_rank == 0 && debug)
      std::cerr << "oracle score: " << objective << std::endl;
    
    if (mpi_rank == 0) {
      if (directory_mode) {
	if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
	  boost::filesystem::remove_all(output_file);
	
	boost::filesystem::create_directories(output_file);
	
	boost::filesystem::directory_iterator iter_end;
	for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
	  boost::filesystem::remove_all(*iter);
	
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
    MPI::COMM_WORLD.Abort(1);
    return 1;
  }
  return 0;
}

enum {
  weights_tag = 1000,
  sentence_tag,
  gradients_tag,
  notify_tag,
  termination_tag,
};

inline
int loop_sleep(bool found, int non_found_iter)
{
  if (! found) {
    boost::thread::yield();
    ++ non_found_iter;
  } else
    non_found_iter = 0;
    
  if (non_found_iter >= 64) {
    struct timespec tm;
    tm.tv_sec = 0;
    tm.tv_nsec = 2000001; // above 2ms
    nanosleep(&tm, NULL);
    
    non_found_iter = 0;
  }
  return non_found_iter;
}

template <typename Generator>
struct TaskOracle
{
  TaskOracle(const hypergraph_set_type&           __graphs,
	     const feature_function_ptr_set_type& __features,
	     const scorer_document_type&          __scorers,
	     score_ptr_set_type&                  __scores,
	     sentence_set_type&                   __sentences,
	     hypergraph_set_type&                 __forests,
	     Generator&                           __generator)
    : graphs(__graphs),
      features(__features),
      scorers(__scorers),
      scores(__scores),
      sentences(__sentences),
      forests(__forests),
      generator(__generator)
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

    typedef std::vector<size_t, std::allocator<size_t> > id_set_type;

    id_set_type ids;
    for (size_t id = 0; id != graphs.size(); ++ id)
      if (graphs[id].is_valid())
	ids.push_back(id);
    
    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);
    
    id_set_type::const_iterator iiter_end = ids.end();
    for (id_set_type::const_iterator iiter = ids.begin(); iiter != iiter_end; ++ iiter) {
      const size_t id = *iiter;
      
      score_ptr_type score_curr;
      if (score_optimum)
	score_curr = score_optimum->clone();
      
      if (scores[id])
	*score_curr -= *scores[id];
      
      cicada::feature::Scorer* __scorer = dynamic_cast<cicada::feature::Scorer*>(features[id].get());
      
      __scorer->assign(score_curr);

      typedef cicada::semiring::Logprob<double> weight_type;
      
      model_type model;
      model.push_back(features[id]);
      
      if (apply_exact)
	cicada::apply_exact(model, graphs[id], graph_oracle);
      else
	cicada::apply_cube_prune(model, graphs[id], graph_oracle, cicada::operation::single_scaled_function<weight_type>(feature_scorer, score_factor), cube_size);
      
      // compute viterbi...
      weight_type weight;
      sentence_type sentence;
      cicada::viterbi(graph_oracle, sentence, weight, cicada::operation::sentence_traversal(), cicada::operation::single_scaled_function<weight_type >(feature_scorer, score_factor));
      
      // compute pruned forest
      hypergraph_type forest;
      cicada::prune_beam(graph_oracle, forest, cicada::operation::single_scaled_function<cicada::semiring::Tropical<double> >(feature_scorer, score_factor), beam_size);
      
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
  
  const hypergraph_set_type&           graphs;
  const feature_function_ptr_set_type& features;
  const scorer_document_type&          scorers;
  score_ptr_set_type&                  scores;
  sentence_set_type&                   sentences;
  hypergraph_set_type&                 forests;
  Generator&                           generator;
};

template <typename Generator>
double compute_oracles(const hypergraph_set_type& graphs,
		       const feature_function_ptr_set_type& features,
		       const scorer_document_type& scorers,
		       sentence_set_type& sentences,
		       hypergraph_set_type& forests,
		       Generator& generator)
{
  typedef TaskOracle<Generator> task_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  score_ptr_set_type scores(graphs.size());
  
  score_ptr_type score_optimum;
  
  double objective_prev = - std::numeric_limits<double>::infinity();
  double objective_best = - std::numeric_limits<double>::infinity();
  
  sentence_set_type   sentences_best;
  hypergraph_set_type forests_best;
  
  const bool error_metric = scorers.error_metric();
  const double score_factor = (error_metric ? - 1.0 : 1.0);
  
  for (int iter = 0; iter < max_iteration; ++ iter) {
    if (debug && mpi_rank == 0)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    task_type(graphs, features, scorers, scores, sentences, forests, generator)();

    bcast_sentences(sentences, forests);
    
    score_optimum.reset();
    for (size_t id = 0; id != sentences.size(); ++ id)
      if (! sentences[id].empty()) {
	scores[id] = scorers[id]->score(sentences[id]);
	
	if (! score_optimum)
	  score_optimum = scores[id]->clone();
	else
	  *score_optimum += *scores[id];
      }
    
    const double objective = score_optimum->score() * score_factor;
    if (mpi_rank == 0 && debug)
      std::cerr << "oracle score: " << objective << std::endl;
    
    // if we found better objective, keep going!
    if (objective > objective_best) {
      objective_best = objective;
      sentences_best = sentences;
      forests_best   = forests;
    }
    
    int terminate = (objective <= objective_prev) && (iter >= min_iteration);
    MPI::COMM_WORLD.Bcast(&terminate, 1, MPI::INT, 0);
    
    if (terminate)
      break;
    
    objective_prev = objective;
  }
  
  sentences.swap(sentences_best);
  forests.swap(forests_best);
  
  return objective_best;
}


void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 feature_function_ptr_set_type& features)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  std::string line;

  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (mpi_rank == 0 && debug)
      std::cerr << "file: " << *titer << std::endl;
      
    if (boost::filesystem::is_directory(*titer)) {

      size_t id;
      hypergraph_type hypergraph;
      
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
	const path_type path = (*titer) / (utils::lexical_cast<std::string>(i) + ".gz");
	
	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	
	if (! std::getline(is, line))
	  throw std::runtime_error("no line in file-no: " + utils::lexical_cast<std::string>(i));
	
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end  = line.end();
	
	if (! parse_id(id, iter, end))
	  throw std::runtime_error("invalid id input: " + path.string());
	if (id != i)
	  throw std::runtime_error("id mismatch: "  + path.string());
	if (static_cast<int>(id % mpi_size) != mpi_rank)
	  throw std::runtime_error("difference it?");
	
	if (! hypergraph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path.string());
	
	graphs[id].unite(hypergraph);
      }
    } else {
      const path_type& path = *titer;

      utils::compress_istream is(path, 1024 * 1024);
      
      size_t id;
      hypergraph_type hypergraph;

      while (std::getline(is, line)) {
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end  = line.end();
	
	if (! parse_id(id, iter, end))
	  throw std::runtime_error("invalid id input: " + path.string());
	if (id >= graphs.size())
	  throw std::runtime_error("tstset size exceeds refset size?" + utils::lexical_cast<std::string>(id) + ": " + titer->string());
	
	if (static_cast<int>(id % mpi_size) != mpi_rank) continue;
	
	if (! hypergraph.assign(iter, end))
	  throw std::runtime_error("invalid graph format" + path.string());
	if (iter != end)
	  throw std::runtime_error("invalid id ||| graph format" + path.string());
	
	graphs[id].unite(hypergraph);
      }
    }
  }

  if (debug && mpi_rank == 0)
    std::cerr << "assign scorer feature" << std::endl;
    
  for (size_t id = 0; id != graphs.size(); ++ id) 
    if (static_cast<int>(id % mpi_size) == mpi_rank) {
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
  
  // collect weights...
  for (int rank = 0; rank < mpi_size; ++ rank) {
    weight_set_type weights;
    weights.allocate();
    
    for (feature_type::id_type id = 0; id != feature_type::allocated(); ++ id)
      if (! feature_type(id).empty())
	weights[feature_type(id)] = 1.0;
    
    bcast_weights(rank, weights);
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

void bcast_sentences(sentence_set_type& sentences, hypergraph_set_type& forests)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    if (rank == mpi_rank) {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::zlib_compressor());
      os.push(utils::mpi_device_bcast_sink(rank, 4096));
      
      for (size_t id = 0; id != sentences.size(); ++ id)
	if (static_cast<int>(id % mpi_size) == mpi_rank)
	  os << id << " ||| " << sentences[id] << " ||| " << forests[id] << '\n';
      
    } else {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::zlib_decompressor());
      is.push(utils::mpi_device_bcast_source(rank, 4096));
      
      std::string line;
      
      int id;
      sentence_type sentence;
      hypergraph_type forest;

      while (std::getline(is, line)) {
	sentence.clear();
	forest.clear();

	std::string::const_iterator iter_end = line.end();
	std::string::const_iterator iter = line.begin();

	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	namespace phoenix = boost::phoenix;

	if (! qi::phrase_parse(iter, iter_end, qi::int_, standard::space, id)) continue;
	if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) continue;
	if (! sentence.assign(iter, iter_end)) continue;
	if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) continue;
	if (! forest.assign(iter, iter_end)) continue;
	
	if (iter != iter_end) continue;
	
	sentences[id].swap(sentence);
	forests[id].swap(forest);
      }
    }
  }
}


void bcast_weights(const int rank, weight_set_type& weights)
{
  typedef std::vector<char, std::allocator<char> > buffer_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  if (mpi_rank == rank) {
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(utils::mpi_device_bcast_sink(rank, 4096));
    
    static const weight_set_type::feature_type __empty;
    
    weight_set_type::const_iterator witer_begin = weights.begin();
    weight_set_type::const_iterator witer_end = weights.end();
    
    for (weight_set_type::const_iterator witer = witer_begin; witer != witer_end; ++ witer)
      if (*witer != 0.0) {
	const weight_set_type::feature_type feature(witer - witer_begin);
	if (feature != __empty) {
	  os << feature << ' ';
	  utils::encode_base64(*witer, std::ostream_iterator<char>(os));
	  os << '\n';
	}
      }
  } else {
    weights.clear();
    weights.allocate();
    
    boost::iostreams::filtering_istream is;
    is.push(boost::iostreams::zlib_decompressor());
    is.push(utils::mpi_device_bcast_source(rank, 4096));
    
    std::string feature;
    std::string value;
    
    while ((is >> feature) && (is >> value))
      weights[feature] = utils::decode_base64<double>(value);
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
    
    ("max-iteration", po::value<int>(&max_iteration), "# of hill-climbing iteration")
    ("min-iteration", po::value<int>(&min_iteration), "# of hill-climbing iteration")
    
    ("apply-exact", po::bool_switch(&apply_exact), "exact application")
    ("cube-size", po::value<int>(&cube_size)->default_value(cube_size),    "cube pruning size")
    ("beam-size", po::value<double>(&beam_size)->default_value(beam_size), "beam pruning size")
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
