//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// kbest variant of forest-oracle computer

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

#include "cicada_text_impl.hpp"
#include "cicada_kbeset_impl.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

struct real_precision20 : boost::spirit::karma::real_policies<double>
{
  static unsigned int precision(double) 
  { 
    return 20;
  }
};

path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";

bool directory_mode = false;

std::string scorer_name = "bleu:order=4,exact=true";

int max_iteration = 10;
int min_iteration = 1;

int debug = 0;

void read_refset(const path_set_type& file,
		 scorer_document_type& scorers);
void read_tstset(const path_set_type& files,
		 hypothesis_map_type& hypotheses);
void initialize_score(hypothesis_map_type& hypotheses,
		      const scorer_document_type& scorers);

template <typename Generator>
void compute_oracles(const scorer_document_type& scorers,
		     const hypothesis_map_type& hypotheses,
		     hypothesis_map_type& oracles,
		     Generator& generator);

void bcast_kbest(hypothesis_map_type& kbests);
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
    
    read_refset(refset_files, scorers);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;
    
    if (mpi_rank == 0 && debug)
      std::cerr << "reading tstset" << std::endl;
    
    hypothesis_map_type hypotheses(scorers.size());
    read_tstset(tstset_files, hypotheses);
    
    initialize_score(hypotheses, scorers);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;

    boost::mt19937 generator;
    generator.seed(time(0) * getpid());
    
    hypothesis_map_type oracles(scorers.size());
    compute_oracles(scorers, hypotheses, oracles, generator);

    boost::spirit::karma::real_generator<double, real_precision20> double20;
    
    if (mpi_rank == 0) {
      if (directory_mode) {
	if (boost::filesystem::exists(output_file) && ! boost::filesystem::is_directory(output_file))
	  boost::filesystem::remove_all(output_file);
	
	boost::filesystem::create_directories(output_file);
	
	boost::filesystem::directory_iterator iter_end;
	for (boost::filesystem::directory_iterator iter(output_file); iter != iter_end; ++ iter)
	  boost::filesystem::remove_all(*iter);
	
	for (size_t id = 0; id != oracles.size(); ++ id)
	  if (! oracles[id].empty()) {
	    namespace karma = boost::spirit::karma;
	    namespace standard = boost::spirit::standard;
	    
	    utils::compress_ostream os(output_file / (utils::lexical_cast<std::string>(id) + ".gz"), 1024 * 1024);
	    os.precision(10);
	    
	    hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	    for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	      const hypothesis_type& hyp(*oiter);
	      
	      os << id << " ||| ";
	    
	      if (! karma::generate(std::ostream_iterator<char>(os), -(standard::string % ' '), hyp.sentence))
		throw std::runtime_error("tokens generation failed...?");
	      os << " ||| ";
	      if (! karma::generate(std::ostream_iterator<char>(os), -((standard::string << '=' << double20) % ' '), hyp.features))
		throw std::runtime_error("tokens generation failed...?");
	      os << '\n';
	    }
	  }
	
      } else {
	utils::compress_ostream os(output_file, 1024 * 1024);
	os.precision(10);
	
	for (size_t id = 0; id != oracles.size(); ++ id)
	  if (! oracles[id].empty()) {
	    namespace karma = boost::spirit::karma;
	    namespace standard = boost::spirit::standard;
	    
	    hypothesis_set_type::const_iterator oiter_end = oracles[id].end();
	    for (hypothesis_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	      const hypothesis_type& hyp(*oiter);
	      
	      os << id << " ||| ";
	      
	      if (! karma::generate(std::ostream_iterator<char>(os), -(standard::string % ' '), hyp.sentence))
		throw std::runtime_error("tokens generation failed...?");
	      os << " ||| ";
	      if (! karma::generate(std::ostream_iterator<char>(os), -((standard::string << '=' << double20) % ' '), hyp.features))
		throw std::runtime_error("tokens generation failed...?");
	      os << '\n';
	    }
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
  typedef std::vector<const hypothesis_type*, std::allocator<const hypothesis_type*> > oracle_set_type;
  typedef std::vector<oracle_set_type, std::allocator<oracle_set_type> > oracle_map_type;
  
  TaskOracle(const scorer_document_type& __scorers,
	     const hypothesis_map_type&  __hypotheses,
	     oracle_map_type&            __oracles,
	     Generator&                  __generator)
    : scorers(__scorers),
      hypotheses(__hypotheses),
      oracles(__oracles),
      generator(__generator)
  {
    score.reset();
    
    for (size_t id = 0; id != oracles.size(); ++ id) 
      if (! oracles[id].empty()) {
	if (! score)
	  score = oracles[id].front()->score->clone();
	else
	  *score += *(oracles[id].front()->score);
      }
  }
  
  void operator()()
  {
    // we will try maximize    
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);

    double objective_optimum = (score
				? score->score() * score_factor
				: - std::numeric_limits<double>::infinity());
    
    for (size_t id = 0; id != hypotheses.size(); ++ id) 
      if (! hypotheses[id].empty()) {
	
	score_ptr_type score_curr = (score ? score->clone() : score_ptr_type());
	
	if (score_curr && ! oracles[id].empty())
	  *score_curr -= *(oracles[id].front()->score);
	
	oracles[id].clear();
	
	hypothesis_set_type::const_iterator hiter_end = hypotheses[id].end();
	for (hypothesis_set_type::const_iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter) {
	
	  score_ptr_type score_sample;

	  if (score_curr) {
	    score_sample = score_curr->clone();
	    *score_sample += *hiter->score;
	  } else
	    score_sample = hiter->score->clone();
	  
	  const double objective = score_sample->score() * score_factor;
	  
	  if (objective > objective_optimum || oracles[id].empty()) {
	    oracles[id].clear();
	    oracles[id].push_back(&(*hiter));
	    
	    objective_optimum = objective;
	    score = score_sample;
	  } else if (objective == objective_optimum)
	    oracles[id].push_back(&(*hiter));
	}
      }
  }
  
  score_ptr_type score;
  
  const scorer_document_type& scorers;
  const hypothesis_map_type& hypotheses;
  oracle_map_type& oracles;
  
  Generator&                           generator;
};

template <typename Generator>
void compute_oracles(const scorer_document_type& scorers,
		     const hypothesis_map_type& hypotheses,
		     hypothesis_map_type& oracles,
		     Generator& generator)
{
  typedef TaskOracle<Generator> task_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  score_ptr_type score_optimum;
  double objective_optimum = - std::numeric_limits<double>::infinity();

  sentence_set_type   sentences_optimum;
  hypergraph_set_type forests_optimum;
  
  const bool error_metric = scorers.error_metric();
  const double score_factor = (error_metric ? - 1.0 : 1.0);
  
  for (int iter = 0; iter < max_iteration; ++ iter) {
    if (debug && mpi_rank == 0)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    sentences_optimum = sentences;
    forests_optimum   = forests;
    
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
    
    int terminate = (objective <= objective_optimum) && (iter >= min_iteration);
    MPI::COMM_WORLD.Bcast(&terminate, 1, MPI::INT, 0);
    
    if (terminate) {
      sentences.swap(sentences_optimum);
      forests.swap(forests_optimum);
      break;
    }
    
    objective_optimum = objective;
  }
  
  for (size_t id = 0; id != graphs.size(); ++ id)
    if (features[id]) {
      if (! scores[id])
	throw std::runtime_error("no scores?");
      
      score_ptr_type score_curr = score_optimum->clone();
      *score_curr -= *scores[id];
      
      cicada::feature::Bleu*       __bleu = dynamic_cast<cicada::feature::Bleu*>(features[id].get());
      cicada::feature::BleuLinear* __bleu_linear = dynamic_cast<cicada::feature::BleuLinear*>(features[id].get());
      
      if (__bleu)
	__bleu->assign(score_curr);
      else
	__bleu_linear->assign(score_curr);
    }
}

void initialize_score(hypothesis_map_type& hypotheses,
		      const scorer_document_type& scorers)
{
  for (size_t id = 0; id != hypotheses.size(); ++ id)
    if (! hypotheses[id].empty()) {
      hypothesis_set_type(hypotheses[id]).swap(hypotheses[id]);
      
      hypothesis_set_type::iterator hiter_end = hypotheses[id].end();
      for (hypothesis_set_type::iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter)
	hiter->score = scorers[id]->score(sentence_type(hiter->sentence.begin(), hiter->sentence.end()));
    }
}

void read_tstset(const path_set_type& files,
		 hypothesis_map_type& hypotheses)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  if (files.empty())
    throw std::runtime_error("no files?");

  parser_type parser;
  kbest_feature_type kbest;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no file: " + fiter->string());

    if (boost::filesystem::is_directory(*fiter)) {
      for (int i = mpi_rank; /**/; i += mpi_size) {
	const path_type path = (*fiter) / (utils::lexical_cast<std::string>(i) + ".gz");

	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iter_type iter(is);
	iter_type iter_end;
	
	while (iter != iter_end) {
	  boost::fusion::get<1>(kbest).clear();
	  boost::fusion::get<2>(kbest).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	    if (iter != iter_end)
	      throw std::runtime_error("kbest parsing failed");
	  
	  const size_t& id = boost::fusion::get<0>(kbest);
	  
	  if (id >= hypotheses.size())
	    throw std::runtime_error("invalid id: " + utils::lexical_cast<std::string>(id));
	  if (id != i)
	    throw std::runtime_error("invalid id: " + utils::lexical_cast<std::string>(id));
	  
	  hypotheses[id].push_back(hypothesis_type(kbest));
	}
      }
    } else {
      utils::compress_istream is(*fiter, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iter_type iter(is);
      iter_type iter_end;
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
	
	const size_t& id = boost::fusion::get<0>(kbest);
	
	if (id >= hypotheses.size())
	  throw std::runtime_error("invalid id: " + utils::lexical_cast<std::string>(id));
	
	if (id % mpi_size == mpi_rank)
	  hypotheses[id].push_back(hypothesis_type(kbest));
      }
    }
  }
}


void read_refset(const path_set_type& files,
		 scorer_document_type& scorers)
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
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      scorers[id]->insert(id_sentence.second);
    }
  }
}

void bcast_kbest(hypothesis_map_type& kbests)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    if (rank == mpi_rank) {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::gzip_compressor());
      os.push(utils::mpi_device_bcast_sink(rank, 4096));
      
      for (size_t id = 0; id != sentences.size(); ++ id)
	if (static_cast<int>(id % mpi_size) == mpi_rank)
	  os << id << " ||| " << sentences[id] << " ||| " << forests[id] << '\n';
      
    } else {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::gzip_decompressor());
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
    os.push(utils::mpi_device_bcast_sink(rank, 1024));
    
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
    is.push(utils::mpi_device_bcast_source(rank, 1024));
    
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
    ("cube-size", po::value<int>(&cube_size), "cube-pruning size")
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
