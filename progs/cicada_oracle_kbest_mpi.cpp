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
#include "utils/sgi_hash_set.hpp"
#include "utils/random_seed.hpp"

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
#include "cicada_kbest_impl.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

struct real_precision_inf : boost::spirit::karma::real_policies<double>
{
  static unsigned int precision(double) 
  { 
    return std::numeric_limits<double>::digits10 + 1;
  }
};

struct real_precision_20 : boost::spirit::karma::real_policies<double>
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
int min_iteration = 5;

int debug = 0;

void read_refset(const path_set_type& file,
		 scorer_document_type& scorers);
void read_tstset(const path_set_type& files,
		 hypothesis_map_type& hypotheses);
void initialize_score(hypothesis_map_type& hypotheses,
		      const scorer_document_type& scorers);

template <typename Generator>
double compute_oracles(const scorer_document_type& scorers,
		       const hypothesis_map_type& hypotheses,
		       hypothesis_map_type& oracles,
		       Generator& generator);

void bcast_kbest(hypothesis_map_type& kbests);
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
    
    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    hypothesis_map_type oracles(scorers.size());
    const double objective = compute_oracles(scorers, hypotheses, oracles, generator);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "oracle score: " << objective << std::endl;

    boost::spirit::karma::real_generator<double, real_precision_20> double20;
    
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
	      if (! hyp.features.empty()) {
		os << " ||| ";
		if (! karma::generate(std::ostream_iterator<char>(os), -((standard::string << '=' << double20) % ' '), hyp.features))
		  throw std::runtime_error("tokens generation failed...?");
	      }
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
	      if (! hyp.features.empty()) {
		os << " ||| ";
		if (! karma::generate(std::ostream_iterator<char>(os), -((standard::string << '=' << double20) % ' '), hyp.features))
		  throw std::runtime_error("tokens generation failed...?");
	      }
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
  TaskOracle(const scorer_document_type& __scorers,
	     const hypothesis_map_type&  __hypotheses,
	     hypothesis_map_type&        __oracles,
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
	  score = oracles[id].front().score->clone();
	else
	  *score += *(oracles[id].front().score);
      }
  }
  
  void operator()()
  {
    typedef std::vector<size_t, std::allocator<size_t> > id_set_type;

    // we will try maximize    
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);

    double objective_optimum = (score
				? score->score() * score_factor
				: - std::numeric_limits<double>::infinity());
    
    id_set_type ids;
    for (size_t id = 0; id != hypotheses.size(); ++ id)
      if (! hypotheses[id].empty())
	ids.push_back(id);
    
    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);
    
    id_set_type::const_iterator iiter_end = ids.end();
    for (id_set_type::const_iterator iiter = ids.begin(); iiter != iiter_end; ++ iiter) {
      const size_t id = *iiter;
      
      score_ptr_type score_curr = (score ? score->clone() : score_ptr_type());
	
      if (score_curr && ! oracles[id].empty())
	*score_curr -= *(oracles[id].front().score);
      
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
	  oracles[id].push_back(*hiter);
	    
	  objective_optimum = objective;
	  score = score_sample;
	} else if (objective == objective_optimum)
	  oracles[id].push_back(*hiter);
      }
    }
  }
  
  score_ptr_type score;
  
  const scorer_document_type& scorers;
  const hypothesis_map_type&  hypotheses;
  hypothesis_map_type&        oracles;
  Generator&                  generator;
};

template <typename Generator>
double compute_oracles(const scorer_document_type& scorers,
		       const hypothesis_map_type& hypotheses,
		       hypothesis_map_type& oracles,
		       Generator& generator)
{
  typedef TaskOracle<Generator> task_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  // TODO!
  //
  // we need sentence+feature bcast to keep oracles consistent
  //
  
  score_ptr_type score_optimum;
  double objective_prev = - std::numeric_limits<double>::infinity();
  double objective_best = - std::numeric_limits<double>::infinity();
  hypothesis_map_type oracles_best(oracles.size());
  
  const bool error_metric = scorers.error_metric();
  const double score_factor = (error_metric ? - 1.0 : 1.0);
  
#if 1
  // initialize...
  for (size_t id = 0; id != hypotheses.size(); ++ id)
    if (! hypotheses[id].empty()) {
      double objective_best = - std::numeric_limits<double>::infinity();
      
      hypothesis_set_type::const_iterator hiter_end = hypotheses[id].end();
      for (hypothesis_set_type::const_iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter) {
	const double score = hiter->score->score() * score_factor;
	
	if (score > objective_best) {
	  oracles[id].clear();
	  oracles[id].push_back(*hiter);
	    
	  objective_best = score;
	} else if (score == objective_best)
	  oracles[id].push_back(*hiter);
      }
    }
    
  bcast_kbest(oracles);
    
  for (size_t id = 0; id != oracles.size(); ++ id) {
    hypothesis_set_type::iterator oiter_end = oracles[id].end();
    for (hypothesis_set_type::iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
      hypothesis_type& hyp = *oiter;
	
      if (! hyp.score)
	hyp.score = scorers[id]->score(sentence_type(hyp.sentence.begin(), hyp.sentence.end()));
    }
  }
#endif
  
  for (int iter = 0; iter < max_iteration; ++ iter) {
    if (debug && mpi_rank == 0)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    task_type(scorers, hypotheses, oracles, generator)();
    
    bcast_kbest(oracles);
    
    score_optimum.reset();
    for (size_t id = 0; id != oracles.size(); ++ id)
      if (! oracles[id].empty()) {
	hypothesis_type& hyp = oracles[id].front();

	if (! hyp.score)
	  hyp.score = scorers[id]->score(sentence_type(hyp.sentence.begin(), hyp.sentence.end()));
	
	if (! score_optimum)
	  score_optimum = hyp.score->clone();
	else
	  *score_optimum += *hyp.score;
      }
    
    const double objective = score_optimum->score() * score_factor;
    if (mpi_rank == 0 && debug)
      std::cerr << "oracle score: " << objective << std::endl;

    if (objective > objective_best) {
      objective_best = objective;
      oracles_best   = oracles;
    }
    
    int terminate = (objective <= objective_prev) && (iter >= min_iteration);
    MPI::COMM_WORLD.Bcast(&terminate, 1, MPI::INT, 0);
    
    if (terminate)
      break;
    
    objective_prev = objective;
  }
  
  oracles.swap(oracles_best);

  return objective_best;
}

void initialize_score(hypothesis_map_type& hypotheses,
		      const scorer_document_type& scorers)
{
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
    std::allocator<hypothesis_type> > hypothesis_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
			std::allocator<hypothesis_type> > hypothesis_unique_type;
#endif
  
  hypothesis_unique_type uniques;

  for (size_t id = 0; id != hypotheses.size(); ++ id)
    if (! hypotheses[id].empty()) {
      uniques.clear();
      uniques.insert(hypotheses[id].begin(), hypotheses[id].end());
      
      hypotheses[id].clear();
      hypothesis_set_type(hypotheses[id]).swap(hypotheses[id]);
      
      hypotheses[id].reserve(uniques.size());
      hypotheses[id].insert(hypotheses[id].end(), uniques.begin(), uniques.end());
      
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
    if (mpi_rank == 0 && debug)
      std::cerr << "file: " << *fiter << std::endl;
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no file: " + fiter->string());

    if (boost::filesystem::is_directory(*fiter)) {
      for (size_t i = mpi_rank; /**/; i += mpi_size) {
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
	
	if (static_cast<int>(id % mpi_size) == mpi_rank)
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
  typedef boost::spirit::istream_iterator  iiter_type;
  typedef std::ostream_iterator<char>      oiter_type;
  
  typedef kbest_feature_parser<iiter_type>    parser_type;
  typedef kbest_feature_generator<oiter_type> generator_type;
  
  namespace qi       = boost::spirit::qi;
  namespace karma    = boost::spirit::karma;
  namespace standard = boost::spirit::standard;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  parser_type    parser;
  generator_type generator;
  
  boost::spirit::karma::real_generator<double, real_precision_inf> double_inf;

  kbest_feature_type kbest;

  // clear non-rank kbests
  for (size_t id = 0; id != kbests.size(); ++ id)
    if (static_cast<int>(id % mpi_size) != mpi_rank)
      kbests[id].clear();
  
  for (int rank = 0; rank < mpi_size; ++ rank) {
    if (rank == mpi_rank) {
      boost::iostreams::filtering_ostream os;
      os.push(boost::iostreams::zlib_compressor());
      os.push(utils::mpi_device_bcast_sink(rank, 4096));
      
      for (size_t id = 0; id != kbests.size(); ++ id)
	if (static_cast<int>(id % mpi_size) == mpi_rank) {
	  
	  hypothesis_set_type::const_iterator oiter_end = kbests[id].end();
	  for (hypothesis_set_type::const_iterator oiter = kbests[id].begin(); oiter != oiter_end; ++ oiter)
	    if (! karma::generate(oiter_type(os),
				  (generator.size
				   << " ||| " << -(standard::string % ' ')
				   << " ||| " << -((standard::string << '=' << double_inf) % ' ')
				   << '\n'),
				  id,
				  oiter->sentence,
				  oiter->features))
	      throw std::runtime_error("error in generating kbest");
	}
      
    } else {
      boost::iostreams::filtering_istream is;
      is.push(boost::iostreams::zlib_decompressor());
      is.push(utils::mpi_device_bcast_source(rank, 4096));
      is.unsetf(std::ios::skipws);
      
      iiter_type iter(is);
      iiter_type iter_end;
      
      while (iter != iter_end) {
	boost::fusion::get<1>(kbest).clear();
	boost::fusion::get<2>(kbest).clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, kbest))
	  if (iter != iter_end)
	    throw std::runtime_error("kbest parsing failed");
	
	const size_t& id = boost::fusion::get<0>(kbest);
	
	kbests[id].push_back(hypothesis_type(kbest));
      }
    }
  }
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s) (in kbest format)")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("directory", po::bool_switch(&directory_mode), "output in directory")
    
    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    
    ("max-iteration", po::value<int>(&max_iteration), "# of hill-climbing iteration")
    ("min-iteration", po::value<int>(&min_iteration), "# of hill-climbing iteration");
  
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
