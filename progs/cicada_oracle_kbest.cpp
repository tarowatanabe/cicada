//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// refset format:
// 0 |||  reference translatin for source sentence 0
// 0 |||  another reference
// 1 |||  reference translation for source sentence 1
//


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

#include "cicada_text_impl.hpp"
#include "cicada_kbest_impl.hpp"

typedef std::vector<const hypothesis_type*, std::allocator<const hypothesis_type*> > oracle_set_type;
typedef std::vector<oracle_set_type, std::allocator<oracle_set_type> > oracle_map_type;

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

int threads = 4;

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
		     oracle_map_type& oracles,
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
    
    read_refset(refset_files, scorers);
    
    if (debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;
    
    // read test set
    if (tstset_files.empty())
      tstset_files.push_back("-");

    if (debug)
      std::cerr << "reading tstset" << std::endl;

    hypothesis_map_type hypotheses(scorers.size());
    read_tstset(tstset_files, hypotheses);

    initialize_score(hypotheses, scorers);
    
    boost::mt19937 generator;
    generator.seed(time(0) * getpid());
    
    oracle_map_type oracles(scorers.size());
    compute_oracles(scorers, hypotheses, oracles, generator);
    
    boost::spirit::karma::real_generator<double, real_precision20> double20;
    
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
	  
	  oracle_set_type::const_iterator oiter_end = oracles[id].end();
	  for (oracle_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	    const hypothesis_type& hyp(*(*oiter));
	    
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
	  
	  oracle_set_type::const_iterator oiter_end = oracles[id].end();
	  for (oracle_set_type::const_iterator oiter = oracles[id].begin(); oiter != oiter_end; ++ oiter) {
	    const hypothesis_type& hyp(*(*oiter));
	    
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
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct TaskOracle
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;

  TaskOracle(queue_type&                 __queue,
	     const scorer_document_type& __scorers,
	     const hypothesis_map_type&  __hypotheses,
	     oracle_map_type&            __oracles,
	     const score_ptr_type&       __score)
    : queue(__queue), scorers(__scorers), hypotheses(__hypotheses), oracles(__oracles), score(__score ? __score->clone() : score_ptr_type()) {}
  
  void operator()()
  {
    const bool error_metric = scorers.error_metric();
    const double score_factor = (error_metric ? - 1.0 : 1.0);
    
    double objective_optimum = (score
				? score->score() * score_factor
				: - std::numeric_limits<double>::infinity());
    for (;;) {
      int id = 0;
      queue.pop(id);
      if (id < 0) break;
      
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
  
  queue_type& queue;
  const scorer_document_type& scorers;
  const hypothesis_map_type& hypotheses;
  oracle_map_type& oracles;
  score_ptr_type score;
};

template <typename Generator>
void compute_oracles(const scorer_document_type& scorers,
		     const hypothesis_map_type& hypotheses,
		     oracle_map_type& oracles,
		     Generator& generator)
{
  typedef TaskOracle task_type;
  typedef task_type::queue_type queue_type;
  
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  typedef std::vector<int, std::allocator<int> > id_set_type;
  
  id_set_type ids;
  for (size_t id = 0; id != scorers.size(); ++ id)
    if (! hypotheses[id].empty())
      ids.push_back(id);
  
  score_ptr_type  score_optimum;
  double          objective_optimum = - std::numeric_limits<double>::infinity();
  oracle_map_type oracles_optimum(oracles.size());
  
  const bool error_metric = scorers.error_metric();
  const double score_factor = (error_metric ? - 1.0 : 1.0);
  
  for (int iter = 0; iter < max_iteration; ++ iter) {
    if (debug)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    oracles_optimum = oracles;
    
    queue_type queue;
    task_set_type tasks(threads, task_type(queue, scorers, hypotheses, oracles, score_optimum));
    
    boost::thread_group workers;
    for (int i = 0; i < threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    id_set_type::const_iterator iiter_end = ids.end();
    for (id_set_type::const_iterator iiter = ids.begin(); iiter != iiter_end; ++ iiter)
      queue.push(*iiter);
    
    for (int i = 0; i < threads; ++ i)
      queue.push(-1);
    
    boost::random_number_generator<Generator> gen(generator);
    std::random_shuffle(ids.begin(), ids.end(), gen);
    
    workers.join_all();
    
    score_optimum.reset();
    oracle_map_type::const_iterator oiter_end = oracles.end();
    for (oracle_map_type::const_iterator oiter = oracles.begin(); oiter != oiter_end; ++ oiter)
      if (! oiter->empty()) {
	if (! score_optimum)
	  score_optimum = oiter->front()->score->clone();
	else
	  *score_optimum = *(oiter->front()->score);
      }
    
    const double objective = score_optimum->score() * score_factor;
    if (debug)
      std::cerr << "oracle score: " << objective << std::endl;
    
    if (objective <= objective_optimum && iter >= min_iteration) {
      oracles.swap(oracles_optimum);
      break;
    }
    
    objective_optimum = objective;
  }
}

struct TaskInit
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
  TaskInit(queue_type&                 __queue,
	   hypothesis_map_type&        __hypotheses,
	   const scorer_document_type& __scorers)
    : queue(__queue), hypotheses(__hypotheses), scorers(__scorers) {}

  void operator()()
  {
    for (;;) {
      int id = 0;
      queue.pop(id);
      if (id < 0) break;
      
      hypothesis_set_type(hypotheses[id]).swap(hypotheses[id]);
      
      hypothesis_set_type::iterator hiter_end = hypotheses[id].end();
      for (hypothesis_set_type::iterator hiter = hypotheses[id].begin(); hiter != hiter_end; ++ hiter)
	hiter->score = scorers[id]->score(sentence_type(hiter->sentence.begin(), hiter->sentence.end()));
    }
  }

  queue_type&                 queue;
  hypothesis_map_type&        hypotheses;
  const scorer_document_type& scorers;
};

void initialize_score(hypothesis_map_type& hypotheses,
		      const scorer_document_type& scorers)
{
  typedef TaskInit task_type;
  typedef task_type::queue_type queue_type;

  queue_type queue;
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, hypotheses, scorers)));
  
  for (int id = 0; id != hypotheses.size(); ++ id)
    queue.push(id);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
}

void read_tstset(const path_set_type& files,
		 hypothesis_map_type& hypotheses)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;
  
  if (files.empty())
    throw std::runtime_error("no files?");

  parser_type parser;
  kbest_feature_type kbest;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    if (debug)
      std::cerr << "file: " << *fiter << std::endl;

    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no file: " + fiter->string());

    if (boost::filesystem::is_directory(*fiter)) {
      for (int i = 0; /**/; ++ i) {
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

void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s) (in kbest format)")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("directory", po::bool_switch(&directory_mode), "output in directory")
        
    ("scorer",    po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    
    ("max-iteration", po::value<int>(&max_iteration), "# of hill-climbing iteration")
    ("min-iteration", po::value<int>(&min_iteration), "# of hill-climbing iteration")
    
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
