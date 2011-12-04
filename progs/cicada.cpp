//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unistd.h>

#include "cicada_impl.hpp"

#include "utils/program_options.hpp"
#include "utils/filesystem.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/filesystem.hpp"
#include "utils/bithack.hpp"

#include <boost/program_options.hpp>
#include <boost/thread.hpp>

typedef std::string op_type;
typedef std::vector<op_type, std::allocator<op_type> > op_set_type;

path_type input_file = "-";

bool input_id_mode = false;
bool input_bitext_mode = false;
bool input_sentence_mode = false;
bool input_lattice_mode = false;
bool input_forest_mode = false;
bool input_span_mode = false;
bool input_alignment_mode = false;
bool input_dependency_mode = false;
bool input_directory_mode = false;

std::string symbol_goal         = vocab_type::S;

grammar_file_set_type grammar_files;
bool grammar_list = false;

grammar_file_set_type tree_grammar_files;
bool tree_grammar_list = false;

feature_parameter_set_type feature_parameters;
bool feature_list = false;

op_set_type ops;
bool op_list = false;

int threads = 1;

int debug = 0;

void cicada_file(const operation_set_type& operations,
		 const model_type& model,
		 const grammar_type& grammar,
		 const tree_grammar_type& tree_grammar,
		 operation_set_type::statistics_type& stats);
void cicada_directory(const operation_set_type& operations,
		      const model_type& model,
		      const grammar_type& grammar,
		      const tree_grammar_type& tree_grammar,
		      operation_set_type::statistics_type& stats);

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (feature_list) {
      std::cout << cicada::FeatureFunction::lists();
      return 0;
    }

    if (op_list) {
      std::cout << operation_set_type::lists();
      return 0;
    }
    
    if (grammar_list) {
      std::cout << grammar_type::lists();
      return 0;
    }

    if (tree_grammar_list) {
      std::cout << tree_grammar_type::lists();
      return 0;
    }
    
    threads = utils::bithack::max(1, threads);

    // read grammars...
    grammar_type grammar(grammar_files.begin(), grammar_files.end());
    if (debug)
      std::cerr << "grammar: " << grammar.size() << std::endl;
    
    tree_grammar_type tree_grammar(tree_grammar_files.begin(), tree_grammar_files.end());
    if (debug)
      std::cerr << "tree grammar: " << tree_grammar.size() << std::endl;
    
    // read features...
    model_type model;
    for (feature_parameter_set_type::const_iterator piter = feature_parameters.begin(); piter != feature_parameters.end(); ++ piter)
      model.push_back(feature_function_type::create(*piter));
    model.initialize();


    operation_set_type operations(ops.begin(), ops.end(),
				  model,
				  grammar,
				  tree_grammar,
				  symbol_goal,
				  input_id_mode || input_directory_mode,
				  input_sentence_mode,
				  input_lattice_mode,
				  input_forest_mode,
				  input_span_mode,
				  input_alignment_mode,
				  input_dependency_mode,
				  input_bitext_mode,
				  false,
				  debug);

    if (! operations.get_output_data().directory.empty()) {
      const path_type& directory = operations.get_output_data().directory;
      
      if (boost::filesystem::exists(directory) && ! boost::filesystem::is_directory(directory))
	utils::filesystem::remove_all(directory);
      
      boost::filesystem::create_directories(directory);
      
      boost::filesystem::directory_iterator iter_end;
      for (boost::filesystem::directory_iterator iter(directory); iter != iter_end; ++ iter)
	utils::filesystem::remove_all(*iter);
      
      ::sync();
    }

    operation_set_type::statistics_type stats;
    
    if (! operations.get_output_data().file.empty())
      cicada_file(operations, model, grammar, tree_grammar, stats);
    else
      cicada_directory(operations, model, grammar, tree_grammar, stats);
    
    if (debug)
      std::cerr << "statistics"<< '\n'
		<< stats;
    
#if 0
    // we will force non directory-input-mode....
    if (input_directory_mode) {
      std::string line;
      
      for (size_t i = 0; /**/; ++ i) {
	const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
	
	const path_type path_input = input_file / file_name;
	
	if (! boost::filesystem::exists(path_input)) break;
	
	utils::compress_istream is(path_input, 1024 * 1024);
	
	if (std::getline(is, line)) {
	  operations(line);
	  operations.clear();
	}
      }
      
    } else {
      utils::compress_istream is(input_file, 1024 * 1024);
      
      std::string line;
      while (std::getline(is, line)) {
	operations(line);
	operations.clear();
      }
    }
    
    if (debug)
      std::cerr << "statistics"<< '\n'
		<< operations.get_statistics();
#endif
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct TaskFile
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;
  
  TaskFile(queue_type&   __queue_is,
	   queue_type&   __queue_os,
	   const model_type& __model,
	   const grammar_type& __grammar,
	   const tree_grammar_type& __tree_grammar)
    : queue_is(__queue_is),
      queue_os(__queue_os),
      _model(__model),
      _grammar(__grammar),
      _tree_grammar(__tree_grammar) {}
  
  void operator()()
  {
    // cloning should be performed in thread... otherwise, strangething may happen
    const model_type        model(_model.clone());
    const grammar_type      grammar(_grammar.clone());
    const tree_grammar_type tree_grammar(_tree_grammar.clone());
    
    operation_set_type operations(ops.begin(), ops.end(),
				  model,
				  grammar,
				  tree_grammar,
				  symbol_goal,
				  true,
				  input_sentence_mode,
				  input_lattice_mode,
				  input_forest_mode,
				  input_span_mode,
				  input_alignment_mode,
				  input_dependency_mode,
				  input_bitext_mode,
				  true,
				  debug);
    
    std::string line;
    while (1) {
      queue_is.pop_swap(line);
      if (line.empty()) break;
      
      if (input_directory_mode) {
	utils::compress_istream is(line, 1024 * 1024);
	
	if (std::getline(is, line) && ! line.empty())
	  operations(line);
      } else
	operations(line);
      
      queue_os.push(utils::lexical_cast<std::string>(operations.get_data().id) + ' ' + operations.get_output_data().buffer);
    }
    
    operations.clear();

    stats = operations.get_statistics();
  }
  
  queue_type&   queue_is;
  queue_type&   queue_os;
  const model_type& _model;
  const grammar_type& _grammar;
  const tree_grammar_type& _tree_grammar;  
  
  operation_set_type::statistics_type stats;
};

struct ReduceFile
{
  typedef TaskFile::queue_type queue_type;
  
  ReduceFile(queue_type& __queue, const path_type& __path)
    : queue(__queue), path(__path) {}
  
  void operator()()
  {
    typedef size_t id_type;
    typedef std::map<id_type, std::string, std::less<id_type>, std::allocator<std::pair<const id_type, std::string> > > buffer_map_type;
    
    buffer_map_type maps;
    std::string buffer;
    
    id_type     id = 0;
    
    utils::compress_ostream os(path, 1024 * 1024);
    
    for (;;) {
      queue.pop_swap(buffer);
      
      if (buffer.empty()) break;

      bool dump = false;
      
      utils::piece buffer_piece(buffer);
      
      utils::piece::const_iterator iter = buffer_piece.begin();
      for (/**/; iter != buffer_piece.end() && ! std::isspace(*iter); ++ iter);
      
      // tokenize here...
      const id_type      buffer_id        = utils::lexical_cast<id_type>(buffer_piece.substr(0, iter - buffer_piece.begin()));
      const utils::piece buffer_tokenized = buffer_piece.substr(iter + 1 - buffer_piece.begin());
      
      if (buffer_id == id) {
	os << buffer_tokenized;
	dump = true;
	
	++ id;
      } else
	maps[buffer_id] = static_cast<std::string>(buffer_tokenized);
      
      for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; ++ id) {
	os << iter->second;
	dump = true;
	maps.erase(iter ++);
      }
      
      if (dump)
	os << std::flush;
    }
    
    for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; ++ id) {
      os << iter->second;
      maps.erase(iter ++);
    }
    
    // we will do twice, in case we have wrap-around for id...!
    if (! maps.empty())
      for (buffer_map_type::iterator iter = maps.find(id); iter != maps.end() && iter->first == id; ++ id) {
	os << iter->second;
	maps.erase(iter ++);
      }
    
    os << std::flush;
    
    if (! maps.empty())
      throw std::runtime_error("id mismatch! expecting: " + utils::lexical_cast<std::string>(id)
			       + " next: " + utils::lexical_cast<std::string>(maps.begin()->first));
  }
  
  queue_type& queue;
  path_type   path;
};

struct TaskDirectory
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;
  
  TaskDirectory(queue_type&   __queue,
		const model_type& __model,
		const grammar_type& __grammar,
		const tree_grammar_type& __tree_grammar)
    : queue(__queue),
      _model(__model),
      _grammar(__grammar),
      _tree_grammar(__tree_grammar) {}
  
  void operator()()
  {
    // cloning should be performed in thread... otherwise, strangething may happen
    const model_type        model(_model.clone());
    const grammar_type      grammar(_grammar.clone());
    const tree_grammar_type tree_grammar(_tree_grammar.clone());
    
    operation_set_type operations(ops.begin(), ops.end(),
				  model,
				  grammar,
				  tree_grammar,
				  symbol_goal,
				  true,
				  input_sentence_mode,
				  input_lattice_mode,
				  input_forest_mode,
				  input_span_mode,
				  input_alignment_mode,
				  input_dependency_mode,
				  input_bitext_mode,
				  true,
				  debug);
    
    std::string line;
    
    while (1) {
      queue.pop_swap(line);
      if (line.empty()) break;

      if (input_directory_mode) {
	utils::compress_istream is(line, 1024 * 1024);
	
	if (std::getline(is, line) && ! line.empty())
	  operations(line);
      } else 
	operations(line);
    }
    
    operations.clear();

    stats = operations.get_statistics();
  }
  
  queue_type&   queue;
  const model_type& _model;
  const grammar_type& _grammar;
  const tree_grammar_type& _tree_grammar;

  operation_set_type::statistics_type stats;
};


void cicada_file(const operation_set_type& operations,
		 const model_type& model,
		 const grammar_type& grammar,
		 const tree_grammar_type& tree_grammar,
		 operation_set_type::statistics_type& stats)
{
  typedef TaskFile   task_type;
  typedef ReduceFile reducer_type;
  
  task_type::queue_type queue_is(threads);
  task_type::queue_type queue_os;
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(queue_os, operations.get_output_data().file)));
  
  boost::thread_group mapper;
  std::vector<task_type, std::allocator<task_type> > tasks(threads, task_type(queue_is, queue_os, model, grammar, tree_grammar));
  
  for (int i = 0; i != threads; ++ i)
    mapper.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  if (input_directory_mode) {
    std::string line;
    
    for (size_t i = 0; /**/; ++ i) {
      const std::string file_name = utils::lexical_cast<std::string>(i) + ".gz";
      
      const path_type path_input = input_file / file_name;
      
      if (! boost::filesystem::exists(path_input)) break;
      
      queue_is.push(path_input.string());
    }
    
  } else {
    utils::compress_istream is(input_file, 1024 * 1024);
    
    size_t id = 0;
    std::string line;
    
    while (std::getline(is, line)) {

      if (! line.empty()) {
	if (input_id_mode)
	  queue_is.push_swap(line);
	else
	  queue_is.push(utils::lexical_cast<std::string>(id) + " ||| " + line);
      }
      
      ++ id;
    }
  }
  
  for (int i = 0; i != threads; ++ i)
    queue_is.push(std::string());
  
  mapper.join_all();
  
  queue_os.push(std::string());
  reducer.join_all();

  for (int i = 0; i != threads; ++ i)
    stats += tasks[i].stats;
}


void cicada_directory(const operation_set_type& operations,
		      const model_type& model,
		      const grammar_type& grammar,
		      const tree_grammar_type& tree_grammar,
		      operation_set_type::statistics_type& stats)
{
  typedef TaskDirectory task_type;
  
  task_type::queue_type queue(threads);
  
  boost::thread_group mapper;
  std::vector<task_type, std::allocator<task_type> > tasks(threads, task_type(queue, model, grammar, tree_grammar));
  
  for (int i = 0; i != threads; ++ i)
    mapper.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  if (input_directory_mode) {
    boost::filesystem::directory_iterator iter_end;
    for (boost::filesystem::directory_iterator iter(input_file); iter != iter_end; ++ iter) {
      const std::string file = path_type(*iter).string();
      
      if (! file.empty())
	queue.push(file);
    }
  } else {
    utils::compress_istream is(input_file, 1024 * 1024);
    
    size_t id = 0;
    std::string line;
    
    while (std::getline(is, line)) {

      if (! line.empty()) {
	if (input_id_mode)
	  queue.push_swap(line);
	else
	  queue.push(utils::lexical_cast<std::string>(id) + " ||| " + line);
      }
      
      ++ id;
    }
  }
  
  for (int i = 0; i != threads; ++ i)
    queue.push(std::string());
  
  mapper.join_all();

  for (int i = 0; i != threads; ++ i)
    stats += tasks[i].stats;
}

struct deprecated
{
  deprecated(const boost::program_options::options_description& __desc)
    : desc(__desc) {}
  
  template <typename Tp>
  void operator()(const Tp& x) const
  {
    std::cout << desc << std::endl;
    exit(1);
  }
  
  const boost::program_options::options_description& desc;
};

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    
    // options for input/output format
    ("input-id",         po::bool_switch(&input_id_mode),         "id-prefixed input")
    ("input-bitext",     po::bool_switch(&input_bitext_mode),     "target sentence prefixed input")
    ("input-sentence",   po::bool_switch(&input_sentence_mode),   "sentence input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    ("input-forest",     po::bool_switch(&input_forest_mode),     "forest input")
    ("input-span",       po::bool_switch(&input_span_mode),       "span input")
    ("input-alignment",  po::bool_switch(&input_alignment_mode),  "alignment input")
    ("input-dependency", po::bool_switch(&input_dependency_mode), "dependency input")
    ("input-directory",  po::bool_switch(&input_directory_mode),  "input in directory")
    
    // grammar
    ("goal",              po::value<std::string>(&symbol_goal)->default_value(symbol_goal),    "goal symbol")
    ("grammar",           po::value<grammar_file_set_type >(&grammar_files)->composing(),      "grammar specification(s)")
    ("grammar-list",      po::bool_switch(&grammar_list),                                      "list of available grammar specifications")
    ("tree-grammar",      po::value<grammar_file_set_type >(&tree_grammar_files)->composing(), "tree grammar specification(s)")
    ("tree-grammar-list", po::bool_switch(&tree_grammar_list),                                 "list of available grammar specifications")
    
    // models...
    ("feature-function",      po::value<feature_parameter_set_type >(&feature_parameters)->composing(), "feature function(s)")
    ("feature-function-list", po::bool_switch(&feature_list),                                           "list of available feature function(s)")

    //operatins...
    ("operation",      po::value<op_set_type>(&ops)->composing(), "operations")
    ("operation-list", po::bool_switch(&op_list),                 "list of available operation(s)");

  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config",  po::value<path_type>(),                    "configuration file")
    ("threads", po::value<int>(&threads),                  "# of threads (highly experimental)")
    ("debug",   po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::options_description opts_deprecated("deprecated options");
  opts_deprecated.add_options()
    ("non-terminal",          po::value<std::string>()->notifier(deprecated(opts_deprecated)), "see --grammar-list")
    ("grammar-static",        po::value<grammar_file_set_type >()->composing()->notifier(deprecated(opts_deprecated)), "use --grammar ")
    ("grammar-glue-straight", po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar glue:straight=[true|false],inverted=[true|false],non-terminal=[x]")
    ("grammar-glue-inverted", po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar glue:straight=[true|false],inverted=[true|false],non-terminal=[x]")
    ("grammar-insertion",     po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar insetion:non-terminal=[x]")
    ("grammar-deletion",      po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --grammar deletion:non-terminal=[x]")
    ("tree-grammar-static",   po::value<grammar_file_set_type >()->composing()->notifier(deprecated(opts_deprecated)),  "use --tree-grammar")
    ("tree-grammar-fallback", po::value<bool>()->notifier(deprecated(opts_deprecated))->zero_tokens(), "use --tree-grammar fallback:non-terminal=[x]");

  po::options_description desc_config;
  po::options_description desc_command;
  po::options_description desc_visible;
  
  desc_config.add(opts_config).add(opts_deprecated);
  desc_command.add(opts_config).add(opts_command).add(opts_deprecated);
  desc_visible.add(opts_config).add(opts_command);
  
  po::variables_map variables;

  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  if (variables.count("config")) {
    const path_type path_config = variables["config"].as<path_type>();
    if (! boost::filesystem::exists(path_config))
      throw std::runtime_error("no config file: " + path_config.string());
    
    utils::compress_istream is(path_config);
    po::store(po::parse_config_file(is, desc_config), variables);
  }
  
  po::notify(variables);

  if (variables.count("help")) {
    
    std::cout << argv[0] << " [options]\n"
	      << desc_visible << std::endl;
    exit(0);
  }
}
