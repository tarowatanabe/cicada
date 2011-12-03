//
// query tree given an input hypergraph
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <unistd.h>
#include <set>
#include <iterator>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>
#include <boost/thread.hpp>


#include "cicada/lattice.hpp"
#include "cicada/sentence.hpp"
#include "cicada/hypergraph.hpp"
#include "cicada/query_cky.hpp"
#include "cicada/grammar.hpp"

#include "utils/lockfree_list_queue.hpp"
#include "utils/json_string_generator.hpp"
#include "utils/compress_stream.hpp"
#include "utils/filesystem.hpp"


typedef boost::filesystem::path path_type;

typedef std::vector<std::string, std::allocator<std::string> > grammar_file_set_type;

typedef cicada::Lattice  lattice_type;
typedef cicada::Sentence sentence_type;

typedef cicada::Grammar         grammar_type;
typedef cicada::Transducer      transducer_type;

typedef transducer_type::feature_set_type   feature_set_type;
typedef transducer_type::attribute_set_type attribute_set_type;

typedef transducer_type::rule_type           rule_type;
typedef transducer_type::rule_ptr_type       rule_ptr_type;
typedef transducer_type::rule_pair_type      rule_pair_type;

struct rule_pair_string_type
{
  std::string lhs;
  std::string source;
  std::string target;
  
  feature_set_type   features;
  attribute_set_type attributes;
  
  rule_pair_string_type() {}

  void clear()
  {
    lhs.clear();
    source.clear();
    target.clear();
    features.clear();
    attributes.clear();
  }
};

struct less_rule_pair
{
  bool operator()(const rule_pair_string_type& x, const rule_pair_string_type& y) const
  {
    return (x.lhs < y.lhs
	    || (!(y.lhs < x.lhs)
		&& (x.source < y.source
		    || (!(y.source < x.source)
			&& (x.target < y.target
			    || (!(y.target < x.target)
				&& (x.features < y.features
				    || (!(y.features < x.features)
					&& x.attributes < y.attributes))))))));
  }
};

typedef std::set<rule_pair_string_type, less_rule_pair, std::allocator<rule_pair_string_type> > rule_pair_unique_type;
typedef std::vector<rule_pair_type, std::allocator<rule_pair_type> >                            rule_pair_set_type;

template <typename Iterator>
struct features_generator : boost::spirit::karma::grammar<Iterator, feature_set_type()>
{
  features_generator() : features_generator::base_type(features)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    features %= " ||| " << ((standard::string << '=' << double20) % ' ');
  }
  
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return std::numeric_limits<double>::digits10 + 1;
    }
  };
  
  boost::spirit::karma::real_generator<double, real_precision> double20;
  boost::spirit::karma::rule<Iterator, feature_set_type()>     features;
};

template <typename Iterator>
struct attributes_generator : boost::spirit::karma::grammar<Iterator, attribute_set_type()>
{
  attributes_generator() : attributes_generator::base_type(attributes)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    data %= int64_ | double10 | string;
    attribute %= standard::string << '=' << data;
    attributes %= " ||| " << (attribute % ' ');
  }
    
  struct real_precision : boost::spirit::karma::real_policies<double>
  {
    static unsigned int precision(double) 
    { 
      return std::numeric_limits<double>::digits10 + 1;
    }
  };
    
  boost::spirit::karma::real_generator<double, real_precision> double10;
  boost::spirit::karma::int_generator<attribute_set_type::int_type, 10, false> int64_;
  utils::json_string_generator<Iterator, true> string;
  
  boost::spirit::karma::rule<Iterator, attribute_set_type::data_type()> data;
  boost::spirit::karma::rule<Iterator, attribute_set_type::value_type()> attribute;
  boost::spirit::karma::rule<Iterator, attribute_set_type()> attributes;
};


path_type input_file = "-";
path_type output_file = "-";

grammar_file_set_type grammar_files;
bool grammar_list = false;

bool input_sentence_mode = false;
bool input_lattice_mode = false;

bool treebank_mode = false;
bool pos_mode = false;

int threads = 1;

int debug = 0;


struct Task
{
  typedef utils::lockfree_list_queue<std::string, std::allocator<std::string> > queue_type;
  
  Task(queue_type& __queue,
       const grammar_type& __grammar) : queue(__queue), grammar(__grammar) {}
  
  void operator()()
  {
    grammar_type      grammar_local(grammar.clone());

    cicada::QueryCKY query(grammar_local, treebank_mode, pos_mode);
    
    rule_pair_set_type    rules;
    rule_pair_string_type rule_string;
    
    boost::iostreams::filtering_ostream os_source;
    boost::iostreams::filtering_ostream os_target;

    lattice_type lattice;
    sentence_type sentence;
    
    std::string line;
    for (;;) {
      queue.pop_swap(line);
      if (line.empty()) break;
      
      if (input_lattice_mode)
	lattice.assign(line);
      else {
	sentence.assign(line);
	lattice = lattice_type(sentence);
      }
      
      rules.clear();
      
      query(lattice, std::back_inserter(rules));

      if (debug)
	std::cerr << "# of rules: " << rules.size() << std::endl;
      
      rule_pair_set_type::iterator riter_end = rules.end();
      for (rule_pair_set_type::iterator riter = rules.begin(); riter != riter_end; ++ riter)
	if (riter->source || riter->target) {
	  rule_string.clear();
	
	  rule_string.lhs = (riter->source ? riter->source->lhs : riter->target->lhs);
	
	  os_source.push(boost::iostreams::back_inserter(rule_string.source));
	  os_target.push(boost::iostreams::back_inserter(rule_string.target));
	
	  if (riter->source)
	    os_source << riter->source->rhs;
	  if (riter->target)
	    os_target << riter->target->rhs;
	  
	  os_source.pop();
	  os_target.pop();
	  
	  rule_string.features.swap(riter->features);
	  rule_string.attributes.swap(riter->attributes);
	
	  rules_unique.insert(rule_string);
	}
    }
    
  }
  
  queue_type&              queue;
  const grammar_type&      grammar;
  
  rule_pair_unique_type      rules_unique;
};

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (grammar_list) {
      std::cout << grammar_type::lists();
      return 0;
    }
    
    if (int(input_lattice_mode) + input_sentence_mode == 0)
      input_sentence_mode = true;
    if (int(input_lattice_mode) + input_sentence_mode > 1)
      throw std::runtime_error("either lattice or sentence input");
    
    threads = utils::bithack::max(1, threads);

    // read grammars...
    grammar_type grammar(grammar_files.begin(), grammar_files.end());
    if (debug)
      std::cerr << "grammar: " << grammar.size() << std::endl;

    typedef Task task_type;
    typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
    
    task_type::queue_type queue(threads);
    task_set_type tasks(threads, task_type(queue, grammar));
    
    boost::thread_group workers;
    for (int i = 0; i != threads; ++ i)
      workers.add_thread(new boost::thread(boost::ref(tasks[i])));
    
    utils::compress_istream is(input_file, 1024 * 1024);
    
    std::string line;
    while (std::getline(is, line))
      if (! line.empty())
	queue.push_swap(line);
    
    for (int i = 0; i != threads; ++ i)
      queue.push(std::string());
    
    workers.join_all();

    rule_pair_unique_type rules_unique;
    
    for (int i = 0; i != threads; ++ i) {
      if (rules_unique.empty())
	rules_unique.swap(tasks[i].rules_unique);
      else
	rules_unique.insert(tasks[i].rules_unique.begin(), tasks[i].rules_unique.end());
      
      tasks[i].rules_unique.clear();
    }
    tasks.clear();
    
    typedef std::ostream_iterator<char> oiter_type;

    features_generator<oiter_type>   generate_features;
    attributes_generator<oiter_type> generate_attributes;
    
    if (! output_file.empty()) {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      utils::compress_ostream os(output_file, 1024 * 1024);
      
      rule_pair_unique_type::const_iterator iter_end = rules_unique.end();
      for (rule_pair_unique_type::const_iterator iter = rules_unique.begin(); iter != iter_end; ++ iter) {
	karma::generate(oiter_type(os),
			standard::string << " ||| " << standard::string << " ||| " << standard::string,
			iter->lhs, iter->source, iter->target);
	
	if (! iter->features.empty())
	  karma::generate(oiter_type(os), generate_features, iter->features);
	
	if (! iter->attributes.empty())
	  karma::generate(oiter_type(os), generate_attributes, iter->attributes);
	os << '\n';
      }
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_config("configuration options");
  opts_config.add_options()
    ("input",  po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("input-sentence",   po::bool_switch(&input_sentence_mode),   "sentence input")
    ("input-lattice",    po::bool_switch(&input_lattice_mode),    "lattice input")
    
    ("treebank", po::bool_switch(&treebank_mode), "assume treebank style grammar")
    ("pos",      po::bool_switch(&pos_mode),      "POS annotated input")
    
    // grammar
    ("grammar",           po::value<grammar_file_set_type >(&grammar_files)->composing(),      "grammar specification(s)")
    ("grammar-list",      po::bool_switch(&grammar_list),                                      "list of available grammar specifications")
    ;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("config",  po::value<path_type>(),                    "configuration file")
    ("threads", po::value<int>(&threads),                  "# of threads (highly experimental)")
    ("debug",   po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_config;
  po::options_description desc_command;
  po::options_description desc_visible;
  
  desc_config.add(opts_config);
  desc_command.add(opts_config).add(opts_command);
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
