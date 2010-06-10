// evaluation tool....
// this will be nothing related to hypergraph!
// (or, do we compute oracle bleu after bleu-composition?)
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "cicada/sentence.hpp"
#include "cicada/eval.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Sentence sentence_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type score_ptr_type;

path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";
std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

int debug = 0;

void read_refset(const path_set_type& files, scorer_document_type& scorers);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  options(argc, argv);

  try {
  
    if (scorer_list) {
      std::cout << cicada::eval::Scorer::lists();
      return 0;
    }
  
    // read reference set
    scorer_document_type scorers(scorer_name);
  
    read_refset(refset_files, scorers);
  
    std::vector<bool, std::allocator<bool> > finished(scorers.size(), false);
    score_ptr_type score;

    int id;
    std::string sep;
    sentence_type sentence;

    if (tstset_files.empty())
      tstset_files.push_back("-");
  
    for (path_set_type::const_iterator fiter = tstset_files.begin(); fiter != tstset_files.end(); ++ fiter) {
      
      if (! boost::filesystem::exists(*fiter) && *fiter != "-")
	throw std::runtime_error("no test set file: " + fiter->file_string());
      
      utils::compress_istream is(*fiter, 1024 * 1024);
      
      while (is >> id >> sep >> sentence) {
	if (sep != "|||") throw std::runtime_error("invalid tstset format...");
	
	if (id >= finished.size())
	  throw std::runtime_error("id exceeds our reference set...");
      
	if (finished[id]) continue;
      
	if (! score)
	  score = scorers[id]->score(sentence);
	else
	  *score += *scorers[id]->score(sentence);
	
	finished[id] = true;
      }

      const std::pair<double, double> value = score->score();
    
      utils::compress_ostream os(output_file);
      os << "score: " << value.first << " penalty: " << value.second << '\n';
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void read_refset(const path_set_type& files, scorer_document_type& scorers)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;

  if (files.empty())
    throw std::runtime_error("no reference files?");
    
  scorers.clear();

  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->file_string());

    utils::compress_istream is(*fiter, 1024 * 1024);
    
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
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;

  po::options_description opts_config("configuration options");
  
  opts_config.add_options()
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s) (in hypergraph format)")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")
    
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
