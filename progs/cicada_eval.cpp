// evaluation tool....
// this will be nothing related to hypergraph!
// (or, do we compute oracle bleu after bleu-composition?)
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include "cicada/sentence.hpp"
#include "cicada/sentence_vector.hpp"
#include "cicada/eval.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Sentence       sentence_type;
typedef cicada::SentenceVector sentence_set_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;

path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";
std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool bootstrap = false;
int  samples = 1000;

int debug = 0;

void read_refset(const path_set_type& files, scorer_document_type& scorers);
void read_tstset(const path_set_type& files, sentence_set_type& sentences);

void options(int argc, char** argv);


int main(int argc, char** argv)
{
  try {

    options(argc, argv);
  
    if (scorer_list) {
      std::cout << cicada::eval::Scorer::lists();
      return 0;
    }
  
    // read reference set
    scorer_document_type scorers(scorer_name);
  
    read_refset(refset_files, scorers);

    if (tstset_files.empty())
      tstset_files.push_back("-");

    sentence_set_type sentences(scorers.size());
    
    read_tstset(tstset_files, sentences);
    
    if (bootstrap) {
      // first, collect BLEU stats...
      
      score_ptr_set_type scores;
      
      for (int seg = 0; seg != scorers.size(); ++ seg) 
	if (scorers[seg]) {
	  if (sentences[seg].empty()) {
	    std::cerr << "WARNING: no translation at: " << seg << std::endl;
	    continue;
	  }
	  
	  scores.push_back(scorers[seg]->score(sentences[seg]));
	}

      if (scores.empty())
	throw std::runtime_error("no error counts?");
      
      boost::mt19937 gen;
      gen.seed(time(0) * getpid());
      boost::random_number_generator<boost::mt19937> generator(gen);
      
      std::vector<double, std::allocator<double> > sampled;
      for (int iter = 0; iter < samples; ++ iter) {
	
	score_ptr_type score(scores.front()->zero());
	for (int i = 0; i != scores.size(); ++ i)
	  *score += *scores[generator(scores.size())];
	
	sampled.push_back(score->score().first);
      }

      std::sort(sampled.begin(), sampled.end());
      
      const int clip_size = sampled.size()* 0.25;
      
      utils::compress_ostream os(output_file);
      os << "mean: " << sampled[sampled.size() / 2]
	 << "95%-interval: " << sampled[clip_size] << " " << sampled[sampled.size() - clip_size - 1]
	 << '\n';
      
    } else {
      score_ptr_type score;
      
      for (int seg = 0; seg != scorers.size(); ++ seg) 
	if (scorers[seg]) {
	  
	  if (sentences[seg].empty()) {
	    std::cerr << "WARNING: no translation at: " << seg << std::endl;
	    continue;
	  }
	  
	  score_ptr_type score_segment = scorers[seg]->score(sentences[seg]);
	  
	  if (debug) {
	    const std::pair<double, double> value = score_segment->score();
	    std::cerr << "segment: " << seg << " score: " << value.first << " penalty: " << value.second << std::endl;
	  }
	  
	  if (! score)
	    score = score_segment;
	  else
	    *score += *score_segment;
	}
      
      if (! score)
	throw std::runtime_error("no statistics to compute error score");
      
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

typedef std::pair<int, sentence_type> id_sentence_type;

template <typename Iterator>
struct sentence_parser : boost::spirit::qi::grammar<Iterator, id_sentence_type(), boost::spirit::standard::space_type>
{
    
  sentence_parser() : sentence_parser::base_type(id_sentence)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    using qi::lexeme;
    using qi::int_;
    
    using standard::char_;
    using standard::space;
    
    word        %= lexeme[+(char_ - space) - "|||"];
    sentence    %= *word;
    id_sentence %= int_ >> "|||" >> sentence;
  }
  
  boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type>      word;
  boost::spirit::qi::rule<Iterator, sentence_type(), boost::spirit::standard::space_type>    sentence;
  boost::spirit::qi::rule<Iterator, id_sentence_type(), boost::spirit::standard::space_type> id_sentence;
};

void read_tstset(const path_set_type& files, sentence_set_type& sentences)
{
  sentence_parser<std::string::const_iterator> parser;
  id_sentence_type id_sentence;
  std::string line;

  std::vector<bool, std::allocator<bool> > finished(sentences.size());
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no test file: " + fiter->file_string());
    
    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      id_sentence.second.clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, end, parser, boost::spirit::standard::space, id_sentence))
	continue;
      
      const int& id = id_sentence.first;
      const sentence_type& sentence = id_sentence.second;
      
      if (id >= sentences.size())
	throw std::runtime_error("id exceeds the reference data");
      
      if (finished[id]) continue;
      
      sentences[id] = sentence;
      finished[id] = true;
    }
  }
}

void read_refset(const path_set_type& files, scorer_document_type& scorers)
{
  if (files.empty())
    throw std::runtime_error("no reference files?");
    
  scorers.clear();

  sentence_parser<std::string::const_iterator> parser;
  id_sentence_type id_sentence;
  std::string line;

  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->file_string());

    utils::compress_istream is(*fiter, 1024 * 1024);
    
    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      id_sentence.second.clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, end, parser, boost::spirit::standard::space, id_sentence))
	continue;

      const int& id = id_sentence.first;
      const sentence_type& sentence = id_sentence.second;
      
      if (id >= scorers.size())
	scorers.resize(id + 1);
      
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      scorers[id]->insert(sentence);
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

    ("bootstrap", po::bool_switch(&bootstrap), "bootstrap resampling")
    ("samples",   po::value<int>(&samples),    "# of samples")
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
