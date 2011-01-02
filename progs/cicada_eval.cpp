//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// evaluation tool....
// this will be nothing related to hypergraph!
// (or, do we compute oracle bleu after bleu-composition?)
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>

#include "cicada/sentence.hpp"
#include "cicada/sentence_vector.hpp"
#include "cicada/eval.hpp"

#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/random.hpp>

#include "cicada_text_impl.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Sentence       sentence_type;
typedef cicada::SentenceVector sentence_set_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef scorer_type::score_ptr_type score_ptr_type;
typedef std::vector<score_ptr_type, std::allocator<score_ptr_type> > score_ptr_set_type;

path_set_type tstset_files;
path_set_type tstset2_files;
path_set_type refset_files;
path_type     output_file = "-";
std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool signtest = false;
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

    if (bootstrap && samples <= 0)
      throw std::runtime_error("invalid sample size for bootstrapping");
    if (bootstrap && signtest)
      throw std::runtime_error("either --bootstrap/--signtest");

    if (! tstset2_files.empty() && ! signtest && ! bootstrap)
      throw std::runtime_error("the second system is used, but what test?");
    
    if (signtest && tstset2_files.empty())
      throw std::runtime_error("signtest without the second system?");
  
    // read reference set
    scorer_document_type scorers(scorer_name);
  
    read_refset(refset_files, scorers);

    if (tstset_files.empty())
      tstset_files.push_back("-");

    sentence_set_type hyps(scorers.size());
    sentence_set_type hyps2;
    
    read_tstset(tstset_files, hyps);

    if (! tstset2_files.empty()) {
      hyps2.reserve(scorers.size());
      hyps2.resize(scorers.size());
      read_tstset(tstset2_files, hyps2);
    }

    if (! hyps2.empty()) {
      const bool error_metric = scorers.error_metric();
      
      const sentence_set_type& hyps1 = hyps;
      score_ptr_set_type scores1;
      score_ptr_set_type scores2;
      
      for (size_t seg = 0; seg != scorers.size(); ++ seg) 
	if (scorers[seg]) {
	  if (hyps1[seg].empty() || hyps2[seg].empty()) {
	    if (hyps1[seg].empty())
	      std::cerr << "WARNING: no translation for system1 at: " << seg << std::endl;
	    if (hyps2[seg].empty())
	      std::cerr << "WARNING: no translation for system2 at: " << seg << std::endl;
	    
	    continue;
	  }
	  
	  scores1.push_back(scorers[seg]->score(hyps1[seg]));
	  scores2.push_back(scorers[seg]->score(hyps2[seg]));
	}

      if (bootstrap) {
      
	boost::mt19937 gen;
	gen.seed(time(0) * getpid());
	boost::random_number_generator<boost::mt19937> generator(gen);
      
	int better1 = 0;
	int better2 = 0;
	for (int iter = 0; iter < samples; ++ iter) {
	  score_ptr_type score1(scores1.front()->zero());
	  score_ptr_type score2(scores2.front()->zero());
	  for (size_t i = 0; i != scores1.size(); ++ i) {
	    const int seg = generator(scores1.size());
	    
	    *score1 += *scores1[seg];
	    *score2 += *scores2[seg];
	  }
	
	  const double eval1 = score1->score();
	  const double eval2 = score2->score();
	
	  if (error_metric) {
	    if (eval1 < eval2)
	      ++ better1;
	    if (eval2 < eval1)
	      ++ better2;
	  } else {
	    if (eval1 > eval2)
	      ++ better1;
	    if (eval2 > eval1)
	      ++ better2;
	  }
	}
      
	utils::compress_ostream os(output_file);
	os << "system1: " << better1 << " system2: " << better2 << std::endl;
      } else {
	int better = 0;
	int worse  = 0;
	
	score_ptr_type score1(scores1.front()->zero());
	score_ptr_type score2(scores2.front()->zero());
	for (size_t i = 0; i != scores1.size(); ++ i) {
	  *score1 += *scores1[i];
	  *score2 += *scores2[i];
	}
	
	const double eval1 = score1->score();
	const double eval2 = score2->score();
	
	score_ptr_type score = score1;
	for (size_t i = 0; i != scores2.size(); ++ i) {
	  *score -= *scores1[i];
	  *score += *scores2[i];
	  
	  const double eval2 = score->score();
	  
	  if (debug >= 2)
	    std::cerr << "system2: " << (*score) << std::endl;
	  
	  if (error_metric) {
	    if (eval2 < eval1)
	      ++ better;
	    else if (eval2 > eval1)
	      ++ worse;
	  } else {
	    if (eval2 > eval1)
	      ++ better;
	    else if (eval2 < eval1)
	      ++ worse;
	  }
	  
	  *score -= *scores2[i];
	  *score += *scores1[i];
	}
	
	const double n = better + worse;
	const double mean = double(better) / n;
	const double se = std::sqrt(mean * (1.0 - mean) / n);
	
	utils::compress_ostream os(output_file);
	os << "system1: " << eval1 << " system2: " << eval2 << '\n'
	   << "system2 better: " << better << " worse: " << worse << '\n'
	   << "Pr(better|different): " << mean << '\n'
	   << "95%-confidence: " << (mean-1.96*se) << ' ' << (mean+1.96*se) << '\n'
	   << "99%-confidence: " << (mean-2.58*se) << ' ' << (mean+2.58*se) << '\n';
	
	if (mean - 2.58*se > 0.5)
	  os << "system2 is significantly better (p < 0.01)";
	else if (mean + 2.58*se < 0.5)
	  os << "system2 is significantly worse (p < 0.01)";
	else if (mean - 1.96*se > 0.5)
	  os << "system2 is significantly better (p < 0.05)";
	else if (mean + 1.96*se < 0.5)
	  os << "system2 is significantly worse (p < 0.05)";
	else
	  os << "no significant difference";
	os << '\n';
      }
    } else if (bootstrap) {
      // first, collect BLEU stats...
      
      score_ptr_set_type scores;
      
      for (size_t seg = 0; seg != scorers.size(); ++ seg) 
	if (scorers[seg]) {
	  if (hyps[seg].empty()) {
	    std::cerr << "WARNING: no translation at: " << seg << std::endl;
	    continue;
	  }
	  
	  scores.push_back(scorers[seg]->score(hyps[seg]));
	}

      if (scores.empty())
	throw std::runtime_error("no error counts?");
      
      boost::mt19937 gen;
      gen.seed(time(0) * getpid());
      boost::random_number_generator<boost::mt19937> generator(gen);
      
      std::vector<double, std::allocator<double> > sampled;
      for (int iter = 0; iter < samples; ++ iter) {
	
	score_ptr_type score(scores.front()->zero());
	for (size_t i = 0; i != scores.size(); ++ i)
	  *score += *scores[generator(scores.size())];
	
	sampled.push_back(score->score());
      }

      std::sort(sampled.begin(), sampled.end());
      
      const int clip_size = sampled.size()* 0.25;
      
      utils::compress_ostream os(output_file);
      os << "mean: " << sampled[sampled.size() / 2]
	 << " 95%-interval: " << sampled[clip_size] << " " << sampled[sampled.size() - clip_size - 1]
	 << '\n';
      
    } else {
      score_ptr_type score;
      
      for (size_t seg = 0; seg != scorers.size(); ++ seg) 
	if (scorers[seg]) {
	  
	  if (hyps[seg].empty()) {
	    std::cerr << "WARNING: no translation at: " << seg << std::endl;
	    continue;
	  }
	  
	  score_ptr_type score_segment = scorers[seg]->score(hyps[seg]);
	  
	  if (debug)
	    std::cerr << "segment: " << seg << ' ' << (*score_segment) << std::endl;
	  
	  if (! score)
	    score = score_segment;
	  else
	    *score += *score_segment;
	}
      
      if (! score)
	throw std::runtime_error("no statistics to compute error score");

      if (debug >= 2)
	std::cerr << score->encode() << std::endl;
      
      utils::compress_ostream os(output_file);
      os << *score << '\n';
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

void read_tstset(const path_set_type& files, sentence_set_type& sentences)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  parser_type parser;
  id_sentence_type id_sentence;

  std::vector<bool, std::allocator<bool> > finished(sentences.size());
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no test file: " + fiter->file_string());
    
    utils::compress_istream is(*fiter, 1024 * 1024);
    is.unsetf(std::ios::skipws);

    iter_type iter(is);
    iter_type iter_end;
    
    while (iter != iter_end) {
      id_sentence.second.clear();
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
	if (iter != iter_end)
	  throw std::runtime_error("tstset parsing failed");
      
      const int& id = id_sentence.first;
      
      if (id >= static_cast<int>(sentences.size()))
	throw std::runtime_error("id exceeds the reference data");
      
      if (finished[id]) continue;
      
      sentences[id] = id_sentence.second;
      finished[id] = true;
    }
  }
  
  for (size_t id = 0; id != finished.size(); ++ id)
    if (! finished[id])
      std::cerr << "WARNING: no tst data for segment:" << id << std::endl;
}

void read_refset(const path_set_type& files, scorer_document_type& scorers)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  parser_type parser;
  id_sentence_type id_sentence;
  
  if (files.empty())
    throw std::runtime_error("no reference files?");
  
  scorers.clear();
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no reference file: " + fiter->file_string());

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
    ("tstset",   po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s)")
    ("tstset2",  po::value<path_set_type>(&tstset2_files)->multitoken(), "test set file(s)")
    ("refset",   po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")

    ("signtest",  po::bool_switch(&signtest),  "sign test")
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
