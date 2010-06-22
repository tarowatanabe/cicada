
//
// refset format:
// 0 |||  reference translatin for source sentence 0
// 0 |||  another reference
// 1 |||  reference translation for source sentence 1
//

#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <stdexcept>
#include <numeric>

#include "cicada/sentence.hpp"
#include "cicada/lattice.hpp"
#include "cicada/hypergraph.hpp"
#include "cicada/inside_outside.hpp"

#include "cicada/feature_function.hpp"
#include "cicada/weight_vector.hpp"
#include "cicada/semiring.hpp"
#include "cicada/viterbi.hpp"

#include "cicada/feature/bleu.hpp"

#include "cicada/eval.hpp"

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

#include "lbfgs.h"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Symbol   symbol_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef rule_type::feature_set_type    feature_set_type;
typedef cicada::WeightVector<double>   weight_set_type;
typedef feature_set_type::feature_type feature_type;

typedef std::vector<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef cicada::feature::Bleu bleu_feature_type;
typedef boost::shared_ptr<bleu_feature_type> bleu_feature_ptr_type;
typedef std::vector<bleu_feature_ptr_type, std::allocator<bleu_feature_ptr_type> > bleu_feature_ptr_set_type;

typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
typedef std::vector<sentence_set_type, std::allocator<sentence_set_type> > sentence_document_type;

struct source_length_function
{
  typedef cicada::semiring::Tropical<int> value_type;
  
  template <typename Edge>
  value_type operator()(const Edge& edge) const
  {
    int length = 0;
    rule_type::symbol_set_type::const_iterator siter_end = edge.rule->source.end();
    for (rule_type::symbol_set_type::const_iterator siter = edge.rule->source.begin(); siter != siter_end; ++ siter)
      length += (*siter != vocab_type::EPSILON && siter->is_terminal());
    
    // since we will "max" at operator+, we will collect negative length
    return cicada::semiring::traits<value_type>::log(- length);
  }
};

path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

int iteration = 100;
bool regularize_l1 = false;
bool regularize_l2 = false;
double C = 1.0;

int threads = 4;

int debug = 0;

void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 bleu_feature_ptr_set_type& features);
void read_refset(const path_set_type& file,
		 scorer_document_type& scorers,
		 sentence_document_type& sentences);
double optimize(const scorer_document_type& scorers,
		const hypergraph_set_type& graphs,
		const bleu_feature_ptr_set_type& features,
		weight_set_type& weights);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);

    if (scorer_list) {
      std::cout << cicada::eval::Scorer::lists();
      return 0;
    }
    
    if (regularize_l1 && regularize_l2)
      throw std::runtime_error("you cannot use both of L1 and L2...");
    
    threads = utils::bithack::max(threads, 1);
    
    // read reference set
    scorer_document_type      scorers(scorer_name);
    sentence_document_type    sentences;
    
    
    read_refset(refset_files, scorers, sentences);
    
    if (debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;

    // read test set
    if (tstset_files.empty())
      tstset_files.push_back("-");

    if (debug)
      std::cerr << "reading hypergraphs" << std::endl;

    hypergraph_set_type       graphs(scorers.size());
    bleu_feature_ptr_set_type features(scorers.size());
    
    read_tstset(tstset_files, graphs, sentences, features);

    if (debug)
      std::cerr << "# of features: " << feature_type::allocated() << std::endl;
    
    weight_set_type weights;
    
    const double objective = optimize(scorers, graphs, features, weights);
    
    if (debug)
      std::cerr << "objective: " << objective << std::endl;

    
    
    utils::compress_ostream os(output_file);
    os.precision(10);
    os << weights;
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


double optimize(const scorer_document_type& scorers,
		const hypergraph_set_type& graphs,
		const bleu_feature_ptr_set_type& features,
		weight_set_type& weights)
{
  
}

void read_tstset(const path_set_type& files,
		 hypergraph_set_type& graphs,
		 const sentence_document_type& sentences,
		 bleu_feature_ptr_set_type& features)
{
  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (debug)
      std::cerr << "file: " << *titer << std::endl;
      
    if (boost::filesystem::is_directory(*titer)) {

      for (int i = 0; /**/; ++ i) {
	const path_type path = (*titer) / (boost::lexical_cast<std::string>(i) + ".gz");

	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	
	int id;
	std::string sep;
	hypergraph_type hypergraph;
      
	weight_set_type origin;
	weight_set_type direction;
      
	while (is >> id >> sep >> hypergraph) {
	
	  if (sep != "|||")
	    throw std::runtime_error("format error?");
	
	  if (id >= graphs.size())
	    throw std::runtime_error("tstset size exceeds refset size?" + boost::lexical_cast<std::string>(id));
	
	  graphs[id].unite(hypergraph);
	}
      }
    } else {
      
      utils::compress_istream is(*titer, 1024 * 1024);
      
      int id;
      std::string sep;
      hypergraph_type hypergraph;
      
      weight_set_type origin;
      weight_set_type direction;
      
      while (is >> id >> sep >> hypergraph) {
	
	if (sep != "|||")
	  throw std::runtime_error("format error?");
	
	if (id >= graphs.size())
	  throw std::runtime_error("tstset size exceeds refset size?" + boost::lexical_cast<std::string>(id));
	
	graphs[id].unite(hypergraph);
      }
    }
  }
  
  for (int id = 0; id < graphs.size(); ++ id) {
    if (graphs[id].goal == hypergraph_type::invalid)
      std::cerr << "invalid graph at: " << id << std::endl;
    else {
      // we will enumerate forest structure... and collect min-size...
      std::vector<source_length_function::value_type, std::allocator<source_length_function::value_type> > lengths(graphs[id].nodes.size());
      
      cicada::inside(graphs[id], lengths, source_length_function());

      cicada::FeatureFunction::feature_function_ptr_type __feature = cicada::FeatureFunction::create(scorer_name + ",exact=true");
      
      cicada::feature::Bleu* __bleu = dynamic_cast<cicada::feature::Bleu*>(__feature.get());
      if (! __bleu)
	throw std::runtime_error("invalid bleu feature... with invalid scorer");
      
      features[id].reset(__bleu);
      
      const int source_length = - log(lengths.back());
      
      sentence_set_type::const_iterator siter_end = sentences[id].end();
      for (sentence_set_type::const_iterator siter = sentences[id].begin(); siter != siter_end; ++ siter)
	features[id]->insert(source_length, *siter);
    }
  }

  
}

void read_refset(const path_set_type& files,
		 scorer_document_type& scorers,
		 sentence_document_type& sentences)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;

  if (files.empty())
    throw std::runtime_error("no reference files?");
  
  scorers.clear();
  sentences.clear();

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
      
      if (id >= sentences.size())
	sentences.resize(id + 1);
      
      if (! scorers[id])
	scorers[id] = scorers.create();
      
      sentences[id].push_back(sentence_type(iter, tokenizer.end()));
      scorers[id]->insert(sentences[id].back());
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
    
    ("iteration",          po::value<int>(&iteration),          "# of mert iteration")
    
    ("regularize-l1", po::bool_switch(&regularize_l1), "regularization via L1")
    ("regularize-l2", po::bool_switch(&regularize_l2), "regularization via L2")
    ("C"            , po::value<double>(&C),           "regularization constant")
    
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
