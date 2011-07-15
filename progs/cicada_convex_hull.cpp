//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// compute convex hull given a direction
//
// output is a set of points with evaluation score (BLEU)
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

#include "cicada/eval.hpp"

#include "cicada/optimize/line_search.hpp"

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

#include "cicada_text_impl.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

typedef cicada::Symbol   symbol_type;
typedef cicada::Vocab    vocab_type;
typedef cicada::Sentence sentence_type;

typedef cicada::HyperGraph hypergraph_type;
typedef cicada::Rule       rule_type;

typedef hypergraph_type::feature_set_type    feature_set_type;
typedef cicada::WeightVector<double>   weight_set_type;
typedef feature_set_type::feature_type feature_type;

typedef std::vector<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef cicada::optimize::LineSearch line_search_type;

typedef line_search_type::segment_type          segment_type;
typedef line_search_type::segment_set_type      segment_set_type;
typedef line_search_type::segment_document_type segment_document_type;


void read_tstset(const path_set_type& files, hypergraph_set_type& graphs);
void read_refset(const path_set_type& file, scorer_document_type& scorers);
void compute_envelope(const scorer_document_type& scorers,
		      const hypergraph_set_type&  graphs,
		      const weight_set_type& origin,
		      const weight_set_type& direction,
		      segment_document_type& segments);

void options(int argc, char** argv);

path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";

double value_lower = -100;
double value_upper =  100;

path_type weights_file;
path_type direction_file;

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool yield_sentence = false;
bool yield_alignment = false;
bool yield_span = false;

int threads = 4;

int debug = 0;

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);

    cicada::optimize::LineSearch::value_min = value_lower;
    cicada::optimize::LineSearch::value_max = value_upper;
    
    if (scorer_list) {
      std::cout << cicada::eval::Scorer::lists();
      return 0;
    }
    
    if (int(yield_sentence) + yield_alignment + yield_span > 1)
      throw std::runtime_error("specify either sentence|alignment|span yield");
    if (int(yield_sentence) + yield_alignment + yield_span == 0)
      yield_sentence = true;

    if (weights_file.empty() || ! boost::filesystem::exists(weights_file))
      throw std::runtime_error("no weight file? " + weights_file.string());
    if (direction_file.empty() || ! boost::filesystem::exists(direction_file))
      throw std::runtime_error("no weight file? " + direction_file.string());
    
    threads = utils::bithack::max(threads, 1);
    
    // read reference set
    scorer_document_type scorers(scorer_name);
    
    read_refset(refset_files, scorers);
    
    if (debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;

    // read test set
    if (tstset_files.empty())
      tstset_files.push_back("-");

    if (debug)
      std::cerr << "reading hypergraphs" << std::endl;

    hypergraph_set_type graphs(scorers.size());
    
    read_tstset(tstset_files, graphs);
    
    weight_set_type weights;
    weight_set_type direction;

    {
      utils::compress_istream is(weights_file, 1024 * 1024);
      is >> weights;
    }
    
    {
      utils::compress_istream is(direction_file, 1024 * 1024);
      is >> direction;
    }
    
    segment_document_type segments(graphs.size());

    compute_envelope(scorers, graphs, weights, direction, segments);
    
    
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct EnvelopeTask
{
  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;

  typedef cicada::semiring::Envelope envelope_type;
  typedef std::vector<envelope_type, std::allocator<envelope_type> >  envelope_set_type;

  typedef utils::lockfree_list_queue<int, std::allocator<int> >  queue_type;

    
  EnvelopeTask(queue_type& __queue,
	       segment_document_type&      __segments,
	       const weight_set_type&      __origin,
	       const weight_set_type&      __direction,
	       const scorer_document_type& __scorers,
	       const hypergraph_set_type&  __graphs)
    : queue(__queue),
      segments(__segments),
      origin(__origin),
      direction(__direction),
      scorers(__scorers),
      graphs(__graphs) {}

  void operator()()
  {
    envelope_set_type envelopes;

    int seg;
    
    while (1) {
      queue.pop(seg);
      if (seg < 0) break;
      
      envelopes.clear();
      envelopes.resize(graphs[seg].nodes.size());

      cicada::inside(graphs[seg], envelopes, cicada::semiring::EnvelopeFunction<weight_set_type>(origin, direction));
      
      const envelope_type& envelope = envelopes[graphs[seg].goal];
      const_cast<envelope_type&>(envelope).sort();
      
      envelope_type::const_iterator eiter_end = envelope.end();
      for (envelope_type::const_iterator eiter = envelope.begin(); eiter != eiter_end; ++ eiter) {
	const envelope_type::line_ptr_type& line = *eiter;
	
	const sentence_type yield = line->yield(cicada::operation::sentence_traversal());
	
	scorer_type::score_ptr_type score = scorers[seg]->score(yield);
	
	if (debug >= 4)
	  std::cerr << "segment: " << seg << " x: " << line->x << std::endl;
	
	segments[seg].push_back(std::make_pair(line->x, score));
      }
    }
  }
  
  queue_type& queue;

  segment_document_type& segments;

  const weight_set_type& origin;
  const weight_set_type& direction;
  
  const scorer_document_type& scorers;
  const hypergraph_set_type&  graphs;
};

void compute_envelope(const scorer_document_type& scorers,
		      const hypergraph_set_type&  graphs,
		      const weight_set_type& origin,
		      const weight_set_type& direction,
		      segment_document_type& segments)
{
  typedef EnvelopeTask task_type;
  typedef task_type::queue_type queue_type;
  
  queue_type queue;
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, segments, origin, direction, scorers, graphs)));
  
  for (size_t seg = 0; seg != graphs.size(); ++ seg)
    if (graphs[seg].goal != hypergraph_type::invalid)
      queue.push(seg);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
}


void read_tstset(const path_set_type& files, hypergraph_set_type& graphs)
{
  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (debug)
      std::cerr << "file: " << *titer << std::endl;
      
    if (boost::filesystem::is_directory(*titer)) {

      for (int i = 0; /**/; ++ i) {
	const path_type path = (*titer) / (utils::lexical_cast<std::string>(i) + ".gz");

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
	
	  if (id >= static_cast<int>(graphs.size()))
	    throw std::runtime_error("tstset size exceeds refset size?" + utils::lexical_cast<std::string>(id));
	
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
	
	if (id >= static_cast<int>(graphs.size()))
	  throw std::runtime_error("tstset size exceeds refset size?" + utils::lexical_cast<std::string>(id));
	
	graphs[id].unite(hypergraph);
      }
    }
  }
  
  for (size_t id = 0; id != graphs.size(); ++ id)
    if (graphs[id].goal == hypergraph_type::invalid)
      std::cerr << "invalid graph at: " << id << std::endl;
}

void read_refset(const path_set_type& files, scorer_document_type& scorers)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;

  if (files.empty())
    throw std::runtime_error("no reference files?");
    
  scorers.clear();

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
    ("tstset",  po::value<path_set_type>(&tstset_files)->multitoken(), "test set file(s) (in hypergraph format)")
    ("refset",  po::value<path_set_type>(&refset_files)->multitoken(), "reference set file(s)")
    
    ("output", po::value<path_type>(&output_file)->default_value(output_file), "output file")
    
    ("value-lower", po::value<double>(&value_lower)->default_value(value_lower), "default lower bounds")
    ("value-upper", po::value<double>(&value_upper)->default_value(value_upper), "default upper_bounds")
    
    ("weights",   po::value<path_type>(&weights_file),   "weights")
    ("direction", po::value<path_type>(&direction_file), "direction")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")

    ("yield-sentence",  po::bool_switch(&yield_sentence),  "optimize wrt sentence yield")
    ("yield-alignment", po::bool_switch(&yield_alignment), "optimize wrt alignment yield")
    ("yield-span",      po::bool_switch(&yield_span),      "optimize wrt span yield")

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
