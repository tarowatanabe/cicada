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
#include "utils/space_separator.hpp"
#include "utils/base64.hpp"
#include "utils/mpi.hpp"
#include "utils/mpi_device.hpp"

#include <boost/tokenizer.hpp>
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
std::string direction_name;

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

bool yield_sentence = false;
bool yield_alignment = false;
bool yield_span = false;

int debug = 0;

struct OutputIterator
{
  OutputIterator(std::ostream& __os, const double& __weight)
    : os(__os), weight(__weight) {}
  
  OutputIterator& operator=(const line_search_type::value_type& value)
  {
    const double feature_lower = weight + value.lower;
    const double feature_upper = weight + value.upper;
    const double average = (value.lower + value.upper) * 0.5;
    
    const double point = (weight + average) * double(! (feature_lower * feature_upper < 0.0));
    
    os << point << ' '  << value.objective << '\n';
    
    return *this;
  }
  
  OutputIterator& operator*() { return *this; }
  OutputIterator& operator++() { return *this; }
  OutputIterator  operator++(int) { return *this; }
 
  
  std::ostream& os;
  double weight;
};

int main(int argc, char ** argv)
{
  utils::mpi_world mpi_world(argc, argv);
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
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
    if (direction_name.empty())
      throw std::runtime_error("no direction?");
    
    // read reference set
    scorer_document_type scorers(scorer_name);
    
    read_refset(refset_files, scorers);
    
    if (mpi_rank == 0 && debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;

    // read test set
    
    if (mpi_rank == 0 && debug)
      std::cerr << "reading hypergraphs" << std::endl;

    hypergraph_set_type graphs(scorers.size());
    
    read_tstset(tstset_files, graphs);
    
    weight_set_type weights;
    {
      utils::compress_istream is(weights_file, 1024 * 1024);
      is >> weights;
    }
    
    weight_set_type direction;
    direction[direction_name] = 1.0;
    
    segment_document_type segments(graphs.size());
    
    compute_envelope(scorers, graphs, weights, direction, segments);

    if (mpi_rank == 0) {
      line_search_type line_search(debug);
      
      utils::compress_ostream os(output_file, 1024 * 1024);
      
      line_search(segments, value_lower, value_upper, scorers.error_metric(), OutputIterator(os, weights[direction_name]));
    }
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

enum {
  envelope_tag = 1000,
};

void compute_envelope(const scorer_document_type& scorers,
		      const hypergraph_set_type&  graphs,
		      const weight_set_type& origin,
		      const weight_set_type& direction,
		      segment_document_type& segments)
{
  typedef utils::mpi_device_source idevice_type;
  typedef utils::mpi_device_sink   odevice_type;
  
  typedef boost::iostreams::filtering_ostream ostream_type;
  typedef boost::iostreams::filtering_istream istream_type;


  typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;

  typedef EnvelopeTask task_type;
  typedef task_type::queue_type queue_type;
  
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  queue_type queue;
  
  boost::thread worker(task_type(queue, segments, origin, direction, scorers, graphs));
  
  for (size_t seg = 0; seg != graphs.size(); ++ seg)
    if (graphs[seg].is_valid())
      queue.push(seg);
  queue.push(-1);
  
  worker.join();
  
  // merge segments into root
  if (mpi_rank == 0) {

    for (int rank = 1; rank != mpi_size; ++ rank) {
      istream_type is;
      is.push(boost::iostreams::zlib_decompressor());
      is.push(idevice_type(rank, envelope_tag, 4096));
      
      std::string line;
      while (std::getline(is, line)) {
	const utils::piece line_piece(line);
	tokenizer_type tokenizer(line_piece);
	
	tokenizer_type::iterator iter = tokenizer.begin();
	if (iter == tokenizer.end()) continue;
	const utils::piece id_str = *iter; 
	
	++ iter;
	if (iter == tokenizer.end() || *iter != "|||") continue;
	
	++ iter;
	if (iter == tokenizer.end()) continue;
	const utils::piece x_str = *iter;
	
	++ iter;
	if (iter == tokenizer.end() || *iter != "|||") continue;
	
	++ iter;
	if (iter == tokenizer.end()) continue;
	const utils::piece score_str = *iter;
	
	const int id = utils::lexical_cast<int>(id_str);
	
	if (id >= static_cast<int>(segments.size()))
	  throw std::runtime_error("invali id?");
	
	segments[id].push_back(std::make_pair(utils::decode_base64<double>(x_str),
					      scorer_type::score_type::decode(score_str)));
      }
    }
    
  } else {
    ostream_type os;
    os.push(boost::iostreams::zlib_compressor());
    os.push(odevice_type(0, envelope_tag, 4096));
    
    for (size_t seg = 0; seg != segments.size(); ++ seg)
      if (! segments[seg].empty()) {
	segment_set_type::const_iterator siter_end = segments[seg].end();
	for (segment_set_type::const_iterator siter = segments[seg].begin(); siter != siter_end; ++ siter) {
	  os << seg << " ||| ";
	  utils::encode_base64(siter->first, std::ostream_iterator<char>(os));
	  os << " ||| " << siter->second->encode()
	     << '\n';
	}
      }
  }
}


void read_tstset(const path_set_type& files, hypergraph_set_type& graphs)
{
  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();

  path_set_type::const_iterator titer_end = tstset_files.end();
  for (path_set_type::const_iterator titer = tstset_files.begin(); titer != titer_end; ++ titer) {
    
    if (debug)
      std::cerr << "file: " << *titer << std::endl;
      
    if (boost::filesystem::is_directory(*titer)) {

      for (int i = mpi_rank; /**/; i += mpi_size) {
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
	
	if (id % mpi_size != mpi_rank) continue;
	
	graphs[id].unite(hypergraph);
      }
    }
  }
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
    
    ("weights",   po::value<path_type>(&weights_file),     "weights")
    ("direction", po::value<std::string>(&direction_name), "direction (feature name)")

    ("scorer",      po::value<std::string>(&scorer_name)->default_value(scorer_name), "error metric")
    ("scorer-list", po::bool_switch(&scorer_list),                                    "list of error metric")

    ("yield-sentence",  po::bool_switch(&yield_sentence),  "optimize wrt sentence yield")
    ("yield-alignment", po::bool_switch(&yield_alignment), "optimize wrt alignment yield")
    ("yield-span",      po::bool_switch(&yield_span),      "optimize wrt span yield")
    
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
