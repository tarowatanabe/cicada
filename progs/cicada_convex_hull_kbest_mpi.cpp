//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// compute convex hull given a direction
//
// output is a set of points with evaluation score (BLEU)
//

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
#include "cicada_kbest_impl.hpp"
#include "cicada_mert_kbest_impl.hpp"

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

typedef std::vector<weight_set_type, std::allocator<weight_set_type> > weight_set_collection_type;

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef cicada::optimize::LineSearch line_search_type;

typedef line_search_type::segment_type          segment_type;
typedef line_search_type::segment_set_type      segment_set_type;
typedef line_search_type::segment_document_type segment_document_type;

path_set_type tstset_files;
path_set_type refset_files;
path_type     output_file = "-";

double value_lower = -100;
double value_upper =  100;

path_type weights_file;
std::string direction_name;

std::string scorer_name = "bleu:order=4";
bool scorer_list = false;

int debug = 0;

void read_tstset(const path_set_type& files, hypothesis_map_type& kbests);
void read_refset(const path_set_type& file, scorer_document_type& scorers);

void initialize_score(hypothesis_map_type& hypotheses,
		      const scorer_document_type& scorers);

void compute_envelope(const scorer_document_type& scorers,
		      const hypothesis_map_type&  kbests,
		      const weight_set_type& origin,
		      const weight_set_type& direction,
		      segment_document_type& segments);


void options(int argc, char** argv);

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
    
    // read reference set
    scorer_document_type scorers(scorer_name);
    read_refset(refset_files, scorers);

    if (mpi_rank == 0 && debug)
      std::cerr << "# of references: " << scorers.size() << std::endl;

    // read test set
    if (mpi_rank == 0 && debug)
      std::cerr << "reading kbests" << std::endl;

    hypothesis_map_type kbests(scorers.size());
    
    read_tstset(tstset_files, kbests);
    
    initialize_score(kbests, scorers);
    
    weight_set_type weights;
    {
      utils::compress_istream is(weights_file, 1024 * 1024);
      is >> weights;
    }
    
    weight_set_type direction;
    direction[direction_name] = 1.0;
    
    segment_document_type segments(scorers.size());
    
    compute_envelope(scorers, kbests, weights, direction, segments);
    
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
	       const hypothesis_map_type&  __kbests)
    : queue(__queue),
      segments(__segments),
      origin(__origin),
      direction(__direction),
      scorers(__scorers),
      kbests(__kbests) {}
  
  void operator()()
  {
    EnvelopeKBest::line_set_type lines;
    int seg;

    EnvelopeKBest envelopes(origin, direction);
    
    while (1) {
      queue.pop(seg);
      if (seg < 0) break;
      
      envelopes(kbests[seg], lines);
      
      EnvelopeKBest::line_set_type::const_iterator liter_end = lines.end();
      for (EnvelopeKBest::line_set_type::const_iterator liter = lines.begin(); liter != liter_end; ++ liter) {
	const EnvelopeKBest::line_type& line = *liter;
	
	if (debug >= 4)
	  std::cerr << "segment: " << seg << " x: " << line.x << std::endl;
	
	segments[seg].push_back(std::make_pair(line.x, line.hypothesis->score));
      }
    }
  }
  
  queue_type& queue;

  segment_document_type& segments;

  const weight_set_type& origin;
  const weight_set_type& direction;
  
  const scorer_document_type& scorers;
  const hypothesis_map_type&  kbests;
};

enum {
  envelope_tag = 1000,
};

void compute_envelope(const scorer_document_type& scorers,
		      const hypothesis_map_type&  kbests,
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
  
  boost::thread worker(task_type(queue, segments, origin, direction, scorers, kbests));
  
  for (size_t seg = 0; seg != kbests.size(); ++ seg)
    if (! kbests[seg].empty())
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
	  utils::encode_base64(siter->x, std::ostream_iterator<char>(os));
	  os << " ||| " << siter->score->encode()
	     << '\n';
	}
      }
  }
}


void read_tstset(const path_set_type& files,
		 hypothesis_map_type& hypotheses)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef kbest_feature_parser<iter_type> parser_type;

  const int mpi_rank = MPI::COMM_WORLD.Get_rank();
  const int mpi_size = MPI::COMM_WORLD.Get_size();
  
  parser_type parser;
  kbest_feature_type kbest;
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
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
	  
	  if (id != i)
	    throw std::runtime_error("id mismatch?");
	  if (i >= hypotheses.size())
	    throw std::runtime_error("invalid id?");
	  
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
	
	if (id % mpi_size != mpi_rank) continue;
	
	if (id >= hypotheses.size())
	  throw std::runtime_error("invalid id?");
	
	hypotheses[id].push_back(hypothesis_type(kbest));
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

struct TaskInit
{
  typedef utils::lockfree_list_queue<int, std::allocator<int> > queue_type;
  
#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
				  std::allocator<hypothesis_type> > hypothesis_unique_type;
#else
  typedef sgi::hash_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
			std::allocator<hypothesis_type> > hypothesis_unique_type;
#endif

  TaskInit(queue_type&                 __queue,
	   hypothesis_map_type&        __hypotheses,
	   const scorer_document_type& __scorers)
    : queue(__queue), hypotheses(__hypotheses), scorers(__scorers) {}

  void operator()()
  {
    hypothesis_unique_type kbests;

    for (;;) {
      int id = 0;
      queue.pop(id);
      if (id < 0) break;
      
      kbests.clear();
      kbests.insert(hypotheses[id].begin(), hypotheses[id].end());
      
      hypotheses[id].clear();
      hypothesis_set_type(hypotheses[id]).swap(hypotheses[id]);
      
      hypotheses[id].reserve(kbests.size());
      hypotheses[id].insert(hypotheses[id].end(), kbests.begin(), kbests.end());
      
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
  
  boost::thread worker(task_type(queue, hypotheses, scorers));
  
  for (size_t id = 0; id != hypotheses.size(); ++ id)
    if (! hypotheses[id].empty())
      queue.push(id);
  queue.push(-1);
  
  worker.join();
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
