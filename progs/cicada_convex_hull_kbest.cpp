//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
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
#include "utils/unordered_map.hpp"

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

int threads = 4;

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
  try {
    options(argc, argv);

    cicada::optimize::LineSearch::value_min = value_lower;
    cicada::optimize::LineSearch::value_max = value_upper;

    if (scorer_list) {
      std::cout << cicada::eval::Scorer::lists();
      return 0;
    }    
    
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
    
    line_search_type line_search(debug);
    
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    line_search(segments, value_lower, value_upper, scorers.error_metric(), OutputIterator(os, weights[direction_name]));
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

void compute_envelope(const scorer_document_type& scorers,
		      const hypothesis_map_type&  kbests,
		      const weight_set_type& origin,
		      const weight_set_type& direction,
		      segment_document_type& segments)
{
  typedef EnvelopeTask task_type;
  typedef task_type::queue_type queue_type;
  
  queue_type queue;
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, segments, origin, direction, scorers, kbests)));
  
  for (size_t seg = 0; seg != kbests.size(); ++ seg)
    if (! kbests[seg].empty())
      queue.push(seg);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
}

struct TaskRead
{
  typedef utils::lockfree_list_queue<path_type, std::allocator<path_type> > queue_type;
  
  TaskRead(queue_type& __queue)
    : queue(__queue) {}
  
  void operator()()
  {
    typedef boost::spirit::istream_iterator iter_type;
    typedef kbest_feature_parser<iter_type> parser_type;
    
    parser_type parser;
    kbest_feature_type kbest;

    for (;;) {
      path_type path;
      queue.pop(path);
      
      if (path.empty()) break;
      
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
	
	if (id >= kbests.size())
	  kbests.resize(id + 1);
	
	kbests[id].push_back(hypothesis_type(kbest));
      }
    }
  }
  
  queue_type& queue;
  hypothesis_map_type kbests;
};

void read_tstset(const path_set_type& files,
		 hypothesis_map_type& hypotheses)
{
  typedef TaskRead task_type;
  typedef task_type::queue_type queue_type;
  
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  queue_type queue(threads);
  task_set_type tasks(threads, task_type(queue));
    
  boost::thread_group workers;
  for (int i = 0; i != threads; ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (path_set_type::const_iterator fiter = files.begin(); fiter != files.end(); ++ fiter) {
    if (! boost::filesystem::exists(*fiter) && *fiter != "-")
      throw std::runtime_error("no file: " + fiter->string());
    
    if (boost::filesystem::is_directory(*fiter)) {
      for (size_t i = 0; /**/; ++ i) {
	const path_type path = (*fiter) / (utils::lexical_cast<std::string>(i) + ".gz");

	if (! boost::filesystem::exists(path)) break;
	
	queue.push(path);
      }
    } else
      queue.push(*fiter);
  }
  
  for (int i = 0; i != threads; ++ i)
    queue.push(path_type());
  
  workers.join_all();
  
  for (int i = 0; i != threads; ++ i) {
    if (tasks[i].kbests.size() > hypotheses.size())
      throw std::runtime_error("invalid kbests");
    
    for (size_t id = 0; id != tasks[i].kbests.size(); ++ id)
      hypotheses[id].insert(hypotheses[id].end(), tasks[i].kbests[id].begin(), tasks[i].kbests[id].end());
    
    tasks[i].kbests.clear();
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
  
  typedef utils::unordered_set<hypothesis_type, boost::hash<hypothesis_type>, std::equal_to<hypothesis_type>,
			       std::allocator<hypothesis_type> >::type hypothesis_unique_type;

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
  
  boost::thread_group workers;
  for (int i = 0; i < threads; ++ i)
    workers.add_thread(new boost::thread(task_type(queue, hypotheses, scorers)));
  
  for (size_t id = 0; id != hypotheses.size(); ++ id)
    queue.push(id);
  
  for (int i = 0; i < threads; ++ i)
    queue.push(-1);
  
  workers.join_all();
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
