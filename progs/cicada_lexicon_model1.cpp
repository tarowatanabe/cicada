//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

#include "cicada_lexicon_impl.hpp"

path_type source_file = "-";
path_type target_file = "-";
path_type span_source_file;
path_type span_target_file;
path_type lexicon_source_target_file;
path_type lexicon_target_source_file;
path_type output_source_target_file;
path_type output_target_source_file;
path_type viterbi_source_target_file;
path_type viterbi_target_source_file;

int iteration = 5;

bool symmetric_mode = false;
bool posterior_mode = false;
bool variational_bayes_mode = false;

bool moses_mode = false;
bool itg_mode = false;
bool max_match_mode = false;

// parameter...
double p0    = 1e-4;
double prior = 0.1;
double smooth = 1e-7;

double threshold = 0.0;

int threads = 2;

int debug = 0;

#include "cicada_lexicon_maximize_impl.hpp"
#include "cicada_lexicon_model1_impl.hpp"

template <typename Learner, typename Maximizer>
void learn(ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source,
	   aligned_type& aligned_source_target,
	   aligned_type& aligned_target_source);

template <typename Aligner>
void viterbi(const ttable_type& ttable_source_target,
	     const ttable_type& ttable_target_source);

void read(const path_type& path, ttable_type& lexicon);
void dump(const path_type& path, const ttable_type& lexicon, const aligned_type& aligned);

void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (itg_mode && max_match_mode)
      throw std::runtime_error("you cannot specify both of ITG and max-match for Viterbi alignment");
    
    threads = utils::bithack::max(threads, 1);
    
    ttable_type ttable_source_target(smooth);
    ttable_type ttable_target_source(smooth);
    
    aligned_type aligned_source_target;
    aligned_type aligned_target_source;
    
    if (! lexicon_source_target_file.empty())
      if (lexicon_source_target_file != "-" && ! boost::filesystem::exists(lexicon_source_target_file))
	throw std::runtime_error("no file: " + lexicon_source_target_file.file_string());

    if (! lexicon_target_source_file.empty())
      if (lexicon_target_source_file != "-" && ! boost::filesystem::exists(lexicon_target_source_file))
	throw std::runtime_error("no file: " + lexicon_target_source_file.file_string());

    boost::thread_group workers_read;
    
    if (! lexicon_source_target_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read, boost::cref(lexicon_source_target_file), boost::ref(ttable_source_target))));
    if (! lexicon_target_source_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read, boost::cref(lexicon_target_source_file), boost::ref(ttable_target_source))));
    
    workers_read.join_all();
    
    if (iteration > 0) {
      if (variational_bayes_mode) {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnModel1SymmetricPosterior, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnModel1Symmetric, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnModel1Posterior, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnModel1, MaximizeBayes>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	}
	
      } else {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnModel1SymmetricPosterior, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnModel1Symmetric, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnModel1Posterior, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	  else
	    learn<LearnModel1, Maximize>(ttable_source_target, ttable_target_source, aligned_source_target, aligned_target_source);
	}
      }
    }
    
    if (! viterbi_source_target_file.empty() || ! viterbi_target_source_file.empty()) {
      if (itg_mode)
	viterbi<ITGModel1>(ttable_source_target, ttable_target_source);
      else if (max_match_mode)
	viterbi<MaxMatchModel1>(ttable_source_target, ttable_target_source);
      else
	viterbi<ViterbiModel1>(ttable_source_target, ttable_target_source);
    }

      
    // final dumping...
    boost::thread_group workers_dump;

    if (! output_source_target_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_source_target_file),
							    boost::cref(ttable_source_target),
							    boost::cref(aligned_source_target))));
    
    if (! output_target_source_file.empty())
      workers_dump.add_thread(new boost::thread(boost::bind(dump,
							    boost::cref(output_target_source_file),
							    boost::cref(ttable_target_source),
							    boost::cref(aligned_target_source))));
    
    workers_dump.join_all();
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}

struct greater_second
{
  template <typename Tp>
  bool operator()(const Tp* x, const Tp* y) const
  {
    return x->second > y->second;
  }
};

void read(const path_type& path, ttable_type& lexicon)
{
  typedef boost::fusion::tuple<std::string, std::string, double > lexicon_parsed_type;
  typedef boost::spirit::istream_iterator iterator_type;

  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  
  qi::rule<iterator_type, std::string(), standard::blank_type>         word;
  qi::rule<iterator_type, lexicon_parsed_type(), standard::blank_type> parser; 
  
  word   %= qi::lexeme[+(standard::char_ - standard::space)];
  parser %= word >> word >> qi::double_ >> (qi::eol | qi::eoi);
  
  lexicon.clear();
  lexicon.smooth = boost::numeric::bounds<double>::highest();
  
  utils::compress_istream is(path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iterator_type iter(is);
  iterator_type iter_end;
  
  lexicon_parsed_type lexicon_parsed;
  
  while (iter != iter_end) {
    boost::fusion::get<0>(lexicon_parsed).clear();
    boost::fusion::get<1>(lexicon_parsed).clear();
    
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, lexicon_parsed))
      if (iter != iter_end)
	throw std::runtime_error("global lexicon parsing failed");
    
    const word_type target(boost::fusion::get<0>(lexicon_parsed));
    const word_type source(boost::fusion::get<1>(lexicon_parsed));
    const double&   prob(boost::fusion::get<2>(lexicon_parsed));
    
    lexicon[source][target] = prob;
    lexicon.smooth = std::min(lexicon.smooth, prob * 0.1);
  }
}

void dump(const path_type& path, const ttable_type& lexicon, const aligned_type& aligned)
{
  typedef ttable_type::count_map_type::value_type value_type;
  typedef std::vector<const value_type*, std::allocator<const value_type*> > sorted_type;

  utils::compress_ostream os(path, 1024 * 1024);
  os.precision(10);

  const aligned_type::aligned_map_type __empty;
  sorted_type sorted;
  
  ttable_type::count_dict_type::const_iterator siter_begin = lexicon.ttable.begin();
  ttable_type::count_dict_type::const_iterator siter_end   = lexicon.ttable.end();
  for (ttable_type::count_dict_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) 
    if (*siter) {
      const word_type source(word_type::id_type(siter - siter_begin));
      const ttable_type::count_map_type& dict = *(*siter);
      
      if (dict.empty()) continue;
      
      sorted.clear();
      sorted.reserve(dict.size());
      
      ttable_type::count_map_type::const_iterator titer_end = dict.end();
      for (ttable_type::count_map_type::const_iterator titer = dict.begin(); titer != titer_end; ++ titer)
	if (titer->second >= 0.0)
	  sorted.push_back(&(*titer));
      
      std::sort(sorted.begin(), sorted.end(), greater_second());
      
      if (threshold > 0.0) {
	const double prob_max       = sorted.front()->second;
	const double prob_threshold = prob_max * threshold;
	
	const aligned_type::aligned_map_type& viterbi = (aligned.exists(source) ? aligned[source] : __empty);
	
	// TODO: extra checking to keep Viterbi alignemnt in the final output!
	
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  if ((*iter)->second >= prob_threshold || viterbi.find((*iter)->first) != viterbi.end())
	    os << (*iter)->first << ' ' << source << ' '  << (*iter)->second << '\n';
      } else {
	sorted_type::const_iterator iter_end = sorted.end();
	for (sorted_type::const_iterator iter = sorted.begin(); iter != iter_end; ++ iter)
	  os << (*iter)->first << ' ' << source << ' '  << (*iter)->second << '\n';
      }
    }
}

template <typename LearnerSet, typename Maximizer>
struct TaskMaximize : public Maximizer
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  TaskMaximize(const LearnerSet& __learners,
	       const int __id,
	       ttable_type& __ttable_source_target,
	       ttable_type& __ttable_target_source,
	       aligned_type& __aligned_source_target,
	       aligned_type& __aligned_target_source)
    : learners(__learners),
      id(__id),
      ttable_source_target(__ttable_source_target),
      ttable_target_source(__ttable_target_source),
      aligned_source_target(__aligned_source_target),
      aligned_target_source(__aligned_target_source) {}
  
  void operator()()
  {
    for (word_type::id_type source_id = id; source_id < ttable_source_target.size(); source_id += learners.size()) {
      for (size_t i = 0; i != learners.size(); ++ i) {
	if (learners[i].counts_source_target.exists(source_id))
	  ttable_source_target[source_id] += learners[i].counts_source_target[source_id];
	
	if (learners[i].aligned_source_target.exists(source_id))
	  aligned_source_target[source_id] += learners[i].aligned_source_target[source_id];
      }
      
      if (ttable_source_target.exists(source_id))
	Maximizer::operator()(ttable_source_target[source_id]);
    }
    
    for (word_type::id_type target_id = id; target_id < ttable_target_source.size(); target_id += learners.size()) {
      for (size_t i = 0; i != learners.size(); ++ i) {
	if (learners[i].counts_target_source.exists(target_id))
	  ttable_target_source[target_id] += learners[i].counts_target_source[target_id];
	
	if (learners[i].aligned_target_source.exists(target_id))
	  aligned_target_source[target_id] += learners[i].aligned_target_source[target_id];
      }
      
      if (ttable_target_source.exists(target_id))
	Maximizer::operator()(ttable_target_source[target_id]);
    }
  }
  
  const LearnerSet& learners;
  const int id;

  ttable_type& ttable_source_target;
  ttable_type& ttable_target_source;
  
  aligned_type& aligned_source_target;
  aligned_type& aligned_target_source;
};

template <typename Learner>
struct TaskLearn : public Learner
{
  struct bitext_type
  {
    bitext_type() : source(), target() {}
    bitext_type(const sentence_type& __source, const sentence_type& __target)
      : source(__source), target(__target) {}
    
    sentence_type source;
    sentence_type target;
  };
  
  typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
  typedef utils::lockfree_list_queue<bitext_set_type, std::allocator<bitext_set_type> > queue_type;
  
  TaskLearn(queue_type& __queue,
	    const ttable_type& ttable_source_target,
	    const ttable_type& ttable_target_source)
    : Learner(ttable_source_target, ttable_target_source), queue(__queue) {}
  
  void operator()()
  {
    Learner::initialize();

    bitext_set_type bitexts;
    
    for (;;) {
      bitexts.clear();
      queue.pop_swap(bitexts);
      if (bitexts.empty()) break;
      
      typename bitext_set_type::const_iterator biter_end = bitexts.end();
      for (typename bitext_set_type::const_iterator biter = bitexts.begin(); biter != biter_end; ++ biter)
	Learner::operator()(biter->source, biter->target);
    }
  }
  
  queue_type& queue;
};

template <typename Learner, typename Maximizer>
void learn(ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source,
	   aligned_type& aligned_source_target,
	   aligned_type& aligned_target_source)
{
  typedef TaskLearn<Learner> learner_type;
  
  typedef typename learner_type::bitext_type     bitext_type;
  typedef typename learner_type::bitext_set_type bitext_set_type;
  typedef typename learner_type::queue_type      queue_type;
  
  typedef std::vector<learner_type, std::allocator<learner_type> > learner_set_type;
  
  typedef TaskMaximize<learner_set_type, Maximizer> maximizer_type;

  queue_type       queue(threads * 64);
  learner_set_type learners(threads, learner_type(queue, ttable_source_target, ttable_target_source));
  
  for (int iter = 0; iter < iteration; ++ iter) {
    if (debug)
      std::cerr << "iteration: " << iter << std::endl;

    utils::resource accumulate_start;
    
    boost::thread_group workers_learn;
    for (size_t i = 0; i != learners.size(); ++ i)
      workers_learn.add_thread(new boost::thread(boost::ref(learners[i])));
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    
    bitext_type     bitext;
    bitext_set_type bitexts;
    
    size_t num_bitext = 0;
    
    for (;;) {
      is_src >> bitext.source;
      is_trg >> bitext.target;
      
      if (! is_src || ! is_trg) break;      
      
      if (bitext.source.empty() || bitext.target.empty()) continue;
      
      bitexts.push_back(bitext);

      ++ num_bitext;
      if (debug) {
	if (num_bitext % 10000 == 0)
	  std::cerr << '.';
	if (num_bitext % 1000000 == 0)
	  std::cerr << '\n';
      }
      
      if (bitexts.size() == 64) {
	queue.push_swap(bitexts);
	bitexts.clear();
      }
    }

    if (! bitexts.empty())
      queue.push_swap(bitexts);
    
    if (debug && num_bitext >= 10000)
      std::cerr << std::endl;
    if (debug)
      std::cerr << "# of bitexts: " << num_bitext << std::endl;

    if (is_src || is_trg)
      throw std::runtime_error("# of samples do not match");
        
    for (size_t i = 0; i != learners.size(); ++ i) {
      bitexts.clear();
      queue.push_swap(bitexts);
    }
    
    workers_learn.join_all();
    
    // merge and normalize...! 
    ttable_source_target.resize(word_type::allocated());
    ttable_target_source.resize(word_type::allocated());
    aligned_source_target.resize(word_type::allocated());
    aligned_target_source.resize(word_type::allocated());
    
    ttable_source_target.initialize();
    ttable_target_source.initialize();
    aligned_source_target.initialize();
    aligned_target_source.initialize();

    double objective_source_target = 0;
    double objective_target_source = 0;
    
    boost::thread_group workers_maximize;
    for (size_t i = 0; i != learners.size(); ++ i) {
      objective_source_target += learners[i].objective_source_target;
      objective_target_source += learners[i].objective_target_source;
      
      workers_maximize.add_thread(new boost::thread(maximizer_type(learners, i,
								   ttable_source_target, ttable_target_source,
								   aligned_source_target, aligned_target_source)));
    }
    
    if (debug)
      std::cerr << "perplexity for P(target | source): " << objective_source_target << '\n'
		<< "perplexity for P(source | target): " << objective_target_source << '\n';
    
    workers_maximize.join_all();
    
    utils::resource accumulate_end;
    
    if (debug)
      std::cerr << "cpu time:  " << accumulate_end.cpu_time() - accumulate_start.cpu_time() << std::endl
		<< "user time: " << accumulate_end.user_time() - accumulate_start.user_time() << std::endl;
    
  }
}

struct ViterbiMapReduce
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  struct bitext_type
  {
    bitext_type() : id(size_type(-1)), source(), target(), alignment() {}
    bitext_type(const size_type& __id, const sentence_type& __source, const sentence_type& __target)
      : id(__id), source(__source), target(__target), alignment() {}
    bitext_type(const size_type& __id, const sentence_type& __source, const sentence_type& __target, const alignment_type& __alignment)
      : id(__id), source(__source), target(__target), alignment(__alignment) {}
    
    bitext_type(const size_type& __id,
		const sentence_type& __source, const sentence_type& __target,
		const span_set_type& __span_source, const span_set_type& __span_target,
		const alignment_type& __alignment)
      : id(__id),
	source(__source), target(__target),
	span_source(__span_source), span_target(__span_target),
	alignment(__alignment) {}
    
    size_type     id;
    sentence_type source;
    sentence_type target;
    span_set_type span_source;
    span_set_type span_target;
    alignment_type alignment;

    void clear()
    {
      id = size_type(-1);
      source.clear();
      target.clear();
      span_source.clear();
      span_target.clear();
      alignment.clear();
    }
    
    void swap(bitext_type& x)
    {
      std::swap(id, x.id);
      source.swap(x.source);
      target.swap(x.target);
      span_source.swap(x.span_source);
      span_target.swap(x.span_target);
      alignment.swap(x.alignment);
    }
  };
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> > queue_type;
};

namespace std
{
  inline
  void swap(ViterbiMapReduce::bitext_type& x, ViterbiMapReduce::bitext_type& y)
  {
    x.swap(y);
  }
};

template <typename Aligner>
struct ViterbiMapper : public ViterbiMapReduce, public Aligner
{
  queue_type& mapper;
  queue_type& reducer_source_target;
  queue_type& reducer_target_source;
  
  ViterbiMapper(const Aligner& __aligner, queue_type& __mapper, queue_type& __reducer_source_target, queue_type& __reducer_target_source)
    : Aligner(__aligner),
      mapper(__mapper),
      reducer_source_target(__reducer_source_target),
      reducer_target_source(__reducer_target_source)
  {}

  void operator()()
  {
    bitext_type bitext;
    alignment_type alignment_source_target;
    alignment_type alignment_target_source;
    
    for (;;) {
      mapper.pop_swap(bitext);
      if (bitext.id == size_type(-1)) break;
      
      alignment_source_target.clear();
      alignment_target_source.clear();
      
      if (! bitext.source.empty() && ! bitext.target.empty())
	Aligner::operator()(bitext.source, bitext.target, bitext.span_source, bitext.span_target, alignment_source_target, alignment_target_source);
      
      reducer_source_target.push(bitext_type(bitext.id, bitext.source, bitext.target, alignment_source_target));
      reducer_target_source.push(bitext_type(bitext.id, bitext.target, bitext.source, alignment_target_source));
    }
  }
};

struct ViterbiReducer : public ViterbiMapReduce
{
  struct less_bitext
  {
    bool operator()(const bitext_type& x, const bitext_type& y) const
    {
      return x.id < y.id;
    }
  };
  typedef std::set<bitext_type, less_bitext, std::allocator<bitext_type> > bitext_set_type;

  path_type   path;
  queue_type& queue;
  
  ViterbiReducer(const path_type& __path, queue_type& __queue) : path(__path), queue(__queue) {}

  typedef int index_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > align_set_type;
  typedef std::set<index_type, std::less<index_type>, std::allocator<index_type> > align_none_type;
  
  align_set_type  aligns;
  align_none_type aligns_none;
  
  void dump(std::ostream& os, const bitext_type& bitext)
  {
    if (moses_mode)
      os << bitext.alignment << '\n';
    else {
      os << "# Sentence pair (" << (bitext.id + 1) << ')'
	 << " source length " << bitext.source.size()
	 << " target length " << bitext.target.size()
	 << " alignment score : " << 0 << '\n';
      os << bitext.target << '\n';
    
      if (bitext.alignment.empty() || bitext.source.empty() || bitext.target.empty()) {
	os << "NULL ({ })";
	sentence_type::const_iterator siter_end = bitext.source.end();
	for (sentence_type::const_iterator siter = bitext.source.begin(); siter != siter_end; ++ siter)
	  os << ' ' << *siter << " ({ })";
	os << '\n';
      } else {
	aligns.clear();
	aligns.resize(bitext.source.size());
      
	aligns_none.clear();
	for (size_type trg = 0; trg != bitext.target.size(); ++ trg)
	  aligns_none.insert(trg + 1);
      
	alignment_type::const_iterator aiter_end = bitext.alignment.end();
	for (alignment_type::const_iterator aiter = bitext.alignment.begin(); aiter != aiter_end; ++ aiter) {
	  aligns[aiter->source].push_back(aiter->target + 1);
	  aligns_none.erase(aiter->target + 1);
	}
      
	os << "NULL";
	os << " ({ ";
	std::copy(aligns_none.begin(), aligns_none.end(), std::ostream_iterator<index_type>(os, " "));
	os << "})";
      
	for (size_type src = 0; src != bitext.source.size(); ++ src) {
	  os << ' ' << bitext.source[src];
	  os << " ({ ";
	  std::copy(aligns[src].begin(), aligns[src].end(), std::ostream_iterator<index_type>(os, " "));
	  os << "})";
	}
	os << '\n';
      }
    }
  }
  
  void operator()()
  {
    if (path.empty()) {
      bitext_type bitext;
      for (;;) {
	queue.pop_swap(bitext);
	if (bitext.id == size_type(-1)) break;
      }
    } else {
      bitext_set_type bitexts;

      const bool flush_output = (path == "-"
			       || (boost::filesystem::exists(path)
				   && ! boost::filesystem::is_regular_file(path)));
      
      utils::compress_ostream os(path, 1024 * 1024 * (! flush_output));
      
      size_type id = 0;
      bitext_type bitext;
      for (;;) {
	queue.pop_swap(bitext);
	if (bitext.id == size_type(-1)) break;
	
	if (bitext.id == id) {
	  dump(os, bitext);
	  ++ id;
	} else
	  bitexts.insert(bitext);
	
	while (! bitexts.empty() && bitexts.begin()->id == id) {
	  dump(os, *bitexts.begin());
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id == id) {
	dump(os, *bitexts.begin());
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while dumping viterbi output?");
    }
  }
};

template <typename Aligner>
void viterbi(const ttable_type& ttable_source_target,
	     const ttable_type& ttable_target_source)
{
  typedef ViterbiReducer         reducer_type;
  typedef ViterbiMapper<Aligner> mapper_type;
  
  typedef reducer_type::bitext_type bitext_type;
  typedef reducer_type::queue_type  queue_type;
  
  typedef std::vector<mapper_type, std::allocator<mapper_type> > mapper_set_type;

  if (debug)
    std::cerr << "Viterbi alignment" << std::endl;
  
  queue_type queue(threads * 4096);
  queue_type queue_source_target(threads * 4096);
  queue_type queue_target_source(threads * 4096);
  
  boost::thread_group mapper;
  for (int i = 0; i != threads; ++ i)
    mapper.add_thread(new boost::thread(mapper_type(Aligner(ttable_source_target, ttable_target_source),
						    queue,
						    queue_source_target,
						    queue_target_source)));

  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(viterbi_source_target_file, queue_source_target)));
  reducer.add_thread(new boost::thread(reducer_type(viterbi_target_source_file, queue_target_source)));
  
  bitext_type bitext;
  bitext.id = 0;
  
  utils::resource viterbi_start;
  
  utils::compress_istream is_src(source_file, 1024 * 1024);
  utils::compress_istream is_trg(target_file, 1024 * 1024);

  std::auto_ptr<std::istream> is_span_src(! span_source_file.empty()
					  ? new utils::compress_istream(span_source_file, 1024 * 1024) : 0);
  std::auto_ptr<std::istream> is_span_trg(! span_target_file.empty()
					  ? new utils::compress_istream(span_target_file, 1024 * 1024) : 0);
  
  for (;;) {
    is_src >> bitext.source;
    is_trg >> bitext.target;
    
    if (is_span_src.get())
      *is_span_src >> bitext.span_source;
    if (is_span_trg.get())
      *is_span_trg >> bitext.span_target;
    
    if (! is_src || ! is_trg || (is_span_src.get() && ! *is_span_src) || (is_span_trg.get() && ! *is_span_trg)) break;
    
    queue.push(bitext);
    
    ++ bitext.id;
    
    if (debug) {
      if (bitext.id % 10000 == 0)
	std::cerr << '.';
      if (bitext.id % 1000000 == 0)
	std::cerr << '\n';
    }
  }

  if (debug && bitext.id >= 10000)
    std::cerr << std::endl;
  if (debug)
    std::cerr << "# of bitexts: " << bitext.id << std::endl;
  
  if (is_src || is_trg || (is_span_src.get() && *is_span_src) || (is_span_trg.get() && *is_span_trg))
    throw std::runtime_error("# of samples do not match");
  
  for (int i = 0; i != threads; ++ i) {
    bitext.clear();
    queue.push_swap(bitext);
  }

  mapper.join_all();

  bitext.clear();
  queue_source_target.push_swap(bitext);
  bitext.clear();
  queue_target_source.push_swap(bitext);
  
  reducer.join_all();

  utils::resource viterbi_end;
  
  if (debug)
    std::cerr << "cpu time:  " << viterbi_end.cpu_time() - viterbi_start.cpu_time() << std::endl
	      << "user time: " << viterbi_end.user_time() - viterbi_start.user_time() << std::endl;
}



void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source", po::value<path_type>(&source_file), "source file")
    ("target", po::value<path_type>(&target_file), "target file")
    
    ("span-source", po::value<path_type>(&span_source_file), "source span file")
    ("span-target", po::value<path_type>(&span_target_file), "target span file")
    
    ("lexicon-source-target", po::value<path_type>(&lexicon_source_target_file), "lexicon model for P(target | source)")
    ("lexicon-target-source", po::value<path_type>(&lexicon_target_source_file), "lexicon model for P(source | target)")
    
    ("output-source-target", po::value<path_type>(&output_source_target_file), "lexicon model output for P(target | source)")
    ("output-target-source", po::value<path_type>(&output_target_source_file), "lexicon model output for P(source | target)")
    
    ("viterbi-source-target", po::value<path_type>(&viterbi_source_target_file), "viterbi for P(target | source)")
    ("viterbi-target-source", po::value<path_type>(&viterbi_target_source_file), "viterbi for P(source | target)")
    
    ("iteration", po::value<int>(&iteration)->default_value(iteration), "max iteration")
    
    ("symmetric",  po::bool_switch(&symmetric_mode),  "symmetric model1 training")
    ("posterior",  po::bool_switch(&posterior_mode),  "posterior constrained model1 training")
    ("variational-bayes", po::bool_switch(&variational_bayes_mode), "variational Bayes estimates")
    
    ("itg",       po::bool_switch(&itg_mode),       "ITG alignment")
    ("max-match", po::bool_switch(&max_match_mode), "maximum matching alignment")
    ("moses",     po::bool_switch(&moses_mode),     "Moses alignment foramt")
    
    ("p0",     po::value<double>(&p0)->default_value(p0),         "parameter for NULL alignment")
    ("prior",  po::value<double>(&prior)->default_value(prior),   "Dirichlet prior for variational Bayes")
    ("smooth", po::value<double>(&smooth)->default_value(smooth), "smoothing parameter for uniform distribution")

    ("threshold", po::value<double>(&threshold)->default_value(threshold), "dump with beam-threshold (<= 0.0 implies no beam)")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
