//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <queue>

#include "cicada_alignment_impl.hpp"

#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/bithack.hpp"
#include "utils/mathop.hpp"
#include "utils/double_base64_parser.hpp"
#include "utils/double_base64_generator.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type source_file = "-";
path_type target_file = "-";
path_type alignment_file;
path_type dependency_source_file;
path_type dependency_target_file;
path_type span_source_file;
path_type span_target_file;
path_type classes_source_file;
path_type classes_target_file;
path_type lexicon_source_target_file;
path_type lexicon_target_source_file;
path_type alignment_source_target_file;
path_type alignment_target_source_file;
double length_source_target = 1.0;
double length_target_source = 1.0;
path_type output_lexicon_source_target_file;
path_type output_lexicon_target_source_file;
path_type output_alignment_source_target_file;
path_type output_alignment_target_source_file;
path_type viterbi_source_target_file;
path_type viterbi_target_source_file;
path_type projected_source_file;
path_type projected_target_file;
path_type posterior_source_target_file;
path_type posterior_target_source_file;

int iteration_model1 = 5;
int iteration_hmm = 5;

bool symmetric_mode = false;
bool posterior_mode = false;
bool variational_bayes_mode = false;
bool pgd_mode = false;

bool moses_mode = false;
bool itg_mode = false;
bool max_match_mode = false;

bool permutation_mode = false;
bool hybrid_mode = false;
bool degree2_mode = false;
bool mst_mode = false;
bool single_root_mode = false;

// parameter...
double p0    = 0.01;
double prior_lexicon = 0.01;
double smooth_lexicon = 1e-100;
double prior_alignment = 0.01;
double smooth_alignment = 1e-100;

double l0_alpha = 100;
double l0_beta = 0.01;

double threshold = 0.0;

int threads = 2;

int debug = 0;

#include "cicada_alignment_maximize_impl.hpp"
#include "cicada_alignment_model1_impl.hpp"
#include "cicada_alignment_hmm_impl.hpp"

template <typename Learner, typename Maximizer>
void learn(const Maximizer& maximier,
	   const int iteration,
	   ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source,
	   atable_type& atable_source_target,
	   atable_type& atalbe_target_source,
	   const classes_type& classes_source,
	   const classes_type& classes_target,
	   aligned_type& aligned_source_target,
	   aligned_type& aligned_target_source);

template <typename Aligner>
void viterbi(const ttable_type& ttable_source_target,
	     const ttable_type& ttable_target_source,
	     const atable_type& atable_source_target,
	     const atable_type& atalbe_target_source,
	     const classes_type& classes_source,
	     const classes_type& classes_target);

template <typename Analyzer>
void project_dependency(const ttable_type& ttable_source_target,
			const ttable_type& ttable_target_source,
			const atable_type& atable_source_target,
			const atable_type& atalbe_target_source,
			const classes_type& classes_source,
			const classes_type& classes_target);

template <typename Infer>
void posterior(const ttable_type& ttable_source_target,
	       const ttable_type& ttable_target_source,
	       const atable_type& atable_source_target,
	       const atable_type& atalbe_target_source,
	       const classes_type& classes_source,
	       const classes_type& classes_target);


void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    if (itg_mode && max_match_mode)
      throw std::runtime_error("you cannot specify both of ITG and max-match for Viterbi alignment");
    
    if (variational_bayes_mode && pgd_mode)
      throw std::runtime_error("either variational-bayes, pgd or none");

    if (! projected_target_file.empty())    
      if (dependency_source_file != "-" && ! boost::filesystem::exists(dependency_source_file))
	throw std::runtime_error("no source side dependency");
    
    if (! projected_source_file.empty())
      if (dependency_target_file != "-" && ! boost::filesystem::exists(dependency_target_file))
	throw std::runtime_error("no target side dependency");
    
    if (int(hybrid_mode) + degree2_mode + mst_mode + permutation_mode > 1)
      throw std::runtime_error("you cannot specify both of Hybrid, Degree2 and MST dependency, permutation parsing");
    
    if (int(hybrid_mode) + degree2_mode + mst_mode + permutation_mode == 0)
      hybrid_mode = true;
    
    threads = utils::bithack::max(threads, 1);
    
    ttable_type ttable_source_target(prior_lexicon, smooth_lexicon);
    ttable_type ttable_target_source(prior_lexicon, smooth_lexicon);
    
    atable_type atable_source_target(prior_alignment, smooth_alignment);
    atable_type atable_target_source(prior_alignment, smooth_alignment);
    
    classes_type classes_source;
    classes_type classes_target;
    
    aligned_type aligned_source_target;
    aligned_type aligned_target_source;
    
    if (! lexicon_source_target_file.empty())
      if (lexicon_source_target_file != "-" && ! boost::filesystem::exists(lexicon_source_target_file))
	throw std::runtime_error("no file: " + lexicon_source_target_file.string());

    if (! lexicon_target_source_file.empty())
      if (lexicon_target_source_file != "-" && ! boost::filesystem::exists(lexicon_target_source_file))
	throw std::runtime_error("no file: " + lexicon_target_source_file.string());

    if (! alignment_source_target_file.empty())
      if (alignment_source_target_file != "-" && ! boost::filesystem::exists(alignment_source_target_file))
	throw std::runtime_error("no file: " + alignment_source_target_file.string());

    if (! alignment_target_source_file.empty())
      if (alignment_target_source_file != "-" && ! boost::filesystem::exists(alignment_target_source_file))
	throw std::runtime_error("no file: " + alignment_target_source_file.string());

    if (! classes_source_file.empty())
      if (classes_source_file != "-" && ! boost::filesystem::exists(classes_source_file))
	throw std::runtime_error("no file: " + classes_source_file.string());

    if (! classes_target_file.empty())
      if (classes_target_file != "-" && ! boost::filesystem::exists(classes_target_file))
	throw std::runtime_error("no file: " + classes_target_file.string());
    
    boost::thread_group workers_read;
    
    // read lexicon
    if (! lexicon_source_target_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read_lexicon, boost::cref(lexicon_source_target_file), boost::ref(ttable_source_target))));
    if (! lexicon_target_source_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read_lexicon, boost::cref(lexicon_target_source_file), boost::ref(ttable_target_source))));
    
    // read alignment
    if (! alignment_source_target_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read_alignment, boost::cref(alignment_source_target_file), boost::ref(atable_source_target))));
    if (! alignment_target_source_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read_alignment, boost::cref(alignment_target_source_file), boost::ref(atable_target_source))));

    // read classes
    if (! classes_source_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read_classes, boost::cref(classes_source_file), boost::ref(classes_source))));
    if (! classes_target_file.empty())
      workers_read.add_thread(new boost::thread(boost::bind(read_classes, boost::cref(classes_target_file), boost::ref(classes_target))));
    
    workers_read.join_all();
    
    if (iteration_model1 > 0) {
      if (debug)
	std::cerr << "start Model1 training" << std::endl;
      
      if (variational_bayes_mode) {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnModel1SymmetricPosterior, MaximizeBayes>(MaximizeBayes(),
								iteration_model1,
								ttable_source_target,
								ttable_target_source,
								atable_source_target,
								atable_target_source,
								classes_source,
								classes_target,
								aligned_source_target,
								aligned_target_source);
	  else
	    learn<LearnModel1Symmetric, MaximizeBayes>(MaximizeBayes(),
						       iteration_model1,
						       ttable_source_target,
						       ttable_target_source,
						       atable_source_target,
						       atable_target_source,
						       classes_source,
						       classes_target,
						       aligned_source_target,
						       aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnModel1Posterior, MaximizeBayes>(MaximizeBayes(),
						       iteration_model1,
						       ttable_source_target,
						       ttable_target_source,
						       atable_source_target,
						       atable_target_source,
						       classes_source,
						       classes_target,
						       aligned_source_target,
						       aligned_target_source);
	  else
	    learn<LearnModel1, MaximizeBayes>(MaximizeBayes(),
					      iteration_model1,
					      ttable_source_target,
					      ttable_target_source,
					      atable_source_target,
					      atable_target_source,
					      classes_source,
					      classes_target,
					      aligned_source_target,
					      aligned_target_source);
	}

      } else if (pgd_mode) {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnModel1SymmetricPosterior, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
							     iteration_model1,
							     ttable_source_target,
							     ttable_target_source,
							     atable_source_target,
							     atable_target_source,
							     classes_source,
							     classes_target,
							     aligned_source_target,
							     aligned_target_source);
	  else
	    learn<LearnModel1Symmetric, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
						    iteration_model1,
						    ttable_source_target,
						    ttable_target_source,
						    atable_source_target,
						    atable_target_source,
						    classes_source,
						    classes_target,
						    aligned_source_target,
						    aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnModel1Posterior, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
						    iteration_model1,
						    ttable_source_target,
						    ttable_target_source,
						    atable_source_target,
						    atable_target_source,
						    classes_source,
						    classes_target,
						    aligned_source_target,
						    aligned_target_source);
	  else
	    learn<LearnModel1, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
					   iteration_model1,
					   ttable_source_target,
					   ttable_target_source,
					   atable_source_target,
					   atable_target_source,
					   classes_source,
					   classes_target,
					   aligned_source_target,
					   aligned_target_source);
	}
	
      } else {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnModel1SymmetricPosterior, Maximize>(Maximize(),
							   iteration_model1,
							   ttable_source_target,
							   ttable_target_source,
							   atable_source_target,
							   atable_target_source,
							   classes_source,
							   classes_target,
							   aligned_source_target,
							   aligned_target_source);
	  else
	    learn<LearnModel1Symmetric, Maximize>(Maximize(),
						  iteration_model1,
						  ttable_source_target,
						  ttable_target_source,
						  atable_source_target,
						  atable_target_source,
						  classes_source,
						  classes_target,
						  aligned_source_target,
						  aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnModel1Posterior, Maximize>(Maximize(),
						  iteration_model1,
						  ttable_source_target,
						  ttable_target_source,
						  atable_source_target,
						  atable_target_source,
						  classes_source,
						  classes_target,
						  aligned_source_target,
						  aligned_target_source);
	  else
	    learn<LearnModel1, Maximize>(Maximize(),
					 iteration_model1,
					 ttable_source_target,
					 ttable_target_source,
					 atable_source_target,
					 atable_target_source,
					 classes_source,
					 classes_target,
					 aligned_source_target,
					 aligned_target_source);
	}
      }
    }
    
    if (iteration_hmm > 0) {
      if (debug)
	std::cerr << "start HMM training" << std::endl;

      if (variational_bayes_mode) {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnHMMSymmetricPosterior, MaximizeBayes>(MaximizeBayes(),
							     iteration_hmm,
							     ttable_source_target,
							     ttable_target_source,
							     atable_source_target,
							     atable_target_source,
							     classes_source,
							     classes_target,
							     aligned_source_target,
							     aligned_target_source);
	  else
	    learn<LearnHMMSymmetric, MaximizeBayes>(MaximizeBayes(),
						    iteration_hmm,
						    ttable_source_target,
						    ttable_target_source,
						    atable_source_target,
						    atable_target_source,
						    classes_source,
						    classes_target,
						    aligned_source_target,
						    aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnHMMPosterior, MaximizeBayes>(MaximizeBayes(),
						    iteration_hmm,
						    ttable_source_target,
						    ttable_target_source,
						    atable_source_target,
						    atable_target_source,
						    classes_source,
						    classes_target,
						    aligned_source_target,
						    aligned_target_source);
	  else
	    learn<LearnHMM, MaximizeBayes>(MaximizeBayes(),
					   iteration_hmm,
					   ttable_source_target,
					   ttable_target_source,
					   atable_source_target,
					   atable_target_source,
					   classes_source,
					   classes_target,
					   aligned_source_target,
					   aligned_target_source);
	}

      } else if (pgd_mode) {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnHMMSymmetricPosterior, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
							  iteration_hmm,
							  ttable_source_target,
							  ttable_target_source,
							  atable_source_target,
							  atable_target_source,
							  classes_source,
							  classes_target,
							  aligned_source_target,
							  aligned_target_source);
	  else
	    learn<LearnHMMSymmetric, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
						 iteration_hmm,
						 ttable_source_target,
						 ttable_target_source,
						 atable_source_target,
						 atable_target_source,
						 classes_source,
						 classes_target,
						 aligned_source_target,
						 aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnHMMPosterior, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
						 iteration_hmm,
						 ttable_source_target,
						 ttable_target_source,
						 atable_source_target,
						 atable_target_source,
						 classes_source,
						 classes_target,
						 aligned_source_target,
						 aligned_target_source);
	  else
	    learn<LearnHMM, MaximizeL0>(MaximizeL0(l0_alpha, l0_beta),
					iteration_hmm,
					ttable_source_target,
					ttable_target_source,
					atable_source_target,
					atable_target_source,
					classes_source,
					classes_target,
					aligned_source_target,
					aligned_target_source);
	}	
	
      } else {
	if (symmetric_mode) {
	  if (posterior_mode)
	    learn<LearnHMMSymmetricPosterior, Maximize>(Maximize(),
							iteration_hmm,
							ttable_source_target,
							ttable_target_source,
							atable_source_target,
							atable_target_source,
							classes_source,
							classes_target,
							aligned_source_target,
							aligned_target_source);
	  else
	    learn<LearnHMMSymmetric, Maximize>(Maximize(),
					       iteration_hmm,
					       ttable_source_target,
					       ttable_target_source,
					       atable_source_target,
					       atable_target_source,
					       classes_source,
					       classes_target,
					       aligned_source_target,
					       aligned_target_source);
	} else {
	  if (posterior_mode)
	    learn<LearnHMMPosterior, Maximize>(Maximize(),
					       iteration_hmm,
					       ttable_source_target,
					       ttable_target_source,
					       atable_source_target,
					       atable_target_source,
					       classes_source,
					       classes_target,
					       aligned_source_target,
					       aligned_target_source);
	  else
	    learn<LearnHMM, Maximize>(Maximize(),
				      iteration_hmm,
				      ttable_source_target,
				      ttable_target_source,
				      atable_source_target,
				      atable_target_source,
				      classes_source,
				      classes_target,
				      aligned_source_target,
				      aligned_target_source);
	}
      }
    }
    
    if (! viterbi_source_target_file.empty() || ! viterbi_target_source_file.empty()) {
      if (itg_mode) {
	if (debug)
	  std::cerr << "ITG alignment" << std::endl;

	viterbi<ITGHMM>(ttable_source_target,
			ttable_target_source,
			atable_source_target,
			atable_target_source,
			classes_source,
			classes_target);
      } else if (max_match_mode) {
	if (debug)
	  std::cerr << "Max-Match alignment" << std::endl;
	
	viterbi<MaxMatchHMM>(ttable_source_target,
			     ttable_target_source,
			     atable_source_target,
			     atable_target_source,
			     classes_source,
			     classes_target);
      } else {
	if (debug)
	  std::cerr << "Viterbi alignment" << std::endl;
	
	viterbi<ViterbiHMM>(ttable_source_target,
			    ttable_target_source,
			    atable_source_target,
			    atable_target_source,
			    classes_source,
			    classes_target);
      }
    }

    // dependency parsing projection
    if (! projected_source_file.empty() || ! projected_target_file.empty()) {
      if (hybrid_mode) {
	if (debug)
	  std::cerr << "hybrid projective dependency" << std::endl;
	
	if (single_root_mode)
	  project_dependency<DependencyHybridSingleRootHMM>(ttable_source_target,
							    ttable_target_source,
							    atable_source_target,
							    atable_target_source,
							    classes_source,
							    classes_target);
	else
	  project_dependency<DependencyHybridHMM>(ttable_source_target,
						  ttable_target_source,
						  atable_source_target,
						  atable_target_source,
						  classes_source,
						  classes_target);

      } else if (degree2_mode) {
	if (debug)
	  std::cerr << "degree2 non-projective dependency" << std::endl;
	
	if (single_root_mode)
	  project_dependency<DependencyDegree2SingleRootHMM>(ttable_source_target,
							     ttable_target_source,
							     atable_source_target,
							     atable_target_source,
							     classes_source,
							     classes_target);
	else
	  project_dependency<DependencyDegree2HMM>(ttable_source_target,
						   ttable_target_source,
						   atable_source_target,
						   atable_target_source,
						   classes_source,
						   classes_target);
      } else if (mst_mode) {
	if (debug)
	  std::cerr << "MST non-projective dependency" << std::endl;
	
	if (single_root_mode)
	  project_dependency<DependencyMSTSingleRootHMM>(ttable_source_target,
							 ttable_target_source,
							 atable_source_target,
							 atable_target_source,
							 classes_source,
							 classes_target);
	else
	  project_dependency<DependencyMSTHMM>(ttable_source_target,
					       ttable_target_source,
					       atable_source_target,
					       atable_target_source,
					       classes_source,
					       classes_target);
      } else if (permutation_mode) {
	if (debug)
	  std::cerr << "permutation" << std::endl;
	
	project_dependency<PermutationHMM>(ttable_source_target,
					   ttable_target_source,
					   atable_source_target,
					   atable_target_source,
					   classes_source,
					   classes_target);
      } else
	throw std::runtime_error("no dependency algorithm?");
    }
    
    if (! posterior_source_target_file.empty() || ! posterior_target_source_file.empty()) {
      if (debug)
	std::cerr << "compute posterior" << std::endl;
      
      posterior<PosteriorHMM>(ttable_source_target,
			      ttable_target_source,
			      atable_source_target,
			      atable_target_source,
			      classes_source,
			      classes_target);
    }
    
    // final writing
    boost::thread_group workers_write;
    
    // write lexicon
    if (! output_lexicon_source_target_file.empty())
      workers_write.add_thread(new boost::thread(boost::bind(write_lexicon,
							     boost::cref(output_lexicon_source_target_file),
							     boost::cref(ttable_source_target),
							     boost::cref(aligned_source_target),
							     threshold)));
    
    if (! output_lexicon_target_source_file.empty())
      workers_write.add_thread(new boost::thread(boost::bind(write_lexicon,
							     boost::cref(output_lexicon_target_source_file),
							     boost::cref(ttable_target_source),
							     boost::cref(aligned_target_source),
							     threshold)));

    // write alignment
    if (! output_alignment_source_target_file.empty())
      workers_write.add_thread(new boost::thread(boost::bind(write_alignment,
							     boost::cref(output_alignment_source_target_file),
							     boost::cref(atable_source_target))));
    
    if (! output_alignment_target_source_file.empty())
      workers_write.add_thread(new boost::thread(boost::bind(write_alignment,
							     boost::cref(output_alignment_target_source_file),
							     boost::cref(atable_target_source))));
    
    
    workers_write.join_all();
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

struct LearnMapReduce
{
  struct bitext_type
  {
    bitext_type() : source(), target(), alignment() {}
    bitext_type(const sentence_type& __source, const sentence_type& __target)
      : source(__source), target(__target), alignment() {}
    bitext_type(const sentence_type& __source, const sentence_type& __target, const alignment_type& __alignment)
      : source(__source), target(__target), alignment(__alignment) {}
    
    sentence_type source;
    sentence_type target;
    alignment_type alignment;
  };
  
  typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
  
  struct ttable_counts_type
  {
    word_type                      word;
    ttable_type::count_map_type    counts;
    aligned_type::aligned_map_type aligned;
    
    ttable_counts_type() : word(), counts(), aligned() {}
    
    void swap(ttable_counts_type& x)
    {
      word.swap(x.word);
      counts.swap(x.counts);
      aligned.swap(x.aligned);
    }
  };

  
  typedef utils::lockfree_list_queue<bitext_set_type, std::allocator<bitext_set_type> >       queue_bitext_type;
  typedef utils::lockfree_list_queue<ttable_counts_type, std::allocator<ttable_counts_type> > queue_ttable_type;
  typedef std::vector<queue_ttable_type, std::allocator<queue_ttable_type> >                  queue_ttable_set_type;
};

namespace std
{
  inline
  void swap(LearnMapReduce::ttable_counts_type& x, LearnMapReduce::ttable_counts_type& y)
  {
    x.swap(y);
  }
};

template <typename Maximizer>
struct LearnReducer : public Maximizer
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef LearnMapReduce map_reduce_type;
  
  typedef map_reduce_type::ttable_counts_type ttable_counts_type;
  typedef map_reduce_type::queue_ttable_type  queue_ttable_type;
  
  LearnReducer(queue_ttable_type& __queue,
	       const ttable_type& __ttable,
	       const aligned_type& __aligned,
	       ttable_type& __ttable_new,
	       aligned_type& __aligned_new,
	       const Maximizer& __base)
    : Maximizer(__base),
      queue(__queue),
      ttable(__ttable),
      aligned(__aligned),
      ttable_new(__ttable_new),
      aligned_new(__aligned_new) {}
  
  void operator()()
  {
    ttable_type  ttable_reduced;
    aligned_type aligned_reduced;
    
    for (;;) {
      ttable_counts_type counts;
      
      queue.pop_swap(counts);
      
      if (counts.counts.empty() && counts.aligned.empty()) break;
      
      if (! counts.counts.empty())
	ttable_reduced[counts.word] += counts.counts;
      if (! counts.aligned.empty())
	aligned_reduced[counts.word] += counts.aligned;
    }
    
    for (word_type::id_type word_id = 0; word_id < aligned_reduced.size(); ++ word_id)
      if (aligned_reduced.exists(word_id)) {
	aligned_new[word_id].swap(aligned_reduced[word_id]);
	
	aligned_reduced.clear(word_id);
      }

    ttable_type::count_map_type ttable_empty;
    
    for (word_type::id_type word_id = 0; word_id < ttable_reduced.size(); ++ word_id)
      if (ttable_reduced.exists(word_id)) { 
	if (ttable.exists(word_id))
	  Maximizer::operator()(ttable_reduced[word_id], ttable[word_id], ttable_new[word_id], ttable.prior, ttable.smooth);
	else
	  Maximizer::operator()(ttable_reduced[word_id], ttable_empty, ttable_new[word_id], ttable.prior, ttable.smooth);
	
	ttable_reduced.clear(word_id);
      }
  }
  
  queue_ttable_type& queue;
  
  const ttable_type&  ttable;
  const aligned_type& aligned;

  ttable_type&  ttable_new;
  aligned_type& aligned_new;
};

template <typename Learner>
struct LearnMapper : public Learner
{
  typedef LearnMapReduce map_reduce_type;
  
  typedef map_reduce_type::bitext_set_type    bitext_set_type;
  typedef map_reduce_type::ttable_counts_type ttable_counts_type;
  
  typedef map_reduce_type::queue_bitext_type      queue_bitext_type;
  typedef map_reduce_type::queue_ttable_type      queue_ttable_type;
  typedef map_reduce_type::queue_ttable_set_type  queue_ttable_set_type;
  
  LearnMapper(queue_bitext_type& __queue_bitext,
	      queue_ttable_set_type& __queue_ttable_source_target,
	      queue_ttable_set_type& __queue_ttable_target_source,
	      const LearnBase& __base)
    : Learner(__base),
      queue_bitext(__queue_bitext),
      queue_ttable_source_target(__queue_ttable_source_target),
      queue_ttable_target_source(__queue_ttable_target_source) {}
  
  void operator()()
  {
    Learner::initialize();

    bitext_set_type    bitexts;
    
    const int iter_mask = (1 << 5) - 1;
    
    for (int iter = 0;; ++ iter) {
      bitexts.clear();
      queue_bitext.pop_swap(bitexts);
      if (bitexts.empty()) break;
      
      typename bitext_set_type::const_iterator biter_end = bitexts.end();
      for (typename bitext_set_type::const_iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {
	if (biter->alignment.empty())
	  Learner::operator()(biter->source, biter->target);
	else
	  Learner::operator()(biter->source, biter->target, biter->alignment);
      }

      if ((iter & iter_mask) == iter_mask) {
	dump();
	Learner::shrink();
      }
    }
    
    dump();
    Learner::shrink();
  }

  void dump()
  {
    const word_type::id_type source_max = utils::bithack::max(Learner::ttable_counts_source_target.size(),
							      Learner::aligned_source_target.size());
    const word_type::id_type target_max = utils::bithack::max(Learner::ttable_counts_target_source.size(),
							      Learner::aligned_target_source.size());
    
    for (word_type::id_type source_id = 0; source_id != source_max; ++ source_id) {
      ttable_counts_type counts;
	  
      if (Learner::ttable_counts_source_target.exists(source_id) && ! Learner::ttable_counts_source_target[source_id].empty())
	counts.counts.swap(Learner::ttable_counts_source_target[source_id]);
	  
      if (Learner::aligned_source_target.exists(source_id) && ! Learner::aligned_source_target[source_id].empty())
	counts.aligned.swap(Learner::aligned_source_target[source_id]);
	  
      if (! counts.counts.empty() || ! counts.aligned.empty()) {
	counts.word = word_type(source_id);
	
	queue_ttable_source_target[source_id % queue_ttable_source_target.size()].push_swap(counts);
      }
    }
	
    for (word_type::id_type target_id = 0; target_id != target_max; ++ target_id) {
      ttable_counts_type counts;
      
      if (Learner::ttable_counts_target_source.exists(target_id) && ! Learner::ttable_counts_target_source[target_id].empty())
	counts.counts.swap(Learner::ttable_counts_target_source[target_id]);
      
      if (Learner::aligned_target_source.exists(target_id) && ! Learner::aligned_target_source[target_id].empty())
	counts.aligned.swap(Learner::aligned_target_source[target_id]);
	  
      if (! counts.counts.empty() || ! counts.aligned.empty()) {
	counts.word = word_type(target_id);
	
	queue_ttable_target_source[target_id % queue_ttable_target_source.size()].push_swap(counts);
      }
    }
	
    Learner::ttable_counts_source_target.clear();
    Learner::ttable_counts_target_source.clear();
    Learner::aligned_source_target.clear();
    Learner::aligned_target_source.clear();
  }
  
  queue_bitext_type& queue_bitext;
  queue_ttable_set_type& queue_ttable_source_target;
  queue_ttable_set_type& queue_ttable_target_source;
};

template <typename TableSet, typename Table>
inline
void merge_tables(TableSet& tables, Table& merged)
{
  merged.clear();

  size_t size_max = 0;
  for (size_t i = 0; i != tables.size(); ++ i)
    size_max = utils::bithack::max(size_max, tables[i].size());

  merged.reserve(size_max);
  merged.resize(size_max);
  
  for (size_t i = 0; i != tables.size(); ++ i)
    for (word_type::id_type word_id = 0; word_id != tables[i].size(); ++ word_id)
      if (tables[i].exists(word_id))
	merged[word_type(word_id)].swap(tables[i][word_type(word_id)]);
}

template <typename Learner, typename Maximizer>
void learn(const Maximizer& maximizer,
	   const int iteration,
	   ttable_type& ttable_source_target,
	   ttable_type& ttable_target_source,
	   atable_type& atable_source_target,
	   atable_type& atable_target_source,
	   const classes_type& classes_source,
	   const classes_type& classes_target,
	   aligned_type& aligned_source_target,
	   aligned_type& aligned_target_source)
{
  typedef LearnMapReduce map_reduce_type;
  typedef LearnMapper<Learner> mapper_type;
  typedef LearnReducer<Maximizer> reducer_type;
  
  typedef map_reduce_type::bitext_type     bitext_type;
  typedef map_reduce_type::bitext_set_type bitext_set_type;
  
  typedef map_reduce_type::queue_bitext_type      queue_bitext_type;
  typedef map_reduce_type::queue_ttable_type      queue_ttable_type;
  typedef map_reduce_type::queue_ttable_set_type  queue_ttable_set_type;
  
  typedef std::vector<mapper_type, std::allocator<mapper_type> > mapper_set_type;
  
  queue_bitext_type queue_bitext(threads * 64);
  queue_ttable_set_type queue_ttable_source_target(utils::bithack::max(1, threads / 2));
  queue_ttable_set_type queue_ttable_target_source(utils::bithack::max(1, threads / 2));
  
  mapper_set_type mappers(threads, mapper_type(queue_bitext,
					       queue_ttable_source_target,
					       queue_ttable_target_source,
					       LearnBase(ttable_source_target, ttable_target_source,
							 atable_source_target, atable_target_source,
							 classes_source, classes_target)));
  
  etable_type etable_source_target(length_source_target);
  etable_type etable_target_source(length_target_source);
  
  for (int iter = 0; iter < iteration; ++ iter) {
    if (debug)
      std::cerr << "iteration: " << (iter + 1) << std::endl;
    
    utils::resource accumulate_start;

    std::vector<ttable_type, std::allocator<ttable_type> > ttable_source_target_new(threads, ttable_type(ttable_source_target.prior,
													 ttable_source_target.smooth));
    std::vector<ttable_type, std::allocator<ttable_type> > ttable_target_source_new(threads, ttable_type(ttable_target_source.prior,
													 ttable_target_source.smooth));
    
    std::vector<aligned_type, std::allocator<aligned_type> > aligned_source_target_new(threads);
    std::vector<aligned_type, std::allocator<aligned_type> > aligned_target_source_new(threads);
    
    boost::thread_group workers_mapper;
    boost::thread_group workers_reducer_source_target;
    boost::thread_group workers_reducer_target_source;
    
    for (size_t i = 0; i != mappers.size(); ++ i)
      workers_mapper.add_thread(new boost::thread(boost::ref(mappers[i])));
    
    for (size_t i = 0; i != queue_ttable_source_target.size(); ++ i)
      workers_reducer_source_target.add_thread(new boost::thread(reducer_type(queue_ttable_source_target[i],
									      ttable_source_target,
									      aligned_source_target,
									      ttable_source_target_new[i],
									      aligned_source_target_new[i],
									      maximizer)));
    for (size_t i = 0; i != queue_ttable_target_source.size(); ++ i)
      workers_reducer_target_source.add_thread(new boost::thread(reducer_type(queue_ttable_target_source[i],
									      ttable_target_source,
									      aligned_target_source,
									      ttable_target_source_new[i],
									      aligned_target_source_new[i],
									      maximizer)));
    
    utils::compress_istream is_src(source_file, 1024 * 1024);
    utils::compress_istream is_trg(target_file, 1024 * 1024);
    std::auto_ptr<std::istream> is_align(! alignment_file.empty()
					 ? new utils::compress_istream(alignment_file, 1024 * 1024) : 0);
    
    bitext_type     bitext;
    bitext_set_type bitexts;
    
    size_t num_bitext = 0;
    size_t length_source = 0;
    size_t length_target = 0;
    double objective_source_target = 0.0;
    double objective_target_source = 0.0;
    
    for (;;) {
      is_src >> bitext.source;
      is_trg >> bitext.target;
      if (is_align.get())
	*is_align >> bitext.alignment;
      
      if (! is_src || ! is_trg || (is_align.get() && ! *is_align)) break;
      
      if (bitext.source.empty() || bitext.target.empty()) continue;
      
      length_source += bitext.source.size();
      length_target += bitext.target.size();
      
      objective_source_target += (utils::mathop::log(etable_source_target(bitext.target.size(), bitext.source.size()))
				  / bitext.target.size());
      objective_target_source += (utils::mathop::log(etable_target_source(bitext.source.size(), bitext.target.size()))
				  / bitext.source.size());

      bitexts.push_back(bitext);

      ++ num_bitext;
      if (debug) {
	if (num_bitext % DEBUG_DOT == 0)
	  std::cerr << '.';
	if (num_bitext % DEBUG_LINE == 0)
	  std::cerr << '\n';
      }
      
      if (bitexts.size() == 64) {
	queue_bitext.push_swap(bitexts);
	bitexts.clear();
      }
    }

    if (! bitexts.empty())
      queue_bitext.push_swap(bitexts);
    
    if (debug && ((num_bitext / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;
    if (debug)
      std::cerr << "# of bitexts: " << num_bitext << std::endl;

    if (is_src || is_trg || (is_align.get() && *is_align))
      throw std::runtime_error("# of samples do not match");
    
    for (size_t i = 0; i != mappers.size(); ++ i) {
      bitexts.clear();
      queue_bitext.push_swap(bitexts);
    }
    
    workers_mapper.join_all();
    
    // send termination to reducer by sending nulls
    for (size_t i = 0; i != queue_ttable_source_target.size(); ++ i)
      queue_ttable_source_target[i].push(queue_ttable_type::value_type());
    
    for (size_t i = 0; i != queue_ttable_target_source.size(); ++ i)
      queue_ttable_target_source[i].push(queue_ttable_type::value_type());
    
    objective_source_target /= num_bitext;
    objective_target_source /= num_bitext;

    for (size_t i = 0; i != mappers.size(); ++ i) {
      objective_source_target += mappers[i].objective_source_target / num_bitext;
      objective_target_source += mappers[i].objective_target_source / num_bitext;
    }
    
    if (debug)
      std::cerr << "P(target | source) log-likelihood: " << objective_source_target << '\n'
		<< "                          entropy: " << std::exp(- objective_source_target) << '\n'
		<< "                       perplexity: " << (- objective_source_target / std::log(2.0)) << '\n'
		<< "P(source | target) log-likelihood: " << objective_target_source << '\n'
		<< "                          entropy: " << std::exp(- objective_target_source) << '\n'
		<< "                       perplexity: " << (- objective_target_source / std::log(2.0)) << '\n';
    
    // merge atable counts... (we will dynamically create probability table!)
    // first, initialize all the threaded alignment
    // we need to initialize cache-local atable, first, since we will clear all the entries again!
    atable_source_target.initialize();
    atable_target_source.initialize();
    
    // second, merge counts
    for (size_t i = 0; i != mappers.size(); ++ i) {
      atable_source_target += mappers[i].atable_counts_source_target;
      atable_target_source += mappers[i].atable_counts_target_source;
    }
    
    // third, estimate unk
    atable_source_target.estimate_unk();
    atable_target_source.estimate_unk();

    // initialize cache
    atable_source_target.initialize_cache();
    atable_target_source.initialize_cache();

    // etable..
    length_source_target = double(length_target) / length_source;
    length_target_source = double(length_source) / length_target;
    
    etable_source_target.assign(length_source_target);
    etable_target_source.assign(length_target_source);
    
    workers_reducer_source_target.join_all();
    workers_reducer_target_source.join_all();
    
    // merge ttable and aligned...
    merge_tables(ttable_source_target_new, ttable_source_target);
    merge_tables(ttable_target_source_new, ttable_target_source);
    merge_tables(aligned_source_target_new, aligned_source_target);
    merge_tables(aligned_target_source_new, aligned_target_source);
    
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
    bitext_type() : id(size_type(-1)), source(), target() {}
    bitext_type(const size_type& __id, const sentence_type& __source, const sentence_type& __target)
      : id(__id), source(__source), target(__target) {}
    
    bitext_type(const size_type& __id,
		const sentence_type& __source, const sentence_type& __target,
		const span_set_type& __span_source, const span_set_type& __span_target)
      : id(__id),
	source(__source), target(__target),
	span_source(__span_source), span_target(__span_target) {}
    
    size_type     id;
    sentence_type source;
    sentence_type target;
    span_set_type span_source;
    span_set_type span_target;

    void clear()
    {
      id = size_type(-1);
      source.clear();
      target.clear();
      span_source.clear();
      span_target.clear();
    }
    
    void swap(bitext_type& x)
    {
      std::swap(id, x.id);
      source.swap(x.source);
      target.swap(x.target);
      span_source.swap(x.span_source);
      span_target.swap(x.span_target);
    }
    
    friend
    bool operator<(const bitext_type& x, const bitext_type& y)
    {
      return x.id < y.id;
    }

    friend
    bool operator>(const bitext_type& x, const bitext_type& y)
    {
      return x.id > y.id;
    }
  };
  
  struct viterbi_type
  {
    size_type id;
    std::string output;
    
    viterbi_type() : id(size_type(-1)), output() {}
    
    void clear()
    {
      id = size_type(-1);
      output.clear();
    }

    void swap(viterbi_type& x)
    {
      std::swap(id, x.id);
      output.swap(x.output);
    }

    friend
    bool operator<(const viterbi_type& x, const viterbi_type& y)
    {
      return x.id < y.id;
    }
    
    friend
    bool operator>(const viterbi_type& x, const viterbi_type& y)
    {
      return x.id > y.id;
    }
  };

  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> > queue_mapper_type;
  typedef utils::lockfree_list_queue<viterbi_type, std::allocator<viterbi_type> > queue_reducer_type;
  
  typedef std::vector<queue_reducer_type, std::allocator<queue_reducer_type> > queue_reducer_set_type;
};

namespace std
{
  inline
  void swap(ViterbiMapReduce::bitext_type& x, ViterbiMapReduce::bitext_type& y)
  {
    x.swap(y);
  }

  inline
  void swap(ViterbiMapReduce::viterbi_type& x, ViterbiMapReduce::viterbi_type& y)
  {
    x.swap(y);
  }
};

template <typename Aligner>
struct ViterbiMapper : public ViterbiMapReduce, public Aligner
{
  queue_mapper_type& mapper;
  queue_reducer_type& reducer_source_target;
  queue_reducer_type& reducer_target_source;
  
  ViterbiMapper(const Aligner& __aligner,
		queue_mapper_type& __mapper,
		queue_reducer_type& __reducer_source_target,
		queue_reducer_type& __reducer_target_source)
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
    
    viterbi_type viterbi_source_target;
    viterbi_type viterbi_target_source;

    const int iter_mask = (1 << 8) - 1;
    
    for (int iter = 0;; ++ iter) {
      mapper.pop_swap(bitext);
      if (bitext.id == size_type(-1)) break;
      
      alignment_source_target.clear();
      alignment_target_source.clear();
      
      if (! bitext.source.empty() && ! bitext.target.empty())
	Aligner::operator()(bitext.source, bitext.target, bitext.span_source, bitext.span_target, alignment_source_target, alignment_target_source);
      
      viterbi_source_target.id = bitext.id;
      viterbi_target_source.id = bitext.id;
      
      write(viterbi_source_target, bitext.id, bitext.source, bitext.target, alignment_source_target);
      write(viterbi_target_source, bitext.id, bitext.target, bitext.source, alignment_target_source);
      
      reducer_source_target.push_swap(viterbi_source_target);
      reducer_target_source.push_swap(viterbi_target_source);

      if ((iter & iter_mask) == iter_mask)
	Aligner::shrink();
    }

    reducer_source_target.push(viterbi_type());
    reducer_target_source.push(viterbi_type());
  }

  typedef int index_type;
  typedef std::vector<index_type, std::allocator<index_type> > index_set_type;
  typedef std::vector<index_set_type, std::allocator<index_set_type> > align_set_type;
  typedef std::set<index_type, std::less<index_type>, std::allocator<index_type> > align_none_type;

  typedef std::vector<char, std::allocator<char> > buffer_type;
  
  align_set_type  aligns;
  align_none_type aligns_none;

  buffer_type buffer;

  void write(viterbi_type& viterbi,
	     const size_type& id,
	     const sentence_type& source,
	     const sentence_type& target,
	     const alignment_type& alignment)
  {
    typedef std::ostream_iterator<char> oiter_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;

    buffer.clear();
    
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::back_inserter(buffer));
    
    if (moses_mode)
      os << alignment << '\n';
    else {
      os << "# Sentence pair (" << (id + 1) << ')'
	 << " source length " << source.size()
	 << " target length " << target.size()
	 << " alignment score : " << 0 << '\n';
      os << target << '\n';
    
      if (source.empty() || target.empty() || alignment.empty()) {
	os << "NULL";
	os << " ({";
	for (size_type trg = 0; trg != target.size(); ++ trg)
	  os << ' ' << (trg + 1);
	os << " })";
	
	for (size_type src = 0; src != source.size(); ++ src) {
	  os << ' ' << source[src];
	  os << " ({ })";
	}
	os << '\n';
      } else {
	aligns.clear();
	aligns.resize(source.size());
      
	aligns_none.clear();
	for (size_type trg = 0; trg != target.size(); ++ trg)
	  aligns_none.insert(trg + 1);
      
	alignment_type::const_iterator aiter_end = alignment.end();
	for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
	  aligns[aiter->source].push_back(aiter->target + 1);
	  aligns_none.erase(aiter->target + 1);
	}
      
	os << "NULL";
	os << " ({ ";
	std::copy(aligns_none.begin(), aligns_none.end(), std::ostream_iterator<index_type>(os, " "));
	os << "})";
      
	for (size_type src = 0; src != source.size(); ++ src) {
	  os << ' ' << source[src];
	  os << " ({ ";
	  std::copy(aligns[src].begin(), aligns[src].end(), std::ostream_iterator<index_type>(os, " "));
	  os << "})";
	}
	os << '\n';
      }
    }
    
    os.reset();
    
    viterbi.output = std::string(buffer.begin(), buffer.end());
  }
};

struct ViterbiReducer : public ViterbiMapReduce
{
  typedef std::pair<viterbi_type, queue_reducer_type*> viterbi_queue_type;
  typedef std::vector<viterbi_queue_type, std::allocator<viterbi_queue_type> > viterbi_queue_set_type;

  struct heap_compare_type
  {
    bool operator()(const viterbi_queue_type* x, const viterbi_queue_type* y) const
    {
      return x->first > y->first;
    }
  };
  
  typedef std::vector<viterbi_queue_type*, std::allocator<viterbi_queue_type*> > heap_base_type;
  typedef std::priority_queue<viterbi_queue_type*, heap_base_type, heap_compare_type> heap_type;
  
  typedef boost::shared_ptr<std::ostream> ostream_ptr_type;
  
  ostream_ptr_type os;
  queue_reducer_set_type& queues;
  bool flush_;
  
  ViterbiReducer(const path_type& path, queue_reducer_set_type& __queues) : os(), queues(__queues), flush_(false)
  {
    if (! path.empty()) {
      flush_ = (path == "-"
		|| (boost::filesystem::exists(path)
		    && ! boost::filesystem::is_regular_file(path)));
      
      os.reset(new utils::compress_ostream(path, 1024 * 1024));
    }
  }
  
  void operator()()
  {
    std::ostream* stream = os.get();

    viterbi_queue_set_type viterbis(queues.size());
    heap_type heap;

    size_type id = 0;
    
    for (size_type shard = 0; shard != queues.size(); ++ shard) {
      viterbis[shard].second = &queues[shard];
      
      queues[shard].pop_swap(viterbis[shard].first);
      
      if (viterbis[shard].first.id != size_type(-1))
	heap.push(&viterbis[shard]);
    }
    
    while (! heap.empty()) {
      viterbi_queue_type* viterbi_queue = heap.top();
      heap.pop();

      if (viterbi_queue->first.id != id)
	throw std::runtime_error("invalid id");
      ++ id;
      
      if (stream) {
	*stream << viterbi_queue->first.output;
	
	if (flush_)
	  *stream << std::flush;
      }
      
      viterbi_queue->second->pop_swap(viterbi_queue->first);
      
      if (viterbi_queue->first.id != size_type(-1))
	heap.push(viterbi_queue);
    }
  }
};

template <typename Aligner>
void viterbi(const ttable_type& ttable_source_target,
	     const ttable_type& ttable_target_source,
	     const atable_type& atable_source_target,
	     const atable_type& atable_target_source,
	     const classes_type& classes_source,
	     const classes_type& classes_target)
{
  typedef ViterbiReducer         reducer_type;
  typedef ViterbiMapper<Aligner> mapper_type;
  
  typedef reducer_type::bitext_type    bitext_type;

  typedef reducer_type::queue_mapper_type      queue_mapper_type;
  typedef reducer_type::queue_reducer_type     queue_reducer_type;
  typedef reducer_type::queue_reducer_set_type queue_reducer_set_type;
  
  typedef std::vector<mapper_type, std::allocator<mapper_type> > mapper_set_type;

  queue_mapper_type     queue(threads * 1024);
  queue_reducer_set_type queue_source_target(threads, queue_reducer_type(1024));
  queue_reducer_set_type queue_target_source(threads, queue_reducer_type(1024));
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(viterbi_source_target_file, queue_source_target)));
  reducer.add_thread(new boost::thread(reducer_type(viterbi_target_source_file, queue_target_source)));

  boost::thread_group mapper;
  for (int i = 0; i != threads; ++ i)
    mapper.add_thread(new boost::thread(mapper_type(Aligner(ttable_source_target, ttable_target_source,
							    atable_source_target, atable_target_source,
							    classes_source, classes_target),
						    queue,
						    queue_source_target[i],
						    queue_target_source[i])));
    
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

    bitext.span_source.clear();
    bitext.span_target.clear();
    
    if (is_span_src.get())
      *is_span_src >> bitext.span_source;
    if (is_span_trg.get())
      *is_span_trg >> bitext.span_target;
    
    if (! is_src || ! is_trg || (is_span_src.get() && ! *is_span_src) || (is_span_trg.get() && ! *is_span_trg)) break;
    
    queue.push(bitext);
    
    ++ bitext.id;
    
    if (debug) {
      if (bitext.id % DEBUG_DOT == 0)
	std::cerr << '.';
      if (bitext.id % DEBUG_LINE == 0)
	std::cerr << '\n';
    }
  }
  
  if (debug && ((bitext.id / DEBUG_DOT) % DEBUG_WRAP))
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
  reducer.join_all();

  utils::resource viterbi_end;
  
  if (debug)
    std::cerr << "cpu time:  " << viterbi_end.cpu_time() - viterbi_start.cpu_time() << std::endl
	      << "user time: " << viterbi_end.user_time() - viterbi_start.user_time() << std::endl;
}

struct ProjectionMapReduce
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  struct bitext_type
  {
    size_type       id;
    sentence_type   source;
    sentence_type   target;
    dependency_type dependency_source;
    dependency_type dependency_target;

    bitext_type() : id(size_type(-1)), source(), target(), dependency_source(), dependency_target() {}
    
    void clear()
    {
      id = size_type(-1);
      source.clear();
      target.clear();
      dependency_source.clear();
      dependency_target.clear();
    }
    
    void swap(bitext_type& x)
    {
      std::swap(id, x.id);
      source.swap(x.source);
      target.swap(x.target);
      dependency_source.swap(x.dependency_source);
      dependency_target.swap(x.dependency_target);
    }
  };
  
  struct projected_type
  {
    size_type id;
    dependency_type dependency;
    
    projected_type() : id(size_type(-1)), dependency() {}
    
    void clear()
    {
      id = size_type(-1);
      dependency.clear();
    }
    
    void swap(projected_type& x)
    {
      std::swap(id, x.id);
      dependency.swap(x.dependency);
    }
  };
  

  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> >       queue_mapper_type;
  typedef utils::lockfree_list_queue<projected_type, std::allocator<projected_type> > queue_reducer_type;
};

namespace std
{
  inline
  void swap(ProjectionMapReduce::bitext_type& x, ProjectionMapReduce::bitext_type& y)
  {
    x.swap(y);
  }
  
  inline
  void swap(ProjectionMapReduce::projected_type& x, ProjectionMapReduce::projected_type& y)
  {
    x.swap(y);
  }
};

template <typename Analyzer>
struct ProjectionMapper : public ProjectionMapReduce, public Analyzer
{
  queue_mapper_type& mapper;
  queue_reducer_type& reducer_source;
  queue_reducer_type& reducer_target;
  
  ProjectionMapper(const Analyzer& __analyzer,
		   queue_mapper_type& __mapper,
		   queue_reducer_type& __reducer_source,
		   queue_reducer_type& __reducer_target)
    : Analyzer(__analyzer),
      mapper(__mapper),
      reducer_source(__reducer_source),
      reducer_target(__reducer_target) {}
  
  void operator()()
  {
    bitext_type bitext;
    projected_type projected_source;
    projected_type projected_target;
    
    const int iter_mask = (1 << 8) - 1;
    
    for (int iter = 0; /**/; ++ iter) {
      mapper.pop_swap(bitext);
      if (bitext.id == size_type(-1)) break;

      projected_source.clear();
      projected_target.clear();
      
      if (! bitext.source.empty() && ! bitext.target.empty())
	Analyzer::operator()(bitext.source,
			     bitext.target,
			     bitext.dependency_source,
			     bitext.dependency_target,
			     projected_source.dependency,
			     projected_target.dependency);
      
      projected_source.id = bitext.id;
      projected_target.id = bitext.id;
      
      reducer_source.push_swap(projected_source);
      reducer_target.push_swap(projected_target);
      
      if ((iter & iter_mask) == iter_mask)
	Analyzer::shrink();
    }
  }
  
};

struct ProjectionReducer : public ProjectionMapReduce
{
  struct less_projected
  {
    bool operator()(const projected_type& x, const projected_type& y) const
    {
      return x.id < y.id;
    }
  };
  typedef std::set<projected_type, less_projected, std::allocator<projected_type> > projected_set_type;
  
  typedef boost::shared_ptr<std::ostream> ostream_ptr_type;
  
  ostream_ptr_type    os;
  queue_reducer_type& queue;
  bool flush_;
  
  ProjectionReducer(const path_type& path, queue_reducer_type& __queue) : os(), queue(__queue), flush_(false)
  {
    if (! path.empty()) {
      flush_ = (path == "-"
		|| (boost::filesystem::exists(path)
		    && ! boost::filesystem::is_regular_file(path)));
      
      os.reset(new utils::compress_ostream(path, 1024 * 1024));
    }
  }
  
  void operator()()
  {
    if (os) {
      projected_type projected;
      for (;;) {
	queue.pop_swap(projected);
	if (projected.id == size_type(-1)) break;
      }
    } else { 
      size_type id = 0;
      projected_type     projected;
      projected_set_type buffer;
      for (;;) {
	queue.pop_swap(projected);
	if (projected.id == size_type(-1)) break;

	bool written = false;
       
	if (projected.id == id) {
	  *os << projected.dependency << '\n';
	  written = true;
	  ++ id;
	} else
	  buffer.insert(projected);
       
	while (! buffer.empty() && buffer.begin()->id == id) {
	  *os << buffer.begin()->dependency << '\n';
	  bool written = false;
	  buffer.erase(buffer.begin());
	  ++ id;
	}
	
	if (written && flush_)
	  *os << std::flush;
      }
     
      while (! buffer.empty() && buffer.begin()->id == id) {
	*os << buffer.begin()->dependency << '\n';
	buffer.erase(buffer.begin());
	++ id;
      }
      
      if (flush_)
	*os << std::flush;
     
      if (! buffer.empty())
	throw std::runtime_error("error while writing dependency output?");
    }
  }
};

template <typename Analyzer>
void project_dependency(const ttable_type& ttable_source_target,
			const ttable_type& ttable_target_source,
			const atable_type& atable_source_target,
			const atable_type& atable_target_source,
			const classes_type& classes_source,
			const classes_type& classes_target)
{
  typedef ProjectionMapper<Analyzer> mapper_type;
  typedef ProjectionReducer          reducer_type;
  
  typedef reducer_type::bitext_type    bitext_type;
  typedef reducer_type::projected_type projected_type;
  
  typedef reducer_type::queue_mapper_type  queue_mapper_type;
  typedef reducer_type::queue_reducer_type queue_reducer_type;

  typedef std::vector<mapper_type, std::allocator<mapper_type> > mapper_set_type;
  
  queue_mapper_type  queue(threads * 4096);
  queue_reducer_type queue_source;
  queue_reducer_type queue_target;
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(projected_source_file, queue_source)));
  reducer.add_thread(new boost::thread(reducer_type(projected_target_file, queue_target)));

  boost::thread_group mapper;
  for (int i = 0; i != threads; ++ i)
    mapper.add_thread(new boost::thread(mapper_type(Analyzer(ttable_source_target, ttable_target_source,
							     atable_source_target, atable_target_source,
							     classes_source, classes_target),
						    queue,
						    queue_source,
						    queue_target)));
    
  bitext_type bitext;
  bitext.id = 0;

  utils::resource projection_start;
  
  utils::compress_istream is_src(source_file, 1024 * 1024);
  utils::compress_istream is_trg(target_file, 1024 * 1024);
  
  std::auto_ptr<std::istream> is_dep_src(! dependency_source_file.empty()
					 ? new utils::compress_istream(dependency_source_file, 1024 * 1024) : 0);
  std::auto_ptr<std::istream> is_dep_trg(! dependency_target_file.empty()
					 ? new utils::compress_istream(dependency_target_file, 1024 * 1024) : 0);
  
  for (;;) {
    is_src >> bitext.source;
    is_trg >> bitext.target;
    
    bitext.dependency_source.clear();
    bitext.dependency_target.clear();

    if (is_dep_src.get())
      *is_dep_src >> bitext.dependency_source;
    if (is_dep_trg.get())
      *is_dep_trg >> bitext.dependency_target;
    
    if (! is_src || ! is_trg || (is_dep_src.get() && ! *is_dep_src) || (is_dep_trg.get() && ! *is_dep_trg)) break;
    
    queue.push(bitext);
    
    ++ bitext.id;
    
    if (debug) {
      if (bitext.id % DEBUG_DOT == 0)
	std::cerr << '.';
      if (bitext.id % DEBUG_LINE == 0)
	std::cerr << '\n';
    }
  }
  
  if (debug && ((bitext.id / DEBUG_DOT) % DEBUG_WRAP))
    std::cerr << std::endl;
  if (debug)
    std::cerr << "# of bitexts: " << bitext.id << std::endl;

  if (is_src || is_trg || (is_dep_src.get() && *is_dep_src) || (is_dep_trg.get() && *is_dep_trg))
    throw std::runtime_error("# of samples do not match");
  
  for (int i = 0; i != threads; ++ i) {
    bitext.clear();
    queue.push_swap(bitext);
  }
  mapper.join_all();
  
  queue_source.push(projected_type());
  queue_target.push(projected_type());
  reducer.join_all();

  utils::resource projection_end;
  
  if (debug)
    std::cerr << "cpu time:  " << projection_end.cpu_time() - projection_start.cpu_time() << std::endl
	      << "user time: " << projection_end.user_time() - projection_start.user_time() << std::endl;
}


struct PosteriorMapReduce
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef utils::vector2<double, std::allocator<double> > matrix_type;
  
  struct bitext_type
  {
    size_type       id;
    sentence_type   source;
    sentence_type   target;

    bitext_type() : id(size_type(-1)), source(), target() {}
    
    void clear()
    {
      id = size_type(-1);
      source.clear();
      target.clear();
    }
    
    void swap(bitext_type& x)
    {
      std::swap(id, x.id);
      source.swap(x.source);
      target.swap(x.target);
    }
  };
  
  struct posterior_type
  {
    size_type id;
    std::string output;
    
    posterior_type() : id(size_type(-1)), output() {}
    
    void clear()
    {
      id = size_type(-1);
      output.clear();
    }
    
    void swap(posterior_type& x)
    {
      std::swap(id, x.id);
      output.swap(x.output);
    }
    
    friend
    bool operator<(const posterior_type& x, const posterior_type& y)
    {
      return x.id < y.id;
    }
    
    friend
    bool operator>(const posterior_type& x, const posterior_type& y)
    {
      return x.id > y.id;
    }
  };
  
  typedef utils::lockfree_list_queue<bitext_type, std::allocator<bitext_type> >       queue_mapper_type;
  typedef utils::lockfree_list_queue<posterior_type, std::allocator<posterior_type> > queue_reducer_type;
  
  typedef std::vector<queue_reducer_type, std::allocator<queue_reducer_type> > queue_reducer_set_type;
};

namespace std
{
  inline
  void swap(PosteriorMapReduce::bitext_type& x, PosteriorMapReduce::bitext_type& y)
  {
    x.swap(y);
  }
  
  inline
  void swap(PosteriorMapReduce::posterior_type& x, PosteriorMapReduce::posterior_type& y)
  {
    x.swap(y);
  }
};

template <typename Infer>
struct PosteriorMapper : public PosteriorMapReduce
{
  Infer infer;

  queue_mapper_type&  mapper;
  queue_reducer_type& reducer_source_target;
  queue_reducer_type& reducer_target_source;
  
  PosteriorMapper(const Infer& __infer,
		  queue_mapper_type& __mapper,
		  queue_reducer_type& __reducer_source_target,
		  queue_reducer_type& __reducer_target_source)
    : infer(__infer),
      mapper(__mapper),
      reducer_source_target(__reducer_source_target),
      reducer_target_source(__reducer_target_source)
  {}
  
  void operator()()
  {
    bitext_type    bitext;

    matrix_type matrix_source_target;
    matrix_type matrix_target_source;

    posterior_type posterior_source_target;
    posterior_type posterior_target_source;
    
    const int iter_mask = (1 << 8) - 1;
    
    for (int iter = 0; /**/; ++ iter) {
      mapper.pop_swap(bitext);
      if (bitext.id == size_type(-1)) break;

      matrix_source_target.clear();
      matrix_target_source.clear();
      
      if (! bitext.source.empty() && ! bitext.target.empty())
	infer(bitext.source, bitext.target, matrix_source_target, matrix_target_source);
      
      posterior_source_target.id = bitext.id;
      posterior_target_source.id = bitext.id;
      
      write(posterior_source_target, matrix_source_target);
      write(posterior_target_source, matrix_target_source);
      
      reducer_source_target.push_swap(posterior_source_target);
      reducer_target_source.push_swap(posterior_target_source);
      
      if ((iter & iter_mask) == iter_mask)
	infer.shrink();
    }
    
    reducer_source_target.push(posterior_type());
    reducer_target_source.push(posterior_type());
  }
  
  struct real_precision : boost::spirit::karma::real_policies<long double>
  {
    static unsigned int precision(long double) 
    { 
      return 20;
    }
  };
  
  typedef std::vector<char, std::allocator<char> > buffer_type;

  buffer_type buffer;

  void write(posterior_type& posterior, const matrix_type& matrix)
  {
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
   
    typedef std::back_insert_iterator<buffer_type> iterator_type;
    
    //karma::real_generator<long double, real_precision> real;
    utils::double_base64_generator<iterator_type> base64;

    buffer.clear();
    iterator_type iter(buffer);

    if (! matrix.empty()) {
      karma::generate(iter, karma::lit('('));
      for (size_type i = 0; i != matrix.size1(); ++ i) {
	if (i)
	  karma::generate(iter, karma::lit(", "));
	
	karma::generate(iter, '(' << (('B' << base64) % ", ") << ')',
			boost::make_iterator_range(matrix.begin(i), matrix.end(i)));
      }
      karma::generate(iter, karma::lit(')'));
    }
    karma::generate(iter, karma::lit('\n'));

    posterior.output = std::string(buffer.begin(), buffer.end());
  }
};

struct PosteriorReducer : public PosteriorMapReduce
{
  typedef std::pair<posterior_type, queue_reducer_type*> posterior_queue_type;
  typedef std::vector<posterior_queue_type, std::allocator<posterior_queue_type> > posterior_queue_set_type;

  struct heap_compare_type
  {
    bool operator()(const posterior_queue_type* x, const posterior_queue_type* y) const
    {
      return x->first > y->first;
    }
  };
  
  typedef std::vector<posterior_queue_type*, std::allocator<posterior_queue_type*> > heap_base_type;
  typedef std::priority_queue<posterior_queue_type*, heap_base_type, heap_compare_type> heap_type;
  
  typedef boost::shared_ptr<std::ostream> ostream_ptr_type;
  
  ostream_ptr_type os;
  queue_reducer_set_type& queues;
  bool flush_;
  
  PosteriorReducer(const path_type& path, queue_reducer_set_type& __queues) : os(), queues(__queues), flush_(false)
  {
    if (! path.empty()) {
      flush_ = (path == "-"
		|| (boost::filesystem::exists(path)
		    && ! boost::filesystem::is_regular_file(path)));
      
      os.reset(new utils::compress_ostream(path, 1024 * 1024));
      os->precision(20);
    }
  }
  
  void operator()() throw()
  {
    std::ostream* stream = os.get();

    posterior_queue_set_type posteriors(queues.size());
    heap_type heap;

    size_type id = 0;
    
    for (size_type shard = 0; shard != queues.size(); ++ shard) {
      posteriors[shard].second = &queues[shard];
      
      queues[shard].pop_swap(posteriors[shard].first);
      
      if (posteriors[shard].first.id != size_type(-1))
	heap.push(&posteriors[shard]);
    }
    
    while (! heap.empty()) {
      posterior_queue_type* posterior_queue = heap.top();
      heap.pop();

      if (posterior_queue->first.id != id)
	throw std::runtime_error("invalid id");
      ++ id;
      
      if (stream) {
	*stream << posterior_queue->first.output;
	
	if (flush_)
	  *stream << std::flush;
      }
      
      posterior_queue->second->pop_swap(posterior_queue->first);
      
      if (posterior_queue->first.id != size_type(-1))
	heap.push(posterior_queue);
    }
  }
};

template <typename Infer>
void posterior(const ttable_type& ttable_source_target,
	       const ttable_type& ttable_target_source,
	       const atable_type& atable_source_target,
	       const atable_type& atable_target_source,
	       const classes_type& classes_source,
	       const classes_type& classes_target)
{
  typedef PosteriorMapper<Infer> mapper_type;
  typedef PosteriorReducer       reducer_type;
  
  typedef reducer_type::bitext_type    bitext_type;
  typedef reducer_type::posterior_type posterior_type;
  
  typedef reducer_type::queue_mapper_type      queue_mapper_type;
  typedef reducer_type::queue_reducer_type     queue_reducer_type;
  typedef reducer_type::queue_reducer_set_type queue_reducer_set_type;
  
  typedef std::vector<mapper_type, std::allocator<mapper_type> > mapper_set_type;
  
  queue_mapper_type  queue(threads * 1024);
  queue_reducer_set_type queue_source_target(threads, queue_reducer_type(1024));
  queue_reducer_set_type queue_target_source(threads, queue_reducer_type(1024));
  
  boost::thread_group reducer;
  reducer.add_thread(new boost::thread(reducer_type(posterior_source_target_file, queue_source_target)));
  reducer.add_thread(new boost::thread(reducer_type(posterior_target_source_file, queue_target_source)));

  boost::thread_group mapper;
  for (int i = 0; i != threads; ++ i)
    mapper.add_thread(new boost::thread(mapper_type(Infer(ttable_source_target, ttable_target_source,
							  atable_source_target, atable_target_source,
							  classes_source, classes_target),
						    queue,
						    queue_source_target[i],
						    queue_target_source[i])));  
  
  bitext_type bitext;
  bitext.id = 0;
  
  utils::resource posterior_start;
  
  utils::compress_istream is_src(source_file, 1024 * 1024);
  utils::compress_istream is_trg(target_file, 1024 * 1024);
  
  for (;;) {
    is_src >> bitext.source;
    is_trg >> bitext.target;
    
    if (! is_src || ! is_trg) break;
    
    queue.push(bitext);
    
    ++ bitext.id;
    
    if (debug) {
      if (bitext.id % DEBUG_DOT == 0)
	std::cerr << '.';
      if (bitext.id % DEBUG_LINE == 0)
	std::cerr << '\n';
    }
  }
  
  if (debug && ((bitext.id / DEBUG_DOT) % DEBUG_WRAP))
    std::cerr << std::endl;
  if (debug)
    std::cerr << "# of bitexts: " << bitext.id << std::endl;

  if (is_src || is_trg)
    throw std::runtime_error("# of samples do not match");
  
  for (int i = 0; i != threads; ++ i) {
    bitext.clear();
    queue.push_swap(bitext);
  }
  mapper.join_all();
  reducer.join_all();

  utils::resource posterior_end;
  
  if (debug)
    std::cerr << "cpu time:  " << posterior_end.cpu_time() - posterior_start.cpu_time() << std::endl
	      << "user time: " << posterior_end.user_time() - posterior_start.user_time() << std::endl;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("source",    po::value<path_type>(&source_file),    "source file")
    ("target",    po::value<path_type>(&target_file),    "target file")
    ("alignment", po::value<path_type>(&alignment_file), "alignment file")
    
    ("dependency-source", po::value<path_type>(&dependency_source_file), "source dependency file")
    ("dependency-target", po::value<path_type>(&dependency_target_file), "target dependency file")
    
    ("span-source", po::value<path_type>(&span_source_file), "source span file")
    ("span-target", po::value<path_type>(&span_target_file), "target span file")

    ("classes-source", po::value<path_type>(&classes_source_file), "source classes file")
    ("classes-target", po::value<path_type>(&classes_target_file), "target classes file")
    
    ("lexicon-source-target", po::value<path_type>(&lexicon_source_target_file), "lexicon model for P(target | source)")
    ("lexicon-target-source", po::value<path_type>(&lexicon_target_source_file), "lexicon model for P(source | target)")
    ("alignment-source-target", po::value<path_type>(&alignment_source_target_file), "alignment model for P(target | source)")
    ("alignment-target-source", po::value<path_type>(&alignment_target_source_file), "alignment model for P(source | target)")
    ("length-source-target",  po::value<double>(&length_source_target)->default_value(length_source_target), "length model for P(target | source)")
    ("length-target-source",  po::value<double>(&length_target_source)->default_value(length_target_source), "length model for P(source | target)")
    
    ("output-lexicon-source-target", po::value<path_type>(&output_lexicon_source_target_file), "lexicon model output for P(target | source)")
    ("output-lexicon-target-source", po::value<path_type>(&output_lexicon_target_source_file), "lexicon model output for P(source | target)")
    ("output-alignment-source-target", po::value<path_type>(&output_alignment_source_target_file), "alignment model output for P(target | source)")
    ("output-alignment-target-source", po::value<path_type>(&output_alignment_target_source_file), "alignment model output for P(source | target)")
    
    ("viterbi-source-target", po::value<path_type>(&viterbi_source_target_file), "viterbi for P(target | source)")
    ("viterbi-target-source", po::value<path_type>(&viterbi_target_source_file), "viterbi for P(source | target)")

    ("projected-source", po::value<path_type>(&projected_source_file), "source dependnecy projected from target")
    ("projected-target", po::value<path_type>(&projected_target_file), "target dependency projected from source")
    
    ("posterior-source-target", po::value<path_type>(&posterior_source_target_file), "posterior for P(target | source)")
    ("posterior-target-source", po::value<path_type>(&posterior_target_source_file), "posterior for P(source | target)")

    ("iteration-model1", po::value<int>(&iteration_model1)->default_value(iteration_model1), "max Model1 iteration")
    ("iteration-hmm", po::value<int>(&iteration_hmm)->default_value(iteration_hmm), "max HMM iteration")
    
    ("symmetric",  po::bool_switch(&symmetric_mode),  "symmetric training")
    ("posterior",  po::bool_switch(&posterior_mode),  "posterior constrained training")
    ("variational-bayes", po::bool_switch(&variational_bayes_mode), "variational Bayes estimates")
    ("pgd",               po::bool_switch(&pgd_mode),               "projected gradient descent")
    
    ("itg",       po::bool_switch(&itg_mode),       "ITG alignment")
    ("max-match", po::bool_switch(&max_match_mode), "maximum matching alignment")
    ("moses",     po::bool_switch(&moses_mode),     "Moses alignment foramt")

    ("permutation", po::bool_switch(&permutation_mode), "permutation")
    ("hybrid",      po::bool_switch(&hybrid_mode),      "hybrid projective dependency parsing")
    ("degree2",     po::bool_switch(&degree2_mode),     "degree2 non-projective dependency parsing")
    ("mst",         po::bool_switch(&mst_mode),         "MST non-projective dependency parsing")
    ("single-root", po::bool_switch(&single_root_mode), "single root dependency")

    ("p0",             po::value<double>(&p0)->default_value(p0),                               "parameter for NULL alignment")
    ("prior-lexicon",  po::value<double>(&prior_lexicon)->default_value(prior_lexicon),         "Dirichlet prior for variational Bayes")
    ("smooth-lexicon", po::value<double>(&smooth_lexicon)->default_value(smooth_lexicon),       "smoothing parameter for uniform distribution")
    ("prior-alignment",  po::value<double>(&prior_alignment)->default_value(prior_alignment),   "Dirichlet prior for variational Bayes")
    ("smooth-alignment", po::value<double>(&smooth_alignment)->default_value(smooth_alignment), "smoothing parameter for uniform distribution")
    
    ("l0-alpha", po::value<double>(&l0_alpha)->default_value(l0_alpha), "L0 regularization")
    ("l0-beta",  po::value<double>(&l0_beta)->default_value(l0_beta),   "L0 regularization")
    
    ("threshold", po::value<double>(&threshold)->default_value(threshold), "write with beam-threshold (<= 0.0 implies no beam)")

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
