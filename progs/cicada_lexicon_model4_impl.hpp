//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_MODEL4_IMPL__HPP__
#define __CICADA_LEXICON_MODEL4_IMPL__HPP__ 1

#include <numeric>
#include <set>

#include "cicada_lexicon_impl.hpp"

#include "cicada/semiring/logprob.hpp"

#include "utils/vector2_aligned.hpp"
#include "utils/vector3_aligned.hpp"
#include "utils/mathop.hpp"
#include "utils/aligned_allocator.hpp"
#include "utils/vector_set.hpp"

#include "kuhn_munkres.hpp"
#include "itg_alignment.hpp"
#include "dependency_hybrid.hpp"
#include "dependency_degree2.hpp"
#include "dependency_mst.hpp"

struct LearnModel4 : public LearnBase
{
  struct AlignmentSet
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef alignment_type::index_type index_type;

    typedef std::pair<index_type, index_type> range_type;
    
    typedef std::vector<index_type, std::allocator<index_type> > aligns_type;
    
    typedef utils::vector_set<index_type, std::less<index_type>, std::allocator<index_type> > aligns_sorted_type;
    typedef std::vector<aligns_sorted_type, std::allocator<aligns_sorted_type> > aligns_map_type;
    typedef std::vector<size_type, std::allocator<size_type> > sum_type;

    AlignmentSet() {}
    AlignmentSet(const sentence_type& source,
		 const sentence_type& target,
		 const alignment_type& alignment)
    { assign(source, target, alignment); }
    
    void alignment(alignment_type& aligns) const
    {
      aligns.clear();
      for (index_type src = 1; src != static_cast<index_type>(mapped.size()); ++ src)
	if (! mapped[src].empty()) {
	  aligns_sorted_type::const_iterator aiter_end = mapped[src].end();
	  for (aligns_sorted_type::const_iterator aiter = mapped[src].begin(); aiter != aiter_end; ++ aiter)
	    aligns.push_back(std::make_pair(src - 1, *aiter - 1));
	}
    }

    void clear()
    {
      aligns.clear();
      mapped.clear();
      sums.clear();
    }

    void assign(const sentence_type& source,
		const sentence_type& target,
		const alignment_type& alignment)
    {
      clear();
      
      aligns.resize(target.size() + 1);
      mapped.resize(source.size() + 1);
      sums.resize(source.size() + 1);
      
      // first, compute one-to-many alignment...
      alignment_type::const_iterator aiter_end = alignment.end();
      for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter)
	aligns[aiter->target + 1] = aiter->source + 1;
      
      // then, compute inverse alignment...
      for (size_type trg = 1; trg != aligns.size(); ++ trg)
	mapped[aligns[trg]].insert(trg);
      
      // summation
      for (size_type src = 0; src != mapped.size(); ++ src)
	sums[src] = std::accumulate(mapped[src].begin(), mapped[src].end(), 0);
    }

    void move(const index_type j, const index_type i_next)
    {
      const index_type i_prev = aligns[j];
      
      sums[i_prev] -= j;
      sums[i_next] += j;
      
      mapped[i_prev].erase(j);
      mapped[i_next].insert(j);
      
      aligns[j] = i_next;
    }

    void swap(const index_type j1, const index_type j2)
    {
      const index_type i1 = aligns[j1];
      const index_type i2 = aligns[j2];
      
      sums[i1] += j2 - j1;
      sums[i2] += j1 - j2;
      
      mapped[i1].erase(j1);
      mapped[i1].insert(j2);
      
      mapped[i2].erase(j2);
      mapped[i2].insert(j1);
      
      std::swap(aligns[j2], aligns[j1]);
    }
    
    size_type fertility(index_type x) const
    {
      return mapped[x].size();
    }
    
    size_type center(index_type x) const
    {
      return (sums[x] ? (sums[x] + mapped[x].size() - 1) / mapped[x].size() : size_type(0));
    }
    
    index_type prev_cept(index_type x) const
    {
      if (! x) return 0;
      
      index_type pos = x - 1;
      while (pos && ! mapped[pos].size())
	-- pos;
      
      return pos;
    }
    
    index_type next_cept(index_type x) const
    {
      size_type pos = x + 1;
      while (pos < mapped.size() && ! mapped[pos].size())
	++ pos;
      
      return pos;
    }

    void shrink()
    {
      clear();
      
      aligns_type(aligns).swap(aligns);
      aligns_map_type(mapped).swap(mapped);
      sum_type(sums).swap(sums);
    }
    
    aligns_type     aligns;
    aligns_map_type mapped;
    sum_type        sums;
  };
  
  typedef AlignmentSet alignment_set_type;

  struct Model4
  {
    Model4(const ttable_type& __ttable,
	   const dtable_type& __dtable,
	   const ntable_type& __ntable,
	   const ptable_type& __ptable,
	   const classes_type& __classes_source,
	   const classes_type& __classes_target)
      : ttable(__ttable),
	dtable(__dtable),
	ntable(__ntable),
	ptable(__ptable),
	classes_source(__classes_source),
	classes_target(__classes_target) {}

    const ttable_type& ttable;
    const dtable_type& dtable;
    const ntable_type& ntable;
    const ptable_type& ptable;
    const classes_type& classes_source;
    const classes_type& classes_target;
  };
  
  typedef Model4 model4_type;

  struct Model4Data
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef alignment_set_type::index_type    index_type;
    typedef std::pair<index_type, index_type> range_type;
    
    typedef cicada::semiring::Logprob<double> logprob_type;
    
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > matrix_swap_type;
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > matrix_move_type;
    
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > ttable_cache_type;
    typedef utils::vector3_aligned<double, utils::aligned_allocator<double> > dtable_head_cache_type;
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > dtable_others_cache_type;
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > ntable_cache_type;

    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > move_score_type;
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > swap_score_type;
    
    typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > posterior_type;
    typedef std::vector<double, std::allocator<double> > posterior_accum_type;
    
    void assign(const sentence_type& source,
		const sentence_type& target,
		const classes_type& classes_source,
		const classes_type& classes_target,
		const alignment_type& alignment)
    {
      aligns.assign(source, target, alignment);
#if 0
      std::cerr << "align: ";
      std::copy(aligns.aligns.begin(), aligns.aligns.end(), std::ostream_iterator<int>(std::cerr, " "));
      std::cerr << std::endl;
#endif
      
      // classes...
      source_class.clear();
      target_class.clear();
      source_class.resize(source.size() + 1,  vocab_type::EPSILON);
      target_class.resize(target.size() + 1, vocab_type::EPSILON);
      
      sentence_type::iterator csiter = source_class.begin() + 1;
      for (sentence_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter, ++ csiter)
	*csiter = classes_source[*siter];
      
      sentence_type::iterator ctiter = target_class.begin() + 1;
      for (sentence_type::const_iterator titer = target.begin(); titer != target.end(); ++ titer, ++ ctiter)
	*ctiter = classes_target[*titer];      
    }

    void assign(const sentence_type& source,
		const sentence_type& target,
		const alignment_type& alignment,
		const model4_type& model4)
    {
      // alignment...
      aligns.assign(source, target, alignment);
      
      // classes...
      source_class.clear();
      target_class.clear();
      source_class.resize(source.size() + 1,  vocab_type::EPSILON);
      target_class.resize(target.size() + 1, vocab_type::EPSILON);
      
      sentence_type::iterator csiter = source_class.begin() + 1;
      for (sentence_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter, ++ csiter)
	*csiter = model4.classes_source[*siter];
      
      sentence_type::iterator ctiter = target_class.begin() + 1;
      for (sentence_type::const_iterator titer = target.begin(); titer != target.end(); ++ titer, ++ ctiter)
	*ctiter = model4.classes_target[*titer];

      //std::cerr << "ttable" << std::endl;
      
      // ttable..
      ttable.clear();
      ttable.reserve(target.size() + 1, source.size() + 1);
      ttable.resize(target.size() + 1, source.size() + 1);
      
      for (size_type trg = 0; trg != target.size(); ++ trg) {
	ttable(trg + 1, 0) = model4.ttable(vocab_type::EPSILON, target[trg]);
	for (size_type src = 0; src != source.size(); ++ src)
	  ttable(trg + 1, src + 1) = model4.ttable(source[src], target[trg]);
      }
      
      //std::cerr << "dtable" << std::endl;
	    
      // dtable...
      dtable_head.clear();
      dtable_head.reserve(source.size() + 1, target.size() + 1, target.size() + 1);
      dtable_head.resize(source.size() + 1, target.size() + 1, target.size() + 1);
      
      for (size_type src = 0; src <= source.size(); ++ src)
	for (int prev = (src != 0); prev <= static_cast<int>(target.size()); ++ prev)
	  for (int next = 1; next <= static_cast<int>(target.size()); ++ next) {
	    dtable_head(src, prev, next) = model4.dtable(source_class[src],
							 target_class[next],
							 source.size(),
							 target.size(),
							 prev,
							 next);
#if 0
	    std::cerr << "head: " << (src ? source[src - 1] : vocab_type::EPSILON)
		      << " prev: " << prev
		      << " next: " << next
		      << " " << dtable_head(src, prev, next)
		      << std::endl;
#endif
	  }
      
      dtable_others.clear();
      dtable_others.reserve(target.size() + 1, target.size() + 1);
      dtable_others.resize(target.size() + 1, target.size() + 1);
      
      for (int prev = 1; prev < static_cast<int>(target.size()); ++ prev)
	for (int next = prev + 1; next <= static_cast<int>(target.size()); ++ next) {
	  dtable_others(prev, next) = model4.dtable(target_class[next],
						    source.size(),
						    target.size(),
						    prev, 
						    next);
#if 0
	  std::cerr << "non-head"
		    << " prev: " << prev
		    << " next: " << next
		    << " " << dtable_others(prev, next)
		    << std::endl;
#endif
	}
      
      //std::cerr << "ntable" << std::endl;
      
      // ntable...
      ntable.clear();
      ntable.reserve(source.size() + 1, target.size() + 1);
      ntable.resize(source.size() + 1, target.size() + 1);
      
      for (size_type src = 0; src != source.size(); ++ src)
	for (int fertility = 0; fertility <= static_cast<int>(target.size()); ++ fertility) {
	  ntable(src + 1, fertility) = model4.ntable(source[src], target.size(), fertility);
	  
	  //std::cerr << "source: " << source[src] << " fert: "  << fertility << " " << ntable(src + 1, fertility) << std::endl;
	}
      
      //std::cerr << "ptable" << std::endl;

      // ptable
      for (int phi0 = 0; phi0 <= static_cast<int>(target.size()); ++ phi0) {
	ntable(0, phi0) = model4.ptable(target.size(), phi0);
	
	//std::cerr << "source: " << vocab_type::EPSILON << " fert: "  << phi0 << " " << ntable(0, phi0) << std::endl;
      }
      
      // compute logprob for the given alignment...
      update();
    }
    
    void update()
    {
      // first, insertion...
      logprob = ntable(0, aligns.fertility(0));

      //std::cerr << "ptable: " << logprob << std::endl;
      
      // second, fertility...
      for (size_type src = 1; src != aligns.mapped.size(); ++ src)
	logprob *= ntable(src, aligns.fertility(src));

      //std::cerr << "ntable: " << logprob << std::endl;
      
      // third, lexicon...
      for (size_type trg = 1; trg != aligns.aligns.size(); ++ trg)
	logprob *= ttable(trg, aligns.aligns[trg]);

      //std::cerr << "ttable: " << logprob << std::endl;
      
      // forth, distortion...
      logprob *= score_distortion<logprob_type>(0, aligns.mapped.size());

      //std::cerr << "logprob: " << logprob << std::endl;

      moves.clear();
      swaps.clear();
      
      moves.reserve(aligns.aligns.size(), aligns.mapped.size());
      swaps.reserve(aligns.aligns.size(), aligns.aligns.size());
      
      moves.resize(aligns.aligns.size(), aligns.mapped.size(), 0.0);
      swaps.resize(aligns.aligns.size(), aligns.aligns.size(), 0.0);
      
      // move/swap matrix...
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	  moves(j, i) = (aligns.aligns[j] != i ? score_move(j, i) : 1.0);
      
      for (index_type j1 = 1; j1 < static_cast<index_type>(aligns.aligns.size()) - 1; ++ j1)
	for (index_type j2 = j1 + 1; j2 < static_cast<index_type>(aligns.aligns.size()); ++ j2)
	  swaps(j1, j2) = (aligns.aligns[j1] != aligns.aligns[j2] ? score_swap(j1, j2) : 1.0);
    }
    
    void update_j(const index_type j)
    {
      for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	moves(j, i) = (aligns.aligns[j] != i ? score_move(j, i) : 1.0);
      
      for (index_type j2 = j + 1; j2 < static_cast<index_type>(aligns.aligns.size()); ++ j2)
	swaps(j, j2) = (aligns.aligns[j] != aligns.aligns[j2] ? score_swap(j, j2) : 1.0);
      
      for (index_type j1 = 1; j1 < j; ++ j1)
	swaps(j1, j) = (aligns.aligns[j1] != aligns.aligns[j] ? score_swap(j1, j) : 1.0);
    }
    
    void update_i(const index_type i)
    {
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	moves(j, i) = (aligns.aligns[j] != i ? score_move(j, i) : 1.0);
    }

    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    const alignment_type& alignment,
		    ttable_type& counts)
    {
      for (size_type trg = 1; trg != aligns.aligns.size(); ++ trg) {
	const word_type word_source(aligns.aligns[trg] ? source[aligns.aligns[trg] - 1] : vocab_type::EPSILON);
	const word_type word_target(target[trg - 1]);
	
	++ counts[word_source][word_target];
      }
    }

    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    const alignment_type& alignment,
		    ntable_type& counts)
    {
      for (size_type src = 1; src != aligns.mapped.size(); ++ src) {
	//std::cerr << "fert: " << src << " " << aligns.fertility(src) << std::endl;
	++ counts(source[src - 1], aligns.fertility(src));
      }
    }
    
    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    const alignment_type& alignmeent,
		    dtable_counts_type& counts)
    {
      size_type cept_prev = 0;
      for (size_type cept = 1; cept != aligns.mapped.size(); ++ cept) 
	if (aligns.fertility(cept)) {
	  //std::cerr << "cept: " << cept << " fettility: " << aligns.fertility(cept) << std::endl;
	  
	  const size_type center_prev = aligns.center(cept_prev);
	  const size_type head = aligns.mapped[cept].front();
#if 0
	  std::cerr << "prev: " << cept_prev
		    << " center: " << center_prev
		    << " head: " << head
		    << std::endl;
#endif
	  
	  ++ counts[std::make_pair(source_class[cept_prev], target_class[head])][head - center_prev];
	  
	  alignment_set_type::aligns_sorted_type::const_iterator jiter_end = aligns.mapped[cept].end();
	  for (alignment_set_type::aligns_sorted_type::const_iterator jiter = aligns.mapped[cept].begin() + 1; jiter != jiter_end; ++ jiter)
	    ++ counts[std::make_pair(vocab_type::NONE, target_class[*jiter])][*jiter - *(jiter - 1)];
	  
	  cept_prev = cept;
	}
    }

    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    ttable_type& counts)
    {
      const int source_size = source.size();
      const int target_size = target.size();
      
      {
	// allocate enough buffer size...
	
	for (int src = 1; src <= source_size; ++ src) {
	  ttable_type::count_map_type& mapped = counts[source[src - 1]];
	  mapped.rehash(mapped.size() + target_size);
	}
	
	ttable_type::count_map_type& mapped = counts[vocab_type::EPSILON];
	mapped.rehash(mapped.size() + target_size);
      }
      
      for (int trg = 1; trg <= target_size; ++ trg) {
	counts[vocab_type::EPSILON][target[trg - 1]] += posterior(trg, 0);
	
	for (int src = 1; src <= source_size; ++ src)
	  counts[source[src - 1]][target[trg - 1]] += posterior(trg, src);
      }
    }
    
    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    ntable_type& counts)
    {
      const int source_size = source.size();
      const int target_size = target.size();
      
      for (int src = 1; src <= source_size; ++ src)
	for (int fertility = 0; fertility <= target_size; ++ fertility)
	  if (posterior_fertility(src, fertility) > 0.0) 
	    counts(source[src - 1], fertility) += posterior_fertility(src, fertility);
    }

    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    dtable_counts_type& counts)
    {
      // we will explicitly enumerate all the alignment and collect counts... very inefficient...
      
      const double factor = 1.0 / total;

      accumulate(source, target, factor, counts);
      
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	  if (aligns.aligns[j] != i) {
	    const index_type i_prev = aligns.aligns[j];
	    aligns.move(j, i);
	    
	    accumulate(source, target, moves(j, i) * factor, counts);
	    
	    aligns.move(j, i_prev);
	  }
      
      for (index_type j1 = 1; j1 < static_cast<index_type>(aligns.aligns.size()) - 1; ++ j1)
	for (index_type j2 = j1 + 1; j2 < static_cast<index_type>(aligns.aligns.size()); ++ j2)
	  if (aligns.aligns[j1] != aligns.aligns[j2]) {
	    aligns.swap(j1, j2);
	    
	    accumulate(source, target, swaps(j1, j2) * factor, counts);
	    
	    aligns.swap(j1, j2);
	  }
    }
    
    void accumulate(const sentence_type& source,
		    const sentence_type& target,
		    const double count,
		    dtable_counts_type& counts)
    {
      size_type cept_prev = 0;
      for (size_type cept = 1; cept != aligns.mapped.size(); ++ cept) 
	if (aligns.fertility(cept)) {
	  const size_type center_prev = aligns.center(cept_prev);
	  const size_type head = aligns.mapped[cept].front();
	  
	  counts[std::make_pair(source_class[cept_prev], target_class[head])][head - center_prev] += count;
	  
	  alignment_set_type::aligns_sorted_type::const_iterator jiter_end = aligns.mapped[cept].end();
	  for (alignment_set_type::aligns_sorted_type::const_iterator jiter = aligns.mapped[cept].begin() + 1; jiter != jiter_end; ++ jiter)
	    counts[std::make_pair(vocab_type::NONE, target_class[*jiter])][*jiter - *(jiter - 1)] += count;
	  
	  cept_prev = cept;
	}
    }

    double objective() const
    {
      return cicada::semiring::log(logprob * total);
    }
    
    void estimate_posterior(const sentence_type& __source,
			    const sentence_type& __target)
    {
      const int source_size = __source.size();
      const int target_size = __target.size();
      
      posterior.clear();
      posterior.reserve(target_size + 1, source_size + 1);
      posterior.resize(target_size + 1, source_size + 1, 0.0);

      posterior_swap.clear();
      posterior_swap.reserve(target_size + 1, source_size + 1);
      posterior_swap.resize(target_size + 1, source_size + 1, 0.0);

      posterior_fertility.clear();
      posterior_fertility.reserve(source_size + 1, target_size + 1);
      posterior_fertility.resize(source_size + 1, target_size + 1, 0.0);
      
      neighbour_move.clear();
      neighbour_swap.clear();
      
      neighbour_move.resize(target_size + 1, 0.0);
      neighbour_swap.resize(target_size + 1, 0.0);
      
      fertility_inc.clear();
      fertility_dec.clear();
      
      fertility_inc.resize(source_size + 1, 0.0);
      fertility_dec.resize(source_size + 1, 0.0);
      
      total = 1.0;
      
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	  if (aligns.aligns[j] != i) {
	    const double score = moves(j, i);
	    
	    total += score;
	    neighbour_move[j] += score;
	    fertility_inc[i] += score;
	    fertility_dec[aligns.aligns[j]] += score;
	  }
      
      for (index_type j1 = 1; j1 < static_cast<index_type>(aligns.aligns.size()) - 1; ++ j1)
	for (index_type j2 = j1 + 1; j2 < static_cast<index_type>(aligns.aligns.size()); ++ j2)
	  if (aligns.aligns[j1] != aligns.aligns[j2]) {
	    const double score = swaps(j1, j2);
	    
	    total += score;
	    neighbour_swap[j1] += score;
	    neighbour_swap[j2] += score;
	    posterior_swap(j1, aligns.aligns[j2]) += score;
	    posterior_swap(j2, aligns.aligns[j1]) += score;
	  }

      const double factor = 1.0 / total;
      
      for (index_type j = 1; j <= target_size; ++ j)
	for (index_type i = 0; i <= source_size; ++ i)
	  posterior(j, i) = (aligns.aligns[j] == i
			     ? total - (neighbour_move[j] + neighbour_swap[j])
			     : moves(j, i) + posterior_swap(j, i)) * factor;
      
      // fertility...
      for (index_type i = 1; i <= source_size; ++ i) {
	posterior_fertility(i, aligns.fertility(i)) += (total - (fertility_inc[i] + fertility_dec[i])) * factor;
	
	if (aligns.fertility(i))
	  posterior_fertility(i, aligns.fertility(i) - 1) += fertility_dec[i] * factor;
	else if (fertility_dec[i] > 0.0)
	  std::cerr << "invalid counts? at fertility_dec: " << fertility_dec[i] << std::endl;
	// otherwise... fail!
	
	if (aligns.fertility(i) + 1 <= static_cast<size_type>(target_size))
	  posterior_fertility(i, aligns.fertility(i) + 1) += fertility_inc[i] * factor;
	else if (fertility_inc[i] > 0.0)
	  std::cerr << "invalid counts? at fertility_inc: " << fertility_inc[i] << std::endl;
	// otherwise.. fail!
      }
    }
    
    bool climb()
    {
      double gain_move = 0.0;
      double gain_swap = 0.0;
      index_type move_j = 0;
      index_type move_i = 0;
      index_type swap_j1 = 0;
      index_type swap_j2 = 0;
      
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	  if (moves(j, i) > gain_move) {
	    gain_move = moves(j, i);
	    move_j = j;
	    move_i = i;
	  }
      
      for (index_type j1 = 1; j1 < static_cast<index_type>(aligns.aligns.size()) - 1; ++ j1)
	for (index_type j2 = j1 + 1; j2 < static_cast<index_type>(aligns.aligns.size()); ++ j2)
	  if (swaps(j1, j2) > gain_swap) {
	    gain_swap = swaps(j1, j2);
	    swap_j1 = j1;
	    swap_j2 = j2;
	  }
      
      if (gain_move <= 1.0 && gain_swap <= 1.0)
	return false;

      if (gain_move >= gain_swap) {
	//std::cerr << "moving: j=" << move_j << " i=" << move_i << " gain: " << gain_move << std::endl;

	move(move_j, move_i);
      } else {
	//std::cerr << "swapping: j1=" << swap_j1 << " j2=" << swap_j2 << " gain: " << gain_swap << std::endl;
	
	swap(swap_j1, swap_j2);
      }
      
      return true;
    }
    
    void move(const index_type j, const index_type i_next)
    {
      const index_type i_prev = aligns.aligns[j];
      
      if (i_prev == i_next) return;
      
      // is this correct...?
      const range_type range1(aligns.prev_cept(i_prev), utils::bithack::min(aligns.next_cept(i_prev),
									    static_cast<index_type>(aligns.mapped.size() - 1)));
      const range_type range2(aligns.prev_cept(i_next), utils::bithack::min(aligns.next_cept(i_next),
									    static_cast<index_type>(aligns.mapped.size() - 1)));
      
      aligns.move(j, i_next);
      logprob *= moves(j, i_next);

      const range_type range1_new(aligns.prev_cept(i_prev), utils::bithack::min(aligns.next_cept(i_prev),
										static_cast<index_type>(aligns.mapped.size() - 1)));
      const range_type range2_new(aligns.prev_cept(i_next), utils::bithack::min(aligns.next_cept(i_next),
										static_cast<index_type>(aligns.mapped.size() - 1)));

      if (range1.second < range2.first || range2.second < range1.first) {
	if (range1_new.first < range1.first || range1_new.second > range1.second)
	  std::cerr << "different prev range: j=" << j << " prev: " << i_prev << " next: " << i_next
		    << " [" << range1.first  << ", " << range1.second << "]"
		    << " [" << range1_new.first  << ", " << range1_new.second << "]"
		    << std::endl;
	if (range2_new.first < range2.first || range2_new.second > range2.second)
	  std::cerr << "different next range: j=" << j << " prev: " << i_prev << " next: " << i_next
		    << " [" << range2.first  << ", " << range2.second << "]"
		    << " [" << range2_new.first  << ", " << range2_new.second << "]"
		    << std::endl;
      } else {
	const range_type range_prev(utils::bithack::min(range1.first, range2.first),
				    utils::bithack::max(range1.second, range2.second));
	const range_type range_new(utils::bithack::min(range1_new.first, range2_new.first),
				   utils::bithack::max(range1_new.second, range2_new.second));
	
	if (range_new.first < range_prev.first || range_new.second > range_prev.second)
	  std::cerr << "different merged range: j=" << j << " prev: " << i_prev << " next: " << i_next
		    << " [" << range_prev.first  << ", " << range_prev.second << "]"
		    << " [" << range_new.first  << ", " << range_new.second << "]"
		    << std::endl;
      }
      
      std::vector<bool, std::allocator<bool> > modified(aligns.mapped.size(), false);
      
      for (index_type i = range1.first; i <= range1.second && i < static_cast<index_type>(aligns.mapped.size()); ++ i)
	modified[i] = true;
      for (index_type i = range2.first; i <= range2.second && i < static_cast<index_type>(aligns.mapped.size()); ++ i)
	modified[i] = true;
      
      for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	if (modified[i])
	  update_i(i);
      
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	if (modified[aligns.aligns[j]])
	  update_j(j);
    }
    
    void swap(const index_type j1, const index_type j2)
    {
      const index_type i1 = aligns.aligns[j1];
      const index_type i2 = aligns.aligns[j2];
      
      if (j1 == j2 || i1 == i2) return;
      
      const range_type range1(aligns.prev_cept(i1), aligns.next_cept(i1));
      const range_type range2(aligns.prev_cept(i2), aligns.next_cept(i2));
      
      aligns.swap(j1, j2);
      logprob *= swaps(utils::bithack::min(j1, j2), utils::bithack::max(j1, j2));
      
      std::vector<bool, std::allocator<bool> > modified(aligns.mapped.size(), false);
      
      for (index_type i = range1.first; i <= range1.second && i < static_cast<index_type>(aligns.mapped.size()); ++ i)
	modified[i] = true;
      for (index_type i = range2.first; i <= range2.second && i < static_cast<index_type>(aligns.mapped.size()); ++ i)
	modified[i] = true;
      
      for (index_type i = 0; i != static_cast<index_type>(aligns.mapped.size()); ++ i)
	if (modified[i])
	  update_i(i);
      
      for (index_type j = 1; j != static_cast<index_type>(aligns.aligns.size()); ++ j)
	if (modified[aligns.aligns[j]])
	  update_j(j);
    }
    
    double score_move(const index_type j, const index_type i_next)
    {
      const index_type i_prev = aligns.aligns[j];
      
      if (i_prev == i_next) return 1.0;
      
      double gain = 1.0;
      
      // ttable...
      gain *= ttable(j, i_next) / ttable(j, i_prev);
      
      // ptable and/or ntable
      gain *= ntable(i_prev, aligns.fertility(i_prev) - 1) / ntable(i_prev, aligns.fertility(i_prev));
      gain *= ntable(i_next, aligns.fertility(i_next) + 1) / ntable(i_next, aligns.fertility(i_next));
      
      // gain in dtable... is this correct...?
      const range_type range1(aligns.prev_cept(i_prev), aligns.next_cept(i_prev));
      const range_type range2(aligns.prev_cept(i_next), aligns.next_cept(i_next));
      
      if (range1.second < range2.first) {
	const double denom1 = score_distortion<double>(range1.first, range1.second);
	const double denom2 = score_distortion<double>(range2.first, range2.second);
	
	aligns.move(j, i_next);
	
	const double numer1 = score_distortion<double>(range1.first, range1.second);
	const double numer2 = score_distortion<double>(range2.first, range2.second);
	
	aligns.move(j, i_prev);
	
	gain *= (numer1 / denom1) * (numer2 / denom2);
      } else if (range2.second < range1.first) {
	const double denom2 = score_distortion<double>(range2.first, range2.second);
	const double denom1 = score_distortion<double>(range1.first, range1.second);
	
	aligns.move(j, i_next);
	
	const double numer2 = score_distortion<double>(range2.first, range2.second);
	const double numer1 = score_distortion<double>(range1.first, range1.second);
	
	aligns.move(j, i_prev);
	
	gain *= (numer1 / denom1) * (numer2 / denom2);
      } else {
	const double denom = score_distortion<double>(utils::bithack::min(range1.first, range2.first), utils::bithack::max(range1.second, range2.second));

	aligns.move(j, i_next);
	
	const double numer = score_distortion<double>(utils::bithack::min(range1.first, range2.first), utils::bithack::max(range1.second, range2.second));
	
	aligns.move(j, i_prev);
	
	gain *= (numer / denom);
      }
      
      return gain;
    }
    
    double score_swap(const index_type j1, const index_type j2)
    {
      const index_type i1 = aligns.aligns[j1];
      const index_type i2 = aligns.aligns[j2];

      if (j1 == j2 || i1 == i2) return 1.0;
      
      double gain = 1.0;
      
      // lexicon model...
      gain *= ttable(j1, i2) / ttable(j1, i1);
      gain *= ttable(j2, i1) / ttable(j2, i2);
      
      // no change in ptable and/or ntable...
      
      // gain in dtable...
      const range_type range1(aligns.prev_cept(i1), aligns.next_cept(i1));
      const range_type range2(aligns.prev_cept(i2), aligns.next_cept(i2));

      if (range1.second < range2.first) {
	const double denom1 = score_distortion<double>(range1.first, range1.second);
	const double denom2 = score_distortion<double>(range2.first, range2.second);
	
	aligns.swap(j1, j2);

	const double numer1 = score_distortion<double>(range1.first, range1.second);
	const double numer2 = score_distortion<double>(range2.first, range2.second);
	
	aligns.swap(j1, j2);
	
	gain *= (numer1 / denom1) * (numer2 / denom2);
      } else if (range2.second < range1.first) {
	const double denom2 = score_distortion<double>(range2.first, range2.second);
	const double denom1 = score_distortion<double>(range1.first, range1.second);
	
	aligns.swap(j1, j2);
	
	const double numer2 = score_distortion<double>(range2.first, range2.second);
	const double numer1 = score_distortion<double>(range1.first, range1.second);
	
	aligns.swap(j1, j2);
	
	gain *= (numer1 / denom1) * (numer2 / denom2);
      } else {
	const double denom = score_distortion<double>(utils::bithack::min(range1.first, range2.first),
						      utils::bithack::max(range1.second, range2.second));
	
	aligns.swap(j1, j2);
	
	const double numer = score_distortion<double>(utils::bithack::min(range1.first, range2.first),
						      utils::bithack::max(range1.second, range2.second));
	
	aligns.swap(j1, j2);
	
	gain *= (numer / denom);
      }
      
      return gain;
    }
    
    template <typename Prob>
    Prob score_distortion(const index_type i1, const index_type i2)
    {
      Prob score = 1.0;
      
      index_type cept_prev = (i1 ? aligns.prev_cept(i1) : 0);
      for (index_type cept = utils::bithack::max(i1, 1); cept <= i2 && cept < static_cast<index_type>(aligns.mapped.size()); ++ cept) 
	if (aligns.fertility(cept)) {
	  const size_type center_prev = aligns.center(cept_prev);
	  
	  score *= dtable_head(cept_prev, center_prev, aligns.mapped[cept].front());
#if 0
	  std::cerr << "cept: " << cept
		    << " prev: " << cept_prev
		    << " center: " << center_prev
		    << " head: " << aligns.mapped[cept].front()
		    << " = " << dtable_head(cept_prev, center_prev, aligns.mapped[cept].front())
		    << std::endl;
#endif
	  
	  alignment_set_type::aligns_sorted_type::const_iterator jiter_end = aligns.mapped[cept].end();
	  for (alignment_set_type::aligns_sorted_type::const_iterator jiter = aligns.mapped[cept].begin() + 1; jiter != jiter_end; ++ jiter) {
	    score *= dtable_others(*(jiter - 1), *jiter);

#if 0	    
	    std::cerr << "prev: " << *(jiter - 1)
		      << " next: " << *jiter
		      << " = " << dtable_others(*(jiter - 1), *jiter)
		      << std::endl;
#endif
	  }
	  
	  cept_prev = cept;
	}
      
      return score;
    }

    void shrink()
    {
      aligns.shrink();
      
      ttable.clear();
      dtable_head.clear();
      dtable_others.clear();
      ntable.clear();
      ttable_cache_type(ttable).swap(ttable);
      dtable_head_cache_type(dtable_head).swap(dtable_head);
      dtable_others_cache_type(dtable_others).swap(dtable_others);
      ntable_cache_type(ntable).swap(ntable);

      moves.clear();
      swaps.clear();
      move_score_type(moves).swap(moves);
      swap_score_type(swaps).swap(swaps);

      posterior.clear();
      posterior_swap.clear();
      posterior_fertility.clear();
      posterior_type(posterior).swap(posterior);
      posterior_type(posterior_swap).swap(posterior_swap);
      posterior_type(posterior_fertility).swap(posterior_fertility);
    }
    
    // alignment and model score
    logprob_type       logprob;
    alignment_set_type aligns;
    
    // caching...
    sentence_type source_class;
    sentence_type target_class;
    
    ttable_cache_type        ttable;
    dtable_head_cache_type   dtable_head;
    dtable_others_cache_type dtable_others;
    ntable_cache_type        ntable;
    
    // move swap matric
    move_score_type moves;
    swap_score_type swaps;

    // posteriors...
    double total;
    posterior_type posterior;
    posterior_type posterior_swap;
    posterior_type posterior_fertility;
    
    posterior_accum_type neighbour_move;
    posterior_accum_type neighbour_swap;
    posterior_accum_type fertility_inc;
    posterior_accum_type fertility_dec;
  };
  

  typedef Model4Data model4_data_type;
  
  LearnModel4(const LearnBase& __base)
    : LearnBase(__base) {}
  
  void sample(const sentence_type& source,
	      const sentence_type& target,
	      alignment_type& alignment,
	      const ttable_type& ttable,
	      const dtable_type& dtable,
	      const ntable_type& ntable,
	      const ptable_type& ptable,
	      const classes_type& classes_source,
	      const classes_type& classes_target,
	      ttable_type& counts_ttable,
	      dtable_counts_type& counts_dtable,
	      ntable_type& counts_ntable,
	      aligned_type& aligned,
	      double& objective)
  {
    model4.assign(source,
		  target,
		  alignment,
		  model4_type(ttable, dtable, ntable, ptable, classes_source, classes_target));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4.climb())
	break;
    
    // update alignment...
    model4.aligns.alignment(alignment);
    
    // udpate objective
    objective += model4.objective() / target.size();

    // compute posterior
    model4.estimate_posterior(source, target);
    
    // accumulate...
    model4.accumulate(source, target, counts_ttable);
    model4.accumulate(source, target, counts_dtable);
    model4.accumulate(source, target, counts_ntable);
  }
  
  // sampling learning... we will compute a set of alignemnt, and perform summation
  void sample(const sentence_type& source,
	      const sentence_type& target,
	      alignment_type& alignment_source_target,
	      alignment_type& alignment_target_source)
  {
    sample(source,
	   target,
	   alignment_source_target,
	   ttable_source_target,
	   dtable_source_target,
	   ntable_source_target,
	   ptable_source_target,
	   classes_source,
	   classes_target,
	   ttable_counts_source_target,
	   dtable_counts_source_target,
	   ntable_counts_source_target,
	   aligned_source_target,
	   objective_source_target);
    
    sample(target,
	   source,
	   alignment_target_source,
	   ttable_target_source,
	   dtable_target_source,
	   ntable_target_source,
	   ptable_target_source,
	   classes_target,
	   classes_source,
	   ttable_counts_target_source,
	   dtable_counts_target_source,
	   ntable_counts_target_source,
	   aligned_target_source,
	   objective_target_source);
  }
  
  // viterbi, or hard learning, and no ttable counts learning!
  void viterbi(const sentence_type& source,
	       const sentence_type& target,
	       const alignment_type& alignment_source_target,
	       const alignment_type& alignment_target_source)
  {
    model4.assign(source, target, classes_source, classes_target, alignment_source_target);
    
    model4.accumulate(source, target, alignment_source_target, ntable_counts_source_target);
    model4.accumulate(source, target, alignment_source_target, dtable_counts_source_target);
    
    model4.assign(target, source, classes_target, classes_source, alignment_target_source);
    
    model4.accumulate(target, source, alignment_target_source, ntable_counts_target_source);
    model4.accumulate(target, source, alignment_target_source, dtable_counts_target_source);
  }

  void shrink()
  {
    model4.shrink();
    
    LearnBase::shrink();
  }
  
  model4_data_type model4;
};

struct LearnModel4Posterior : public LearnBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  typedef std::vector<double, std::allocator<double> > phi_set_type;
  
  LearnModel4Posterior(const LearnBase& __base)
    : LearnBase(__base) {}
  
  void sample(const sentence_type& source,
	      const sentence_type& target,
	      alignment_type& alignment,
	      const ttable_type& ttable,
	      const dtable_type& dtable,
	      const ntable_type& ntable,
	      const ptable_type& ptable,
	      const classes_type& classes_source,
	      const classes_type& classes_target,
	      ttable_type& counts_ttable,
	      dtable_counts_type& counts_dtable,
	      ntable_type& counts_ntable,
	      aligned_type& aligned,
	      double& objective)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    model4.assign(source,
		  target,
		  alignment,
		  model4_type(ttable, dtable, ntable, ptable, classes_source, classes_target));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4.climb())
	break;
    
    // udpate objective
    objective += model4.objective() / target.size();
    
    phi.clear();
    exp_phi_old.clear();
    
    phi.reserve(source_size + 1);
    exp_phi_old.reserve(source_size + 1);
    
    phi.resize(source_size + 1, 0.0);
    exp_phi_old.resize(source_size + 1, 1.0);

    for (int iter = 0; iter < 5; ++ iter) {
      // compute posterior
      model4.estimate_posterior(source, target);
      
      exp_phi.clear();
      exp_phi.reserve(source_size + 1);
      exp_phi.resize(source_size + 1, 1.0);
      
      size_type count_zero = 0;
      for (size_type src = 1; src <= source_size; ++ src) {
	double sum_posterior = 0.0;
	for (size_type trg = 1; trg <= target_size; ++ trg)
	  sum_posterior += model4.posterior(trg, src);
	
	phi[src] += 1.0 - sum_posterior;
	if (phi[src] > 0.0)
	  phi[src] = 0.0;
	
	count_zero += (phi[src] == 0.0);
	exp_phi[src] = utils::mathop::exp(phi[src]);
      }
      
      if (count_zero == source_size) break;
      
      // rescale emission table...
      for (size_type trg = 1; trg <= target_size; ++ trg)
	for (size_type src = 1; src <= source_size; ++ src)
	  model4.ttable(trg, src) *= exp_phi[src] / exp_phi_old[src];
      
      // swap...
      exp_phi_old.swap(exp_phi);

      model4.update();
      
      // maximum 30 iterations...
      for (int iter = 0; iter != 30; ++ iter)
	if (! model4.climb())
	  break;
    }

    model4.aligns.alignment(alignment);
    
    // accumulate...
    model4.accumulate(source, target, counts_ttable);
    model4.accumulate(source, target, counts_dtable);
    model4.accumulate(source, target, counts_ntable);
  }

  // sampling learning... we will compute a set of alignemnt, and perform summation
  void sample(const sentence_type& source,
	      const sentence_type& target,
	      alignment_type& alignment_source_target,
	      alignment_type& alignment_target_source)
  {
    sample(source,
	   target,
	   alignment_source_target,
	   ttable_source_target,
	   dtable_source_target,
	   ntable_source_target,
	   ptable_source_target,
	   classes_source,
	   classes_target,
	   ttable_counts_source_target,
	   dtable_counts_source_target,
	   ntable_counts_source_target,
	   aligned_source_target,
	   objective_source_target);
    
    sample(target,
	   source,
	   alignment_target_source,
	   ttable_target_source,
	   dtable_target_source,
	   ntable_target_source,
	   ptable_target_source,
	   classes_target,
	   classes_source,
	   ttable_counts_target_source,
	   dtable_counts_target_source,
	   ntable_counts_target_source,
	   aligned_target_source,
	   objective_target_source);
  }

  // viterbi, or hard learning, and no ttable counts learning!
  void viterbi(const sentence_type& source,
	       const sentence_type& target,
	       const alignment_type& alignment_source_target,
	       const alignment_type& alignment_target_source)
  {
    model4.assign(source, target, classes_source, classes_target, alignment_source_target);
    
    model4.accumulate(source, target, alignment_source_target, ntable_counts_source_target);
    model4.accumulate(source, target, alignment_source_target, dtable_counts_source_target);
    
    model4.assign(target, source, classes_target, classes_source, alignment_target_source);
    
    model4.accumulate(target, source, alignment_target_source, ntable_counts_target_source);
    model4.accumulate(target, source, alignment_target_source, dtable_counts_target_source);
  }

  void shrink()
  {
    phi.clear();
    exp_phi.clear();
    exp_phi_old.clear();
    
    phi_set_type(phi).swap(phi);
    phi_set_type(exp_phi).swap(exp_phi);
    phi_set_type(exp_phi_old).swap(exp_phi_old);
    
    model4.shrink();
    
    LearnBase::shrink();
  }
  
  model4_data_type model4;

  phi_set_type phi;
  phi_set_type exp_phi;
  phi_set_type exp_phi_old;
};

struct LearnModel4Symmetric : public LearnBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  LearnModel4Symmetric(const LearnBase& __base)
    : LearnBase(__base) {}
  
  // sampling learning... we will compute a set of alignemnt, and perform summation
  void sample(const sentence_type& source,
	      const sentence_type& target,
	      alignment_type& alignment_source_target,
	      alignment_type& alignment_target_source)
  {
    model4_source_target.assign(source,
				target,
				alignment_source_target,
				model4_type(ttable_source_target,
					    dtable_source_target,
					    ntable_source_target,
					    ptable_source_target,
					    classes_source,
					    classes_target));
    
    model4_target_source.assign(target,
				source,
				alignment_target_source,
				model4_type(ttable_target_source,
					    dtable_target_source,
					    ntable_target_source,
					    ptable_target_source,
					    classes_target,
					    classes_source));
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_source_target.climb())
	break;

    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_target_source.climb())
	break;

    // update alignment...
    model4_source_target.aligns.alignment(alignment_source_target);
    model4_target_source.aligns.alignment(alignment_target_source);

    objective_source_target += model4_source_target.objective() / target.size();
    objective_target_source += model4_target_source.objective() / source.size();
    
    // compute posterior
    model4_source_target.estimate_posterior(source, target);
    model4_target_source.estimate_posterior(target, source);

    // ttable update...
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = 0; trg <= target_size; ++ trg) {
	double count = (! trg ? 1.0 : model4_source_target.posterior(trg, src)) * (! src ? 1.0 : model4_target_source.posterior(src, trg));
	
	if (src && trg)
	  count = utils::mathop::sqrt(count);
	
	const word_type word_source = (! src ? vocab_type::EPSILON : source[src - 1]);
	const word_type word_target = (! trg ? vocab_type::EPSILON : target[trg - 1]);
	
	if (trg)
	  ttable_counts_source_target[word_source][word_target] += count;
	
	if (src)
	  ttable_counts_target_source[word_target][word_source] += count;
      }
    
    // accumulate...
    model4_source_target.accumulate(source, target, dtable_counts_source_target);
    model4_source_target.accumulate(source, target, ntable_counts_source_target);
    
    model4_target_source.accumulate(target, source, dtable_counts_target_source);
    model4_target_source.accumulate(target, source, ntable_counts_target_source);
  }

  // viterbi, or hard learning, and no ttable counts learning!
  void viterbi(const sentence_type& source,
	       const sentence_type& target,
	       const alignment_type& alignment_source_target,
	       const alignment_type& alignment_target_source)
  {
    model4_data_type& model4 = model4_source_target;
    
    model4.assign(source, target, classes_source, classes_target, alignment_source_target);
    
    model4.accumulate(source, target, alignment_source_target, ntable_counts_source_target);
    model4.accumulate(source, target, alignment_source_target, dtable_counts_source_target);
    
    model4.assign(target, source, classes_target, classes_source, alignment_target_source);
    
    model4.accumulate(target, source, alignment_target_source, ntable_counts_target_source);
    model4.accumulate(target, source, alignment_target_source, dtable_counts_target_source);
  }

  void shrink()
  {
    model4_source_target.shrink();
    model4_target_source.shrink();
    
    LearnBase::shrink();
  }
  
  model4_data_type model4_source_target;
  model4_data_type model4_target_source;  
};

struct LearnModel4SymmetricPosterior : public LearnBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  LearnModel4SymmetricPosterior(const LearnBase& __base)
    : LearnBase(__base) {}

  typedef utils::vector2_aligned<double, utils::aligned_allocator<double> > phi_set_type;

  // sampling learning... we will compute a set of alignemnt, and perform summation
  void sample(const sentence_type& source,
	      const sentence_type& target,
	      alignment_type& alignment_source_target,
	      alignment_type& alignment_target_source)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    model4_source_target.assign(source,
				target,
				alignment_source_target,
				model4_type(ttable_source_target,
					    dtable_source_target,
					    ntable_source_target,
					    ptable_source_target,
					    classes_source,
					    classes_target));
    
    model4_target_source.assign(target,
				source,
				alignment_target_source,
				model4_type(ttable_target_source,
					    dtable_target_source,
					    ntable_target_source,
					    ptable_target_source,
					    classes_target,
					    classes_source));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_source_target.climb())
	break;

    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_target_source.climb())
	break;
    
    objective_source_target += model4_source_target.objective() / target.size();
    objective_target_source += model4_target_source.objective() / source.size();

    phi.clear();
    exp_phi.clear();
    exp_phi_old.clear();
    
    phi.reserve(target_size + 1, source_size + 1);
    exp_phi.reserve(target_size + 1, source_size + 1);
    exp_phi_old.reserve(target_size + 1, source_size + 1);
    
    phi.resize(target_size + 1, source_size + 1, 0.0);
    exp_phi.resize(target_size + 1, source_size + 1, 1.0);
    exp_phi_old.resize(target_size + 1, source_size + 1, 1.0);
    
    for (int iter = 0; iter != 5; ++ iter) {
      model4_source_target.estimate_posterior(source, target);
      model4_target_source.estimate_posterior(target, source);
      
      bool updated = false;
      
      // update phi...
      for (size_type src = 1; src <= source_size; ++ src)
	for (size_type trg = 1; trg <= target_size; ++ trg) {
	  const double epsi = model4_source_target.posterior(trg, src) - model4_target_source.posterior(src, trg);
	  const double update = - epsi;
	  
	  phi(trg, src) += update;
	  
	  updated |= (phi(trg, src) != 0.0);
	  exp_phi(trg, src) = utils::mathop::exp(phi(trg, src));
	}
      
      if (! updated) break;
      
      // rescale  ttable...
      for (size_type trg = 1; trg <= target_size; ++ trg)
	for (size_type src = 1; src <= source_size; ++ src)
	  model4_source_target.ttable(trg, src) *= exp_phi(trg, src) / exp_phi_old(trg, src);
      
      for (size_type src = 1; src <= source_size; ++ src)
	for (size_type trg = 1; trg <= target_size; ++ trg)
	  model4_target_source.ttable(src, trg) *= exp_phi_old(trg, src) / exp_phi(trg, src);
      
      // swap...
      exp_phi_old.swap(exp_phi);
      
      // update matrix...
      model4_source_target.update();
      model4_target_source.update();

      // maximum 30 iterations...
      for (int iter = 0; iter != 30; ++ iter)
	if (! model4_source_target.climb())
	  break;
      
      for (int iter = 0; iter != 30; ++ iter)
	if (! model4_target_source.climb())
	  break;
    }
    
    // update alignment...
    model4_source_target.aligns.alignment(alignment_source_target);
    model4_target_source.aligns.alignment(alignment_target_source);

    // accumulate...
    model4_source_target.accumulate(source, target, ttable_counts_source_target);
    model4_source_target.accumulate(source, target, dtable_counts_source_target);
    model4_source_target.accumulate(source, target, ntable_counts_source_target);
    
    model4_target_source.accumulate(target, source, ttable_counts_target_source);
    model4_target_source.accumulate(target, source, dtable_counts_target_source);
    model4_target_source.accumulate(target, source, ntable_counts_target_source);
  }

  // viterbi, or hard learning, and no ttable counts learning!
  void viterbi(const sentence_type& source,
	       const sentence_type& target,
	       const alignment_type& alignment_source_target,
	       const alignment_type& alignment_target_source)
  {
    model4_data_type& model4 = model4_source_target;
    
    model4.assign(source, target, classes_source, classes_target, alignment_source_target);
    
    model4.accumulate(source, target, alignment_source_target, ntable_counts_source_target);
    model4.accumulate(source, target, alignment_source_target, dtable_counts_source_target);
    
    model4.assign(target, source, classes_target, classes_source, alignment_target_source);
    
    model4.accumulate(target, source, alignment_target_source, ntable_counts_target_source);
    model4.accumulate(target, source, alignment_target_source, dtable_counts_target_source);
  }

  void shrink()
  {
    phi.clear();
    exp_phi.clear();
    exp_phi_old.clear();
    
    phi_set_type(phi).swap(phi);
    phi_set_type(exp_phi).swap(exp_phi);
    phi_set_type(exp_phi_old).swap(exp_phi_old);

    model4_source_target.shrink();
    model4_target_source.shrink();
    
    LearnBase::shrink();
  }

  model4_data_type model4_source_target;
  model4_data_type model4_target_source;

  phi_set_type phi;
  phi_set_type exp_phi;
  phi_set_type exp_phi_old;
};

struct ViterbiModel4 : public ViterbiBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  ViterbiModel4(const ttable_type& __ttable_source_target,
		const ttable_type& __ttable_target_source,
		const atable_type& __atable_source_target,
		const atable_type& __atable_target_source,
		const dtable_type& __dtable_source_target,
		const dtable_type& __dtable_target_source,
		const ntable_type& __ntable_source_target,
		const ntable_type& __ntable_target_source,
		const ptable_type& __ptable_source_target,
		const ptable_type& __ptable_target_source,
		const classes_type& __classes_source,
		const classes_type& __classes_target)
  : ViterbiBase(__ttable_source_target, __ttable_target_source,
		__atable_source_target, __atable_target_source,
		__dtable_source_target, __dtable_target_source,
		__ntable_source_target, __ntable_target_source,
		__ptable_source_target, __ptable_target_source,
		__classes_source, __classes_target) {}
  
  
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    model4.assign(source,
		  target,
		  alignment_source_target,
		  model4_type(ttable_source_target,
			      dtable_source_target,
			      ntable_source_target,
			      ptable_source_target,
			      classes_source,
			      classes_target));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4.climb())
	break;
    
    model4.aligns.alignment(alignment_source_target);
    
    model4.assign(target,
		  source,
		  alignment_target_source,
		  model4_type(ttable_target_source,
			      dtable_target_source,
			      ntable_target_source,
			      ptable_target_source,
			      classes_target,
			      classes_source));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4.climb())
	break;
    
    model4.aligns.alignment(alignment_target_source);
    
  }
  
  void shrink()
  {
    model4.shrink();
    
    ViterbiBase::shrink();
  }

  model4_data_type model4;
};

struct PosteriorModel4 : public ViterbiBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  PosteriorModel4(const ttable_type& __ttable_source_target,
		  const ttable_type& __ttable_target_source,
		  const atable_type& __atable_source_target,
		  const atable_type& __atable_target_source,
		  const dtable_type& __dtable_source_target,
		  const dtable_type& __dtable_target_source,
		  const ntable_type& __ntable_source_target,
		  const ntable_type& __ntable_target_source,
		  const ptable_type& __ptable_source_target,
		  const ptable_type& __ptable_target_source,
		  const classes_type& __classes_source,
		  const classes_type& __classes_target)
  : ViterbiBase(__ttable_source_target, __ttable_target_source,
		__atable_source_target, __atable_target_source,
		__dtable_source_target, __dtable_target_source,
		__ntable_source_target, __ntable_target_source,
		__ptable_source_target, __ptable_target_source,
		__classes_source, __classes_target) {}

  template <typename Matrix>
  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source,
		  Matrix& posterior_source_target,
		  Matrix& posterior_target_source)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    model4.assign(source,
		  target,
		  alignment_source_target,
		  model4_type(ttable_source_target,
			      dtable_source_target,
			      ntable_source_target,
			      ptable_source_target,
			      classes_source,
			      classes_target));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4.climb())
	break;
    
    for (size_type trg = 0; trg <= target_size; ++ trg)
      std::copy(model4.posterior.begin(trg), model4.posterior.end(trg), posterior_source_target.begin(trg));
    
    model4.assign(target,
		  source,
		  alignment_target_source,
		  model4_type(ttable_target_source,
			      dtable_target_source,
			      ntable_target_source,
			      ptable_target_source,
			      classes_target,
			      classes_source));
    
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4.climb())
	break;
    
    for (size_type src = 0; src <= source_size; ++ src)
      std::copy(model4.posterior.begin(src), model4.posterior.end(src), posterior_target_source.begin(src));
  } 
  
  void shrink()
  {
    model4.shrink();
    
    ViterbiBase::shrink();
  }

  model4_data_type model4;
};

struct ITGModel4 : public ViterbiBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  typedef utils::vector2<double, std::allocator<double> > matrix_type;

  ITGModel4(const ttable_type& __ttable_source_target,
	    const ttable_type& __ttable_target_source,
	    const atable_type& __atable_source_target,
	    const atable_type& __atable_target_source,
	    const dtable_type& __dtable_source_target,
	    const dtable_type& __dtable_target_source,
	    const ntable_type& __ntable_source_target,
	    const ntable_type& __ntable_target_source,
	    const ptable_type& __ptable_source_target,
	    const ptable_type& __ptable_target_source,
	    const classes_type& __classes_source,
	    const classes_type& __classes_target)
  : ViterbiBase(__ttable_source_target, __ttable_target_source,
		__atable_source_target, __atable_target_source,
		__dtable_source_target, __dtable_target_source,
		__ntable_source_target, __ntable_target_source,
		__ptable_source_target, __ptable_target_source,
		__classes_source, __classes_target) {}


  class insert_align
  {
    alignment_type& alignment_source_target;
    alignment_type& alignment_target_source;
    
  public:
    insert_align(alignment_type& __alignment_source_target,
		 alignment_type& __alignment_target_source)
      : alignment_source_target(__alignment_source_target),
	alignment_target_source(__alignment_target_source) {}
    
    template <typename Edge>
    insert_align& operator=(const Edge& edge)
    {	
      alignment_source_target.push_back(edge);
      alignment_target_source.push_back(std::make_pair(edge.second, edge.first));
      
      return *this;
    }
    
    insert_align& operator*() { return *this; }
    insert_align& operator++() { return *this; }
    insert_align operator++(int) { return *this; }
  };

  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    model4_source_target.assign(source,
				target,
				alignment_source_target,
				model4_type(ttable_source_target,
					    dtable_source_target,
					    ntable_source_target,
					    ptable_source_target,
					    classes_source,
					    classes_target));
    
    model4_target_source.assign(target,
				source,
				alignment_target_source,
				model4_type(ttable_target_source,
					    dtable_target_source,
					    ntable_target_source,
					    ptable_target_source,
					    classes_target,
					    classes_source));
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_source_target.climb())
	break;

    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_target_source.climb())
	break;
    
    // compute posterior
    model4_source_target.estimate_posterior(source, target);
    model4_target_source.estimate_posterior(target, source);
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    costs.clear();
    costs.reserve(source_size + 1, target_size + 1);
    costs.resize(source_size + 1, target_size + 1, boost::numeric::bounds<double>::lowest());
    
    for (size_type src = 1; src <= source_size; ++ src)
      for (size_type trg = 1; trg <= target_size; ++ trg)
	costs(src, trg) = 0.5 * (utils::mathop::log(model4_source_target.posterior(trg, src)) 
				 + utils::mathop::log(model4_target_source.posterior(src, trg)));
    
    for (size_type trg = 1; trg <= target_size; ++ trg)
      costs(0, trg) = utils::mathop::log(model4_source_target.posterior(trg, 0));
    
    for (size_type src = 1; src <= source_size; ++ src)
      costs(src, 0) = utils::mathop::log(model4_target_source.posterior(src, 0));

    alignment_source_target.clear();
    alignment_target_source.clear();
    
    aligner(costs, insert_align(alignment_source_target, alignment_target_source));
    
    std::sort(alignment_source_target.begin(), alignment_source_target.end());
    std::sort(alignment_target_source.begin(), alignment_target_source.end());
  }

  void shrink()
  {
    costs.clear();
    matrix_type(costs).swap(costs);
    
    model4_source_target.shrink();
    model4_target_source.shrink();

    aligner.shrink();
    
    ViterbiBase::shrink();
  }


  matrix_type costs;
  
  model4_data_type model4_source_target;
  model4_data_type model4_target_source;
  
  detail::ITGAlignment aligner;
};

struct MaxMatchModel4 : public ViterbiBase
{
  typedef LearnModel4::model4_type      model4_type;
  typedef LearnModel4::model4_data_type model4_data_type;

  typedef model4_data_type::size_type       size_type;
  typedef model4_data_type::difference_type difference_type;

  typedef utils::vector2<double, std::allocator<double> > matrix_type;

  MaxMatchModel4(const ttable_type& __ttable_source_target,
		 const ttable_type& __ttable_target_source,
		 const atable_type& __atable_source_target,
		 const atable_type& __atable_target_source,
		 const dtable_type& __dtable_source_target,
		 const dtable_type& __dtable_target_source,
		 const ntable_type& __ntable_source_target,
		 const ntable_type& __ntable_target_source,
		 const ptable_type& __ptable_source_target,
		 const ptable_type& __ptable_target_source,
		 const classes_type& __classes_source,
		 const classes_type& __classes_target)
  : ViterbiBase(__ttable_source_target, __ttable_target_source,
		__atable_source_target, __atable_target_source,
		__dtable_source_target, __dtable_target_source,
		__ntable_source_target, __ntable_target_source,
		__ptable_source_target, __ptable_target_source,
		__classes_source, __classes_target) {}
  
  class insert_align
  {
    int source_size;
    int target_size;
    
    alignment_type& alignment_source_target;
    alignment_type& alignment_target_source;
    
  public:
    insert_align(const int& _source_size,
		 const int& _target_size,
		 alignment_type& __alignment_source_target,
		 alignment_type& __alignment_target_source)
      : source_size(_source_size), target_size(_target_size),
	alignment_source_target(__alignment_source_target),
	alignment_target_source(__alignment_target_source) {}
    
    template <typename Edge>
    insert_align& operator=(const Edge& edge)
    {	
      if (edge.first < source_size && edge.second < target_size) {
	alignment_source_target.push_back(edge);
	alignment_target_source.push_back(std::make_pair(edge.second, edge.first));
      }
      
      return *this;
    }
    
    insert_align& operator*() { return *this; }
    insert_align& operator++() { return *this; }
    insert_align operator++(int) { return *this; }
  };

  void operator()(const sentence_type& source,
		  const sentence_type& target,
		  alignment_type& alignment_source_target,
		  alignment_type& alignment_target_source)
  {
    model4_source_target.assign(source,
				target,
				alignment_source_target,
				model4_type(ttable_source_target,
					    dtable_source_target,
					    ntable_source_target,
					    ptable_source_target,
					    classes_source,
					    classes_target));
    
    model4_target_source.assign(target,
				source,
				alignment_target_source,
				model4_type(ttable_target_source,
					    dtable_target_source,
					    ntable_target_source,
					    ptable_target_source,
					    classes_target,
					    classes_source));
    // maximum 30 iterations...
    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_source_target.climb())
	break;

    for (int iter = 0; iter != 30; ++ iter)
      if (! model4_target_source.climb())
	break;
    
    // compute posterior
    model4_source_target.estimate_posterior(source, target);
    model4_target_source.estimate_posterior(target, source);
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    costs.clear();
    costs.reserve(source_size + target_size, target_size + source_size);
    costs.resize(source_size + target_size, target_size + source_size, 0.0);
    
    for (size_type src = 0; src != source_size; ++ src)
      for (size_type trg = 0; trg != target_size; ++ trg) {
	costs(src, trg) = 0.5 * (utils::mathop::log(model4_source_target.posterior(trg + 1, src + 1))
				 + utils::mathop::log(model4_target_source.posterior(src + 1, trg + 1)));
	
	costs(src, trg + source_size) = utils::mathop::log(model4_target_source.posterior(src + 1, 0));
	costs(src + target_size, trg) = utils::mathop::log(model4_source_target.posterior(trg + 1, 0));
      }

    alignment_source_target.clear();
    alignment_target_source.clear();
    
    kuhn_munkres_assignment(costs, insert_align(source_size, target_size, alignment_source_target, alignment_target_source));
    
    std::sort(alignment_source_target.begin(), alignment_source_target.end());
    std::sort(alignment_target_source.begin(), alignment_target_source.end());
  }
  
  void shrink()
  {
    costs.clear();
    matrix_type(costs).swap(costs);    
    
    model4_source_target.shrink();
    model4_target_source.shrink();

    ViterbiBase::shrink();
  }
  
  matrix_type costs;

  model4_data_type model4_source_target;
  model4_data_type model4_target_source;
};


#endif