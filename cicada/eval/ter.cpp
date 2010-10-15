 
#include <algorithm>
#include <set>
#include <iterator>

#include "ter.hpp"

#include <google/dense_hash_set>

#include <boost/functional/hash.hpp>

#include <utils/vector2.hpp>
#include <utils/sgi_hash_map.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  namespace eval
  {

    struct TERScorerConstant
    {
      struct COSTS
      {
	static const double insertion;
	static const double deletion;
	static const double substitution;
	static const double shift;
      };
      
      struct TRANSITION
      {
	enum transition_type {
	  match,
	  substitution,
	  insertion,
	  deletion,
	};
      };
      
      static const int max_shift_size = 10;
      static const int max_shift_dist = 50;
    };
    
    const double TERScorerConstant::COSTS::insertion    = 1.0;
    const double TERScorerConstant::COSTS::deletion     = 1.0;
    const double TERScorerConstant::COSTS::substitution = 1.0;
    const double TERScorerConstant::COSTS::shift        = 1.0;
    
    // Do we really implement this...???
    class TERScorerImpl : public TERScorerConstant
    {
    private:
      friend class TERScorer;
      
    public:
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;
      
      typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;

      struct Score
      {
	double score;
	int match;
	int insertion;
	int deletion;
	int substitution;
	int shift;
	
	Score() : score(), match(0), insertion(0), deletion(0), substitution(0), shift(0) {}
      };

      struct Shift
      {
	Shift() : begin(0), end(0), moveto(0), reloc(0) {}
	Shift(const int __begin, const int __end, const int __moveto, const int __reloc)
	  : begin(__begin), end(__end), moveto(__moveto), reloc(__reloc) {}
	
	int begin;
	int end;
	int moveto;
	int reloc;
      };

      typedef TRANSITION::transition_type transition_type;
      
      typedef Score value_type;
      typedef Shift shift_type;

      typedef std::vector<transition_type, std::allocator<transition_type> > path_type;

      typedef std::vector<shift_type, std::allocator<shift_type> > shift_set_type;
      typedef std::vector<shift_set_type, std::allocator<shift_set_type> > shift_matrix_type;

      typedef std::vector<bool, std::allocator<bool> > error_set_type;
      typedef std::vector<int, std::allocator<int> > align_set_type;

      typedef sentence_type ngram_type;
      typedef std::set<int, std::less<int>, std::allocator<int> > index_set_type;
      
#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<ngram_type, index_set_type, boost::hash<ngram_type>, std::equal_to<ngram_type>,
				      std::allocator<std::pair<const ngram_type, index_set_type> > > ngram_index_map_type;
#else
      typedef sgi::hash_map<ngram_type, index_set_type, boost::hash<ngram_type>, std::equal_to<ngram_type>,
			    std::allocator<std::pair<const ngram_type, index_set_type> > > ngram_index_map_type;
#endif
      
      
      TERScorerImpl(const sentence_type& __ref)
	: ref(__ref) {  }
      
      
      value_type operator()(const sentence_type& sentence) const
      {
	value_type value;
	
	calculate_shifts(sentence, ref, value);
	
	return value;
      }

    private:

      double calculate_shifts(const sentence_type& hyp_orig, const sentence_type& ref, value_type& value) const
      {
	value = value_type();
	
	sentence_type hyp = hyp_orig;
	path_type     path;
	double        cost = minimum_edit_distance(hyp, ref, path);

	//std::cerr << "initial cost: " << cost << std::endl;
	
	sentence_type hyp_new;
	path_type     path_new;
	double        cost_new;
	
	ngram_index_map_type ngram_index;
	build_ngram_matches(hyp, ngram_index);
	
	while (1) {
	  hyp_new.clear();
	  path_new.clear();
	  cost_new = 0;

	  if (! calculate_best_shift(hyp, path, cost, hyp_new, path_new, cost_new, ngram_index))
	    break;
	  
	  value.score += COSTS::shift;
	  ++ value.shift;
	  
	  hyp.swap(hyp_new);
	  path.swap(path_new);
	  cost = cost_new;

	  //std::cerr << "new cost: " << cost << " shift: " << value.shift << std::endl;
	}
	
	path_type::const_iterator piter_end = path.end();
	for (path_type::const_iterator piter = path.begin(); piter != piter_end; ++ piter) {
	  switch (*piter) {
	  case TRANSITION::match:        ++ value.match; break;
	  case TRANSITION::substitution: ++ value.substitution; break;
	  case TRANSITION::insertion:    ++ value.insertion; break;
	  case TRANSITION::deletion:     ++ value.deletion; break;
	  }
	}
	
	value.score += cost;

#if 0
	std::cerr << "score: " << value.score
		  << " match: " << value.match
		  << " substitution: " << value.substitution
		  << " insertion: " << value.insertion
		  << " deletion: " << value.deletion
		  << " shift: " << value.shift
		  << std::endl;
#endif
	
	return value.score;
      }

      void build_ngram_matches(const sentence_type& hyp, ngram_index_map_type& ngram_index) const
      {
	ngram_index.clear();
	
	word_set_type words_intersect;
	words_intersect.set_empty_key(word_type());
	words_intersect.insert(hyp.begin(), hyp.end());

	word_set_type::const_iterator iiter_end = words_intersect.end();
	
	ngram_type ngram;
	
	for (int start = 0; start != ref.size(); ++ start) {
	  ngram.clear();
	  const int max_length = utils::bithack::min(max_shift_size, static_cast<int>(ref.size() - start));
	  for (int length = 0; length != max_length && words_intersect.find(ref[start + length]) != iiter_end; ++ length) {
	    ngram.push_back(ref[start + length]);
	    
	    ngram_index[ngram].insert(start);
	  }
	}
      }
      
      void find_alignment_error(const path_type& path,
				error_set_type& herr,
				error_set_type& rerr,
				align_set_type& ralign) const
      {
	int hpos = -1;
	int rpos = -1;
	
	//std::cerr << "edit distance: ";
	for (int i = 0; i < path.size(); ++ i) {
	  switch (path[i]) {
	  case TRANSITION::match:
	    //std::cerr << " M";
	    ++ hpos;
	    ++ rpos;
	    herr.push_back(false);
	    rerr.push_back(false);
	    ralign.push_back(hpos);
	    break;
	  case TRANSITION::substitution:
	    //std::cerr << " S";
	    ++ hpos;
	    ++ rpos;
	    herr.push_back(true);
	    rerr.push_back(true);
	    ralign.push_back(hpos);
	    break;
	  case TRANSITION::insertion:
	    //std::cerr << " I";
	    ++ hpos;
	    herr.push_back(true);
	    break;
	  case TRANSITION::deletion:
	    //std::cerr << " D";
	    ++ rpos;
	    rerr.push_back(true);
	    ralign.push_back(hpos);
	    break;
	  }
	}
	//std::cerr << std::endl;
      }
      
      bool calculate_best_shift(const sentence_type& hyp,
				const path_type& path,
				const double cost,
				sentence_type& hyp_best,
				path_type& path_best,
				double& cost_best,
				const ngram_index_map_type& ngram_index) const
      {
	//std::cerr << "hyp: " << hyp << std::endl;
	//std::cerr << "ref: " << ref << std::endl;

	error_set_type herr;
	error_set_type rerr;
	align_set_type ralign;
	
	herr.reserve(path.size());
	rerr.reserve(path.size());
	ralign.reserve(path.size());
	
	find_alignment_error(path, herr, rerr, ralign);
	
#if 0
	std::cerr << "ref-align: ";
	std::copy(ralign.begin(), ralign.end(), std::ostream_iterator<int>(std::cerr, " "));
	std::cerr << std::endl;
	
	std::cerr << "ref-err: ";
	std::copy(rerr.begin(), rerr.end(), std::ostream_iterator<bool>(std::cerr, " "));
	std::cerr << std::endl;

	std::cerr << "hyp-err: ";
	std::copy(herr.begin(), herr.end(), std::ostream_iterator<bool>(std::cerr, " "));
	std::cerr << std::endl;
#endif
	
	shift_matrix_type shifts(max_shift_size + 1);
	
	gather_all_possible_shifts(hyp, ralign, herr, rerr, ngram_index, shifts);
	
	double cost_shift_best = 0;
	cost_best = cost;
	bool found = false;
	
	// enumerate from max-shifts
	sentence_type hyp_shifted(hyp.size());
	path_type     path_shifted;
	
	for (int i = shifts.size() - 1; i >= 0; -- i) 
	  if (! shifts[i].empty()) {
	    const double curfix = cost - (cost_shift_best + cost_best);
	    const double maxfix = 2.0 * (i + 1);
	    
	    if (curfix > maxfix || (cost_shift_best != 0 && curfix == maxfix)) break;
	  
	    for (int j = 0; j < shifts[i].size(); ++ j) {
	      const shift_type& shift = shifts[i][j];
	      
	      const double curfix = cost - (cost_shift_best + cost_best);
	      
	      if (curfix > maxfix || (cost_shift_best != 0 && curfix == maxfix)) break;
	    
	      //std::cerr << "candidate shift: [" << shift.begin << ", " << shift.end << "]: " << shift.reloc << std::endl;
	      
	      perform_shift(hyp, shift, hyp_shifted);
	    
	      const double cost_shifted = minimum_edit_distance(hyp_shifted, ref, path_shifted);
	      const double gain = (cost_best + cost_shift_best) - (cost_shifted + COSTS::shift);
	      
	      //std::cerr << "hyp original: " << hyp << std::endl;
	      //std::cerr << "hyp shifted:  " << hyp_shifted << std::endl;
	      //std::cerr << "cost shifted: " << cost_shifted << " gain: " << gain << std::endl;
	    
	      if (gain > 0 || (cost_shift_best == 0 && gain == 0)) {
		cost_best       = cost_shifted;
		cost_shift_best = COSTS::shift;
		
		path_best.swap(path_shifted);
		hyp_best.swap(hyp_shifted);
		found = true;
		
		//std::cerr << "better shift: [" << shift.begin << ", " << shift.end << "]: " << shift.reloc << std::endl;
	      }
	    }
	  }

	return found;
      }

      void perform_shift(const sentence_type& sentence,
			 const shift_type& shift,
			 sentence_type& shifted) const
      {
	shifted.clear();
	
	std::back_insert_iterator<sentence_type> oiter(shifted);
	
	sentence_type::const_iterator siter_begin = sentence.begin();
	sentence_type::const_iterator siter_end = sentence.end();

	if (shift.reloc == -1) {
	  std::copy(siter_begin + shift.begin, siter_begin + shift.end + 1, oiter);
	  std::copy(siter_begin, siter_begin + shift.begin, oiter);
	  std::copy(siter_begin + shift.end + 1, siter_end, oiter);
	} else if (shift.reloc < shift.begin) {
	  std::copy(siter_begin, siter_begin + shift.reloc + 1, oiter);
	  std::copy(siter_begin + shift.begin, siter_begin + shift.end + 1, oiter);
	  std::copy(siter_begin + shift.reloc + 1, siter_begin + shift.begin, oiter);
	  std::copy(siter_begin + shift.end + 1, siter_end, oiter);
	} else if (shift.end < shift.reloc) {
	  std::copy(siter_begin, siter_begin + shift.begin, oiter);
	  std::copy(siter_begin + shift.end + 1, siter_begin + shift.reloc + 1, oiter);
	  std::copy(siter_begin + shift.begin, siter_begin + shift.end + 1, oiter);
	  std::copy(siter_begin + shift.reloc + 1, siter_end, oiter);
	} else {
	  std::copy(siter_begin, siter_begin + shift.begin, oiter);
	  std::copy(siter_begin + shift.end + 1, std::min(siter_end, siter_begin + shift.end + shift.reloc - shift.begin + 1), oiter);
	  std::copy(siter_begin + shift.begin, siter_begin + shift.end + 1, oiter);
	  std::copy(siter_begin + shift.end + shift.reloc - shift.begin + 1, siter_end, oiter);
	}

	if (sentence.size() != shifted.size())
	  throw std::runtime_error(std::string("size do not match:")
				   + " original: " + boost::lexical_cast<std::string>(sentence.size())
				   + " shifted: "  + boost::lexical_cast<std::string>(shifted.size()));
      }

      void gather_all_possible_shifts(const sentence_type& hyp,
				      const align_set_type& ralign,
				      const error_set_type& herr,
				      const error_set_type& rerr,
				      const ngram_index_map_type& ngram_index,
				      shift_matrix_type& shifts) const
      {
	ngram_type ngram;
	
	for (int start = 0; start != hyp.size(); ++ start) {
	  ngram.clear();
	  ngram.push_back(hyp[start]);
	  
	  ngram_index_map_type::const_iterator niter = ngram_index.find(ngram);
	  if (niter == ngram_index.end()) continue;
	  
	  bool found = false;
	  index_set_type::const_iterator iiter_end = niter->second.end();
	  for (index_set_type::const_iterator iiter = niter->second.begin(); iiter != iiter_end && ! found; ++ iiter) {
	    const int moveto = *iiter;
	    found = (start != ralign[moveto] && (ralign[moveto] - start <= max_shift_dist) && (start - ralign[moveto] - 1 <= max_shift_dist));
	  }
	  
	  if (! found) continue;
	  
	  ngram.clear();
	  const int last = utils::bithack::min(start + max_shift_size, static_cast<int>(hyp.size()));
	  for (int end = start; found && end != last; ++ end) {
	    ngram.push_back(hyp[end]);

	    //std::cerr << "range: [" << start << ", " << end << "]" << std::endl;
	    
	    found = false;
	    
	    ngram_index_map_type::const_iterator niter = ngram_index.find(ngram);
	    if (niter == ngram_index.end()) break;

	    //std::cerr << "found ngram: " << niter->first << std::endl;
	    
	    error_set_type::const_iterator hiter_begin = herr.begin() + start;
	    error_set_type::const_iterator hiter_end   = herr.begin() + end + 1;
	    if (std::find(hiter_begin, hiter_end, true) == hiter_end) {
	      found = true;
	      continue;
	    }
	    
	    index_set_type::const_iterator iiter_end = niter->second.end();
	    for (index_set_type::const_iterator iiter = niter->second.begin(); iiter != iiter_end; ++ iiter) {
	      const int moveto = *iiter;
	      
	      if (ralign[moveto] != start
		  && (ralign[moveto] < start || end < ralign[moveto])
		  && ralign[moveto] - start <= max_shift_dist
		  && start - ralign[moveto] - 1 <= max_shift_dist) {
		
		found = true;
		
		error_set_type::const_iterator riter_begin = rerr.begin() + moveto;
		error_set_type::const_iterator riter_end   = rerr.begin() + end - start + moveto + 1;
		
		if (std::find(riter_begin, riter_end, true) == riter_end) continue;
		
		shift_set_type& sshifts = shifts[end - start];
		
		for (int roff = -1; roff <= end - start; ++ roff) {
		  
		  if (roff == -1 && moveto == 0)
		    sshifts.push_back(shift_type(start, end, -1, -1));
		  else if (start != ralign[moveto + roff] && (roff == 0 || ralign[moveto + roff] != ralign[moveto]))
		    sshifts.push_back(shift_type(start, end, moveto + roff, ralign[moveto + roff]));
		}
	      }
	    }
	  }
	}
      }
      
      double minimum_edit_distance(const sentence_type& hyp, const sentence_type& ref, path_type& path) const
      {
	typedef utils::vector2<transition_type, std::allocator<transition_type> > matrix_transition_type;
	typedef utils::vector2<double, std::allocator<double> > matrix_cost_type;

	matrix_transition_type trans(hyp.size() + 1, ref.size() + 1, TRANSITION::match);
	matrix_cost_type       costs(hyp.size() + 1, ref.size() + 1, 0.0);
	
	for (int i = 0; i <= hyp.size(); ++ i)
	  costs(i, 0) = i * COSTS::insertion;
	for (int j = 0; j <= ref.size(); ++ j)
	  costs(0, j) = j * COSTS::deletion;
	
	for (int i = 1; i <= hyp.size(); ++ i)
	  for (int j = 1; j <= ref.size(); ++ j) {
	    double&          cur_cost = costs(i, j);
	    transition_type& cur_tran = trans(i, j);
	    
	    if (hyp[i - 1] == ref[j - 1]) {
	      cur_cost = costs(i - 1, j - 1);
	      cur_tran = TRANSITION::match;
	    } else {
	      cur_cost = costs(i - 1, j - 1) + COSTS::substitution;
	      cur_tran = TRANSITION::substitution;
	    }
	    
	    const double ins = costs(i - 1, j) + COSTS::insertion;
	    if (cur_cost > ins) {
	      cur_cost = ins;
	      cur_tran = TRANSITION::insertion;
	    }
	    const double del = costs(i, j - 1) + COSTS::deletion;
	    if (cur_cost > del) {
	      cur_cost = del;
	      cur_tran = TRANSITION::deletion;
	    }
	  }

	path.clear();
	int i = hyp.size();
	int j = ref.size();
	while (i > 0 || j > 0) {
	  if (j == 0) {
	    -- i;
	    path.push_back(TRANSITION::insertion);
	  } else if (i == 0) {
	    -- j;
	    path.push_back(TRANSITION::deletion);
	  } else {
	    const transition_type& t = trans(i, j);
	    path.push_back(t);
	    switch (t) {
	    case TRANSITION::substitution:
	    case TRANSITION::match:        -- i; -- j; break;
	    case TRANSITION::insertion:    -- i; break;
	    case TRANSITION::deletion:     -- j; break;
	    }
	  }
	}
	
	std::reverse(path.begin(), path.end());
	
	return costs(hyp.size(), ref.size());
      }
      
    private:
      sentence_type ref;
    };
   
    TERScorer::TERScorer(const TERScorer& x)
      : Scorer(static_cast<const Scorer&>(*this))
    {
      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
    }
    
    TERScorer::~TERScorer()
    {
      clear();
    }
    
    TERScorer& TERScorer::operator=(const TERScorer& x)
    {
      clear();
      
      static_cast<Scorer&>(*this) = static_cast<const Scorer&>(x);

      for (impl_set_type::const_iterator iter = x.impl.begin(); iter != x.impl.end(); ++ iter)
	impl.push_back(new impl_type(*(*iter)));
      
      return *this;
    }
    
    void TERScorer::clear()
    {
      for (impl_set_type::iterator iter = impl.begin(); iter != impl.end(); ++ iter)
	delete *iter;
      
      impl.clear();
    }
    
    void TERScorer::insert(const sentence_type& __sentence)
    {
      sentence_type sentence_split;
      sentence_type sentence_lower;
      
      if (split)
	split_non_ascii_characters(__sentence, sentence_split);
      const sentence_type& __sentence_split = (split ? sentence_split : __sentence);
      
      if (lower)
	lower_case(__sentence_split, sentence_lower);
      const sentence_type& sentence = (lower ? sentence_lower : __sentence_split);

      impl.push_back(new impl_type(sentence));
    }
    
    TERScorer::score_ptr_type TERScorer::score(const sentence_type& __sentence) const
    {
      sentence_type sentence_split;
      sentence_type sentence_lower;
      
      if (split)
	split_non_ascii_characters(__sentence, sentence_split);
      const sentence_type& __sentence_split = (split ? sentence_split : __sentence);
      
      if (lower)
	lower_case(__sentence_split, sentence_lower);
      const sentence_type& sentence = (lower ? sentence_lower : __sentence_split);
      
      double score_best = std::numeric_limits<double>::infinity();

      std::auto_ptr<TER> ter(new TER());
      
      for (impl_set_type::const_iterator iter = impl.begin(); iter != impl.end(); ++ iter) {
	impl_type& evaluator = const_cast<impl_type&>(*(*iter));
	
	const impl_type::value_type value = evaluator(sentence);
	
	if (value.score < score_best) {
	  score_best = value.score;
	  
	  ter->insertion    = value.insertion;
	  ter->deletion     = value.deletion;
	  ter->substitution = value.substitution;
	  ter->shift        = value.shift;
	}

	ter->references += evaluator.ref.size();
      }
      
      if (! impl.empty())
	ter->references /= impl.size();
      
#if 0
      std::cerr << "final: references: " << ter->references
		<< " insertion: " << ter->insertion
		<< " deletion: " << ter->deletion
		<< " substitution: " << ter->substitution
		<< " shift: " << ter->shift
		  << std::endl;
#endif

      return score_ptr_type(ter.release());
    }
  };
};

