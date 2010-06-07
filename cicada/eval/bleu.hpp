// -*- mode: c++ -*-

#ifndef __CICADA__EVAL__BLEU__HPP__
#define __CICADA__EVAL__BLEU__HPP__ 1

#include <cmath>
#include <cfloat>

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <utility>
#include <map>

#include <cicada/eval/score.hpp>

#include <utils/compact_trie.hpp>
#include <utils/bithack.hpp>

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace eval
  {
    class Bleu : public Score
    {
    public:
      typedef double count_type;
      
      typedef std::vector<count_type, std::allocator<count_type> > ngram_counts_type;

    public:      
      Bleu(const int order)
	: ngrams_reference(order, 0),  ngrams_hypothesis(order, 0),
	  length_reference(0), length_hypothesis(0) {}
      
      std::pair<double, double> score() const
      {
	const double factor = 1.0 / ngrams_hypothesis.size();
	const double penalty = std::min(1.0 - length_reference / length_hypothesis, 0.0);
	
	double smooth = 0.5;
	double score = 0.0;
	
	for (int n = 0; n < ngrams_hypothesis.size(); ++ n) {
	  const double p = (ngrams_reference[n] > 0
			    ? (ngrams_hypothesis[n] > 0 ? ngrams_hypothesis[n] : smooth) / ngrams_reference[n]
			    : 0.0);
	  
	  score += p > 0.0 ? std::log(p) : 0.0;
	  smooth *= 0.5;
	}
	
	score /= ngrams_hypothesis.size();
	score += penalty;
	
	return std::make_pair(std::exp(score), std::exp(penalty));
      }

      void assign(const score_type& score)
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	ngrams_reference  = rhs->ngrams_reference;
	ngrams_hypothesis = rhs->ngrams_hypothesis;
	
	length_reference  = rhs->length_reference;
	length_hypothesis = rhs->length_hypothesis;
      }
      
      void plus_equal(const score_type& score)
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	std::transform(ngrams_reference.begin(),  ngrams_reference.end(),  rhs->ngrams_reference.begin(),  ngrams_reference.begin(),  std::plus<count_type>());
	std::transform(ngrams_hypothesis.begin(), ngrams_hypothesis.end(), rhs->ngrams_hypothesis.begin(), ngrams_hypothesis.begin(), std::plus<count_type>());
	
	length_reference  += rhs->length_reference;
	length_hypothesis += rhs->length_hypothesis;
      }
      
      void minus_equal(const score_type& score)
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	std::transform(ngrams_reference.begin(),  ngrams_reference.end(),  rhs->ngrams_reference.begin(),  ngrams_reference.begin(),  std::minus<count_type>());
	std::transform(ngrams_hypothesis.begin(), ngrams_hypothesis.end(), rhs->ngrams_hypothesis.begin(), ngrams_hypothesis.begin(), std::minus<count_type>());
	
	length_reference  -= rhs->length_reference;
	length_hypothesis -= rhs->length_hypothesis;
      }
      
      void multiplies_equal(const double& scale)
      {
	std::transform(ngrams_reference.begin(),  ngrams_reference.end(),  ngrams_reference.begin(),  std::bind2nd(std::multiplies<count_type>(), scale));
	std::transform(ngrams_hypothesis.begin(), ngrams_hypothesis.end(), ngrams_hypothesis.begin(), std::bind2nd(std::multiplies<count_type>(), scale));
	
	length_reference  *= scale;
	length_hypothesis *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	std::transform(ngrams_reference.begin(),  ngrams_reference.end(),  ngrams_reference.begin(),  std::bind2nd(std::divides<count_type>(), scale));
	std::transform(ngrams_hypothesis.begin(), ngrams_hypothesis.end(), ngrams_hypothesis.begin(), std::bind2nd(std::divides<count_type>(), scale));
	
	length_reference  /= scale;
	length_hypothesis /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new Bleu(ngrams_hypothesis.size()));
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new Bleu(*this));
      }


      ngram_counts_type ngrams_reference;
      ngram_counts_type ngrams_hypothesis;
      count_type        length_reference;
      count_type        length_hypothesis;
    };
    

    class BleuScorer : public Scorer
    {
    private:
      typedef double count_type;
      
      typedef std::allocator<std::pair<const word_type, count_type> > ngram_allocator_type;
      typedef utils::compact_trie<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>, ngram_allocator_type> ngram_set_type;
      
      typedef std::vector<int, std::allocator<int> > size_set_type;

    public:
      BleuScorer(int __order = 4) : order(__order) {}

      bool error_metric() const { return false; }

      void clear()
      {
	ngrams.clear();
	sizes.clear();
	order = 0;
      }
      
      void insert(const sentence_type& sentence)
      {
	typedef ngram_set_type::id_type id_type;
	typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
	counts_type counts;
	
	
	sentence_type::const_iterator siter_end = sentence.end();
	for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	  ngram_set_type::id_type id = ngrams.root();
	  
	  for (sentence_type::const_iterator iter = siter; iter != std::min(siter + order, siter_end); ++ iter) {
	    id = ngrams.insert(id, *iter);
	    ++ counts[id];
	  }
	}
	
	counts_type::const_iterator citer_end = counts.end();
	for (counts_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer)
	  ngrams[citer->first] = std::max(ngrams[citer->first], citer->second);
	
	sizes.push_back(sentence.size());
	
	std::sort(sizes.begin(), sizes.end());
      }

      score_ptr_type score(const sentence_type& sentence) const
      {
	typedef ngram_set_type::id_type id_type;
	typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
	typedef std::vector<counts_type, std::allocator<counts_type> > counts_set_type;
	
	
	std::auto_ptr<Bleu> bleu(new Bleu(order));
	counts_set_type counts(order);
	
	const int hypothesis_size = sentence.size();
	
	int reference_size = 0;
	int min_diff = boost::numeric::bounds<int>::highest();
	
	for (size_set_type::const_iterator siter = sizes.begin(); siter != sizes.end(); ++ siter) {
	  const int diff = utils::bithack::abs(*siter - hypothesis_size);
	  
	  if (diff < min_diff) {
	    min_diff = diff;
	    reference_size = *siter;
	  } else if (diff == min_diff)
	    reference_size = utils::bithack::min(reference_size, *siter);
	}
	
	bleu->length_hypothesis += hypothesis_size;
	bleu->length_reference  += reference_size;
	
	// collect total counts...
	for (int n = 0; n < utils::bithack::min(order, hypothesis_size); ++ n)
	  bleu->ngrams_reference[n] += hypothesis_size - n;
	
	// collect ngrams matched with references
	sentence_type::const_iterator siter_end = sentence.end();
	for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
	
	  int n = 0;
	  id_type id = ngrams.root();
	  for (sentence_type::const_iterator iter = siter; iter != std::min(siter + order, siter_end); ++ iter, ++ n) {
	    id = ngrams.find(id, *iter);
	    
	    if (ngrams.is_root(id)) break;
	    
	    // ngram at [siter, iter] with id
	    ++ counts[n][id];
	  }
	}

	// collect clip counts...
	for (int n = 0; n < order; ++ n) {
	  counts_type::const_iterator citer_end = counts[n].end();
	  for (counts_type::const_iterator citer = counts[n].begin(); citer != citer_end; ++ citer)
	    bleu->ngrams_hypothesis[n] += std::min(citer->second, ngrams[citer->first]);
	}
	
	return score_ptr_type(bleu.release());
      }

    private:
      ngram_set_type ngrams;
      size_set_type  sizes;
      int            order;
    };
  };
};

#endif
