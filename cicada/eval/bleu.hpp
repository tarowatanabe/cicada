// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__EVAL__BLEU__HPP__
#define __CICADA__EVAL__BLEU__HPP__ 1

#include <cmath>
#include <cfloat>

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <utility>

#include <cicada/eval/score.hpp>

#include <utils/compact_trie_dense.hpp>

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
      
      double score() const
      {
	const double penalty = std::min(1.0 - length_reference / length_hypothesis, 0.0);
	
	double smooth = 0.5;
	double score = 0.0;
	
	for (size_t n = 0; n < ngrams_hypothesis.size(); ++ n) {
	  const double p = (ngrams_reference[n] > 0
			    ? (ngrams_hypothesis[n] > 0 ? ngrams_hypothesis[n] : smooth) / ngrams_reference[n]
			    : 0.0);
	  
	  score += p > 0.0 ? std::log(p) : 0.0;
	  smooth *= 0.5;
	}
	
	score /= ngrams_hypothesis.size();
	score += penalty;
	
	return std::exp(score);
      }

      bool equal(const score_type& score) const
      {
	const Bleu* rhs = dynamic_cast<const Bleu*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEU score");
	
	return (ngrams_reference == rhs->ngrams_reference
		&& ngrams_hypothesis == rhs->ngrams_hypothesis
		&& length_reference == rhs->length_reference
		&& length_hypothesis == rhs->length_hypothesis);
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

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);

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
      typedef utils::compact_trie_dense<word_type, count_type, boost::hash<word_type>, std::equal_to<word_type>, ngram_allocator_type> ngram_set_type;
      
      typedef std::vector<int, std::allocator<int> > size_set_type;

    public:
      BleuScorer(int __order = 4) : ngrams(word_type()), sizes(), order(__order) { }

      bool error_metric() const { return false; }

      scorer_ptr_type clone() const { return scorer_ptr_type(new BleuScorer(*this)); }

      void clear()
      {
	ngrams.clear();
	sizes.clear();
      }
      
      void insert(const sentence_type& __sentence);

      score_ptr_type score(const sentence_type& __sentence) const;

    private:
      ngram_set_type ngrams;
      size_set_type  sizes;
      
      int            order;
    };
  };
};

#endif
