// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//
// BleuS scorer

#ifndef __CICADA__EVAL__BLEUS__HPP__
#define __CICADA__EVAL__BLEUS__HPP__ 1

#include <cmath>
#include <cfloat>

#include <stdexcept>
#include <algorithm>
#include <vector>
#include <utility>

#include <cicada/eval/score.hpp>
#include <cicada/symbol_vector.hpp>

#include <utils/trie_compact.hpp>
#include <utils/simple_vector.hpp>

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace eval
  {
    class BleuS : public Score
    {
    public:
      typedef double count_type;

    public:
      BleuS() : bleu(0), penalty(0) {}
      
      double score() const
      {
	return (penalty == 0.0 ? 0.0 : bleu / penalty);
      }

      double loss() const { return 1.0 - score(); }

      double reward() const { return score(); }
      
      bool error_metric() const { return false; }

      bool equal(const score_type& score) const
      {
	const BleuS* rhs = dynamic_cast<const BleuS*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEUS score");
	
	return (bleu == rhs->bleu && penalty == rhs->penalty);
      }

      void assign(const score_type& score)
      {
	const BleuS* rhs = dynamic_cast<const BleuS*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEUS score");
	
	bleu    = rhs->bleu;
	penalty = rhs->penalty;
      }
      
      void plus_equal(const score_type& score)
      {
	const BleuS* rhs = dynamic_cast<const BleuS*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEUS score");

	bleu    += rhs->bleu;
	penalty += rhs->penalty;
      }
      
      void minus_equal(const score_type& score)
      {
	const BleuS* rhs = dynamic_cast<const BleuS*>(&score);
	if (! rhs)
	  throw std::runtime_error("invalid BLEUS score");

	bleu    -= rhs->bleu;
	penalty -= rhs->penalty;
      }
      
      void multiplies_equal(const double& scale)
      {
	bleu    *= scale;
	penalty *= scale;
      }
      
      void divides_equal(const double& scale)
      {
	bleu    /= scale;
	penalty /= scale;
      }
      
      score_ptr_type zero() const
      {
	return score_ptr_type(new BleuS());
      }

      score_ptr_type clone() const
      {
	return score_ptr_type(new BleuS(*this));
      }

      std::string description() const;
      std::string encode() const;

      static score_ptr_type decode(std::string::const_iterator& iter, std::string::const_iterator end);
      static score_ptr_type decode(const utils::piece& encoded);

      count_type bleu;
      count_type penalty;
    };
    

    class BleuSScorer : public Scorer
    {
    private:
      typedef double count_type;
      
      typedef std::allocator<std::pair<const word_type, count_type> > ngram_allocator_type;
      typedef utils::trie_compact<word_type, count_type,
				  utils::unassigned<word_type>, 
				  boost::hash<word_type>, std::equal_to<word_type>,
				  ngram_allocator_type> ngram_set_type;
      
      typedef std::vector<int, std::allocator<int> > size_set_type;

    public:
      BleuSScorer(int __order = 4) : ngrams(), sizes(), order(__order) { }

      bool error_metric() const { return false; }

      scorer_ptr_type clone() const { return scorer_ptr_type(new BleuSScorer(*this)); }

      void clear()
      {
	ngrams.clear();
	sizes.clear();
      }
      
      void insert(const sentence_type& __sentence);
      
      score_ptr_type score(const sentence_type& __sentence) const;
      
      // bleus specific interface...
      count_type reference_length(const double& length) const;
      count_type find(const SymbolVector& ngram) const;

    private:
      ngram_set_type ngrams;
      size_set_type  sizes;
      
      int            order;
    };
  };
};

#endif
