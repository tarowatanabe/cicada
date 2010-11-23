// -*- mode: c++ -*-

#ifndef __CICADA__EVAL__SCORE__HPP__
#define __CICADA__EVAL__SCORE__HPP__ 1

// base-class for scoring...
#include <vector>
#include <utility>

#include <cicada/sentence.hpp>
#include <cicada/sentence_vector.hpp>
#include <cicada/tokenizer.hpp>

#include <boost/shared_ptr.hpp>

namespace cicada
{
  namespace eval
  {
    class Score
    {
    public:
      typedef Score score_type;
      typedef boost::shared_ptr<score_type>  score_ptr_type;
      
    public:
      virtual ~Score() {}
      
      
      std::pair<double, double> operator()() const { return score(); }
      virtual std::pair<double, double> score() const = 0;
      
      virtual void assign(const score_type& score) = 0;
      
      virtual void plus_equal(const score_type& score) = 0;
      virtual void minus_equal(const score_type& score) = 0;
      
      virtual void multiplies_equal(const double& scale) = 0;
      virtual void divides_equal(const double& scale) = 0;

      virtual score_ptr_type zero() const = 0;
      virtual score_ptr_type clone() const = 0;
      
    public:
      
      score_type& operator=(const score_type& score)
      {
	this->assign(score);
	return *this;
      }
      
      score_type& operator+=(const score_type& score)
      {
	this->plus_equal(score);
	return *this;
      }
      
      score_type& operator-=(const score_type& score)
      {
	this->minus_equal(score);
	return *this;
      }
      
      score_type& operator*=(const double& scale)
      {
	this->multiplies_equal(scale);
	return *this;
      }
      
      score_type& operator/=(const double& scale)
      {
	this->divides_equal(scale);
	return *this;
      }
      
    };
    


    class Scorer
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol         word_type;
      typedef cicada::Sentence       sentence_type;
      typedef cicada::SentenceVector sentence_set_type;
      
      typedef Score score_type;
      typedef Scorer scorer_type;
      
      typedef boost::shared_ptr<score_type>  score_ptr_type;
      typedef boost::shared_ptr<scorer_type> scorer_ptr_type;

      typedef cicada::Tokenizer tokenizer_type;
      
    public:
      Scorer() : tokenizer(0) {}
      Scorer(const Scorer& x)
	: tokenizer(0)
      {
	if (x.tokenizer)
	  tokenizer = &tokenizer_type::create(x.tokenizer->algorithm());
      }
      
      virtual ~Scorer() {}

      Scorer& operator=(const Scorer& x)
      {
	tokenizer = 0;
	if (x.tokenizer)
	  tokenizer = &tokenizer_type::create(x.tokenizer->algorithm());
	return *this;
      }
      
      // insert a sentence for scoring
      virtual void clear() = 0;
      virtual void insert(const sentence_type& sentence) = 0;

      score_ptr_type operator()(const sentence_type& sentence) const { return score(sentence); }
      virtual score_ptr_type score(const sentence_type& sentence) const = 0;
      virtual bool error_metric() const = 0;
      virtual scorer_ptr_type clone() const = 0;
      
      static const char*     lists();
      static scorer_ptr_type create(const std::string& parameter);

      void tokenize(const sentence_type& source, sentence_type& tokenized) const
      {
	if (tokenizer)
	  tokenizer->operator()(source, tokenized);
	else
	  tokenized = source;
      }
      
    protected:
      const tokenizer_type* tokenizer;
    };    
    
    class ScorerDocument
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Sentence sentence_type;
      typedef cicada::Symbol   word_type;

      typedef Scorer scorer_type;

      typedef scorer_type::score_type      score_type;
      typedef scorer_type::score_ptr_type  score_ptr_type;
      typedef scorer_type::scorer_ptr_type scorer_ptr_type;

      
    private:
      typedef std::vector<scorer_ptr_type, std::allocator<scorer_ptr_type> > scorer_set_type;

    public:
      typedef scorer_set_type::const_iterator const_iterator;
      typedef scorer_set_type::iterator       iterator;
      
    public:
      ScorerDocument(const std::string& __parameter)
	: parameter(__parameter) {}
      
      inline       scorer_ptr_type& operator[](size_type pos)       { return scorers[pos]; }
      inline const scorer_ptr_type& operator[](size_type pos) const { return scorers[pos]; }

      inline const_iterator begin() const { return scorers.begin(); }
      inline       iterator begin()       { return scorers.begin(); }
      
      inline const_iterator end() const { return scorers.end(); }
      inline       iterator end()       { return scorers.end(); }
      
      bool error_metric() const
      {
	return scorer_type::create(parameter)->error_metric();
      }
      
      void resize(size_type x) { scorers.resize(x); }
      void clear() { scorers.clear(); }

      size_type size() const { return scorers.size(); }
      bool empty() const { return scorers.empty(); }
      
      scorer_ptr_type create() { return scorer_type::create(parameter); }
      
    private:
      scorer_set_type scorers;
      std::string     parameter;
    };
    
  };
};

#endif
