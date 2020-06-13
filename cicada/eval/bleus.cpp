//
//  Copyright(C) 2012-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"
#include "bleus.hpp"

#include <sstream>
#include <iterator>
#include <memory>

#include "utils/bithack.hpp"
#include "utils/unordered_map.hpp"

#include <boost/functional/hash/hash.hpp>

namespace cicada
{
  namespace eval
  {
    std::string BleuS::description() const
    {
      std::ostringstream stream;
      
      stream << "bleus: " << score() << " bleu: " << bleu << " norm: " << norm;
      
      return stream.str();
    }
    
    
    std::string BleuS::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"bleus\",";
      stream << "\"bleu\":" << escaper(bleu) 
	     << ','
	     << "\"norm\":" << escaper(norm)
	     << '}';
      
      return stream.str();
    }

    typedef std::pair<double, double> bleus_parsed_pair_type;

    template <typename Iterator>
    struct bleus_parser : boost::spirit::qi::grammar<Iterator, bleus_parsed_pair_type(), boost::spirit::standard::space_type>
    {
      bleus_parser() : bleus_parser::base_type(bleus_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	bleus_parsed %= (qi::lit('{')
			 >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"bleus\"") >> qi::lit(',')
			 >> qi::lit("\"bleu\"") >> qi::lit(':') >> double_value >> qi::lit(',')
			 >> qi::lit("\"norm\"") >> qi::lit(':') >> double_value
			 >> qi::lit('}'));
	
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, bleus_parsed_pair_type(), space_type> bleus_parsed;
    };

    Score::score_ptr_type BleuS::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      const char* citer_begin = &(*iter);
      const char* citer       = &(*iter);
      const char* citer_end   = &(*end);
      
      score_ptr_type result = decode(citer, citer_end);
      
      iter += citer - citer_begin;
      
      return result;
    }

    Score::score_ptr_type BleuS::decode(utils::piece::const_iterator& iter, utils::piece::const_iterator end)
    {
      typedef utils::piece::const_iterator iterator_type;

      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      bleus_parser<iterator_type> parser;
      bleus_parsed_pair_type bleus_parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, bleus_parsed);
      if (! result)
	return score_ptr_type();
      
      std::unique_ptr<BleuS> bleus(new BleuS());
      
      bleus->bleu = bleus_parsed.first;
      bleus->norm = bleus_parsed.second;
      
      return score_ptr_type(bleus.release());
    }
    
    Score::score_ptr_type BleuS::decode(const utils::piece& encoded)
    {
      utils::piece::const_iterator iter(encoded.begin());
      utils::piece::const_iterator iter_end (encoded.end());
      
      return decode(iter, iter_end);
    }
    
    void BleuSScorer::insert(const sentence_type& __sentence)
    {
      typedef ngram_set_type::id_type id_type;
      typedef utils::unordered_map<id_type, count_type, boost::hash<id_type>, std::equal_to<id_type>, std::allocator<std::pair<const id_type, count_type> > >::type counts_type;

      sentence_type sentence;
      counts_type counts;
	
      tokenize(__sentence, sentence);
	
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
    
    BleuSScorer::score_ptr_type BleuSScorer::score(const sentence_type& __sentence) const
    {
      typedef ngram_set_type::id_type id_type;
      typedef utils::unordered_map<id_type, count_type, boost::hash<id_type>, std::equal_to<id_type>, std::allocator<std::pair<const id_type, count_type> > >::type counts_type;
	
      typedef std::vector<counts_type, std::allocator<counts_type> > counts_set_type;
      typedef utils::simple_vector<count_type, std::allocator<count_type> > ngram_counts_type;
	
      sentence_type sentence;
      tokenize(__sentence, sentence);
      
      std::unique_ptr<BleuS> bleus(new BleuS());
      counts_set_type counts(order);
      
      ngram_counts_type ngrams_hypothesis(order, 0);
      ngram_counts_type ngrams_matched(order, 0);
      count_type        length_reference(0);
      count_type        length_hypothesis(0);
      
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
	
      length_hypothesis += hypothesis_size;
      length_reference  += reference_size;
	
      // collect total counts...
      for (int n = 0; n < utils::bithack::min(order, hypothesis_size); ++ n)
	ngrams_hypothesis[n] += hypothesis_size - n;
	
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
	  ngrams_matched[n] += std::min(citer->second, ngrams[citer->first]);
      }
      
      bleus->norm = 1;
      if (ngrams_matched[0] == 0.0)
	bleus->bleu = 0.0;
      else {
	const double penalty = std::min(1.0 - length_reference / length_hypothesis, 0.0);
	const int norm = utils::bithack::min(order, hypothesis_size);
	
	double score = 0.0;
	for (size_t n = 0; n < ngrams_hypothesis.size(); ++ n)
	  score += std::log(ngrams_matched[n] + (n != 0)) - std::log(ngrams_hypothesis[n] + (n != 0));
	
	bleus->bleu = std::exp(score / norm + penalty);
      }
      
      return score_ptr_type(bleus.release());
    }
    
    BleuSScorer::count_type BleuSScorer::reference_length(const double& length) const
    {
      int reference_size = 0;
      double min_diff = boost::numeric::bounds<double>::highest();
      
      for (size_set_type::const_iterator siter = sizes.begin(); siter != sizes.end(); ++ siter) {
	const double diff = ::fabs(*siter - length);
	
	if (diff < min_diff) {
	  min_diff = diff;
	  reference_size = *siter;
	} else if (diff == min_diff)
	  reference_size = utils::bithack::min(reference_size, *siter);
      }
      
      return reference_size;
    }
    
    BleuSScorer::count_type BleuSScorer::find(const SymbolVector& ngram) const
    {
      typedef SymbolVector ngram_type;
      
      ngram_set_type::id_type id = ngrams.root();
      ngram_type::const_iterator iter_end = ngram.end();
      for (ngram_type::const_iterator iter = ngram.begin(); iter != iter_end; ++ iter) {
	id = ngrams.find(id, *iter);
	
	if (ngrams.is_root(id)) break;
      }
      
      return (ngrams.is_root(id) ? 0.0 : ngrams[id]);
    }
    
  };
};
