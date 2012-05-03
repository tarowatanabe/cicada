//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"
#include "bleus.hpp"

#include <map>
#include <sstream>
#include <iterator>
#include <memory>

#include "utils/bithack.hpp"

namespace cicada
{
  namespace eval
  {
    std::string BleuS::description() const
    {
      std::vector<double, std::allocator<double> > precisions;
      
      const double penalty = std::min(1.0 - length_reference / length_hypothesis, 0.0);
      
      double smooth = 0.5;
      double score = 0.0;
      
      for (size_t n = 0; n < ngrams_hypothesis.size(); ++ n) {
	const double p = (ngrams_reference[n] > 0
			  ? (ngrams_hypothesis[n] > 0 ? ngrams_hypothesis[n] : smooth) / ngrams_reference[n]
			  : 0.0);
	
	precisions.push_back(p);
	
	score += p > 0.0 ? std::log(p) : 0.0;
	smooth *= 0.5;
      }
      
      score /= ngrams_hypothesis.size();
      score += penalty;
      
      std::ostringstream stream;
      stream << "bleus: " << std::exp(score) << ' ';
      if (! precisions.empty()) {
	std::copy(precisions.begin(), precisions.end() - 1, std::ostream_iterator<double>(stream, "|"));
	stream << precisions.back();
      }
      stream << " penalty: " << std::exp(penalty);
      
      return stream.str();
    }
    
    
    std::string BleuS::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"bleus\",";
      stream << "\"reference\":[";
      stream << escaper(length_reference);
      for (size_t i = 0; i != ngrams_reference.size(); ++ i)
	stream << ',' << escaper(ngrams_reference[i]);
      stream << "],";
      stream << "\"hypothesis\":[";
      stream << escaper(length_hypothesis);
      for (size_t i = 0; i != ngrams_hypothesis.size(); ++ i)
	stream << ',' << escaper(ngrams_hypothesis[i]);
      stream << "]";
      stream << '}';
      
      return stream.str();
    }

    typedef std::vector<double, std::allocator<double> > bleus_parsed_type;
    typedef std::pair<bleus_parsed_type, bleus_parsed_type> bleus_parsed_pair_type;

    template <typename Iterator>
    struct bleus_parser : boost::spirit::qi::grammar<Iterator, bleus_parsed_pair_type(), boost::spirit::standard::space_type>
    {
      bleus_parser() : bleus_parser::base_type(bleus_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	double_values %= double_value % ',';
	
	bleus_parsed %= (qi::lit('{')
			 >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"bleus\"") >> qi::lit(',')
			 >> qi::lit("\"reference\"") >> qi::lit(':')
			 >> qi::lit('[') >> double_values >> qi::lit(']') >> qi::lit(',')
			 >> qi::lit("\"hypothesis\"") >> qi::lit(':')
			 >> qi::lit('[') >> double_values >> qi::lit(']')
			 >> qi::lit('}'));
	  
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, bleus_parsed_type(), space_type> double_values;
      boost::spirit::qi::rule<Iterator, bleus_parsed_pair_type(), space_type> bleus_parsed;
    };

    Score::score_ptr_type BleuS::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;

      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      bleus_parser<iterator_type> parser;
      bleus_parsed_pair_type bleus_parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, bleus_parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<BleuS> bleus(new BleuS(0));
      
      bleus->length_reference = bleus_parsed.first.front();
      bleus->ngrams_reference.insert(bleus->ngrams_reference.end(), bleus_parsed.first.begin() + 1, bleus_parsed.first.end());
      bleus->length_hypothesis = bleus_parsed.second.front();
      bleus->ngrams_hypothesis.insert(bleus->ngrams_hypothesis.end(), bleus_parsed.second.begin() + 1, bleus_parsed.second.end());
      
      return score_ptr_type(bleus.release());
    }
    
    Score::score_ptr_type BleuS::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end (encoded.end());
      
      return decode(iter, iter_end);
    }
    
    void BleuSScorer::insert(const sentence_type& __sentence)
    {
      typedef ngram_set_type::id_type id_type;
      typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;

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
      typedef std::map<id_type, count_type, std::less<id_type>, std::allocator<std::pair<const id_type, count_type> > > counts_type;
	
      typedef std::vector<counts_type, std::allocator<counts_type> > counts_set_type;
	
      sentence_type sentence;
      tokenize(__sentence, sentence);
	
      std::auto_ptr<BleuS> bleus(new BleuS(order));
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
	
      bleus->length_hypothesis += hypothesis_size;
      bleus->length_reference  += reference_size;
	
      // collect total counts...
      for (int n = 0; n < utils::bithack::min(order, hypothesis_size); ++ n)
	bleus->ngrams_reference[n] += hypothesis_size - n;
	
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
	  bleus->ngrams_hypothesis[n] += std::min(citer->second, ngrams[citer->first]);
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
