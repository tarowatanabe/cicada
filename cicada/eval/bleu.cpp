//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"
#include "bleu.hpp"

#include <sstream>
#include <iterator>
#include <memory>


namespace cicada
{
  namespace eval
  {
    std::string Bleu::description() const
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
      stream << "bleu: " << std::exp(score) << ' ';
      if (! precisions.empty()) {
	std::copy(precisions.begin(), precisions.end() - 1, std::ostream_iterator<double>(stream, "|"));
	stream << precisions.back();
      }
      stream << " penalty: " << std::exp(penalty);
      
      return stream.str();
    }
    
    
    std::string Bleu::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"bleu\",";
      stream << "\"reference\":[";
      stream << "\"" << escaper(length_reference) << "\"";
      for (size_t i = 0; i != ngrams_reference.size(); ++ i)
	stream << ",\"" << escaper(ngrams_reference[i]) << "\"";
      stream << "],";
      stream << "\"hypothesis\":[";
      stream << "\"" << escaper(length_hypothesis) << "\"";
      for (size_t i = 0; i != ngrams_hypothesis.size(); ++ i)
	stream << ",\"" << escaper(ngrams_hypothesis[i]) << "\"";
      stream << "]";
      stream << '}';
      
      return stream.str();
    }

    typedef std::vector<double, std::allocator<double> > bleu_parsed_type;
    typedef std::pair<bleu_parsed_type, bleu_parsed_type> bleu_parsed_pair_type;

    template <typename Iterator>
    struct bleu_parser : boost::spirit::qi::grammar<Iterator, bleu_parsed_pair_type(), boost::spirit::standard::space_type>
    {
      bleu_parser() : bleu_parser::base_type(bleu_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	double_values %= double_value % ',';
	
	bleu_parsed %= (qi::lit('{')
			>> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"bleu\"") >> qi::lit(',')
			>> qi::lit("\"reference\"") >> qi::lit(':')
			>> qi::lit('[') >> double_values >> qi::lit(']') >> qi::lit(',')
			>> qi::lit("\"hypothesis\"") >> qi::lit(':')
			>> qi::lit('[') >> double_values >> qi::lit(']')
			>> qi::lit('}'));
	  
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, bleu_parsed_type(), space_type> double_values;
      boost::spirit::qi::rule<Iterator, bleu_parsed_pair_type(), space_type> bleu_parsed;
    };

    Score::score_ptr_type Bleu::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;

      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      bleu_parser<iterator_type> parser;
      bleu_parsed_pair_type bleu_parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, bleu_parsed);
      if (! result)
	return score_ptr_type();
      
      std::auto_ptr<Bleu> bleu(new Bleu(0));
      
      bleu->length_reference = bleu_parsed.first.front();
      bleu->ngrams_reference.insert(bleu->ngrams_reference.end(), bleu_parsed.first.begin() + 1, bleu_parsed.first.end());
      bleu->length_hypothesis = bleu_parsed.second.front();
      bleu->ngrams_hypothesis.insert(bleu->ngrams_hypothesis.end(), bleu_parsed.second.begin() + 1, bleu_parsed.second.end());
      
      return score_ptr_type(bleu.release());
    }
    
    Score::score_ptr_type Bleu::decode(const std::string& encoded)
    {
      std::string::const_iterator iter     = encoded.begin();
      std::string::const_iterator iter_end = encoded.end();
      
      return decode(iter, iter_end);
    }
  };
};
