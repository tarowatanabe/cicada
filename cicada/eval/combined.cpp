//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>

#include "combined.hpp"

namespace cicada
{
  namespace eval
  {
    std::string Combined::description() const
    {
      std::ostringstream stream;
      stream << "combined: " << score() << ' ';
      
      stream << '{';
      if (! scores.empty()) {
	score_ptr_set_type::const_iterator siter_end = scores.end() - 1;
	for (score_ptr_set_type::const_iterator siter = scores.begin(); siter != siter_end; ++ siter)
	  stream << (*siter)->description() << ", ";
	stream << scores.back()->description();
      }
      stream << '}';
      
      return stream.str();
    }
    
    std::string Combined::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"combined\",";
      stream << "\"score\":[";
      if (! scores.empty()) {
	for (size_t i = 0; i != scores.size() - 1; ++ i)
	  stream << scores[i]->encode() << ',';
	stream << scores.back()->encode();
      }
      stream << "],";
      stream << "\"weight\":[";
      if (! weights.empty()) {
	for (size_t i = 0; i != weights.size() - 1; ++ i)
	  stream << '\"' << escaper(weights[i]) << "\",";
	stream << '\"' << escaper(weights.back()) << '\"';
      }
      stream << "]";
      stream << '}';
      
      return stream.str();
    }    
    
    Score::score_ptr_type Combined::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;

      std::auto_ptr<Combined> combined(new Combined());

      if (! qi::phrase_parse(iter, end, qi::lit('{') >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit("\"combined\"") >> qi::lit(','), standard::space))
	return score_ptr_type();
      if (! qi::phrase_parse(iter, end, qi::lit("\"score\"") >> qi::lit(':') >> qi::lit('['), standard::space))
	return score_ptr_type();
      
      while (iter != end) {
	score_ptr_type score = Score::decode(iter, end);

	if (! score)
	  std::cerr << "remaining: " << std::string(iter, end) << std::endl;
	
	if (! score) break;
	
	combined->scores.push_back(score);
	
	if (! qi::phrase_parse(iter, end, qi::lit(','), standard::space))
	  break;
      }
      
      if (! qi::phrase_parse(iter, end, qi::lit(']') >> qi::lit(','), standard::space))
	return score_ptr_type();
      
      double_base64_parser<iterator_type> double_value;
      if (! qi::phrase_parse(iter, end, (qi::lit("\"weight\"") >> qi::lit(':') >> qi::lit('[') >> (double_value % ',') >> qi::lit(']') >> qi::lit('}')),
			     standard::space,
			     combined->weights))
	return score_ptr_type();
      
      return score_ptr_type(combined.release());
    }
    
    Score::score_ptr_type Combined::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }

  };
};
