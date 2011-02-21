//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "decode.hpp"
#include "encode.hpp"

#include <sstream>
#include <iterator>

#include "f.hpp"
#include "sk.hpp"
#include "sb.hpp"
#include "wlcs.hpp"

namespace cicada
{
  namespace eval
  {
    std::string F::description() const
    {
      std::ostringstream stream;
      stream << __description() << ": " << score()
	     << ' ' << (norm_hyp != 0.0 ? match_hyp / norm_hyp : 0.0)
	     << '|' << (norm_ref != 0.0 ? match_ref / norm_ref : 0.0);
      
      return stream.str();
    }
    
    std::string F::encode() const
    {
      std::ostringstream stream;
      stream << '{' << "\"eval\":\"" << __description() << "\",";
      stream << "\"reference\":[";
      stream << "\"" << escaper(match_ref) << "\",";
      stream << "\"" << escaper(norm_ref) << "\"";
      stream << "],";
      stream << "\"hypothesis\":[";
      stream << "\"" << escaper(match_hyp) << "\",";
      stream << "\"" << escaper(norm_hyp) << "\"";
      stream << "]";
      stream << '}';
      return stream.str();
    }
    
    typedef boost::fusion::tuple<std::string, double, double, double, double> f_parsed_type;
    
    
    template <typename Iterator>
    struct f_parser : boost::spirit::qi::grammar<Iterator, f_parsed_type(), boost::spirit::standard::space_type>
    {
      f_parser() : f_parser::base_type(f_parsed)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	f_parsed %= (qi::lit('{')
		     >> qi::lit("\"eval\"") >> qi::lit(':') >> qi::lit('\"') >> qi::lexeme[+(~standard::char_('\"'))] >> qi::lit('\"') >> qi::lit(',')
		     >> qi::lit("\"reference\"") >> qi::lit(':')
		     >> qi::lit('[') >> double_value >> qi::lit(',') >> double_value >> qi::lit(']') >> qi::lit(',')
		     >> qi::lit("\"hypothesis\"") >> qi::lit(':')
		     >> qi::lit('[') >> double_value >> qi::lit(',') >> double_value >> qi::lit(']')
		     >> qi::lit('}'));
	
      }
      
      typedef boost::spirit::standard::space_type space_type;
      
      double_base64_parser<Iterator> double_value;
      boost::spirit::qi::rule<Iterator, f_parsed_type(), space_type> f_parsed;
    };
    
    Score::score_ptr_type F::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      typedef std::string::const_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      f_parser<iterator_type> parser;
      f_parsed_type           parsed;
      
      const bool result = qi::phrase_parse(iter, end, parser, standard::space, parsed);
      if (! result)
	return score_ptr_type();

      if (boost::fusion::get<0>(parsed) == "sk") {
	
	std::auto_ptr<SK> sk(new SK());
	sk->match_ref = boost::fusion::get<1>(parsed);
	sk->norm_ref  = boost::fusion::get<2>(parsed);
	sk->match_hyp = boost::fusion::get<3>(parsed);
	sk->norm_hyp  = boost::fusion::get<4>(parsed);
	      
	return score_ptr_type(sk.release());
      } else if (boost::fusion::get<0>(parsed) == "sb") {
	std::auto_ptr<SB> sb(new SB());
	sb->match_ref = boost::fusion::get<1>(parsed);
	sb->norm_ref  = boost::fusion::get<2>(parsed);
	sb->match_hyp = boost::fusion::get<3>(parsed);
	sb->norm_hyp  = boost::fusion::get<4>(parsed);
	      
	return score_ptr_type(sb.release());
      } else if (boost::fusion::get<0>(parsed) == "wlcs") {
	std::auto_ptr<WLCS> wlcs(new WLCS());
	wlcs->match_ref = boost::fusion::get<1>(parsed);
	wlcs->norm_ref  = boost::fusion::get<2>(parsed);
	wlcs->match_hyp = boost::fusion::get<3>(parsed);
	wlcs->norm_hyp  = boost::fusion::get<4>(parsed);
	
	return score_ptr_type(wlcs.release());
      } else
	return score_ptr_type();
    }
    
    Score::score_ptr_type F::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }
    
  };  
};
