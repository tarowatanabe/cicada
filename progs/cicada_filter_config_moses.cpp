//
//  Copyright(C) 2011-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// a simple moses.ini filter (coupled with weights file...)
//
// read weights files, and replace moses.ini weights (if none, use zero-weight)
//

#include <string>
#include <vector>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/fusion/tuple.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/tokenizer.hpp>
#include <boost/xpressive/xpressive.hpp>

#include "utils/lexical_cast.hpp"
#include "utils/compress_stream.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/hashmurmur.hpp"

typedef boost::filesystem::path path_type;

struct feature_hash : public utils::hashmurmur<size_t>
{
  typedef utils::hashmurmur<size_t> hasher_type;
  
  size_t operator()(const std::string& x) const
  {
    return hasher_type::operator()(x.begin(), x.end(), 0);
  }
};

typedef boost::fusion::tuple<std::string, int, double> feature_weight_type;

typedef std::vector<double, std::allocator<double> > weight_set_type;

typedef utils::unordered_map<std::string, weight_set_type, feature_hash, std::equal_to<std::string>,
			     std::allocator<std::pair<const std::string, weight_set_type> > >::type feature_weight_set_type;
typedef utils::unordered_map<std::string, std::string, feature_hash, std::equal_to<std::string>,
			     std::allocator<std::pair<const std::string, std::string> > >::type feature_map_type;

typedef utils::unordered_set<std::string, feature_hash, std::equal_to<std::string>, std::allocator<std::string> >::type feature_consumed_type;

template <typename Iterator>
struct feature_weight_parser : boost::spirit::qi::grammar<Iterator, feature_weight_type()>
{
  feature_weight_parser() : feature_weight_parser::base_type(feature_weight)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    feature_weight %= (qi::omit[*standard::blank] >> qi::lexeme[+(standard::char_ - standard::space - ':')] >> ':' >> qi::int_
		       >> qi::omit[+standard::blank] >> qi::double_
		       >> qi::omit[*standard::blank] >> (qi::eol | qi::eoi));
  }
  
  boost::spirit::qi::rule<Iterator, feature_weight_type()>  feature_weight;
};

typedef std::vector<std::string, std::allocator<std::string> > feature_set_type;

path_type input_file = "-";
path_type output_file = "-";

path_type weights_file;

feature_set_type features_biased;
double           weight_biased = -1;

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);
    
    if (weights_file.empty()) {
      // we will extract feature value!
      
      feature_map_type feature_map;
      feature_map["weight-d"]          = "d";
      feature_map["weight-lr"]         = "lr";
      feature_map["weight-generation"] = "g";
      feature_map["weight-i"]          = "I";
      feature_map["weight-l"]          = "lm";
      feature_map["weight-lex"]        = "l";
      feature_map["weight-t"]          = "tm";
      feature_map["weight-w"]          = "w";
      feature_map["weight-u"]          = "u";
      feature_map["weight-e"]          = "e";
      feature_map["weight-slm"]        = "slm";
      
      utils::compress_istream is(input_file);
      utils::compress_ostream os(output_file);
      os.precision(10);
    
      namespace xpressive = boost::xpressive;

      xpressive::sregex pattern_comment = *xpressive::_s >> '#' >> *xpressive::_;
      xpressive::sregex pattern_empty = *xpressive::_s;
      
      xpressive::sregex pattern_section_weight = (*xpressive::_s
						  >> '[' >> (xpressive::s1= "weight-" >> +(~xpressive::_s)) >> ']'
						  >> *xpressive::_s);
      
      xpressive::sregex pattern_section = (*xpressive::_s
					   >> '[' >> (xpressive::s1= +(~xpressive::_s)) >> ']'
					   >> *xpressive::_s);
      
      std::string weight_prefix;
      std::string line;
      bool consumed = std::getline(is, line);
      while (consumed) {
	xpressive::smatch what;
      
	if (! xpressive::regex_match(line, what, pattern_section_weight)) {
	  //os << line << '\n';
	  consumed = std::getline(is, line);
	  continue;
	}
	
	// this is the section for "weight"
	//os << line << '\n';
	
	const std::string section = what[1];
	
	feature_map_type::const_iterator fiter = feature_map.find(section);
	
	weight_prefix = (fiter == feature_map.end() ? std::string("") : fiter->second);
	
	int pos = 0;
	do {
	  consumed = std::getline(is, line);
	  if (! consumed || xpressive::regex_match(line, what, pattern_section)) break;
	  
	  if (xpressive::regex_match(line, what, pattern_comment) || xpressive::regex_match(line, what, pattern_empty)) {
	    //os << line << '\n';
	  } else {
	    if (! weight_prefix.empty())
	      os << weight_prefix << ':' << utils::lexical_cast<std::string>(pos) << ' ' << line << '\n';
	    
	    ++ pos;
	  }
	} while (consumed);
      }
      
    } else {
      if (weights_file != "-" && ! boost::filesystem::exists(weights_file))
	throw std::runtime_error("no weight file? " + weights_file.string());

      feature_map_type feature_map;
      feature_map["d"]   = "weight-d";
      feature_map["lr"]  = "weight-lr";
      feature_map["g"]   = "weight-generation";
      feature_map["I"]   = "weight-i";
      feature_map["lm"]  = "weight-l";
      feature_map["l"]   = "weight-lex";
      feature_map["tm"]  = "weight-t";
      feature_map["w"]   = "weight-w";
      feature_map["u"]   = "weight-u";
      feature_map["e"]   = "weight-e";
      feature_map["slm"] = "weight-slm";

      feature_weight_set_type feature_weights;
      feature_weight_type feature_weight;
      
      {
	utils::compress_istream is(weights_file, 1024 * 1024);
	
	std::string line;
	feature_weight_parser<std::string::const_iterator> parser;
      
	while (std::getline(is, line)) {
	  boost::fusion::get<0>(feature_weight).clear();
	  
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator iter_end  = line.end();
	
	  if (! boost::spirit::qi::parse(iter, iter_end, parser, feature_weight))
	    if (iter != iter_end) {
	      std::cerr << "WARNING: ignoring: " << line << std::endl;
	      continue;
	    }
	
	  const std::string& feature = boost::fusion::get<0>(feature_weight);
	  const int&         pos     = boost::fusion::get<1>(feature_weight);
	  const double&      value   = boost::fusion::get<2>(feature_weight);

	  if (feature.empty()) continue;

	  std::string feat = feature;
	  feature_map_type::const_iterator miter = feature_map.find(feature);
	  if (miter != feature_map.end())
	    feat = miter->second;
	
	  weight_set_type& weights = feature_weights[feat];
	  if (pos >= weights.size())
	    weights.resize(pos + 1);
	  weights[pos] = value;
	}
      }
    
      if (! features_biased.empty()) {
	feature_weight_parser<std::string::const_iterator> parser;

	for (feature_set_type::const_iterator fiter = features_biased.begin(); fiter != features_biased.end(); ++ fiter) {
	  const std::string line = *fiter + " " + utils::lexical_cast<std::string>(weight_biased);
	
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator end  = line.end();

	  boost::fusion::get<0>(feature_weight).clear();

	  // if failed, continue!
	  if (! boost::spirit::qi::parse(iter, end, parser, feature_weight))
	    if (iter != end) {
	      std::cerr << "WARNING: ignoring: " << line << std::endl;
	      continue;
	    }
	
	  const std::string& feature = boost::fusion::get<0>(feature_weight);
	  const int&         pos     = boost::fusion::get<1>(feature_weight);
	  const double&      value   = boost::fusion::get<2>(feature_weight);
	
	  if (feature.empty()) continue;
	
	  std::string feat = feature;
	  feature_map_type::const_iterator miter = feature_map.find(feature);
	  if (miter != feature_map.end())
	    feat = miter->second;
	
	  weight_set_type& weights = feature_weights[feat];
	  if (pos >= weights.size())
	    weights.resize(pos + 1);
	  weights[pos] = value;
	}
      }
    
    
      utils::compress_istream is(input_file);
      utils::compress_ostream os(output_file);
      os.precision(10);

      feature_consumed_type features_consumed;
    
      namespace xpressive = boost::xpressive;

      xpressive::sregex pattern_comment = *xpressive::_s >> '#' >> *xpressive::_;
      xpressive::sregex pattern_empty = *xpressive::_s;
    
      xpressive::sregex pattern_section_weight = (*xpressive::_s
						  >> '[' >> (xpressive::s1= "weight-" >> +(~xpressive::_s)) >> ']'
						  >> *xpressive::_s);
    
      xpressive::sregex pattern_section = (*xpressive::_s
					   >> '[' >> (xpressive::s1= +(~xpressive::_s)) >> ']'
					   >> *xpressive::_s);
    
      std::string line;
      bool consumed = std::getline(is, line);
      while (consumed) {
	xpressive::smatch what;
      
	if (! xpressive::regex_match(line, what, pattern_section_weight)) {
	  os << line << '\n';
	  consumed = std::getline(is, line);
	  continue;
	}
      
	// this is the section for "weight"
	os << line << '\n';
      
	const std::string section = what[1];
	feature_weight_set_type::const_iterator witer = feature_weights.find(section);
      
	static const weight_set_type weights_empty;
      
	const weight_set_type& weights = (witer != feature_weights.end() ? witer->second : weights_empty);
      
	features_consumed.insert(section);
      
	int pos = 0;
	do {
	  consumed = std::getline(is, line);
	  if (! consumed || xpressive::regex_match(line, what, pattern_section)) break;
	
	  if (xpressive::regex_match(line, what, pattern_comment) || xpressive::regex_match(line, what, pattern_empty))
	    os << line << '\n';
	  else {
	    if (pos < static_cast<int>(weights.size()))
	      os << weights[pos] << '\n';
	    else
	      os << 0.0 << '\n';
	    ++ pos;
	  }
	} while (consumed);
      }
    
      // write all additional parameters...
    
      feature_weight_set_type::const_iterator fiter_end = feature_weights.end();
      for (feature_weight_set_type::const_iterator fiter = feature_weights.begin(); fiter != fiter_end; ++ fiter)
	if (features_consumed.find(fiter->first) == features_consumed.end()) {
	  os << '\n' << '[' << fiter->first << ']' << '\n';
	  for (size_t pos = 0; pos != fiter->second.size(); ++ pos)
	    os << fiter->second[pos] << '\n';
	  os << '\n';
	}
    }
  }
  catch (const std::exception& err) {
    std::cerr << "error: " << err.what() << std::endl;
    return 1;
  }
  return 0;
}


void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",     po::value<path_type>(&input_file)->default_value(input_file),   "input file")
    ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output file")
    ("weights",   po::value<path_type>(&weights_file),                            "weights file")

    ("bias-features", po::value<feature_set_type>(&features_biased)->multitoken(),     "biased features")
    ("bias-weight",   po::value<double>(&weight_biased)->default_value(weight_biased), "bias weight")
    
    ("help", "help message");
  
  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }
}

