//
// weights management...
//
// read weights and
// --average to perform averaging
// --sum to perform summing
// --sort to perform sorting
// --sort-abs to performn sorting by absolute values
//
// --normalize-l1 to perform normalization via l1
// --normalize-l2 to perform normalization via l2
//

#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <algorithm>
#include <iterator>

#include "cicada/weight_vector.hpp"
#include "cicada/dot_product.hpp"
#include "cicada/parameter.hpp"

#include <boost/spirit/include/qi.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include "utils/lexical_cast.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/piece.hpp"

typedef boost::filesystem::path path_type;
typedef std::vector<std::string, std::allocator<std::string> > path_set_type;

typedef cicada::Parameter parameter_type;
typedef cicada::WeightVector<double> weight_set_type;

path_set_type input_files;
path_type     output_file = "-";

path_type bound_lower_file;
path_type bound_upper_file;

bool norm_mode = false;
bool distance_mode = false;
bool average_mode = false;
bool sum_mode = false;
bool sort_mode = false;
bool sort_abs_mode = false;
bool normalize_l1 = false;
bool normalize_l2 = false;
bool distance_l1 = false;
bool distance_l2 = false;
bool norm_l1 = false;
bool norm_l2 = false;

int debug = 0;

template <typename Tp>
struct greater_second
{
  bool operator()(const Tp& x, const Tp& y) const
  {
    return x.second > y.second;
  }
};

template <typename Tp>
struct greater_fabs_second
{
  bool operator()(const Tp& x, const Tp& y) const
  {
    return std::fabs(x.second) > std::fabs(y.second);
  }
};

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (int(sum_mode) + average_mode + distance_mode + norm_mode > 1)
      throw std::runtime_error("You cannnot perform both sum, average and dist");
    if (normalize_l1 && normalize_l2)
      throw std::runtime_error("You cannnot perform both normalize-l1 and normalize-l2");
    if (sort_mode && sort_abs_mode)
      throw std::runtime_error("You cannnot perform both sort and sort-abs");

    if (int(sum_mode) + average_mode + distance_mode + norm_mode == 0)
      sum_mode = true;
    
    weight_set_type weights;
    size_t          num = 0;
    
    if (input_files.empty())
      input_files.push_back("-");
    
    if (norm_mode) {
      if (input_files.size() > 2)
	throw std::runtime_error("we need a single weights");
    
      if (int(norm_l1) + norm_l2 > 1)
	throw std::runtime_error("you cannot compute both l1 and l2 norm");
      if (int(norm_l1) + norm_l2 == 0)
	norm_l2 = true;

      weight_set_type weights;
      
      utils::compress_istream is(input_files.front(), 1024 * 1024);
      
      is >> weights;
      
      if (norm_l1) {
	double norm = 0.0;
	for (size_t i = 0; i != weights.size(); ++ i)
	  norm += std::fabs(weights[i]);
	
	utils::compress_ostream os(output_file);
	os << norm << '\n';
      } else {
	double norm = 0.0;
	for (size_t i = 0; i != weights.size(); ++ i)
	  norm += weights[i] * weights[i];
	
	utils::compress_ostream os(output_file);
	os << norm << '\n';
      }
      
    } else if (distance_mode) {
      if (input_files.size() != 2)
	throw std::runtime_error("we need at least two weights for comparison");

      if (int(distance_l1) + distance_l2 > 1)
	throw std::runtime_error("you cannot compute both l1 and l2 distance");
      if (int(distance_l1) + distance_l2 == 0)
	distance_l2 = true;
      
      weight_set_type weights1;
      weight_set_type weights2;
      
      utils::compress_istream is1(input_files.front(), 1024 * 1024);
      utils::compress_istream is2(input_files.back(),  1024 * 1024);
      
      is1 >> weights1;
      is2 >> weights2;

      weights1.allocate();
      weights2.allocate();

      if (distance_l1) {
	double diff = 0.0;
	for (size_t i = 0; i != weights1.size(); ++ i)
	  diff += std::fabs(weights1[i] - weights2[i]);
	
	utils::compress_ostream os(output_file);
	os << diff << '\n';
      } else {
	double diff = 0.0;
	for (size_t i = 0; i != weights1.size(); ++ i)
	  diff += (weights1[i] - weights2[i]) * (weights1[i] - weights2[i]);
	if (diff != 0.0)
	  diff = std::sqrt(diff);
	
	utils::compress_ostream os(output_file);
	os << diff << '\n';
      }
      
      return 0;
    }
    
    weight_set_type weights_input;
    for (path_set_type::const_iterator fiter = input_files.begin(); fiter != input_files.end(); ++ fiter) {
      const parameter_type param(*fiter);

      if (param.name() != "-" && ! boost::filesystem::exists(param.name()))
	throw std::runtime_error("no file: " + param.name());

      double scale = 1.0;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "scale" || utils::ipiece(piter->first) == "weight")
	  scale = utils::lexical_cast<double>(piter->second);
	else
	  throw std::runtime_error("WARNING: unsupported parameter: " + piter->first + "=" + piter->second);
      }
      
      utils::compress_istream is(param.name(), 1024 * 1024);
      
      weights_input.clear();
      is >> weights_input;
      
      if (scale != 1.0)
	weights_input *= scale;
      
      if (weights.empty())
	weights.swap(weights_input);
      else
	weights += weights_input;
      
      ++ num;
    }
    
    if (average_mode && num)
      weights *= (1.0 / num);
    
    if (normalize_l1) {
      double sum = 0.0;
      weight_set_type::iterator witer_end = weights.end();
      for (weight_set_type::iterator witer = weights.begin(); witer != witer_end; ++ witer)
	sum += std::fabs(*witer);
      
      if (sum != 0.0)
	weights *= 1.0 / std::sqrt(sum);
      
    } else if (normalize_l2) {
      const double sum = cicada::dot_product(weights, weights);
      
      if (sum != 0.0)
	weights *= 1.0 / std::sqrt(sum);
    }
    
    if (! bound_lower_file.empty()) {
      typedef boost::spirit::istream_iterator iter_type;

      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      if (bound_lower_file != "-" && ! boost::filesystem::exists(bound_lower_file))
	throw std::runtime_error("no file? " + bound_lower_file.string());
      
      utils::compress_istream is(bound_lower_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iter_type iter(is);
      iter_type iter_end;
      
      std::string feature;
      double      value;
      
      while (iter != iter_end) {
	feature.clear();
	value = 0.0;
	
	if (! qi::phrase_parse(iter, iter_end,
			       qi::lexeme[+(standard::char_ - standard::space)] >> qi::double_ >> (qi::eol | qi::eoi),
			       standard::blank, feature, value))
	  if (iter != iter_end)
	    throw std::runtime_error("weights file parsing failed");
	
	weights[feature] = std::max(weights[feature], value);
      }
    }

    if (! bound_upper_file.empty()) {
      typedef boost::spirit::istream_iterator iter_type;

      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      if (bound_upper_file != "-" && ! boost::filesystem::exists(bound_upper_file))
	throw std::runtime_error("no file? " + bound_upper_file.string());
      
      utils::compress_istream is(bound_upper_file, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iter_type iter(is);
      iter_type iter_end;
      
      std::string feature;
      double      value;
      
      while (iter != iter_end) {
	feature.clear();
	value = 0.0;
	
	if (! qi::phrase_parse(iter, iter_end,
			       qi::lexeme[+(standard::char_ - standard::space)] >> qi::double_ >> (qi::eol | qi::eoi),
			       standard::blank, feature, value))
	  if (iter != iter_end)
	    throw std::runtime_error("weights file parsing failed");
	
	weights[feature] = std::min(weights[feature], value);
      }
    }

    
    if (sort_mode || sort_abs_mode) {
      typedef std::pair<weight_set_type::feature_type, double> value_type;
      typedef std::vector<value_type, std::allocator<value_type> > value_set_type;
     
      value_set_type values;
      
      weight_set_type::feature_type::id_type id = 0;
      weight_set_type::iterator witer_end = weights.end();
      for (weight_set_type::iterator witer = weights.begin(); witer != witer_end; ++ witer, ++ id)
	if (*witer != 0.0)
	  values.push_back(value_type(weight_set_type::feature_type(id), *witer));

      if (sort_mode)
	std::sort(values.begin(), values.end(), greater_second<value_type>());
      else
	std::sort(values.begin(), values.end(), greater_fabs_second<value_type>());
      
      utils::compress_ostream os(output_file, 1024 * 1024);
      os.precision(20);

      value_set_type::const_iterator fiter_end = values.end();
      for (value_set_type::const_iterator fiter = values.begin(); fiter != fiter_end; ++ fiter)
	if (! fiter->first.empty())
	  os << fiter->first << ' ' << fiter->second << '\n';
      
    } else {
      utils::compress_ostream os(output_file, 1024 * 1024);
      os.precision(20);
      
      os << weights;
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
    ("input",     po::value<path_set_type>(&input_files)->multitoken(),           "input files")
    ("output",    po::value<path_type>(&output_file)->default_value(output_file), "output")

    ("bound-lower", po::value<path_type>(&bound_lower_file), "lower bound file")
    ("bound-upper", po::value<path_type>(&bound_upper_file), "upper bound file")

    ("distance",     po::bool_switch(&distance_mode), "compute distance")
    ("norm",         po::bool_switch(&norm_mode),     "compute norm")
    ("average",      po::bool_switch(&average_mode),  "average weights")
    ("sum",          po::bool_switch(&sum_mode),      "sum weights")
    
    ("sort",         po::bool_switch(&sort_mode),     "sort weights")
    ("sort-abs",     po::bool_switch(&sort_abs_mode), "sort weights wrt absolute value")
    ("normalize-l1", po::bool_switch(&normalize_l1),  "weight normalization by L1")
    ("normalize-l2", po::bool_switch(&normalize_l2),  "weight normalization by L2")
    ("distance-l1",  po::bool_switch(&distance_l1),  "weight distance by L1")
    ("distance-l2",  po::bool_switch(&distance_l2),  "weight distance by L2")
    ("norm-l1",      po::bool_switch(&norm_l1),      "weight norm by L1")
    ("norm-l2",      po::bool_switch(&norm_l2),      "weight norm by L2")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    
    ("help", "help message");
  
  po::positional_options_description pos;
  pos.add("input", -1); // all the files

  po::command_line_parser parser(argc, argv);
  parser.style(po::command_line_style::unix_style & (~po::command_line_style::allow_guessing));
  parser.options(desc);
  parser.positional(pos);
  
  po::variables_map vm;
  po::store(parser.run(), vm);
  po::notify(vm);
  
  if (vm.count("help")) {
    std::cout << argv[0] << " [options]" << '\n' << desc << '\n';
    exit(0);
  }

}
