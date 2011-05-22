//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <set>

#include <boost/spirit/include/qi.hpp>

#include "cicada_impl.hpp"
#include "cicada/remove_epsilon.hpp"
#include "cicada/unite.hpp"

#include "utils/sgi_hash_map.hpp"
#include "utils/program_options.hpp"

#include <boost/program_options.hpp>

#include <google/dense_hash_set>

typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
typedef std::vector<feature_type, std::allocator<feature_type> > feature_list_type;

path_set_type input_files;
path_type output_file = "-";
path_type confidence_feature_file;
path_type count_feature_file;

std::string confidence;
std::string count;
double count_weight = 1.0;

bool remove_epsilon = false;
bool multiple_mode = false;

int debug = 0;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);

    feature_list_type features_confidence;
    feature_list_type features_count;

    if (! confidence_feature_file.empty()) {
      if (confidence_feature_file != "-" && ! boost::filesystem::exists(confidence_feature_file))
	throw std::runtime_error("no confidence feature file? " + confidence_feature_file.string());
      
      utils::compress_istream is(confidence_feature_file);
      std::string feature;
      while (is >> feature)
	features_confidence.push_back(feature);
    }
    
    if (! count_feature_file.empty()) {
      if (count_feature_file != "-" && ! boost::filesystem::exists(count_feature_file))
	throw std::runtime_error("no count feature file? " + count_feature_file.string());
      
      utils::compress_istream is(count_feature_file);
      std::string feature;
      while (is >> feature)
	features_count.push_back(feature);
    }

    const bool flush_output = (output_file == "-"
                               || (boost::filesystem::exists(output_file)
                                   && ! boost::filesystem::is_regular_file(output_file)));
    
    cicada::Feature feature_confidence(confidence);
    cicada::Feature feature_count(count);
    
    lattice_type merged;
    lattice_type lattice;

    if (input_files.empty())
      input_files.push_back("-");

    if (multiple_mode) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      // lattice ||| lattice ||| lattice
      
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
      
      for (path_set_type::const_iterator iter = input_files.begin(); iter != input_files.end(); ++ iter) {
	utils::compress_istream is(input_files.front(), 1024 * 1024);
	std::string line;
	
	while (std::getline(is, line)) {
	  int rank = 1;
	  int id = 0;
	  
	  merged.clear();
	  lattice.clear();
	  
	  std::string::const_iterator iter = line.begin();
	  std::string::const_iterator end = line.end();
	  
	  for (/**/; iter != end; ++ id, ++ rank) {
	    if (id != 0)
	      if (! qi::phrase_parse(iter, end, "|||", standard::space))
		break;
	    
	    if (! lattice.assign(iter, end))
	      throw std::runtime_error("invalid lattice format");
	    
	    if (lattice.empty()) continue;
	    
	    const double conf = 1.0 / (1.0 + rank);
	    
	    feature_set_type features;
	    
	    if (! features_confidence.empty()) {
	      if (id >= static_cast<int>(features_confidence.size()))
		throw std::runtime_error("# of confidence features do not match");
	      features[features_confidence[id]] = conf;
	    }
	    if (! features_count.empty()) {
	      if (id >= static_cast<int>(features_count.size()))
		throw std::runtime_error("# of count features do not match");
	      features[features_count[id]] = count_weight;
	    }
	    if (! feature_confidence.empty())
	      features[feature_confidence] = conf;
	    if (! feature_count.empty())
	      features[feature_count] = count_weight;
	    
	    if (! features.empty()) {
	      lattice_type::iterator liter_end = lattice.end();
	      for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
		lattice_type::arc_set_type::iterator aiter_end = liter->end();
		for (lattice_type::arc_set_type::iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
		  aiter->features += features;
	      }
	    }

	    merged.unite(lattice);
	  }

	  if (remove_epsilon)
	    cicada::remove_epsilon(merged);
	  
	  os << merged << '\n';
	}
      }
    } else if (input_files.size() == 1) {
      utils::compress_istream is(input_files.front(), 1024 * 1024);
      std::string line;
    
      int rank = 1;
      int id = 0;
      for (/**/; std::getline(is, line); ++ id, ++ rank) {
	std::string::const_iterator iter = line.begin();
	std::string::const_iterator end = line.end();
      
	if (! lattice.assign(iter, end))
	  throw std::runtime_error("invalid lattice format");

	if (lattice.empty()) continue;
      
	const double conf = 1.0 / (1.0 + rank);
      
	feature_set_type features;
	
	if (! features_confidence.empty()) {
	  if (id >= static_cast<int>(features_confidence.size()))
	    throw std::runtime_error("# of confidence features do not match");
	  features[features_confidence[id]] = conf;
	}
	if (! features_count.empty()) {
	  if (id >= static_cast<int>(features_count.size()))
	    throw std::runtime_error("# of count features do not match");
	  features[features_count[id]] = count_weight;
	}
	if (! feature_confidence.empty())
	  features[feature_confidence] = conf;
	if (! feature_count.empty())
	  features[feature_count] = count_weight;
      
	if (! features.empty()) {
	  lattice_type::iterator liter_end = lattice.end();
	  for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	    lattice_type::arc_set_type::iterator aiter_end = liter->end();
	    for (lattice_type::arc_set_type::iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
	      aiter->features += features;
	  }
	}
	
	merged.unite(lattice);
      }
    
      if (remove_epsilon)
	cicada::remove_epsilon(merged);
    
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
    
      os << merged << '\n';
    } else {
      // we will handle multiple files!
      
      if (! features_confidence.empty())
	if (input_files.size() != features_confidence.size())
	  throw std::runtime_error("input file do not match with # of confidence feature");

      if (! features_count.empty())
	if (input_files.size() != features_count.size())
	  throw std::runtime_error("input file do not match with # of count feature");
      
      typedef std::vector<std::istream*, std::allocator<std::istream*> > istream_set_type;
      
      istream_set_type istreams(input_files.size());
      for (size_t i = 0; i != input_files.size(); ++ i)
	istreams[i] = new utils::compress_istream(input_files[i], 1024 * 1024);
      
      utils::compress_ostream os(output_file, 1024 * 1024 * (! flush_output));
      
      std::string line;
      
      for (;;) {
	int rank = 1;
	
	merged.clear();
	lattice.clear();
	
	size_t num_failed = 0;
	for (size_t id = 0; id != istreams.size(); ++ id, ++ rank) {
	  if (std::getline(*istreams[id], line)) {
	    std::string::const_iterator iter = line.begin();
	    std::string::const_iterator end = line.end();
	    
	    if (! lattice.assign(iter, end))
	      throw std::runtime_error("invalid lattice format");
	    
	    if (lattice.empty()) continue;
	    
	    const double conf = 1.0 / (1.0 + rank);
	    
	    feature_set_type features;
	    
	    if (! features_confidence.empty()) {
	      if (id >= static_cast<int>(features_confidence.size()))
		throw std::runtime_error("# of confidence features do not match");
	      features[features_confidence[id]] = conf;
	    }
	    if (! features_count.empty()) {
	      if (id >= static_cast<int>(features_count.size()))
		throw std::runtime_error("# of count features do not match");
	      features[features_count[id]] = count_weight;
	    }
	    if (! feature_confidence.empty())
	      features[feature_confidence] = conf;
	    if (! feature_count.empty())
	      features[feature_count] = count_weight;
						       
	    if (! features.empty()) {
	      lattice_type::iterator liter_end = lattice.end();
	      for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
		lattice_type::arc_set_type::iterator aiter_end = liter->end();
		for (lattice_type::arc_set_type::iterator aiter = liter->begin(); aiter != aiter_end; ++ aiter)
		  aiter->features += features;
	      }
	    }
	    
	    merged.unite(lattice);
	  } else
	    ++ num_failed;
	}
	
	if (num_failed) {
	  if (num_failed != istreams.size())
	    throw std::runtime_error("# of lines do not match");
	  break;
	}
	
	if (remove_epsilon)
	  cicada::remove_epsilon(merged);
	
	os << merged << '\n';
      }
      
      for (size_t i = 0; i != istreams.size(); ++ i)
	delete istreams[i];
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

  po::variables_map variables;
  
  po::options_description desc("options");
  desc.add_options()
    ("input",  po::value<path_set_type>(&input_files)->multitoken(),   "input lattices")
    ("output", po::value<path_type>(&output_file)->default_value("-"),  "output merged lattice")
    
    ("confidence-feature-file", po::value<path_type>(&confidence_feature_file), "confidence feature file")
    ("count-feature-file",      po::value<path_type>(&count_feature_file),      "count feature file")

    ("confidence",   po::value<std::string>(&confidence),    "add confidence weight feature name")
    ("count",        po::value<std::string>(&count),         "add count weight feature name")
    ("count-weight", po::value<double>(&count_weight),       "count weight")
    
    ("remove-epsilon", po::bool_switch(&remove_epsilon), "remvoe epsilon")
    ("multiple", po::bool_switch(&multiple_mode), "multiple forest in one line")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");

  po::positional_options_description pos;
  pos.add("input", -1); // all the files

  po::command_line_parser parser(argc, argv);
  parser.style(po::command_line_style::unix_style & (~po::command_line_style::allow_guessing));
  parser.options(desc);
  parser.positional(pos);
  
  po::store(parser.run(), variables);
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
