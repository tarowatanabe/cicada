
#include <iostream>
#include <vector>
#include <string>
#include <stdexcept>
#include <set>

#include "cicada_impl.hpp"

#include "utils/sgi_hash_map.hpp"
#include "utils/program_options.hpp"

#include <boost/program_options.hpp>

#include <google/dense_hash_set>

struct TER
{
  struct COSTS
  {
    static const double insertion;
    static const double deletion;
    static const double substitution;
    static const double shift;
  };
  
  struct TRANSITION
  {
    enum transition_type {
      epsilon,
      match,
      substitution,
      insertion,
      deletion,
    };
  };
  
  struct Operation
  {
    int pos;
    TRANSITION::transition_type op;
    
    Operation() : pos(-1), op(TRANSITION::epsilon) {}
    Operation(const TRANSITION::transition_type __op)
      : pos(-1), op(__op) {}
    Operation(const int __pos, const TRANSITION::transition_type __op)
      : pos(__pos), op(__op) {}
  };
  
  static const int max_shift_size = 10;
  static const int max_shift_dist = 50;
};

const double TER::COSTS::insertion    = 1.0;
const double TER::COSTS::deletion     = 1.0;
const double TER::COSTS::substitution = 1.0;
const double TER::COSTS::shift        = 1.0;


struct MinimumEditDistance
{
  typedef TER::Operation operation_type;
  typedef std::vector<operation_type, std::allocator<operation_type> > operation_set_type;


  typedef utils::vector2<operation_type, std::allocator<operation_type> > matrix_operation_type;
  typedef utils::vector2<double, std::allocator<double> > matrix_cost_type;

  double operator()(const sentence_type& hyp, const lattice_type& ref, operation_set_type& operations)
  {
    // vocab_type::EPSILON is treated differently...
    
    matrix_cost_type      costs(hyp.size() + 1, ref.size() + 1, 0.0);
    matrix_operation_type ops(hyp.size() + 1, ref.size() + 1);
    
    for (int i = 0; i <= hyp.size(); ++ i)
      costs(i, 0) = i * TER::COSTS::insertion;
    for (int j = 0; j <= ref.size(); ++ j)
      costs(0, j) = j * TER::COSTS::deletion;
    
    for (int i = 1; i <= hyp.size(); ++ i)
      for (int j = 1; j <= ref.size(); ++ j) {
	double&         cur_cost = costs(i, j);
	operation_type& cur_op   = ops(i, j);
	
	lattice_type::arc_set_type::const_iterator aiter_begin   = ref[j - 1].begin();
	lattice_type::arc_set_type::const_iterator aiter_end     = ref[j - 1].end();
	lattice_type::arc_set_type::const_iterator aiter_match   = aiter_end;
	lattice_type::arc_set_type::const_iterator aiter_epsilon = aiter_end;
	
	for (lattice_type::arc_set_type::const_iterator aiter = aiter_begin; aiter != aiter_end; ++ aiter) {
	  if (aiter->label == hyp[i - 1])
	    aiter_match = aiter;
	  else if (aiter->label == vocab_type::EPSILON)
	    aiter_epsilon = aiter;
	}
	
	if (aiter_match != aiter_end) {
	  cur_cost = costs(i - 1, j - 1);
	  cur_op   = operation_type(aiter_match - aiter_begin, TER::TRANSITION::match);
	} else {
	  cur_cost = costs(i - 1, j - 1) + TER::COSTS::substitution;
	  cur_op   = operation_type(aiter_match - aiter_begin, TER::TRANSITION::substitution);
	}
	
	if (aiter_epsilon != aiter_end) {
	  const double eps = costs(i, j - 1);
	  if (cur_cost > eps) {
	    cur_cost = eps;
	    cur_op   = operation_type(aiter_epsilon - aiter_begin, TER::TRANSITION::epsilon);
	  }
	}
	
	const double ins = costs(i - 1, j) + TER::COSTS::insertion;
	if (cur_cost > ins) {
	  cur_cost = ins;
	  cur_op   = operation_type(TER::TRANSITION::insertion);
	}
	
	const double del = costs(i, j - 1) + TER::COSTS::deletion;
	if (cur_cost > del) {
	  cur_cost = del;
	  cur_op   = operation_type(TER::TRANSITION::deletion);
	}
      }
    
    
    operations.clear();
    int i = hyp.size();
    int j = ref.size();
    while (i > 0 || j > 0) {
      if (j == 0) {
	-- i;
	operations.push_back(operation_type(TER::TRANSITION::insertion));
      } else if (i == 0) {
	-- j;
	operations.push_back(operation_type(TER::TRANSITION::deletion));
      } else {
	const operation_type& op = ops(i, j);
	operations.push_back(op);
	
	switch (op.op) {
	case TER::TRANSITION::substitution:
	case TER::TRANSITION::match:        -- i; -- j; break;
	case TER::TRANSITION::insertion:    -- i; break;
	case TER::TRANSITION::epsilon:
	case TER::TRANSITION::deletion:     -- j; break;
	}
      }
    }
    
    return costs(hyp.size(), ref.size());
  }
};

struct TranslationErrorRate : public TER
{
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  
  typedef sentence_type ngram_type;

  typedef std::pair<int, int> cn_index_type;
  typedef std::set<cn_index_type, std::less<cn_index_type>, std::allocator<cn_index_type> > index_set_type;
  
#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<ngram_type, index_set_type, boost::hash<ngram_type>, std::equal_to<ngram_type>,
				  std::allocator<std::pair<const ngram_type, index_set_type> > > ngram_index_map_type;
#else
  typedef sgi::hash_map<ngram_type, index_set_type, boost::hash<ngram_type>, std::equal_to<ngram_type>,
			std::allocator<std::pair<const ngram_type, index_set_type> > > ngram_index_map_type;
#endif

  typedef google::dense_hash_set<word_type, boost::hash<word_type>, std::equal_to<word_type> > word_set_type;

  typedef MinimumEditDistance minimum_edit_distance_type;
  typedef minimum_edit_distance_type::operation_type     operation_type;
  typedef minimum_edit_distance_type::operation_set_type operation_set_type;
  typedef operation_set_type path_type;

  typedef std::vector<bool, std::allocator<bool> > error_set_type;
  typedef std::vector<cn_index_type, std::allocator<cn_index_type> > align_set_type;
  

  struct Score
  {
    double score;
    int match;
    int insertion;
    int deletion;
    int substitution;
    int shift;
	
    Score() : score(), match(0), insertion(0), deletion(0), substitution(0), shift(0) {}
  };
  
  struct Shift
  {
    Shift() : begin(0), end(0), moveto(0), reloc(0) {}
    Shift(const int __begin, const int __end, const int __moveto, const int __reloc)
      : begin(__begin), end(__end), moveto(__moveto), reloc(__reloc) {}
	
    int begin;
    int end;
    int moveto;
    int reloc;
  };
  
  typedef Score value_type;
  typedef Shift shift_type;
  
  typedef std::vector<shift_type, std::allocator<shift_type> > shift_set_type;
  typedef std::vector<shift_set_type, std::allocator<shift_set_type> > shift_matrix_type;
  
  void find_alignment_error(const path_type& path,
			    error_set_type& herr,
			    error_set_type& rerr,
			    align_set_type& ralign) const
  {
    int hpos = -1;
    int rpos = -1;
    
    //std::cerr << "edit distance: ";
    for (int i = 0; i < path.size(); ++ i) {
      switch (path[i].op) {
      case TRANSITION::match:
	//std::cerr << " M";
	++ hpos;
	++ rpos;
	herr.push_back(false);
	rerr.push_back(false);
	ralign.push_back(cn_index_type(hpos, path[i].pos));
	break;
      case TRANSITION::substitution:
	//std::cerr << " S";
	++ hpos;
	++ rpos;
	herr.push_back(true);
	rerr.push_back(true);
	ralign.push_back(cn_index_type(hpos, path[i].pos));
	break;
      case TRANSITION::insertion:
	//std::cerr << " I";
	++ hpos;
	herr.push_back(true);
	break;
      case TRANSITION::deletion:
	//std::cerr << " D";
	++ rpos;
	rerr.push_back(true);
	ralign.push_back(cn_index_type(hpos, -1));
	break;
      case TRANSITION::epsilon:
	//std::cerr << " E";
	++ rpos;
	rerr.push_back(true);
	ralign.push_back(cn_index_type(hpos, path[i].pos));
	break;
      }
    }
    //std::cerr << std::endl;
  }
  
  bool calculate_best_shift(const lattice_type& ref,
			    const sentence_type& hyp,
			    const path_type& path,
			    const double cost,
			    sentence_type& hyp_best,
			    path_type& path_best,
			    double& cost_best,
			    const ngram_index_map_type& ngram_index)
  {
    
    error_set_type herr;
    error_set_type rerr;
    align_set_type ralign;
    
    herr.reserve(path.size());
    rerr.reserve(path.size());
    ralign.reserve(path.size());
    
    find_alignment_error(path, herr, rerr, ralign);
    

    shift_matrix_type shifts(max_shift_size + 1);
    
    
    
  }

  void gather_all_possible_shifts(const lattice_type& ref,
				  const sentence_type& hyp,
				  const align_set_type& ralign,
				  const error_set_type& herr,
				  const error_set_type& rerr,
				  const ngram_index_map_type& ngram_index,
				  shift_matrix_type& shifts) const
  {
    ngram_type ngram;
	
    for (int start = 0; start != hyp.size(); ++ start) {
      ngram.clear();
      //ngram.push_back(ref[start].);
      
      ngram_index_map_type::const_iterator niter = ngram_index.find(ngram);
      if (niter == ngram_index.end()) continue;
      
      bool found = false;
      index_set_type::const_iterator iiter_end = niter->second.end();
      for (index_set_type::const_iterator iiter = niter->second.begin(); iiter != iiter_end && ! found; ++ iiter) {
	const int moveto = iiter->first;
	found = (start != ralign[moveto].first && (ralign[moveto].first - start <= max_shift_dist) && (start - ralign[moveto].first - 1 <= max_shift_dist));
      }
      
      if (! found) continue;
      
      
    }
  }

  void build_ngram_matches(const lattice_type& ref, const sentence_type& hyp, ngram_index_map_type& ngram_index) const
  {
    typedef std::pair<ngram_type, cn_index_type> intersected_type;
    typedef std::vector<intersected_type, std::allocator<intersected_type> > intersected_set_type;

    // we need to expand lattice....
      ngram_index.clear();
    
    word_set_type words_intersect;
    words_intersect.set_empty_key(word_type());
    words_intersect.insert(hyp.begin(), hyp.end());
    
    word_set_type::const_iterator iiter_end = words_intersect.end();
    
    ngram_type ngram;

    intersected_set_type intersected;
    intersected_set_type intersected_next;
    
    for (int start = 0; start != ref.size(); ++ start) {
      ngram.clear();
      const int max_length = utils::bithack::min(max_shift_size, static_cast<int>(ref.size() - start));
      
      // we will expand in tree structure...
      
      for (int length = 0; length != max_length; ++ length) {
	
	if (length == 0) {
	  intersected.clear();
	  for (int pos = 0; pos < ref[start].size(); ++ pos)
	    if (ref[start][pos].label == vocab_type::EPSILON)
	      intersected.push_back(intersected_type(ngram_type(), cn_index_type(start, pos)));
	    else if (words_intersect.find(ref[start][pos].label) != iiter_end) {
	      intersected.push_back(intersected_type(ngram_type(1, ref[start][pos].label), cn_index_type(start, pos)));
	      ngram_index[intersected.back().first].insert(intersected.back().second);
	    }
	  
	} else {
	  if (intersected.empty()) break;
	  
	  intersected_next.clear();
	  for (int pos = 0; pos < ref[start + length].size(); ++ pos) {
	    if (ref[start][pos].label == vocab_type::EPSILON)
	      intersected_next.insert(intersected_next.end(), intersected.begin(), intersected.end());
	    else if (words_intersect.find(ref[start + length][pos].label) != iiter_end) {
	      intersected_set_type::const_iterator iiter_end = intersected.end();
	      for (intersected_set_type::const_iterator iiter = intersected.begin(); iiter != iiter_end; ++ iiter) {
		ngram_type ngram(iiter->first);
		ngram.push_back(ref[start + length][pos].label);
		
		intersected.push_back(intersected_type(ngram, iiter->second));
		ngram_index[intersected.back().first].insert(intersected.back().second);
	      }
	    }
	  }
	  
	  intersected.swap(intersected_next);
	}
      }
    }
  }
};



path_type input_file = "-";
path_type output_file = "-";

std::string confidence;
std::string count;
double count_weight = 1.0;

bool individual = false;

// input mode... use of one-line lattice input or sentence input?
void options(int argc, char** argv);

int main(int argc, char ** argv)
{
  try {
    options(argc, argv);
    
    utils::compress_istream is(input_file, 1024 * 1024);
    utils::compress_ostream os(output_file, 1024 * 1024);
    
    lattice_type merged;
    lattice_type lattice;

    cicada::Feature feature_confidence(confidence);
    cicada::Feature feature_count(count);
    
    int rank = 1;
    std::string line;
    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();
      
      if (! lattice.assign(iter, end))
	throw std::runtime_error("invalid hypergraph format");
      
      if (! feature_confidence.empty()) {
	const double conf = 1.0 / (1.0 + rank);
	
	lattice_type::iterator liter_end = lattice.end();
	for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type& arcs = *liter;
	  
	  for (lattice_type::arc_set_type::iterator aiter = arcs.begin(); aiter != arcs.end(); ++ aiter)
	    aiter->features[feature_confidence] = conf;
	}
      } 
      
      if (! feature_count.empty()) {
	lattice_type::iterator liter_end = lattice.end();
	for (lattice_type::iterator liter = lattice.begin(); liter != liter_end; ++ liter) {
	  lattice_type::arc_set_type& arcs = *liter;
	  
	  for (lattice_type::arc_set_type::iterator aiter = arcs.begin(); aiter != arcs.end(); ++ aiter)
	    aiter->features[feature_count] = count_weight;
	}
      }
      
      if (individual)
	os << lattice << '\n';
      else {
	// perform merging...
	
      }
	
      
      ++ rank;
    }
    
    if (! individual)
      os << merged << '\n';
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
    ("input",  po::value<path_type>(&input_file)->default_value("-"),   "input lattices")
    ("output", po::value<path_type>(&output_file)->default_value("-"),  "output merged lattice")
    
    ("confidence",   po::value<std::string>(&confidence),    "add confidence weight feature name")
    ("count",        po::value<std::string>(&count),         "add count weight feature name")
    ("count-weight", po::value<double>(&count_weight),       "count weight")
    
    ("individual", po::bool_switch(&individual), "no merging")
    
    ("help", "help message");

  po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);
  
  if (variables.count("help")) {
    std::cout << argv[0] << " [options]\n"
	      << desc << std::endl;
    exit(0);
  }
}
