#include <stdexcept>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <iostream>
#include <numeric>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/tokenizer.hpp>

#include <utils/space_separator.hpp>
#include <utils/base64.hpp>

#include <utils/arc_list.hpp>

#include <utils/lockfree_list_queue.hpp>

inline
path_type add_suffix(const path_type& path, const std::string& suffix)
{
  bool has_suffix_gz  = false;
  bool has_suffix_bz2 = false;
  
  path_type path_added = path;
  
  if (path.extension() == ".gz") {
    path_added = path.parent_path() / path.stem();
    has_suffix_gz = true;
  } else if (path.extension() == ".bz2") {
    path_added = path.parent_path() / path.stem();
    has_suffix_bz2 = true;
  }
  
  path_added = path_added.file_string() + suffix;
  
  if (has_suffix_gz)
    path_added = path_added.file_string() + ".gz";
  else if (has_suffix_bz2)
    path_added = path_added.file_string() + ".bz2";
  
  return path_added;
}

inline
void encode_support_vectors(std::string& data,
			    const size_t& id,
			    const int& source_length,
			    const sentence_type& viterbi,
			    const hypergraph_type& oracle,
			    const hypergraph_type& violated,
			    const sentence_set_type& targets)
{
  data.clear();
  
  boost::iostreams::filtering_ostream os;
  os.push(boost::iostreams::back_inserter(data));
  
  os << id << ' ' << source_length << " |||";
  os << ' ' << viterbi << " |||";
  os << ' ' << oracle << " |||";
  os << ' ' << violated;
  
  sentence_set_type::const_iterator titer_end = targets.end();
  for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
    os << " ||| " << *titer;
  
  os.pop();
}

inline
void decode_support_vectors(const std::string& data,
			    size_t& id,
			    int& source_length,
			    sentence_type& viterbi,
			    hypergraph_type& oracle,
			    hypergraph_type& violated,
			    sentence_set_type& targets)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  namespace phoenix = boost::phoenix;
  
  using qi::phrase_parse;
  using qi::_1;
  using qi::ulong_;
  using qi::int_;
  using standard::space;
  
  using phoenix::ref;
  
  std::string::const_iterator iter = data.begin();
  std::string::const_iterator end = data.end();
  
  if (! phrase_parse(iter, end, ulong_ [ref(id) = _1] >> int_ [ref(source_length) = _1] >> "|||", space))
    throw std::runtime_error("invalid id and source-length");
  
  if (! viterbi.assign(iter, end))
    throw std::runtime_error("invalid sentence");
  
  if (! phrase_parse(iter, end, "|||", space))
    throw std::runtime_error("expected separator");
  
  if (! oracle.assign(iter, end))
    throw std::runtime_error("invalid oracle hypergraph");
  
  if (! phrase_parse(iter, end, "|||", space))
    throw std::runtime_error("expected separator");
  
  if (! violated.assign(iter, end))
    throw std::runtime_error("invalid violated hypergraph");

  targets.clear();
  while (phrase_parse(iter, end, "|||", space)) {
    targets.push_back(sentence_type());
    if (! targets.back().assign(iter, end))
      throw std::runtime_error("invalid sentence format");
  }
  
  if (iter != end)
    throw std::runtime_error("still data remain?");
}


struct OptimizeCP
{
  // cutting-plane optimizer...
  
};

struct OptimizeMIRA
{
  OptimizeMIRA(const weight_set_type& __weights,
	       const double& __lambda,
	       const int __debug=0)
    : lambda(__lambda),
      C(1.0 / __lambda),
      weights(__weights),
      accumulated(),
      updated(1),
      debug(__debug) {}

  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    gradient_type;

  void finalize()
  {
    accumulated *= - 1.0 / updated;
    accumulated += weights;
    accumulated *= updated;
  }
  

  template <typename LabelSet, typename FeatureSet>
  struct h_matrix : public utils::hashmurmur<size_t>
  {
    h_matrix(const LabelSet& __labels, const FeatureSet& __features)
      : labels(__labels), features(__features) {}
    
    typedef utils::hashmurmur<size_t> hasher_type;

    double operator()(int i, int j) const
    {
      return labels[i] * labels[j] * features[i].dot(features[j]);
    }
    
    const LabelSet& labels;
    const FeatureSet& features;
  };
  
  template <typename LabelSet, typename MarginSet, typename FeatureSet>
  void operator()(const LabelSet&   labels,
		  const MarginSet&  margins,
		  const FeatureSet& features)
  {
    // QP solver used on OCAS
    h_matrix<LabelSet, FeatureSet>  H(labels, features);

    alpha.clear();
    gradient.clear();

    alpha.reserve(labels.size());
    gradient.reserve(labels.size());
    
    alpha.resize(labels.size(), 0.0);
    gradient.resize(labels.size(), 0.0);
    
    double alpha_neq = C;
    
    for (int i = 0; i < labels.size(); ++ i)
      gradient[i] = margins[i] - labels[i] * features[i].dot(weights);
    
    if (debug) {
      double obj_primal = 0.0;
      double obj_dual = 0.0;
      for (int i = 0; i < labels.size(); ++ i) {
	obj_primal += (std::max(gradient[i], 0.0) * C) / labels.size();
	obj_dual += gradient[i] * alpha[i];
      }
      
      std::cerr << "initial primal: " << obj_primal << " dual: " << obj_dual << std::endl; 
    }
    
    
    for (int iter = 0; iter != 1000; ++ iter) {
      
      // eq (23) and (26) to compute u and delta
      int u = -1;
      double max_obj = - std::numeric_limits<double>::infinity();
      double delta = 0.0;
      
      for (int i = 0; i < labels.size(); ++ i) {
	delta -= alpha[i] * gradient[i];
	
	if (gradient[i] > max_obj)  {
	  max_obj = gradient[i];
	  u = i;
	}
      }

      // if we relax the summation constraints, we select different u...
      // is this correct?
      //if (gradient[u] < 0.0)
      //  u = -1;
      //else
      delta += C * gradient[u];
      
      // tolerance
      if (delta <= 1e-4) break;
      
      // select v (26)
      if (u >= 0) {
	int v = -1;
	double max_improvement = - std::numeric_limits<double>::infinity();;
	double tau = 1.0;
	
	for (int i = 0; i < labels.size(); ++ i)
	  if (i != u && alpha[i] > 0.0) {
	    // compute (25)
	    const double numer = alpha[i] * (gradient[u] - gradient[i]);
	    const double denom = alpha[i] * alpha[i] * (H(u, u) - 2.0 * H(u, i) + H(i, i));
	    
	    if (denom > 0.0) {
	      const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
	      
	      if (improvement > max_improvement) {
		max_improvement = improvement;
		tau = std::max(0.0, std::min(1.0, numer / denom));
		v = i;
	      }
	    }
	  }
	
#if 0
	// check if virtual variable can be used for update...
	if (alpha_neq > 0.0) {
	  const double numer = alpha_neq * gradient[u];
	  const double denom = alpha_neq * alpha_neq * H(u, u);
	  
	  if (denom > 0.0) {
	    const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
	    
	    if (improvement > max_improvement) {
	      max_improvement = improvement;
	      tau = std::max(0.0, std::min(1.0, numer / denom));
	      v = -1;
	    }
	  }
	}
#endif
	
	if (v >= 0) {
	  // maximize objective, u and v
	  double update = alpha[v] * tau;
	  if (alpha[u] + update < 0.0)
	    update = - alpha[u];
	  
	  // clipping...
	  alpha[u] += update;
	  alpha[v] -= update;
	  
	  // update g...
	  for (int i = 0; i < labels.size(); ++ i)
	    gradient[i] += update * (H(i, v) - H(i, u));
	  
	} else if (alpha_neq > 0.0) {
	  double update = alpha_neq * tau;
	  if (alpha[u] + update < 0.0)
	    update = - alpha[u];
	  
	  alpha[u] += update;
	  alpha_neq -= update;
	  
	  // update g..
	  for (int i = 0; i < labels.size(); ++ i)
	    gradient[i] -= update * H(i, u);
	}
	
      } else {
	int v = -1;
	double max_improvement = - std::numeric_limits<double>::infinity();
	double tau = 1.0;
	
	for (int i = 0; i < labels.size(); ++ i) 
	  if (alpha[i] > 0.0) {
	    
	    const double numer = alpha[i] * gradient[i];
	    const double denom = alpha[i] * alpha[i] * H(i, i);
	    
	    if (denom > 0.0) {
	      const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
	      
	      if (improvement > max_improvement) {
		max_improvement = improvement;
		tau = std::max(0.0, std::min(1.0, numer / denom));
		v = i;
	      }
	    }
	  }
	
	// safe to check this...
	if (v >= 0) {
	  double update = alpha[v] * tau;
	  if (alpha_neq + update < 0.0)
	    update = - alpha_neq;
	  
	  alpha_neq += update;
	  alpha[v] -= update;
	  
	  for (int i = 0; i < labels.size(); ++ i)
	    gradient[i] += update * H(i, v);
	}
      }
    }
    
    bool perform_update = false;
    for (int i = 0; i < labels.size(); ++ i) 
      if (alpha[i] > 0.0) {
	typename FeatureSet::value_type::const_iterator fiter_end = features[i].end();
	for (typename FeatureSet::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	  weights[fiter->first] += alpha[i] * labels[i] * fiter->second;
	  accumulated[fiter->first] += alpha[i] * labels[i] * fiter->second * updated;
	}
	
	perform_update = true;
      }
    updated += perform_update;
    
    if (debug) {
      double obj_primal = 0.0;
      double obj_dual = 0.0;
      for (int i = 0; i < labels.size(); ++ i) {
	obj_primal += (std::max(gradient[i], 0.0) * C) / labels.size();
	obj_dual += gradient[i] * alpha[i];
      }
      
      std::cerr << "final primal: " << obj_primal << " dual: " << obj_dual << std::endl;
    }
  }

  alpha_type    alpha;
  gradient_type gradient;
  
  double lambda;
  double C;
  
  weight_set_type weights;
  weight_set_type accumulated;
  size_t          updated;

  int debug;
};

struct Dumper
{
  typedef std::pair<path_type, weight_set_type > value_type;
  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;
  
  
  Dumper(queue_type& __queue)
    : queue(__queue) {}
  
  void operator()()
  {
    value_type value;
    
    while (1) {
      queue.pop_swap(value);
      if (value.first.empty()) break;
      
      utils::compress_ostream os(value.first, 1024 * 1024);
      os.precision(20);
      os << value.second;
    }
  }

  queue_type& queue;
};
