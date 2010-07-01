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

inline
void encode_support_vectors(std::string& data,
			    const size_t& id,
			    const int& source_length,
			    const sentence_type& viterbi,
			    const hypergraph_type& oracle,
			    const hypergraph_type& violated)
{
  data.clear();
  
  boost::iostreams::filtering_ostream os;
  os.push(boost::iostreams::back_inserter(data));
  
  os << id << ' ' << source_length << " |||";
  os << ' ' << viterbi << " |||";
  os << ' ' << oracle << " |||";
  os << ' ' << violated;
  
  os.pop();
}

inline
void decode_support_vectors(const std::string& data,
			    size_t& id,
			    int& source_length,
			    sentence_type& viterbi,
			    hypergraph_type& oracle,
			    hypergraph_type& violated)
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
  
  if (iter != end)
    throw std::runtime_error("still data remain?");
}

inline
void encode_feature_vectors(std::string& data, const size_t& id, const int& source_length, const double& loss, const sentence_type& viterbi, const feature_set_type& oracle, const feature_set_type& violated)
{
  data.clear();
  boost::iostreams::filtering_ostream os;
  os.push(boost::iostreams::back_inserter(data));
  
  os << id << ' ' << source_length << ' ' << utils::encode_base64(loss) << " |||";
  
  os << ' ' << viterbi << " |||";
  
  feature_set_type::const_iterator oiter_end = oracle.end();
  for (feature_set_type::const_iterator oiter = oracle.begin(); oiter != oiter_end; ++ oiter)
    if (! oiter->first.empty() && oiter->second != feature_set_type::mapped_type())
      os << ' ' << oiter->first << ' ' << utils::encode_base64(oiter->second);
  os << " |||";
  
  feature_set_type::const_iterator viter_end = violated.end();
  for (feature_set_type::const_iterator viter = violated.begin(); viter != viter_end; ++ viter)
    if (! viter->first.empty() && viter->second != feature_set_type::mapped_type())
      os << ' ' << viter->first << ' ' << utils::encode_base64(viter->second);
  os << " |||";

  os.pop();
}

inline
void decode_feature_vectors(const std::string& data, size_t& id, int& source_length, double& loss, sentence_type& viterbi, feature_set_type& oracle, feature_set_type& violated)
{
  typedef boost::tokenizer<utils::space_separator> tokenizer_type;
  
  viterbi.clear();
  oracle.clear();
  violated.clear();
  
  tokenizer_type tokenizer(data);
  
  tokenizer_type::iterator iter = tokenizer.begin();
  tokenizer_type::iterator end = tokenizer.end();
  
  if (iter == end) return;

  id = atoll((*iter).c_str());
  ++ iter;

  if (iter == end)
    throw std::runtime_error("invalid format: invaild end...?");
  
  source_length = atoi((*iter).c_str());
  ++ iter;
  
  if (iter == end)
    throw std::runtime_error("invalid format: invaild end...?");
  
  loss = utils::decode_base64<double>(*iter);
  ++ iter;
  
  if (iter == end || *iter != "|||")
    throw std::runtime_error("invalid format: expected separator");
  ++ iter;
  
  while (iter != end) {
    const std::string token = *iter;
    if (token == "|||") break;
    
    viterbi.push_back(token);
    ++ iter;
  }
  
  if (iter == end || *iter != "|||")
    throw std::runtime_error("invalid format: expected separator after sentence");
  ++ iter;
  
  while (iter != end) {
    const std::string token = *iter;
    if (token == "|||") break;
    
    ++ iter;
    if (iter == end)
      throw std::runtime_error("invalid format: expected base64 encoded data");
    
    oracle[token] = utils::decode_base64<feature_set_type::mapped_type>(*iter);
    ++ iter;
  }

  if (iter == end || *iter != "|||")
    throw std::runtime_error("invalid format: expected separator");
  ++ iter;
  
  while (iter != end) {
    const std::string token = *iter;
    if (token == "|||") break;
    
    ++ iter;
    if (iter == end)
      throw std::runtime_error("invalid format: expected base64 encoded data");
    
    violated[token] = utils::decode_base64<feature_set_type::mapped_type>(*iter);
    ++ iter;
  }
  
}

struct OptimizeSMO
{
  OptimizeSMO(const weight_set_type& __weights, const double& __lambda, const int __debug=0)
    : lambda(__lambda),
      C(1.0 / __lambda),
      weights(__weights),
      accumulated(),
      updated(1),
      debug(__debug) {}

  typedef utils::vector2<double, std::allocator<double> > hessian_type;
  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    gradient_type;
  typedef std::vector<double, std::allocator<double> >    kkt_type;

  void finalize()
  {
    accumulated *= - 1.0 / updated;
    accumulated += weights;
    accumulated *= updated;
  }
  
  
  template <typename LabelSet, typename MarginSet, typename FeatureSet>
  void operator()(const LabelSet&   labels,
		  const MarginSet&  margins,
		  const FeatureSet& features)
  {
    
    // QP solver used on OCAS
    hessian.clear();
    alpha.clear();
    gradient.clear();
    
    hessian.resize(labels.size(), labels.size(), 0.0);
    alpha.resize(labels.size(), 0.0);
    gradient.resize(labels.size(), 0.0);
    
    double alpha_neq = C;
    
    for (int i = 0; i < labels.size(); ++ i) {
      gradient[i] = margins[i] - labels[i] * features[i].dot(weights);
      for (int j = 0; j < labels.size(); ++ j)
	hessian(i, j) = labels[i] * labels[j] * features[i].dot(features[j]);
    }
    
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
	    const double denom = alpha[i] * alpha[i] * (hessian(u, u) - 2.0 * hessian(u, i) + hessian(i, i));
	    
	    if (denom > 0.0) {
	      const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
	      
	      if (improvement > max_improvement) {
		max_improvement = improvement;
		tau = std::max(0.0, std::min(1.0, numer / denom));
		v = i;
	      }
	    }
	  }
	
	// check if virtual variable can be used for update...
#if 0
	if (alpha_neq > 0.0) {
	  const double numer = alpha_neq * gradient[u];
	  const double denom = alpha_neq * alpha_neq * hessian(u, u);
	  
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
	    gradient[i] += update * (hessian(i, v) - hessian(i, u));
	  
	} else if (alpha_neq > 0.0) {
	  double update = alpha_neq * tau;
	  if (alpha[u] + update < 0.0)
	    update = - alpha[u];
	  
	  alpha[u] += update;
	  alpha_neq -= update;
	  
	  // update g..
	  for (int i = 0; i < labels.size(); ++ i)
	    gradient[i] -= update * hessian(i, u);
	}
	
      } else {
	int v = -1;
	double max_improvement = - std::numeric_limits<double>::infinity();;
	double tau = 1.0;
	
	for (int i = 0; i < labels.size(); ++ i) 
	  if (alpha[i] > 0.0) {
	    
	    const double numer = alpha[i] * gradient[i];
	    const double denom = alpha[i] * alpha[i] * hessian(i, i);
	    
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
	    gradient[i] += update * hessian(i, v);
	}
      }
    }
    
    for (int i = 0; i < labels.size(); ++ i) 
      if (alpha[i] > 0.0) {
	typename FeatureSet::value_type::const_iterator fiter_end = features[i].end();
	for (typename FeatureSet::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	  weights[fiter->first] += alpha[i] * labels[i] * fiter->second;
	  accumulated[fiter->first] += alpha[i] * labels[i] * fiter->second * updated;
	}
      }
    ++ updated;
    
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
  
  hessian_type  hessian;
  alpha_type    alpha;
  gradient_type gradient;
  kkt_type      kkt;
  
  double lambda;
  double C;
  
  weight_set_type weights;
  weight_set_type accumulated;
  size_t          updated;
  
  int debug;
};

struct OptimizeSGDL2
{
  OptimizeSGDL2(const double& __C, const bool __bundle=false, const int __debug=0)
    : C(__C),
      weight_norm(0.0),
      weight_scale(0.0),
      epoch(0),
      weights(),
      accumulated(),
      updated(1),
      bundle(__bundle),
      debug(__debug) {}
  
  void rescale(const double& alpha)
  {
    weight_norm *= alpha * alpha;
    if (alpha != 0.0)
      weight_scale *= alpha;
    else {
      weight_scale = 1.0;
      weights.clear();
    }
  }
  
  void operator()(const double loss, const feature_set_type& oracle, const feature_set_type& violated)
  {
    if (weight_scale == 0.0) {
      weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
      weight_scale = 1.0;
    }
    
    if (loss <= 0.0) return;
    
    const double& lambda = C;
    const double eta = 1.0 / (lambda * (epoch + 2));
    ++ epoch;
    
    // perform rescaling...
    rescale(1.0 - eta * lambda);
    
    // update amount: eta * loss...
    
    
    
    // perform final rescaling...
    if (weight_norm > 1.0 / lambda)
      rescale(std::sqrt(1.0 / (lambda * weight_norm)));
  }
  
  void finalize()
  {
    if (weight_scale != 0.0 && weight_scale != 1.0)
      weights *= weight_scale;
    
    weight_norm = std::inner_product(weights.begin(), weights.end(), weights.begin(), 0.0);
    weight_scale = 1.0;
  }
  
  double C;

  double weight_norm;
  double weight_scale;
  size_t epoch;
  
  weight_set_type weights;
  weight_set_type accumulated;
  size_t          updated;
  
  bool bundle;
  
  int debug;

};

struct OptimizeMIRA
{
  OptimizeMIRA(const double& __C, const int __debug=0) : C(__C), weights(), accumulated(), updated(1), debug(__debug) {}

  struct diff_norm
  {
    double operator()(const double& x, const double& y) const{
      return (x - y) * (x - y);
    }
  };
  
  void operator()(const double loss, const feature_set_type& oracle, const feature_set_type& violated)
  {
    const double margin = oracle.dot(weights) - violated.dot(weights);
    const double norm   = oracle.dot(violated, diff_norm());
    
    const double alpha = std::max(0.0, std::min(1.0 / C, (loss - margin) / norm));
    
    if (loss - margin > 0.0) {
      if (debug)
	std::cerr << "loss: " << loss << " margin: " << margin << " norm: " << norm << " alpha: " << alpha << std::endl;
      
      // update...
      feature_set_type::const_iterator riter_end = oracle.end();
      for (feature_set_type::const_iterator riter = oracle.begin(); riter != riter_end; ++ riter) {
	weights[riter->first] += riter->second * alpha;
	accumulated[riter->first] += riter->second * alpha * updated;
      }
      
      feature_set_type::const_iterator piter_end = violated.end();
      for (feature_set_type::const_iterator piter = violated.begin(); piter != piter_end; ++ piter) {
	weights[piter->first] -= piter->second * alpha;
	accumulated[piter->first] -= piter->second * alpha * updated;
      }
      
      ++ updated;
    }
  }

  void finalize()
  {
    accumulated *= - 1.0 / updated;
    accumulated += weights;
    accumulated *= updated;
  }

  double C;
  
  weight_set_type weights;
  weight_set_type accumulated;
  size_t          updated;
  
  int debug;
};
