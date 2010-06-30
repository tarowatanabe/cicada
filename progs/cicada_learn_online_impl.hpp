#include <iostream>
#include <numeric>
#include <stdexcept>

#include <boost/tuple/tuple.hpp>
#include <boost/tokenizer.hpp>

#include <utils/space_separator.hpp>
#include <utils/base64.hpp>

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
    
  }

  double C;
  
  weight_set_type weights;
  weight_set_type accumulated;
  size_t          updated;
  
  int debug;
};
