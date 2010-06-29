#include <numeric>

#include <boost/tuple/tuple.hpp>

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
