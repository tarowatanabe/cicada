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
#include <map>

#include <boost/tuple/tuple.hpp>
#include <boost/tokenizer.hpp>

#include <utils/space_separator.hpp>
#include <utils/base64.hpp>

#include <utils/arc_list.hpp>

#include <utils/lockfree_list_queue.hpp>
#include <utils/vector2.hpp>

#include <cicada/optimize/line_search.hpp>

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

struct LineSearch
{
  typedef cicada::eval::Scorer         scorer_type;
  typedef cicada::eval::ScorerDocument scorer_document_type;
  
  typedef scorer_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_type::score_ptr_type  score_ptr_type;

  typedef cicada::optimize::LineSearch line_search_type;
  
  typedef line_search_type::segment_type          segment_type;
  typedef line_search_type::segment_set_type      segment_set_type;
  typedef line_search_type::segment_document_type segment_document_type;

  typedef line_search_type::value_type optimum_type;
  
  typedef cicada::semiring::Envelope envelope_type;
  typedef std::vector<envelope_type, std::allocator<envelope_type> >  envelope_set_type;

  typedef cicada::Sentence sentence_type;

  LineSearch(const int __debug=0) : line_search(__debug), debug(__debug) {}

  segment_document_type segments;
  envelope_set_type     envelopes;
  sentence_type         yield;
  
  line_search_type line_search;
  int debug;
  
  template <typename HypergraphSet, typename ScorerSet>
  double operator()(const HypergraphSet& graphs, const ScorerSet& scorers, const weight_set_type& origin, const weight_set_type& direction)
  {
    segments.clear();
    segments.resize(graphs.size());
    
    for (int seg = 0; seg < graphs.size(); ++ seg) 
      if (graphs[seg].is_valid()) {
	
	if (debug >= 4)
	  std::cerr << "line-search segment: " << seg << std::endl;

	envelopes.clear();
	envelopes.resize(graphs[seg].nodes.size());
	
	cicada::inside(graphs[seg], envelopes, cicada::semiring::EnvelopeFunction<weight_set_type>(origin, direction));
	
	const envelope_type& envelope = envelopes[graphs[seg].goal];
	const_cast<envelope_type&>(envelope).sort();
	
	envelope_type::const_iterator eiter_end = envelope.end();
	for (envelope_type::const_iterator eiter = envelope.begin(); eiter != eiter_end; ++ eiter) {
	  const envelope_type::line_ptr_type& line = *eiter;
	  
	  line->yield(yield);
	  
	  scorer_type::score_ptr_type score = scorers[seg]->score(yield);
	  
	  if (debug >= 4)
	    std::cerr << "segment: " << seg << " x: " << line->x << std::endl;
	  
	  segments[seg].push_back(std::make_pair(line->x, score));
	}
      }
    
    // traverse segments...
    const optimum_type optimum = line_search(segments, 0.0, 1.0, false);
    
    return (optimum.lower + optimum.upper) * 0.5;
  }
  
};


struct OptimizeCP
{
  // cutting-plane optimizer...

  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::eval::Scorer         scorer_type;
  typedef cicada::eval::ScorerDocument scorer_document_type;
  
  typedef scorer_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_type::score_ptr_type  score_ptr_type;
  
  typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;
  typedef std::deque<scorer_ptr_type, std::allocator<scorer_ptr_type> > scorer_set_type;

  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    gradient_type;
  typedef std::vector<double, std::allocator<double> >    error_bound_type;
  typedef std::vector<double, std::allocator<double> >    objective_type;
  
  typedef std::vector<int, std::allocator<int> > pos_set_type;
  typedef std::vector<pos_set_type, std::allocator<pos_set_type> >  pos_map_type;
  
  typedef std::vector<int, std::allocator<int> > ids_type;
  typedef std::vector<double, std::allocator<double> > labels_type;
  typedef std::vector<double, std::allocator<double> > margins_type;
  typedef std::vector<feature_set_type, std::allocator<feature_set_type> > features_type;
  
  typedef std::vector<int, std::allocator<int> >          timestamp_type;

  typedef utils::vector2<double, std::allocator<double> > matrix_type;
  
  OptimizeCP(const weight_set_type& __weights,
	     const double& __lambda,
	     const double& __tolerance,
	     const int __debug=0)
    : line_search(__debug),
      lambda(__lambda),
      C(1.0 / __lambda),
      tolerance(__tolerance),
      objective_max(- std::numeric_limits<double>::infinity()),
      objective_min(  std::numeric_limits<double>::infinity()),
      weights(__weights),
      accumulated(),
      updated(1),
      debug(__debug) {}
  
  void finalize()
  {
    
  }
  
  void initialize()
  {
    objective_max = - std::numeric_limits<double>::infinity();
    objective_min =   std::numeric_limits<double>::infinity();
    
    accumulated = weights;
    updated = 1;
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

  template <typename IDSet, typename LabelSet, typename MarginSet, typename FeatureSet>
  void operator()(const IDSet&      __ids,
		  const LabelSet&   __labels,
		  const MarginSet&  __margins,
		  const FeatureSet& __features,
		  const bool optimized)
  {
    // QP solver used on OCAS
    //h_matrix<LabelSet, FeatureSet>  H(labels, features);
    
    size_t num_instance = 0;
    for (int i = 0; i < __labels.size(); ++ i) {
      if (__margins[i] <= 0) continue;
      
      const double margin = __labels[i] * __features[i].dot(weights);
      const double grad = __margins[i] - margin;
      
      // we have reached better error bound
      if (__ids[i] < error_bound.size() && error_bound[__ids[i]] >= grad) continue;
      
      // we have achieved some tolerance
      if (__ids[i] < objective.size() && C * grad - objective[__ids[i]] <= tolerance) continue;
      
      ids.push_back(__ids[i]);
      labels.push_back(__labels[i]);
      margins.push_back(__margins[i]);
      features.push_back(__features[i]);
      
      objective_max = std::max(objective_max, grad);
      objective_min = std::min(objective_min, grad);
      
      ++ num_instance;
    }
    
    if (! num_instance) return;
    
    const size_t model_size = labels.size();
    const size_t model_size_prev = model_size - num_instance;
    
    H_new.resize(model_size, model_size, 0.0);
    for (int i = 0; i < model_size_prev; ++ i)
      std::copy(H.begin(i), H.end(i), H_new.begin(i));
    H.swap(H_new);

    // re-compute new H for  i >= size_prev || j >= size_prev
    // 
    for (int i = 0; i < model_size; ++ i)
      for (int j = (i >= model_size_prev ? size_type(0) : model_size_prev); j < model_size; ++ j)
	H(i, j) = labels[i] * labels[j] * features[i].dot(features[j]);
    
    alpha.resize(model_size, 0.0);
    
    gradient.clear();
    gradient.assign(margins.begin(), margins.end());
    
    pos_map.clear();
    alpha_neq.clear();
    for (int i = 0; i < model_size; ++ i) {
      const int pos = ids[i];
      const size_t size = std::max(pos_map.size(), size_t(pos + 1));
      
      pos_map.resize(size);
      alpha_neq.resize(size, C);
      
      pos_map[pos].push_back(i);
      alpha_neq[pos] -= alpha[i];
    }
    
    for (int i = 0; i < model_size; ++ i)
      for (int j = 0; j < model_size; ++ j)
	gradient[i] -= H(i, j) * alpha[j];
    
    double obj_primal = 0.0;
    double obj_dual   = 0.0;
    
    objective.clear();
    objective.resize(pos_map.size(), 0.0);
    
    size_type num_samples = 0;
    for (int k = 0; k < pos_map.size(); ++ k) 
      if (! pos_map[k].empty()) {
	++ num_samples;
	
	double obj_primal_local = 0.0;
	objective[k] = 0.0;
	for (pos_set_type::const_iterator kiter = pos_map[k].begin(); kiter != pos_map[k].end(); ++ kiter) {
	  obj_primal_local = std::max(obj_primal_local, gradient[*kiter]);
	  objective[k] += gradient[*kiter] * alpha[*kiter];
	}
	obj_primal += C * obj_primal_local;
	obj_dual += objective[k];
      }
    
    if (debug)
      std::cerr << "initial primal: " << obj_primal << " dual: " << obj_dual << std::endl;
    
    bool perform_update = false;
    for (int iter = 0; iter != 1000; ++ iter) {

      bool perform_update_local = false;
      
      for (int k = 0; k < pos_map.size(); ++ k) 
	if (! pos_map[k].empty()) {
	  int u = -1;
	  double max_obj = - std::numeric_limits<double>::infinity();
	  double delta = 0.0;
	  
	  for (pos_set_type::const_iterator kiter = pos_map[k].begin(); kiter != pos_map[k].end(); ++ kiter) {
	    delta -= alpha[*kiter] * gradient[*kiter];
	    
	    if (gradient[*kiter] > max_obj) {
	      max_obj = gradient[*kiter];
	      u = *kiter;
	    }
	  }
	  
	  if (gradient[u] < 0.0)
	    u = -1;
	  else
	    delta += C * gradient[u];
	  
	  if (delta <= tolerance) continue;
	  
	  if (u >= 0) {
	    int v = -1;
	    double max_improvement = - std::numeric_limits<double>::infinity();
	    double tau = 1.0;
	    
	    for (pos_set_type::const_iterator kiter = pos_map[k].begin(); kiter != pos_map[k].end(); ++ kiter) 
	      if (*kiter != u && alpha[*kiter] > 0.0) {
		
		const double numer = alpha[*kiter] * (gradient[u] - gradient[*kiter]);
		const double denom = alpha[*kiter] * alpha[*kiter] * (H(u, u) - 2.0 * H(u, *kiter) + H(*kiter, *kiter));
		
		if (denom > 0.0) {
		  const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		  
		  if (improvement > max_improvement) {
		    max_improvement = improvement;
		    tau = std::max(0.0, std::min(1.0, numer / denom));
		    v = *kiter;
		  }
		}
	      }
	    
	    if (alpha_neq[k] > 0.0) {
	      const double numer = alpha_neq[k] * gradient[u];
	      const double denom = alpha_neq[k] * alpha_neq[k] * H(u, u);
	      
	      if (denom > 0.0) {
		const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		
		if (improvement > max_improvement) {
		  max_improvement = improvement;
		  tau = std::max(0.0, std::min(1.0, numer / denom));
		  v = -1;
		}
	      }
	    }
	    
	    if (v >= 0) {
	      // maximize objective, u and v
	      double update = alpha[v] * tau;
	      if (alpha[u] + update < 0.0)
		update = - alpha[u];
	      
	      if (update != 0.0) {
		perform_update_local = true;
		
		alpha[u] += update;
		alpha[v] -= update;
		for (int i = 0; i < model_size; ++ i)
		  gradient[i] += update * (H(i, v) - H(i, u));
	      }
	    } else {
	      double update = alpha_neq[k] * tau;
	      if (alpha[u] + update < 0.0)
		update = - alpha[u];
	      
	      if (update != 0.0) {
		perform_update_local = true;
		
		alpha[u] += update;
		alpha_neq[k] -= update;
		for (int i = 0; i < model_size; ++ i)
		  gradient[i] -= update * H(i, u);
	      }
	    }
	    
	  } else {
	    int v = -1;
	    double max_improvement = - std::numeric_limits<double>::infinity();
	    double tau = 1.0;
	    
	    for (pos_set_type::const_iterator kiter = pos_map[k].begin(); kiter != pos_map[k].end(); ++ kiter) 
	      if (alpha[*kiter] > 0.0) {
		
		const double numer = alpha[*kiter] * gradient[*kiter];
		const double denom = alpha[*kiter] * alpha[*kiter] * H(*kiter, *kiter);
		
		if (denom > 0.0) {
		  const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		  
		  if (improvement > max_improvement) {
		    max_improvement = improvement;
		    tau = std::max(0.0, std::min(1.0, numer / denom));
		    v = *kiter;
		  }
		}
	      }
	    
	    if (v >= 0) {
	      double update = alpha[v] * tau;
	      if (alpha_neq[k] + update < 0.0)
		update = - alpha_neq[k];
	      
	      if (update != 0.0) {
		perform_update_local = true;
		
		alpha_neq[k] += update;
		alpha[v] -= update;
		for (int i = 0; i < model_size; ++ i)
		  gradient[i] += update * H(i, v);
	      }
	    }
	  }
	}
      
      perform_update |= perform_update_local;
      
      if (! perform_update_local) break;
      
      // compute primal/dual
      obj_primal = 0.0;
      obj_dual   = 0.0;
      for (int k = 0; k < pos_map.size(); ++ k) 
	if (! pos_map[k].empty()) {
	  ++ num_samples;
	  
	  double obj_primal_local = 0.0;
	  objective[k] = 0.0;
	  for (pos_set_type::const_iterator kiter = pos_map[k].begin(); kiter != pos_map[k].end(); ++ kiter) {
	    obj_primal_local = std::max(obj_primal_local, gradient[*kiter]);
	    objective[k] += gradient[*kiter] * alpha[*kiter];
	  }
	  obj_primal += C * obj_primal_local;
	  obj_dual += objective[k];
	}
      
      if (obj_primal - obj_dual <= tolerance * num_samples)
	break;
    }
    
    if (debug)
      std::cerr << "final primal: " << obj_primal << " dual: " << obj_dual << std::endl;

    if (! perform_update) return;

    if (optimized) {
      timestamp.resize(model_size, 0);
      weights_new.clear();
      
      for (int i = 0; i < labels.size(); ++ i) {
	if (alpha[i] > 0.0) {
	  typename FeatureSet::value_type::const_iterator fiter_end = features[i].end();
	  for (typename FeatureSet::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter)
	    weights_new[fiter->first] += alpha[i] * labels[i] * fiter->second;
	}
	
	timestamp[i] = ((int(alpha[i] != 0.0) - 1) & (timestamp[i] + 1));
      }
      
      direction = weights_new;
      direction -= weights;
      
      if (! direction.empty()) {
	const double update = line_search(hypergraphs, scorers, weights, direction);
	if (update != 0.0)
	  direction *= update;
	
	if (debug)
	  std::cerr << "optimized update: " << update << std::endl;
		
	weights += direction;
	accumulated += weights;

	++ updated;
      }
    } else {
      timestamp.resize(model_size, 0);
      weights.clear();
      for (int i = 0; i < labels.size(); ++ i) {
	if (alpha[i] > 0.0) {
	  typename FeatureSet::value_type::const_iterator fiter_end = features[i].end();
	  for (typename FeatureSet::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	    const double value = alpha[i] * labels[i] * fiter->second;
	    
	    weights[fiter->first]     += value;
	    accumulated[fiter->first] += value;
	  }
	}
	
	timestamp[i] = ((int(alpha[i] != 0.0) - 1) & (timestamp[i] + 1));
      }
      
      ++ updated;
    }
    
    // if the alpha value for a support vector is zero for a fixed amount of time, 50, then, prune...
    int pos_last = model_size;
    for (int pos = 0; pos < pos_last; /**/) 
      if (timestamp[pos] >= 50) {
	if (pos_last - 1 != pos) {
	  const int pos_swap1 = pos;
	  const int pos_swap2 = pos_last - 1;
	  
	  std::swap(timestamp[pos_swap1], timestamp[pos_swap2]);
	  std::swap(alpha[pos_swap1],     alpha[pos_swap2]);
	  
	  std::swap(ids[pos_swap1],      ids[pos_swap2]);
	  std::swap(labels[pos_swap1],   labels[pos_swap2]);
	  std::swap(margins[pos_swap1],  margins[pos_swap2]);
	  std::swap(features[pos_swap1], features[pos_swap2]);
	  
	  H.swap(pos_swap1, pos_swap2);
	}
	-- pos_last;
      } else
	++ pos;
    
    if (pos_last != model_size) {
      timestamp.resize(pos_last);
      alpha.resize(pos_last);
      
      ids.resize(pos_last);
      labels.resize(pos_last);
      margins.resize(pos_last);
      features.resize(pos_last);
      
      H_new.resize(pos_last, pos_last, 0.0);
      for (int i = 0; i < pos_last; ++ i)
	std::copy(H.begin(i), H.begin(i) + pos_last, H_new.begin(i));
      H.swap(H_new);

      if (debug)
	std::cerr << "support vector size: " << model_size << " pruned: " << pos_last << std::endl;
    }

    error_bound.clear();
    error_bound.resize(objective.size(), - std::numeric_limits<double>::infinity());
    for (int i = 0; i < labels.size(); ++ i)
      error_bound[ids[i]] = std::max(error_bound[ids[i]], margins[i] - labels[i] * features[i].dot(weights));
  }

  LineSearch line_search;

  hypergraph_set_type hypergraphs;
  scorer_set_type     scorers;
  
  ids_type      ids;
  labels_type   labels;
  margins_type  margins;
  features_type features;

  matrix_type H;
  matrix_type H_new;
  
  alpha_type       alpha;
  alpha_type       alpha_neq;
  gradient_type    gradient;
  error_bound_type error_bound;
  objective_type   objective;
  pos_map_type     pos_map;

  timestamp_type   timestamp;
  
  double lambda;
  double C;
  double tolerance;

  double objective_max;
  double objective_min;
  
  weight_set_type weights;
  weight_set_type weights_new;
  weight_set_type accumulated;
  weight_set_type direction;
  size_t          updated;

  int debug;
};

struct OptimizeMIRA
{
  OptimizeMIRA(const weight_set_type& __weights,
	       const double& __lambda,
	       const double& __tolerance,
	       const int __debug=0)
    : line_search(__debug),
      lambda(__lambda),
      C(1.0 / __lambda),
      tolerance(__tolerance),
      objective_max(- std::numeric_limits<double>::infinity()),
      objective_min(  std::numeric_limits<double>::infinity()),
      weights(__weights),
      accumulated(),
      updated(1),
      debug(__debug) {}

  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef cicada::eval::Scorer         scorer_type;
  typedef cicada::eval::ScorerDocument scorer_document_type;
  
  typedef scorer_type::scorer_ptr_type scorer_ptr_type;
  typedef scorer_type::score_ptr_type  score_ptr_type;
  
  typedef std::deque<hypergraph_type, std::allocator<hypergraph_type> > hypergraph_set_type;
  typedef std::deque<scorer_ptr_type, std::allocator<scorer_ptr_type> > scorer_set_type;

  typedef std::vector<double, std::allocator<double> >    alpha_type;
  typedef std::vector<double, std::allocator<double> >    gradient_type;
  typedef std::vector<bool, std::allocator<bool> >        skipped_type;

  typedef std::vector<int, std::allocator<int> > pos_set_type;
  typedef std::pair<pos_set_type, double> id_map_type;
  typedef std::map<int, id_map_type, std::less<int>, std::allocator<std::pair<const int, id_map_type> > >  pos_map_type;

  void finalize()
  {
    accumulated *= - 1.0 / updated;
    accumulated += weights;
    accumulated *= updated;
  }

  void initialize()
  {
    objective_max = - std::numeric_limits<double>::infinity();
    objective_min =   std::numeric_limits<double>::infinity();

    accumulated.clear();
    updated = 1;
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
  
  template <typename IDSet, typename LabelSet, typename MarginSet, typename FeatureSet>
  void operator()(const IDSet&      ids,
		  const LabelSet&   labels,
		  const MarginSet&  margins,
		  const FeatureSet& features,
		  const bool optimized)
  {
    // QP solver used on OCAS
    h_matrix<LabelSet, FeatureSet>  H(labels, features);
    
    alpha.clear();
    gradient.clear();
    skipped.clear();

    alpha.reserve(labels.size());
    gradient.reserve(labels.size());
    skipped.reserve(labels.size());
    
    alpha.resize(labels.size(), 0.0);
    gradient.resize(labels.size(), 0.0);
    skipped.resize(labels.size(), false);

    pos_map.clear();
    
    size_t num_instance = 0;
    for (int i = 0; i < labels.size(); ++ i) {
      gradient[i] = margins[i] - labels[i] * features[i].dot(weights);
      
      const bool skipping = (gradient[i] <= 0 || margins[i] <= 0);
      
      skipped[i] = skipping;
      num_instance += ! skipping;
      
      if (! skipping) {
	pos_map[ids[i]].first.push_back(i);
	pos_map[ids[i]].second = C;
	objective_max = std::max(objective_max, gradient[i]);
	objective_min = std::min(objective_min, gradient[i]);
      }
    }
        
    if (! num_instance) return;

    const size_t model_size = num_instance;

    double obj_primal = 0.0;
    double obj_dual   = 0.0;

    for (pos_map_type::iterator miter = pos_map.begin(); miter != pos_map.end(); ++ miter) {
      const pos_set_type& pos_set = miter->second.first;
      
      double obj_primal_local = 0.0;
      double objective = 0.0;
      
      for (pos_set_type::const_iterator kiter = pos_set.begin(); kiter != pos_set.end(); ++ kiter) {
	obj_primal_local = std::max(obj_primal_local, gradient[*kiter]);
	objective += gradient[*kiter] * alpha[*kiter];
      }
      
      obj_primal += C * obj_primal_local;
      obj_dual += objective;
    }
    
    if (debug)
      std::cerr << "initial primal: " << obj_primal << " dual: " << obj_dual << std::endl; 
    
    bool perform_update = false;
    for (int iter = 0; iter != 1000; ++ iter) {
      
      bool perform_update_local = false;
      
      for (pos_map_type::iterator miter = pos_map.begin(); miter != pos_map.end(); ++ miter) {
	const int k = miter->first;
	const pos_set_type& pos_set = miter->second.first;
	double& alpha_neq = miter->second.second;
	
	int u = -1;
	double max_obj = - std::numeric_limits<double>::infinity();
	double delta = 0.0;
	
	for (pos_set_type::const_iterator kiter = pos_set.begin(); kiter != pos_set.end(); ++ kiter) {
	  delta -= alpha[*kiter] * gradient[*kiter];
	  
	  if (gradient[*kiter] > max_obj) {
	    max_obj = gradient[*kiter];
	    u = *kiter;
	  }
	}
	
	if (gradient[u] < 0.0)
	  u = -1;
	else
	  delta += C * gradient[u];
	
	if (delta <= tolerance) continue;
	
	if (u >= 0) {
	  int v = -1;
	  double max_improvement = - std::numeric_limits<double>::infinity();
	  double tau = 1.0;
	  
	  for (pos_set_type::const_iterator kiter = pos_set.begin(); kiter != pos_set.end(); ++ kiter) 
	    if (*kiter != u && alpha[*kiter] > 0.0) {
	      const double numer = alpha[*kiter] * (gradient[u] - gradient[*kiter]);
	      const double denom = alpha[*kiter] * alpha[*kiter] * (H(u, u) - 2.0 * H(u, *kiter) + H(*kiter, *kiter));
	      
	      if (denom > 0.0) {
		const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		
		if (improvement > max_improvement) {
		  max_improvement = improvement;
		  tau = std::max(0.0, std::min(1.0, numer / denom));
		  v = *kiter;
		}
	      }
	    }
	  
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
	  
	  if (v >= 0) {
	    // maximize objective, u and v
	    double update = alpha[v] * tau;
	    if (alpha[u] + update < 0.0)
	      update = - alpha[u];
	      
	    if (update != 0.0) {
	      perform_update_local = true;
	      
	      alpha[u] += update;
	      alpha[v] -= update;
	      for (int i = 0; i < model_size; ++ i)
		if (! skipped[i])
		  gradient[i] += update * (H(i, v) - H(i, u));
	    }
	  } else {
	    double update = alpha_neq * tau;
	    if (alpha[u] + update < 0.0)
	      update = - alpha[u];
	      
	    if (update != 0.0) {
	      perform_update_local = true;
		
	      alpha[u] += update;
	      alpha_neq -= update;
	      for (int i = 0; i < model_size; ++ i)
		if (! skipped[i])
		  gradient[i] -= update * H(i, u);
	    }
	  }
	    
	} else {
	  int v = -1;
	  double max_improvement = - std::numeric_limits<double>::infinity();
	  double tau = 1.0;
	    
	  for (pos_set_type::const_iterator kiter = pos_set.begin(); kiter != pos_set.end(); ++ kiter) 
	    if (alpha[*kiter] > 0.0) {
	      
	      const double numer = alpha[*kiter] * gradient[*kiter];
	      const double denom = alpha[*kiter] * alpha[*kiter] * H(*kiter, *kiter);
		
	      if (denom > 0.0) {
		const double improvement = (numer < denom ? (numer * numer) / denom : numer - 0.5 * denom);
		  
		if (improvement > max_improvement) {
		  max_improvement = improvement;
		  tau = std::max(0.0, std::min(1.0, numer / denom));
		  v = *kiter;
		}
	      }
	    }
	    
	  if (v >= 0) {
	    double update = alpha[v] * tau;
	    if (alpha_neq + update < 0.0)
	      update = - alpha_neq;
	    
	    if (update != 0.0) {
	      perform_update_local = true;
		
	      alpha_neq += update;
	      alpha[v] -= update;
	      for (int i = 0; i < model_size; ++ i)
		if (! skipped[i])
		  gradient[i] += update * H(i, v);
	    }
	  }
	}
      }
      
      perform_update |= perform_update_local;
      
      if (! perform_update_local) break;
      
      // compute primal/dual
      obj_primal = 0.0;
      obj_dual   = 0.0;
      for (pos_map_type::iterator miter = pos_map.begin(); miter != pos_map.end(); ++ miter) {
	const pos_set_type& pos_set = miter->second.first;
      
	double obj_primal_local = 0.0;
	double objective = 0.0;
      
	for (pos_set_type::const_iterator kiter = pos_set.begin(); kiter != pos_set.end(); ++ kiter) {
	  obj_primal_local = std::max(obj_primal_local, gradient[*kiter]);
	  objective += gradient[*kiter] * alpha[*kiter];
	}
      
	obj_primal += C * obj_primal_local;
	obj_dual += objective;
      }
      
      if (obj_primal - obj_dual <= tolerance * num_instance)
	break;
    }
    
    if (debug)
      std::cerr << "final primal: " << obj_primal << " dual: " << obj_dual << std::endl;
    
    if (! perform_update) return;
    
    if (optimized) {
      
      direction.clear();
      for (int i = 0; i < labels.size(); ++ i) 
	if (! skipped[i] && alpha[i] > 0.0) {
	  typename FeatureSet::value_type::const_iterator fiter_end = features[i].end();
	  for (typename FeatureSet::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter)
	    direction[fiter->first] += alpha[i] * labels[i] * fiter->second;
	}
      
      if (! direction.empty()) {
	// starting from weights, we will move toward direction
	// we will adjust amount of move!
	
	const double update = line_search(hypergraphs, scorers, weights, direction);
	if (update != 0.0)
	  direction *= update;

	if (debug)
	  std::cerr << "optimized update: " << update << std::endl;
	
	weights += direction;
	
	direction *= updated;
	accumulated += direction;

	++ updated;
      }
    } else {
      bool perform_update = false;
      for (int i = 0; i < labels.size(); ++ i) 
	if (! skipped[i] && alpha[i] > 0.0) {
	  typename FeatureSet::value_type::const_iterator fiter_end = features[i].end();
	  for (typename FeatureSet::value_type::const_iterator fiter = features[i].begin(); fiter != fiter_end; ++ fiter) {
	    weights[fiter->first] += alpha[i] * labels[i] * fiter->second;
	    accumulated[fiter->first] += alpha[i] * labels[i] * fiter->second * updated;
	  }
	  
	  perform_update = true;
	}
      updated += perform_update;
    }
  }

  LineSearch line_search;

  hypergraph_set_type hypergraphs;
  scorer_set_type     scorers;

  alpha_type    alpha;
  gradient_type gradient;
  skipped_type  skipped;
  pos_map_type  pos_map;
  
  double lambda;
  double C;
  double tolerance;

  double objective_max;
  double objective_min;
  
  weight_set_type weights;
  weight_set_type accumulated;
  weight_set_type direction;
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
