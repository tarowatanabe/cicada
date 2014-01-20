//
//  Copyright(C) 2013-2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEARN_BITREE_RNN__IMPL__HPP__
#define __CICADA_LEARN_BITREE_RNN__IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/tokenizer.hpp>

#include <sstream>
#include <vector>
#include <deque>

#include "cicada_kbest_impl.hpp"
#include "cicada_mert_kbest_impl.hpp"

#include "cicada/semiring.hpp"
#include "cicada/eval.hpp"

#include "cicada/kbest.hpp"
#include "cicada/kbest_diverse.hpp"
#include "cicada/operation/traversal.hpp"
#include "cicada/operation/functional.hpp"
#include "cicada/bitree_rnn.hpp"

#include "utils/unordered_set.hpp"
#include "utils/base64.hpp"
#include "utils/space_separator.hpp"
#include "utils/piece.hpp"
#include "utils/config.hpp"
#include "utils/mathop.hpp"
#include "utils/hashmurmur3.hpp"
#include "utils/getline.hpp"
#include "utils/compact_map.hpp"

#include <boost/tokenizer.hpp>

#include <codec/lz4.hpp>

typedef cicada::eval::Scorer         scorer_type;
typedef cicada::eval::ScorerDocument scorer_document_type;

typedef cicada::BiTreeRNN tree_rnn_type;

struct LearnBase : public utils::hashmurmur3<size_t>
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef utils::hashmurmur3<size_t> hasher_type;

  typedef cicada::Symbol word_type;
  typedef cicada::Vocab  vocab_type;

  typedef cicada::BiTreeRNN tree_rnn_type;

  typedef tree_rnn_type::tensor_type tensor_type;

  typedef cicada::HyperGraph hypergraph_type;

  typedef hypergraph_type::rule_type rule_type;

  typedef hypergraph_type::feature_set_type   feature_set_type;
  typedef hypergraph_type::attribute_set_type attribute_set_type;
  
  typedef feature_set_type::feature_type     feature_type;
  typedef attribute_set_type::attribute_type attribute_type;
  
  struct Candidate
  {
    Candidate() : hypothesis_(), graph_() {}
    Candidate(const hypothesis_type& hypothesis)
      : hypothesis_(hypothesis), graph_() {}

    Candidate(const hypothesis_type& hypothesis,
	      const hypergraph_type& graph)
      : hypothesis_(hypothesis), graph_(graph) {}
    
    void swap(Candidate& x)
    {
      std::swap(hypothesis_, x.hypothesis_);
      std::swap(graph_,      x.graph_);
    }
    
    hypothesis_type hypothesis_;
    hypergraph_type graph_;
  };
  
  typedef Candidate candidate_type;
  typedef std::deque<candidate_type, std::allocator<candidate_type> > candidate_set_type;
  typedef std::deque<candidate_set_type, std::allocator<candidate_set_type> > candidate_map_type;
  
  struct Gradient
  {
    feature_set_type weights_;
    tree_rnn_type    theta_;

    size_type        count_;
    
    Gradient()
      : weights_(), theta_(), count_(0) {}
    Gradient(const size_type& hidden, const size_type& embedding)
      : weights_(), theta_(hidden, embedding), count_(0) {}
    Gradient(const feature_set_type& weights, const tree_rnn_type& theta)
      : weights_(weights), theta_(theta), count_(0) {}
    
    void clear()
    {
      weights_.clear();
      theta_.clear();
      
      count_ = 0;
    }
  };
  
  typedef Gradient gradient_type;

  typedef cicada::feature::FrontierBiTreeRNN::feature_name_set_type feature_name_set_type;
  
  typedef std::vector<double, std::allocator<double> > loss_set_type;
  typedef std::vector<double, std::allocator<double> > margin_set_type;
  
  typedef std::vector<size_type, std::allocator<size_type> >  node_map_type;

  typedef std::vector<word_type, std::allocator<word_type> > word_set_type;

  struct State
  {
    State() : source_(vocab_type::EMPTY), target_(vocab_type::EMPTY), node_(size_type(-1)), final_(false) {}
    
    bool final() const { return final_; }
    bool initial() const { return source_ == vocab_type::EMPTY && target_ == vocab_type::EMPTY && node_ == size_type(-1); }
    bool terminal() const { return ! initial() && (source_ != vocab_type::EMPTY || target_ != vocab_type::EMPTY); }
    bool pair() const { return ! initial() && source_ != vocab_type::EMPTY && target_ != vocab_type::EMPTY; }
    bool source() const { return ! initial() && source_ != vocab_type::EMPTY; }
    bool target() const { return ! initial() && target_ != vocab_type::EMPTY; }
    bool non_terminal() const { return ! initial() && node_ != size_type(-1); }
    
    feature_set_type features_;
    tensor_type      layer_;
    tensor_type      delta_;
    
    word_type source_;
    word_type target_;
    size_type node_;
    bool      final_;
  };
  
  typedef State state_type;
  typedef utils::chunk_vector<state_type, 1024 * 1024 / sizeof(state_type), std::allocator<state_type> > state_set_type;

  typedef rule_type::symbol_set_type phrase_type;

  struct cache_phrase_type
  {
    std::string frontier;
    phrase_type phrase;
    
    cache_phrase_type() : frontier(), phrase() {}
  };
  typedef utils::array_power2<cache_phrase_type, 1024 * 4, std::allocator<cache_phrase_type> > cache_phrase_set_type;
    
  LearnBase() 
    : attr_frontier_source_("frontier-source"),
      attr_frontier_target_("frontier-target") {}
  
  bool no_bos_eos_;
  bool skip_sgml_tag_;
  feature_name_set_type names_;
  
  margin_set_type margin_kbests_;
  margin_set_type margin_oracles_;
  
  loss_set_type loss_kbests_;
  loss_set_type loss_oracles_;
  
  node_map_type  node_map_;
  state_set_type states_;
  word_set_type  words_;
  
  cache_phrase_set_type cache_source_;
  cache_phrase_set_type cache_target_;
  phrase_type   phrase_tmp_;

  attribute_type attr_frontier_source_;
  attribute_type attr_frontier_target_;
  
  void initialize(const feature_name_set_type& names,
		  const bool no_bos_eos,
		  const bool skip_sgml_tag,
		  weight_set_type& weights,
		  tensor_type& W,
		  tree_rnn_type& theta)
  {
    names_         = names;
    no_bos_eos_    = no_bos_eos;
    skip_sgml_tag_ = skip_sgml_tag;
    
    finalize(weights, W, theta);
  }
  
  void finalize(weight_set_type& weights,
		tensor_type& W,
		tree_rnn_type& theta) const
  {
    if (names_.size() != theta.hidden_)
      throw std::runtime_error("invalid feature names");
    
    // initialize W from the weights...
    W = tensor_type::Zero(1, theta.hidden_);
    for (size_type i = 0; i != theta.hidden_; ++ i)
      W(0, i) = weights[names_[i]];
  }

  double accumulate(const size_type id,
		    const candidate_set_type& kbests,
		    const candidate_set_type& oracles,
		    const weight_set_type& weights,
		    const tensor_type& W,
		    const tree_rnn_type& theta,
		    gradient_type& gradient)
  {
    // first, compuate pairs...
    margin_kbests_.clear();
    margin_oracles_.clear();
    
    candidate_set_type::const_iterator kiter_begin = kbests.begin();
    candidate_set_type::const_iterator kiter_end   = kbests.end();

    candidate_set_type::const_iterator oiter_begin = oracles.begin();
    candidate_set_type::const_iterator oiter_end   = oracles.end();
    
    for (candidate_set_type::const_iterator kiter = kiter_begin; kiter != kiter_end; ++ kiter)
      margin_kbests_.push_back(cicada::dot_product(weights,
						   kiter->hypothesis_.features.begin(),
						   kiter->hypothesis_.features.end(),
						   0.0));
    
    
    for (candidate_set_type::const_iterator oiter = oiter_begin; oiter != oiter_end; ++ oiter)
      margin_oracles_.push_back(cicada::dot_product(weights,
						    oiter->hypothesis_.features.begin(),
						    oiter->hypothesis_.features.end(),
						    0.0));
    
    size_type num_loss = 0;
    for (size_type k = 0; k != margin_kbests_.size(); ++ k)
      for (size_type o = 0; o != margin_oracles_.size(); ++ o)
	num_loss += (1.0 - (margin_oracles_[o] - margin_kbests_[k])) > 0.0;
    
    // if no errors suffered, we will simply return...
    if (! num_loss)
      return 0.0;

    const double error_factor = 1.0 / num_loss;
    
    loss_kbests_.clear();
    loss_oracles_.clear();
    
    loss_kbests_.resize(margin_kbests_.size());
    loss_oracles_.resize(margin_oracles_.size());
    
    double loss = 0.0;
    
    for (size_type k = 0; k != margin_kbests_.size(); ++ k)
      for (size_type o = 0; o != margin_oracles_.size(); ++ o) {
	const double error = std::max(1.0 - (margin_oracles_[o] - margin_kbests_[k]), 0.0);
	
	if (error == 0.0) continue;
	
	loss_oracles_[o] -= error_factor;
	loss_kbests_[k]  += error_factor;
	
	loss += error;
      }

    ++ gradient.count_;
    
    for (size_type k = 0; k != loss_kbests_.size(); ++ k)
      accumulate(loss_kbests_[k], kbests[k], weights, W, theta, gradient);
    
    for (size_type o = 0; o != loss_oracles_.size(); ++ o)
      accumulate(loss_oracles_[o], oracles[o], weights, W, theta, gradient);
    
    return loss * error_factor;
  }

  struct skipper_epsilon
  {
    bool operator()(const word_type& word) const
    {
      return word == vocab_type::EPSILON;
    }
  };
  
  struct skipper_sgml
  {
    bool operator()(const word_type& word) const
    {
      return word == vocab_type::EPSILON || (word != vocab_type::BOS && word != vocab_type::EOS && word.is_sgml_tag());
    }
  };
  
  void accumulate(const double& loss,
		  const candidate_type& cand,
		  const weight_set_type& weights,
		  const tensor_type& W,
		  const tree_rnn_type& theta,
		  gradient_type& gradient)
  {
    forward(cand, weights, W, theta);
    
    backward(loss, cand, weights, W, theta, gradient);
  }
  
  void forward(const candidate_type& cand,
	       const weight_set_type& weights,
	       const tensor_type& W,
	       const tree_rnn_type& theta)
  {
    node_map_.clear();
    states_.clear();

    const tensor_type init = theta.Bi_.array().unaryExpr(tree_rnn_type::shtanh());
    
    hypergraph_type::node_set_type::const_iterator niter_end = cand.graph_.nodes.end();
    for (hypergraph_type::node_set_type::const_iterator niter = cand.graph_.nodes.begin(); niter != niter_end; ++ niter) {
      typedef hypergraph_type::node_type node_type;
      
      const node_type& node = *niter;
      
      if (node.edges.size() != 1)
	throw std::runtime_error("invalid 1-best hypergraph?");

      const bool final = (node.id == cand.graph_.goal);
      
      node_type::edge_set_type::const_iterator eiter_end = node.edges.end();
      for (node_type::edge_set_type::const_iterator eiter = node.edges.begin(); eiter != eiter_end; ++ eiter) {
	typedef hypergraph_type::edge_type edge_type;
	typedef hypergraph_type::rule_type rule_type;
	
	typedef rule_type::symbol_set_type symbol_set_type;
	
	const edge_type& edge = cand.graph_.edges[*eiter];
	
	states_.push_back(state_type());
	
	state_type& state = states_.back();
	state.layer_ = init;
	state.delta_ = tensor_type::Zero(theta.hidden_, 1);

	if (skip_sgml_tag_)
	  forward(edge, edge.rule->rhs.begin(), edge.rule->rhs.end(), weights, W, theta, skipper_sgml());
	else
	  forward(edge, edge.rule->rhs.begin(), edge.rule->rhs.end(), weights, W, theta, skipper_epsilon());
      }
    }
  }
  
  struct __attribute_string : public boost::static_visitor<const cicada::AttributeVector::string_type&>
  {
    typedef cicada::AttributeVector attribute_set_type;
    
    static 
    const std::string& empty()
    {
      static std::string __empty;
      return __empty;
    }
    
    const attribute_set_type::string_type& operator()(const attribute_set_type::int_type& x) const { return empty(); }
    const attribute_set_type::string_type& operator()(const attribute_set_type::float_type& x) const { return empty(); }
    const attribute_set_type::string_type& operator()(const attribute_set_type::string_type& x) const { return x; }
  };

  template <typename Skipper>
  const phrase_type& cache_phrase(const attribute_set_type& attributes,
				  const attribute_type& attribute,
				  cache_phrase_set_type& caches,
				  Skipper skipper)
  {
    attribute_set_type::const_iterator aiter = attributes.find(attribute);
    if (aiter == attributes.end())
      return phrase_tmp_;
    
    const std::string& frontier = boost::apply_visitor(__attribute_string(), aiter->second);
    
    if (frontier.empty())
      return phrase_tmp_;
    
    const size_type cache_pos = hasher_type::operator()(frontier.begin(), frontier.end(), 0) & (caches.size() - 1);
    
    cache_phrase_type& cache = const_cast<cache_phrase_type&>(caches[cache_pos]);
    
    if (cache.frontier != frontier) {
      typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
      
      utils::piece frontier_piece(frontier);
      tokenizer_type tokenizer(frontier_piece);
      
      cache.frontier = frontier;
      
      cache.phrase.clear();
      tokenizer_type::iterator titer_end = tokenizer.end();
      for (tokenizer_type::iterator titer = tokenizer.begin(); titer != titer_end; ++ titer) {
	const word_type word = *titer;
	
	if (! word.is_non_terminal() && ! skipper(word))
	  cache.phrase.push_back(word);
      }
    }
    
    return cache.phrase;
  }
  
  template <typename Iterator, typename Skipper>
  void forward(const hypergraph_type::edge_type& edge,
	       Iterator first, Iterator last,
	       const weight_set_type& weights,
	       const tensor_type& W,
	       const tree_rnn_type& theta,
	       Skipper skipper)
  {
    const size_type offset1 = 0;
    const size_type offset2 = theta.hidden_;
    const size_type offset_source = theta.hidden_;
    const size_type offset_target = theta.hidden_ + theta.embedding_;
    
    int non_terminal_pos = 0;
    for (/**/; first != last; ++ first)
      if (first->is_non_terminal()) {
	const int __non_terminal_index = first->non_terminal_index();
	const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0,
							    non_terminal_pos,
							    __non_terminal_index - 1);
	++ non_terminal_pos;

	states_.push_back(state_type());
	/**/  state_type& state_next = states_.back();
	const state_type& state_prev = states_[states_.size() - 2];
	const state_type& state_ante = states_[node_map_[edge.tails[antecedent_index]]];
	
	state_next.layer_ = (theta.Bn_
			     + theta.Wn_.block(0, offset1, theta.hidden_, theta.hidden_) * state_prev.layer_
			     + theta.Wn_.block(0, offset2, theta.hidden_, theta.hidden_) * state_ante.layer_
			     ).array().unaryExpr(tree_rnn_type::shtanh());
	state_next.delta_ = tensor_type::Zero(theta.hidden_, 1);
	state_next.node_ = node_map_[edge.tails[antecedent_index]];
      } 
    
    const phrase_type& phrase_source = cache_phrase(edge.attributes, attr_frontier_source_, cache_source_, skipper);
    const phrase_type& phrase_target = cache_phrase(edge.attributes, attr_frontier_target_, cache_target_, skipper);
    
    if (! phrase_source.empty() || ! phrase_target.empty()) {
      
      if (phrase_source.empty()) {
	// only target...
	phrase_type::const_iterator titer_begin = phrase_target.begin();
	phrase_type::const_iterator titer_end   = phrase_target.end();
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	  states_.push_back(state_type());
	  /**/  state_type& state_next = states_.back();
	  const state_type& state_prev = states_[states_.size() - 2];
	  
	  state_next.layer_ = (theta.Bt_
			       + theta.Wt_.block(0, offset1, theta.hidden_, theta.hidden_)    * state_prev.layer_
			       + theta.Wt_.block(0, offset2, theta.hidden_, theta.embedding_) * theta.target_->operator()(*titer)
			       ).array().unaryExpr(tree_rnn_type::shtanh());
	  state_next.delta_ = tensor_type::Zero(theta.hidden_, 1);
	  state_next.target_ = *titer;
	}
      } else if (phrase_target.empty()) {
	// only source...
	phrase_type::const_iterator siter_begin = phrase_source.begin();
	phrase_type::const_iterator siter_end   = phrase_source.end();
	for (phrase_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	  states_.push_back(state_type());
	  /**/  state_type& state_next = states_.back();
	  const state_type& state_prev = states_[states_.size() - 2];
	  
	  state_next.layer_ = (theta.Bs_
			       + theta.Ws_.block(0, offset1, theta.hidden_, theta.hidden_)    * state_prev.layer_
			       + theta.Ws_.block(0, offset2, theta.hidden_, theta.embedding_) * theta.source_->operator()(*siter)
			       ).array().unaryExpr(tree_rnn_type::shtanh());
	  state_next.delta_ = tensor_type::Zero(theta.hidden_, 1);
	  state_next.source_ = *siter;
	}
      } else {
	// pairs...
	phrase_type::const_iterator siter_begin = phrase_source.begin();
	phrase_type::const_iterator siter_end   = phrase_source.end();
	phrase_type::const_iterator titer_begin = phrase_target.begin();
	phrase_type::const_iterator titer_end   = phrase_target.end();
	
	for (phrase_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	  for (phrase_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	    states_.push_back(state_type());
	    /**/  state_type& state_next = states_.back();
	    const state_type& state_prev = states_[states_.size() - 2];
	    
	    state_next.layer_ = (theta.Bp_
				 + theta.Wp_.block(0, offset1, theta.hidden_, theta.hidden_) * state_prev.layer_
				 + (theta.Wp_.block(0, offset_source, theta.hidden_, theta.embedding_)
				    * theta.source_->operator()(*siter))
				 + (theta.Wp_.block(0, offset_target, theta.hidden_, theta.embedding_)
				    * theta.target_->operator()(*titer))
				 ).array().unaryExpr(tree_rnn_type::shtanh());
	    state_next.delta_ = tensor_type::Zero(theta.hidden_, 1);
	    state_next.source_ = *siter;
	    state_next.target_ = *titer;
	  }
      }
    }
    
    // set final!
    state_type& state = states_.back();
    
    state.final_    = true;
    state.features_ = edge.features;
    
    // TODO: check whether the values are similar!
    for (size_type i = 0; i != theta.hidden_; ++ i)
      state.features_[names_[i]] = state.layer_(i, 0);
    
    // set up node-map
    if (edge.head >= node_map_.size())
      node_map_.resize(edge.head + 1, size_type(-1));
    
    node_map_[edge.head] = states_.size() - 1;
  }
  
  void backward(const double& loss,
		const candidate_type& cand,
		const weight_set_type& weights,
		const tensor_type& W,
		const tree_rnn_type& theta,
		gradient_type& gradient)
  {
    const size_type offset1 = 0;
    const size_type offset2 = theta.hidden_;
    const size_type offset_source = theta.hidden_;
    const size_type offset_target = theta.hidden_ + theta.embedding_;
    
    state_set_type::reverse_iterator siter_end = states_.rend();
    for (state_set_type::reverse_iterator siter = states_.rbegin(); siter != siter_end; ++ siter) {
      state_type& state = *siter;
      
      if (state.final()) {
	feature_set_type::const_iterator fiter_end = state.features_.end();
	for (feature_set_type::const_iterator fiter = state.features_.begin(); fiter != fiter_end; ++ fiter)
	  gradient.weights_[fiter->first] += loss * fiter->second;
	
	state.delta_.array() += (state.layer_.array().unaryExpr(tree_rnn_type::dshtanh())
				 * W.transpose().array() * loss);
      }
      
      if (state.initial()) // initial
	gradient.theta_.Bi_ += state.delta_;
      else {
	state_type& state_prev = *(siter + 1);
	
	if (state.terminal()) { // terminal
	  if (state.pair()) {
	    gradient.theta_.Wp_.block(0, offset1, theta.hidden_, theta.hidden_)
	      += state.delta_ * state_prev.layer_.transpose();
	    gradient.theta_.Wp_.block(0, offset_source, theta.hidden_, theta.embedding_)
	      += state.delta_ * theta.source_->operator()(state.source_).transpose();
	    gradient.theta_.Wp_.block(0, offset_target, theta.hidden_, theta.embedding_)
	      += state.delta_ * theta.target_->operator()(state.target_).transpose();
	    gradient.theta_.Bp_ += state.delta_;
	    
	    // propagate back to previous state
	    state_prev.delta_.array() += (state_prev.layer_.array().unaryExpr(tree_rnn_type::dshtanh())
					  * (theta.Wp_.block(0, offset1, theta.hidden_, theta.hidden_).transpose()
					     * state.delta_).array());
	  } else if (state.source()) {
	    gradient.theta_.Ws_.block(0, offset1, theta.hidden_, theta.hidden_)
	      += state.delta_ * state_prev.layer_.transpose();
	    gradient.theta_.Ws_.block(0, offset2, theta.hidden_, theta.embedding_)
	      += state.delta_ * theta.source_->operator()(state.source_).transpose();
	    gradient.theta_.Bs_ += state.delta_;
	    
	    // propagate back to previous state
	    state_prev.delta_.array() += (state_prev.layer_.array().unaryExpr(tree_rnn_type::dshtanh())
					  * (theta.Ws_.block(0, offset1, theta.hidden_, theta.hidden_).transpose()
					     * state.delta_).array());
	    
	  } else if (state.target()) {
	    gradient.theta_.Wt_.block(0, offset1, theta.hidden_, theta.hidden_)
	      += state.delta_ * state_prev.layer_.transpose();
	    gradient.theta_.Wt_.block(0, offset2, theta.hidden_, theta.embedding_)
	      += state.delta_ * theta.target_->operator()(state.target_).transpose();
	    gradient.theta_.Bt_ += state.delta_;
	    
	    // propagate back to previous state
	    state_prev.delta_.array() += (state_prev.layer_.array().unaryExpr(tree_rnn_type::dshtanh())
					  * (theta.Wt_.block(0, offset1, theta.hidden_, theta.hidden_).transpose()
					     * state.delta_).array());
	  } else
	    throw std::runtime_error("invalid state");
	  
	} else { // non-temrinal
	  state_type& state_ante = states_[state.node_];
	  
	  gradient.theta_.Wn_.block(0, offset1, theta.hidden_, theta.hidden_)
	    += state.delta_ * state_prev.layer_.transpose();
	  gradient.theta_.Wn_.block(0, offset2, theta.hidden_, theta.hidden_)
	    += state.delta_ * state_ante.layer_.transpose();
	  gradient.theta_.Bn_ += state.delta_;
	  
	  // propagate back to previous state and antecedent
	  state_prev.delta_.array() += (state_prev.layer_.array().unaryExpr(tree_rnn_type::dshtanh())
					* (theta.Wn_.block(0, offset1, theta.hidden_, theta.hidden_).transpose()
					   * state.delta_).array());
	  
	  state_ante.delta_.array() += (state_ante.layer_.array().unaryExpr(tree_rnn_type::dshtanh())
					* (theta.Wn_.block(0, offset2, theta.hidden_, theta.hidden_).transpose()
					   * state.delta_).array());
	}
      }
    }
  }
};


struct LearnAdaGrad : public LearnBase
{
  LearnAdaGrad(const tree_rnn_type& theta,
	       const double lambda,
	       const double eta0)
    : lambda_(lambda),
      eta0_(eta0),
      G_weights_(),
      G_theta_(theta)
  {
    G_weights_.clear();
    G_theta_.clear();
  }
  
  void learn(weight_set_type& weights,
	     tensor_type& W,
	     tree_rnn_type& theta,
	     const gradient_type& gradient)
  {
    if (! gradient.count_) return;
    
    const double scale = 1.0 / gradient.count_;
    
    // update weights
    feature_set_type::const_iterator giter_end = gradient.weights_.end();
    for (feature_set_type::const_iterator giter = gradient.weights_.begin(); giter != giter_end; ++ giter) 
      if (giter->second != 0.0) {
	double& G = G_weights_[giter->first];
	double& x = weights[giter->first];
	
	G += giter->second * giter->second * scale * scale;
	
	const double rate = eta0_ / std::sqrt(double(1.0) + G);
	const double f = x - rate * scale * giter->second;
	
	x = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
      }
    
    // update theta
    update(theta.Wp_, G_theta_.Wp_, gradient.theta_.Wp_, scale, lambda_ != 0.0);
    update(theta.Bp_, G_theta_.Bp_, gradient.theta_.Bp_, scale, false);

    update(theta.Ws_, G_theta_.Ws_, gradient.theta_.Ws_, scale, lambda_ != 0.0);
    update(theta.Bs_, G_theta_.Bs_, gradient.theta_.Bs_, scale, false);

    update(theta.Wt_, G_theta_.Wt_, gradient.theta_.Wt_, scale, lambda_ != 0.0);
    update(theta.Bt_, G_theta_.Bt_, gradient.theta_.Bt_, scale, false);
    
    update(theta.Wn_, G_theta_.Wn_, gradient.theta_.Wn_, scale, lambda_ != 0.0);
    update(theta.Bn_, G_theta_.Bn_, gradient.theta_.Bn_, scale, false);
    
    update(theta.Bi_, G_theta_.Bi_, gradient.theta_.Bi_, scale, false);
    
    finalize(weights, W, theta);
  }

  template <typename Theta, typename GradVar, typename Grad>
  struct update_visitor_regularize
  {
    update_visitor_regularize(Eigen::MatrixBase<Theta>& theta,
			      Eigen::MatrixBase<GradVar>& G,
			      const Eigen::MatrixBase<Grad>& g,
			      const double& scale,
			      const double& lambda,
			      const double& eta0)
      : theta_(theta), G_(G), g_(g), scale_(scale), lambda_(lambda), eta0_(eta0) {}
    
    void init(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      operator()(value, i, j);
    }
    
    void operator()(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      if (g_(i, j) == 0) return;
      
      G_(i, j) += g_(i, j) * g_(i, j) * scale_ * scale_;
      
      const double rate = eta0_ / std::sqrt(double(1.0) + G_(i, j));
      const double f = theta_(i, j) - rate * scale_ * g_(i, j);

      theta_(i, j) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
    }
    
    Eigen::MatrixBase<Theta>&      theta_;
    Eigen::MatrixBase<GradVar>&    G_;
    const Eigen::MatrixBase<Grad>& g_;
    
    const double scale_;
    const double lambda_;
    const double eta0_;
  };

  struct learning_rate
  {
    learning_rate(const double& eta0) : eta0_(eta0) {}
    
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return (x == 0.0 ? 0.0 : eta0_ / std::sqrt(double(1.0) + x));
    }
    
    const double& eta0_;
  };

  template <typename Theta, typename GradVar, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool regularize)
  {
    if (regularize) {
      update_visitor_regularize<Theta, GradVar, Grad> visitor(theta, G, g, scale, lambda_, eta0_);
      
      theta.visit(visitor);
    } else {
      G.array() += g.array().square() * scale * scale;
      theta.array() -= scale * g.array() * G.array().unaryExpr(learning_rate(eta0_));
    }
  }
  
  double lambda_;
  double eta0_;
  
  weight_set_type G_weights_;
  tree_rnn_type   G_theta_;
};

struct LearnSGD : public LearnBase
{
  LearnSGD(const tree_rnn_type& theta,
	   const double lambda,
	   const double eta0)
    : lambda_(lambda),
      eta0_(eta0)
  { }
  
  void learn(weight_set_type& weights,
	     tensor_type& W,
	     tree_rnn_type& theta,
	     const gradient_type& gradient)
  {
    if (! gradient.count_) return;

    const double scale = 1.0 / gradient.count_;
    
    // update weights
    if (lambda_ != 0.0) {
      feature_set_type::const_iterator giter_end = gradient.weights_.end();
      for (feature_set_type::const_iterator giter = gradient.weights_.begin(); giter != giter_end; ++ giter) 
	if (giter->second != 0.0) {
	  double& x = weights[giter->first];
	  
	  x = x * (1.0 - eta0_ * lambda_) - eta0_ * scale * giter->second;
	}
    } else {
      feature_set_type::const_iterator giter_end = gradient.weights_.end();
      for (feature_set_type::const_iterator giter = gradient.weights_.begin(); giter != giter_end; ++ giter) 
	if (giter->second != 0.0)
	  weights[giter->first] -= eta0_ * scale * giter->second;
    }
    
    // update theta
    update(theta.Wp_, gradient.theta_.Wp_, scale, lambda_ != 0.0);
    update(theta.Bp_, gradient.theta_.Bp_, scale, false);

    update(theta.Ws_, gradient.theta_.Ws_, scale, lambda_ != 0.0);
    update(theta.Bs_, gradient.theta_.Bs_, scale, false);

    update(theta.Wt_, gradient.theta_.Wt_, scale, lambda_ != 0.0);
    update(theta.Bt_, gradient.theta_.Bt_, scale, false);
    
    update(theta.Wn_, gradient.theta_.Wn_, scale, lambda_ != 0.0);
    update(theta.Bn_, gradient.theta_.Bn_, scale, false);
    
    update(theta.Bi_, gradient.theta_.Bi_, scale, false);
    
    finalize(weights, W, theta);
  }

  template <typename Theta, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      const Eigen::MatrixBase<Grad>& g,
	      const double scale,
	      const bool regularize=true) const
  {
    if (regularize)
      theta *= 1.0 - eta0_ * lambda_;
    
    theta.noalias() -= (eta0_ * scale) * g;
  }
  
  double lambda_;
  double eta0_;  
};

struct KBestSentence
{
  typedef LearnBase::candidate_type     candidate_type;
  typedef LearnBase::candidate_set_type candidate_set_type;
  typedef LearnBase::candidate_map_type candidate_map_type;
  
  typedef cicada::semiring::Logprob<double>               weight_type;
  typedef cicada::operation::weight_function<weight_type> function_type;
  
  typedef hypergraph_type::id_type   id_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > edge_set_type;

  template <typename Tp>
  struct unassigned_id
  {
    Tp operator()() const { return Tp(-1); }
  };
  
  typedef utils::compact_map<id_type, id_type,
			     unassigned_id<id_type>, unassigned_id<id_type>,
			     boost::hash<id_type>, std::equal_to<id_type>,
			     std::allocator<std::pair<const id_type, id_type> > > node_map_type;
  
  typedef std::vector<id_type, std::allocator<id_type> > head_set_type;

  struct traversal_type
  {
    typedef cicada::HyperGraph hypergraph_type;
    typedef cicada::Rule       rule_type;
    typedef cicada::Sentence   sentence_type;
    typedef cicada::Vocab      vocab_type;
      
    typedef hypergraph_type::feature_set_type   feature_set_type;
    
    typedef boost::tuple<sentence_type, edge_set_type, feature_set_type> value_type;

    traversal_type() {}

    template <typename Edge, typename Iterator>
    void operator()(const Edge& edge, value_type& yield, Iterator first, Iterator last) const
    {
      boost::get<0>(yield).clear();
      boost::get<1>(yield).clear();
      boost::get<1>(yield).push_back(edge.id);
      boost::get<2>(yield) = edge.features;
      
      int non_terminal_pos = 0;
      rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
      for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	if (titer->is_non_terminal()) {
	  const int __non_terminal_index = titer->non_terminal_index();
	  const int pos = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	  ++ non_terminal_pos;
	  
	  boost::get<0>(yield).insert(boost::get<0>(yield).end(),
				      boost::get<0>(*(first + pos)).begin(),
				      boost::get<0>(*(first + pos)).end());
	} else if (*titer != vocab_type::EPSILON && *titer != vocab_type::BOS && *titer != vocab_type::EOS)
	  boost::get<0>(yield).push_back(*titer);
      
      // collect edge id and features...
      for (/**/; first != last; ++ first) {
	boost::get<1>(yield).insert(boost::get<1>(yield).end(), boost::get<1>(*first).begin(), boost::get<1>(*first).end());
	boost::get<2>(yield) += boost::get<2>(*first);
      }
    }
  };

  typedef cicada::operation::kbest_sentence_filter_unique filter_unique_type;
  typedef cicada::operation::kbest_sentence_filter        filter_type;
  
  typedef traversal_type::value_type derivation_type;
  
  void operator()(operation_set_type& operations, const std::string& input, const weight_set_type& weights, candidate_set_type& kbests)
  {    
    kbests.clear();
    
    if (input.empty()) return;
    
    // assign weights
    operations.assign(weights);
    
    // perform translation
    operations(input);
    
    // generate kbests...
    const hypergraph_type& graph = operations.get_data().hypergraph;

    if (kbest_unique_mode) {
      if (kbest_diversity != 0.0) {
	cicada::KBestDiverse<traversal_type, function_type, filter_unique_type> derivations(graph,
											    kbest_size,
											    traversal_type(),
											    function_type(weights),
											    filter_unique_type(graph),
											    kbest_diversity);
	
	kbest_derivations(graph, derivations, kbests);
      } else {
	cicada::KBest<traversal_type, function_type, filter_unique_type> derivations(graph,
										     kbest_size,
										     traversal_type(),
										     function_type(weights),
										     filter_unique_type(graph));
	
	kbest_derivations(graph, derivations, kbests);
      }
    } else {
      if (kbest_diversity != 0.0) {
	cicada::KBestDiverse<traversal_type, function_type, filter_type> derivations(graph,
										     kbest_size,
										     traversal_type(),
										     function_type(weights),
										     filter_type(),
										     kbest_diversity);
	
	kbest_derivations(graph, derivations, kbests);
      } else {
	cicada::KBest<traversal_type, function_type, filter_type> derivations(graph,
									      kbest_size,
									      traversal_type(),
									      function_type(weights),
									      filter_type());
	
	kbest_derivations(graph, derivations, kbests);
      }
    }
  }

  node_map_type node_maps_;
  head_set_type heads_;
  edge_set_type tails_;
  
  template <typename Derivations>
  void kbest_derivations(const hypergraph_type& graph, const Derivations& derivations, candidate_set_type& kbests)
  {
    typedef Derivations derivation_set_type;
    
    typename derivation_set_type::const_iterator diter_end = derivations.end();
    for (typename derivation_set_type::const_iterator diter = derivations.begin(); diter != diter_end; ++ diter) {
      kbests.push_back(candidate_type(hypothesis_type(boost::get<0>(diter->second).begin(), boost::get<0>(diter->second).end(),
						      boost::get<2>(diter->second).begin(), boost::get<2>(diter->second).end())));
      
      
      const edge_set_type& edges = boost::get<1>(diter->second);
      hypergraph_type& graph_kbest = kbests.back().graph_;
      
      node_maps_.clear();
      
      heads_.clear();
      heads_.reserve(edges.size());

      graph_kbest.clear();
      
      id_type node_id = 0;
      edge_set_type::const_iterator eiter_end = edges.end();
      for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter) {
	std::pair<node_map_type::iterator, bool> result = node_maps_.insert(std::make_pair(graph.edges[*eiter].head, node_id));
	  
	heads_.push_back(result.first->second);
	node_id += result.second;
      }
      
      for (id_type node = 0; node != node_id; ++ node)
	graph_kbest.add_node();
      
      id_type edge_id = 0;
      for (edge_set_type::const_iterator eiter = edges.begin(); eiter != eiter_end; ++ eiter, ++ edge_id) {
	const hypergraph_type::edge_type& edge = graph.edges[*eiter];
	
	tails_.clear();
	hypergraph_type::edge_type::node_set_type::const_iterator titer_end = edge.tails.end();
	for (hypergraph_type::edge_type::node_set_type::const_iterator titer = edge.tails.begin(); titer != titer_end; ++ titer) {
	  node_map_type::const_iterator niter = node_maps_.find(*titer);
	  if (niter == node_maps_.end())
	    throw std::runtime_error("no node?");
	  
	  tails_.push_back(niter->second);
	}
	
	hypergraph_type::edge_type& edge_kbest = graph_kbest.add_edge(tails_.begin(), tails_.end());
	edge_kbest.rule       = edge.rule;
	edge_kbest.features   = edge.features;
	edge_kbest.attributes = edge.attributes;
	
	graph_kbest.connect_edge(edge_kbest.id, heads_[edge_id]);
      }
      
      node_map_type::const_iterator niter = node_maps_.find(graph.goal);
      if (niter == node_maps_.end())
	throw std::runtime_error("did not reach goal?");
      
      graph_kbest.goal = niter->second;
      
      graph_kbest.topologically_sort();
    }
  }
};


struct Oracle
{
  typedef LearnBase::candidate_type     candidate_type;
  typedef LearnBase::candidate_set_type candidate_set_type;
  typedef LearnBase::candidate_map_type candidate_map_type;
  
  typedef std::vector<const candidate_type*, std::allocator<const candidate_type*> > oracle_set_type;
  typedef std::vector<oracle_set_type, std::allocator<oracle_set_type> > oracle_map_type;
  
  template <typename Generator>
  std::pair<score_ptr_type, score_ptr_type>
  operator()(const candidate_map_type& kbests,
	     const scorer_document_type& scorers,
	     candidate_map_type& oracles,
	     Generator& generator)
  {
    typedef std::vector<size_t, std::allocator<size_t> > id_set_type;
    
    score_ptr_type score_1best;    
    
    id_set_type ids;
    boost::random_number_generator<boost::mt19937> gen(generator);
    
    // initialization...
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	ids.push_back(id);
	
	candidate_set_type::const_iterator hiter_end = kbests[id].end();
	for (candidate_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(hiter->hypothesis_);
	  
	  if (! hyp.score)
	    hyp.score = scorers[id]->score(sentence_type(hyp.sentence.begin(), hyp.sentence.end()));
	}
	
	if (! score_1best)
	  score_1best = kbests[id].front().hypothesis_.score->clone();
	else
	  *score_1best += *(kbests[id].front().hypothesis_.score);
      }
    
    if (! score_1best)
      throw std::runtime_error("no evaluation score?");
    
    // assign loss
    for (size_t id = 0; id != kbests.size(); ++ id) 
      if (! kbests[id].empty()) {
	score_ptr_type score_segment = score_1best->clone();
	*score_segment -= *kbests[id].front().hypothesis_.score;
	
	candidate_set_type::const_iterator hiter_end = kbests[id].end();
	for (candidate_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  hypothesis_type& hyp = const_cast<hypothesis_type&>(hiter->hypothesis_);
	  
	  *score_segment += *hyp.score;
	  
	  hyp.loss = score_segment->loss();
	  
	  *score_segment -= *hyp.score;
	}
      }
    
    score_ptr_type score_best;
    score_ptr_type score_curr;
    score_ptr_type score_next;
    
    oracle_map_type oracles_best(kbests.size());
    oracle_map_type oracles_curr(kbests.size());
    oracle_map_type oracles_next(kbests.size());
    
    double objective_best = - std::numeric_limits<double>::infinity();
    double objective_curr = - std::numeric_limits<double>::infinity();
    double objective_next = - std::numeric_limits<double>::infinity();
    
    // 
    // 10 iteration will be fine
    //
    const int max_iter = 10;
    const int min_iter = 4;
    for (int i = 0; i != max_iter; ++ i) {
      
      for (id_set_type::const_iterator iiter = ids.begin(); iiter != ids.end(); ++ iiter) {
	const size_t id = *iiter;
	
	score_ptr_type score_removed = (score_next ? score_next->clone() : score_ptr_type());
	
	if (score_removed && ! oracles_curr[id].empty())
	  *score_removed -= *(oracles_curr[id].front()->hypothesis_.score);
	
	oracles_next[id].clear();
	
	candidate_set_type::const_iterator hiter_end = kbests[id].end();
	for (candidate_set_type::const_iterator hiter = kbests[id].begin(); hiter != hiter_end; ++ hiter) {
	  score_ptr_type score_sample;
	  
	  if (score_removed) {
	    score_sample = score_removed->clone();
	    *score_sample += *(hiter->hypothesis_.score);
	  } else
	    score_sample = hiter->hypothesis_.score->clone();
	  
	  const double objective_sample = score_sample->reward();
	  
	  if (objective_sample > objective_next || oracles_next[id].empty()) {
	    oracles_next[id].clear();
	    oracles_next[id].push_back(&(*hiter));
	    
	    objective_next = objective_sample;
	    score_next     = score_sample;
	  } else if (objective_sample == objective_next) {
	    if (hiter->hypothesis_.sentence.size() > oracles_next[id].front()->hypothesis_.sentence.size()) {
	      oracles_next[id].insert(oracles_next[id].begin(), &(*hiter));
	      
	      score_next = score_sample;
	    } else
	      oracles_next[id].push_back(&(*hiter));
	  }
	}
      }
      
      if (objective_next > objective_best) {
	score_best     = score_next->clone();
	objective_best = objective_next;
	oracles_best   = oracles_next;
      }
      
      if (i >= min_iter && objective_next <= objective_curr) break;
      
      score_curr     = score_next->clone();
      objective_curr = objective_next;
      oracles_curr   = oracles_next;
      
      std::random_shuffle(ids.begin(), ids.end(), gen);
    }
    
    oracles.clear();
    oracles.resize(kbests.size());
    for (size_t id = 0; id != kbests.size(); ++ id)
      for (size_t i = 0; i != oracles_best[id].size(); ++ i)
	oracles[id].push_back(*oracles_best[id][i]);
    
    return std::make_pair(score_1best, score_best);
  }
};

inline
void read_refset(const path_type& refset_path,
		 scorer_document_type& scorers,
		 const size_t shard_rank=0,
		 const size_t shard_size=1)
{
  typedef boost::spirit::istream_iterator iter_type;
  typedef cicada_sentence_parser<iter_type> parser_type;
  
  parser_type parser;
  id_sentence_type id_sentence;

  utils::compress_istream is(refset_path, 1024 * 1024);
  is.unsetf(std::ios::skipws);
  
  iter_type iter(is);
  iter_type iter_end;
  
  while (iter != iter_end) {
    id_sentence.second.clear();
    if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::blank, id_sentence))
      if (iter != iter_end)
	throw std::runtime_error("refset parsing failed");
    
    const size_t id = id_sentence.first;
    
    if (shard_size && (id % shard_size != shard_rank)) continue;
    
    const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
    
    if (id_rank >= scorers.size())
      scorers.resize(id_rank + 1);
    
    if (! scorers[id_rank])
      scorers[id_rank] = scorers.create();
    
    scorers[id_rank]->insert(id_sentence.second);
  }
}

class Event
{
public:
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
public:
  Event() : buffer(), size(0) {}
  Event(const Event& x) : buffer(x.buffer), size(x.size) {}
  Event(const std::string& data) { encode(data); }
  Event& operator=(const std::string& data)
  {
    encode(data);
    return *this;
  }
  
  Event& operator=(const Event& x)
  {
    buffer = x.buffer;
    size   = x.size;
    return *this;
  }

  operator std::string() const { return decode(); }

  bool empty() const { return buffer.empty(); }
  void swap(Event& x)
  {
    buffer.swap(x.buffer);
    std::swap(size, x.size);
  }
  
  void encode(const std::string& data)
  {
    buffer.clear();
    buffer.reserve(data.size() >> 2);
    
    boost::iostreams::filtering_ostream os;
    os.push(codec::lz4_compressor());
    os.push(boost::iostreams::back_inserter(buffer));
    os.write(data.c_str(), data.size());
    os.reset();
    
    buffer_type(buffer).swap(buffer);
    
    size = data.size();
  }
  
  std::string decode() const
  {
    buffer_type output(size);
    
    boost::iostreams::filtering_istream is;
    is.push(codec::lz4_decompressor());
    is.push(boost::iostreams::array_source(&(*buffer.begin()), buffer.size()));
    is.read(&(*output.begin()), size);
    
    return std::string(output.begin(), output.end());
  }
  
private:
  typedef std::vector<char, std::allocator<char> > buffer_type;

private:
  buffer_type buffer;
  size_type   size;
};

namespace std
{
  inline
  void swap(Event& x, Event& y)
  {
    x.swap(y);
  }
};

typedef Event event_type;
typedef std::vector<event_type, std::allocator<event_type> > event_set_type;

inline
void read_events(const path_type& input_path,
		 event_set_type& events,
		 const bool directory_mode,
		 const bool id_mode,
		 const size_t shard_rank,
		 const size_t shard_size)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;

  if (directory_mode) {
    if (! boost::filesystem::is_directory(input_path))
      throw std::runtime_error("input is not directory! " + input_path.string());

    boost::spirit::qi::uint_parser<size_t, 10, 1, -1> id_parser;
    
    size_t id;
    std::string line;
    for (size_t i = 0; /**/; ++ i)
      if (shard_size <= 0 || i % shard_size == shard_rank) {
	const path_type path = input_path / (utils::lexical_cast<std::string>(i) + ".gz");
	
	if (! boost::filesystem::exists(path)) break;
	
	utils::compress_istream is(path, 1024 * 1024);
	utils::getline(is, line);
	
	if (line.empty()) continue;
	
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (id != i)
	  throw std::runtime_error("id doest not match!");
	
	const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	
	if (id_rank >= events.size())
	  events.resize(id_rank + 1);
	
	events[id_rank] = line;
      }
  } else if (id_mode) {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    boost::spirit::qi::uint_parser<size_t, 10, 1, -1> id_parser;
    
    size_t id;
    std::string line;
    while (utils::getline(is, line)) 
      if (! line.empty()) {
	std::string::const_iterator iter     = line.begin();
	std::string::const_iterator iter_end = line.end();
	
	if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::blank, id))
	  throw std::runtime_error("id prefixed input format error");
	
	if (shard_size == 0 || id % shard_size == shard_rank) {
	  const size_t id_rank = (shard_size == 0 ? id : id / shard_size);
	  
	  if (id_rank >= events.size())
	    events.resize(id_rank + 1);
	  
	  events[id_rank] = line;
	}
      }
  } else {
    utils::compress_istream is(input_path, 1024 * 1024);
    
    std::string line;
    for (size_t id = 0; utils::getline(is, line); ++ id) 
      if (shard_size == 0 || id % shard_size == shard_rank) {
	if (! line.empty())
	  events.push_back(utils::lexical_cast<std::string>(id) + " ||| " + line);
	else
	  events.push_back(std::string());
      }
  }
}


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
  
  path_added = path_added.string() + suffix;
  
  if (has_suffix_gz)
    path_added = path_added.string() + ".gz";
  else if (has_suffix_bz2)
    path_added = path_added.string() + ".bz2";
  
  return path_added;
}

#endif
