//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// 
// ITG-based autoencoder
//
// X -> f/e : [vec_f, vec_e]
// X -> [X1, X2] : W_s [vec_x1_f, vec_x2_f, vec_x1_e, vec_x2_e] + B_s
// X -> <X1, X2> : W_i [vec_x1_f, vec_x2_f, vec_x2_e, vec_x1_e] + B_i
//
// - Use ITG beam search of Saers et al., (2009)
// - Learning is performed by autoencoding: recovering representation of children vectors
//
// word vector format: word dim1 dim2 dim3 ...
// tensor: [[dim1, dim2, ...], [dim1, dim2, ...], ... ] (make it compatible with neuron package...?)
// 

// currently, we learn and output the last derivation in a format similar to pialign:
//
// < [ "word" "word" ] < "word" "word" > >
//
// [ ] indicates straight and < > indicates inversion.
// All the words are c-style double quoted, including SGML tags.
//

// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by Pegasos (default)
// SGD with L2 regularizer inspired by AdaGrad
// SGD with L2/L2 regularizer from RDA
//
// + batch algorithm using LBFGS
//


#include <cmath>
#include <climit>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"

#include "utils/bichart.hpp"
#include "utils/bithack.hpp"
#include "utils/unordered_map.hpp"
#include "utils/random_seed.hpp"
#include "utils/repository.hpp"

struct Bitext
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef cicada::Symbol word_type;
  typedef cicada::Sentence sentence_type;
  typedef cicada::Vocab vocab_type;
  
  Bitext() : source_(), target_() {}
  Bitext(const sentence_type& source, const sentence_type& target) : source_(source), target_(target) {}

  void clear()
  {
    source_.clear();
    target_.clear();
  }

  void swap(Bitext& x)
  {
    source_.swap(x.source_);
    target_.swap(x.target_);
  }
  
  sentence_type source_;
  sentence_type target_;
};

namespace std
{
  inline
  void swap(Bitext& x, Bitet& y)
  {
    x.swap(y);
  }
};

struct Model
{
  // model parameters
  
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef Bitext bitext_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  typedef utils::unordered_map<word_type, tensor_type,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, tensor_type> > >::type embedding_type;

  typedef boost::filesystem::path path_type;
  
  Model() : dimension_(0), lambda_(0) {}
  Model(const size_type& dimension) 
    : dimension_(dimension), lambda_(0) { initialize(dimension); }
  Model(const size_type& dimension, const double& lambdam) 
    : dimension_(dimension), lambda_(lambda) { initialize(dimension); }
  
  
  Model& operator+=(const Model& x)
  {
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = source_.find(siter->first);
      
      if (eiter == source_.end())
	eiter = source_.insert(std::make_pair(siter->first, tensor_type(dimension_, 1)));
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = siter->second.col(0)[i] * x.scale_source_;
	
	norm_source_ += 2.0 * eiter->second.col(0)[i] * scale_source_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / scale_source_;
      }
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      embedding_type::iterator eiter = target_.find(titer->first);
      
      if (eiter == target_.end())
	eiter = source_.insert(std::make_pair(siter->first, tensor_type(dimension_, 1)));
      
      for (difference_type i = 0; i != titer->second.rows(); ++ i) {
	const double amount = titer->second.col(0)[i] * x.scale_target_;
	
	norm_target_ += 2.0 * eiter->second.col(0)[i] * scale_source_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / scale_target_;
      }
    }

    Ws1_ += x.Ws1_;
    bs1_ += x.bs1_;
    Wi1_ += x.Wi1_;
    bi1_ += x.bi1_;

    Ws2_ += x.Ws2_;
    bs2_ += x.bs2_;
    Wi2_ += x.Wi2_;
    bi2_ += x.bi2_;
  }
  
  void clear()
  {
    // embedding
    source_.clear();
    target_.clear();
    scale_source_ = 1.0;
    scale_target_ = 1.0;
    norm_source_ = 0.0;
    norm_target_ = 0.0;
    
    // matrices
    Ws1_.setZero();
    bs1_.setZero();
    Wi1_.setZero();
    bi1_.setZero();

    Ws2_.setZero();
    bs2_.setZero();
    Wi2_.setZero();
    bi2_.setZero();
  }

  void finalize()
  {
    if (scale_source_ == 0.0) {
      source_.clear();
      scale_source_ = 1.0;
      norm_source_ = 0.0;
    } else if (scale_source_ != 1.0) {
      norm_source_ = 0.0;
      
      embedding_type::iterator siter_end = source_.end();
      for (embedding_type::iterator siter = source_.begin(); siter != siter_end; ++ siter) {
	siter->second.array() *= scale_source_;
	norm_source_ += siter->second.squaredNorm();
      }
      
      scale_source_ = 1.0;
    }
    
    if (scale_target_ == 0.0) {
      target_.clear();
      scale_target_ = 1.0;
      norm_target_ = 0.0;      
    } else if  (scale_target_ != 1.0) {
      norm_target_ = 0.0;
      
      embedding_type::iterator titer_end = target_.end();
      for (embedding_type::iterator titer = target_.begin(); titer != titer_end; ++ titer) {
	titer->second.array() *= scale_target_;
	norm_target_ += titer->second.squaredNorm();
      }
      
      scale_target_ = 1.0;
    }
    
  }

  void rescale(const double scaling, const bool ignore_bias)
  {
    if (ignore_bias) {
      if (scaling == 0.0) {
	source_.clear();
	target_.clear();
	scale_source_ = 1.0;
	scale_target_ = 1.0;
	norm_source_ = 0.0;
	norm_target_ = 0.0;
	
	Ws1_.setZero();
	Wi1_.setZero();
	
	Ws2_.setZero();
	Wi2_.setZero();
      } else {
	scale_source_ *= scaling;
	scale_target_ *= scaling;
	norm_source_ *= scaling * scaling;
	norm_target_ *= scaling * scaling;
	
	Ws1_.array() *= scaling;
	Wi1_.array() *= scaling;
	
	Ws2_.array() *= scaling;
	Wi2_.array() *= scaling;
      }
    } else {
      if (scaling == 0.0)
	clear();
      else {
	scale_source_ *= scaling;
	scale_target_ *= scaling;
	norm_source_ *= scaling * scaling;
	norm_target_ *= scaling * scaling;
	
	Ws1_.array() *= scaling;
	bs1_.array() *= scaling;
	Wi1_.array() *= scaling;
	bi1_.array() *= scaling;
	
	Ws2_.array() *= scaling;
	bs2_.array() *= scaling;
	Wi2_.array() *= scaling;
	bi2_.array() *= scaling;
      }
    }
  }

  double squared_norm(bool ignore_bias) const
  {
    double norm = norm_source_ + norm_target;

    norm += Ws1_.squaredNorm();
    norm += Wi1_.squaredNorm();
    
    norm += Ws2_.squaredNorm();
    norm += Wi2_.squaredNorm();
    
    if (! ignore_bian) {
      norm += bs1_.squaredNorm();
      norm += bi1_.squaredNorm();
      
      norm += bs2_.squaredNorm();
      norm += bi2_.squaredNorm();
    }
    
    return norm;
  }
  
  void initialize(const size_type dimension)
  {
    // intialize randomly...
    dimension_ = dimension;
    
    // embedding
    source_.clear();
    target_.clear();
    scale_source_ = 1.0;
    scale_target_ = 1.0;
    norm_source_ = 0.0;
    norm_target_ = 0.0;
    
    // score
    Ws1_ = tensor_type::Random(dimension * 2, dimension * 4);
    bs1_ = tensor_type::Random(dimension * 2, 1);
    Wi1_ = tensor_type::Random(dimension * 2, dimension * 4);
    bi1_ = tensor_type::Random(dimension * 2, 1);
    
    // reconstruction
    Ws2_ = tensor_type::Random(dimension * 4, dimension * 2);
    bs2_ = tensor_type::Random(dimension * 4, 1);
    Wi2_ = tensor_type::Random(dimension * 4, dimension * 2);
    bi2_ = tensor_type::Random(dimension * 4, 1);
  }
  
  template <typename Iterator>
  void embeddig(Iterator first, Iterator last)
  {
    source_.clear();
    target_.clear();
    scale_source_ = 1.0;
    scale_target_ = 1.0;
    
    source_[vocab_type::EPSILON].setRandom();
    target_[vocab_type::EPSILON].setRandom();
    
    norm_source_ = source_[vocab_type::EPSILON].squaredNorm();
    norm_target_ = target_[vocab_type::EPSILON].squaredNorm();
    
    for (/**/; first != last; ++ first) {
      norm_source_ += embedding(first->source, source_);
      norm_target_ += embedding(first->target, target_);
    }
  }
  
  double embedding(const sentence_type& sentence, embedding_type& embedding)
  {
    double norm = 0.0;
    
    sentence_type::const_iterator siter_end = sentence.end();
    for (sentence_type::const_iterator siter = sentence.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = embedding.find(*siter);
      
      if (eiter == embedding.end()) {
	tensor_type& we = embedding[*siter];
	
	we = tensor_type::Random(dimension, 1);
	norm += we.squaredNorm();
      }
    }
    
    return norm;
  }
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    
    
    
  }
  
  // dimension...
  size_type dimension_;
  
  // hyperparameter
  double lambda_;
  
  // Embedding
  embedding_type source_;
  embedding_type target_;
  double scale_source_;
  double scale_target_;
  double norm_source_;
  double norm_target_;
  
  // W{s,i}1 and b{s,i}1 for encoding
  tensor_type Ws1_;
  tensor_type bs1_;
  tensor_type Wi1_;
  tensor_type bi1_;
  
  // W{s,i}2 and b{s,i}2 for reconstruction
  tensor_type Ws2_;
  tensor_type bs2_;
  tensor_type Wi2_;
  tensor_type bi2_;
};


struct ITGTree
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;
  
  typedef model_type::tensor_type tensor_type;
  
  typedef Bitext bitext_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  struct Span
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef uint32_t index_type;
    typedef index_type value_type;
    
    Span(const index_type& first, const index_type& last) : first_(first), last_(last) {}
    template <typename _Index>
    Span(const std::pair<_Index, _Index>& x) : first_(x.first), last_(x.second) {}
    Span() : first_(0), last_(0) {}
    
    bool empty() const { return first_ == last_; }
    size_type size() const { return last_ - first_; }
    
    index_type first_;
    index_type last_;
  };
  
  typedef Span span_type;
  
  struct SpanPair
  {
    typedef Span span_type;
    typedef span_type value_type;
    
    SpanPair(const span_type& source, const span_type& target) : source_(source), target_(target) {}
    tempalte <typename _Span>
    SpanPair(const std::pair<_Span, _Span>& x) : source(x.first), target(x.second) {}
    SpanPair() : source_(), target_() {}
        
    span_type source_;
    span_type target_;
  };
  
  typedef SpanPair span_pair_type;
  
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  typedef std::vector<span_pair_set_type, std::allocator<span_pair_set_type> > agenda_type;

  typedef std::pair<span_pair_type, span_pair_type> span_pair_stack_type;

  typedef std::vector<span_pair_stack_type, std::allocator<span_pair_stack_type> > stack_type;

  typedef std::pair<logprob_type, span_pair_type> score_span_pair_type;
  typedef std::vector<score_span_pair_type, std::allocator<score_span_pair_type> > heap_type;

  typedef std::pair<span_pair_type, span_pair_type> tail_set_type;
  
  struct tail_set_unassigned
  {
    tail_set_type operator()() const
    {
      return tail_set_type(span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)),
			   span_pair_type(size_type(-1), size_type(-1), size_type(-1), size_type(-1)));
    }
  };
  
  typedef utils::compact_set<tail_set_type,
			     tail_set_unassigned, tail_set_unassigned,
			     utils::hashmurmur3<size_t>, std::equal_to<tail_set_type>,
			     std::allocator<tail_set_type> > tail_set_unique_type;
  
  struct Node
  {
    typedef uint32_t index_type;
    typedef std::vector<index_type, std::allocator<index_type> > edge_set_type;
    
    Node() : score_(std::numeric_limits<double>::infinity()), total_(0.0) {}

    bool terminal() const { return tails_.first.empty() && tail_.second.empty(); }
    
    double      score_;
    double      total_;
    
    tensor_type output_;
    tensor_type output_norm_;
    tensor_type delta_;
    
    tensor_type reconstruction_;
    tensor_type delta_reconstruction_;
    
    tail_set_type tails_;
  };
  
  typedef Node node_type;
  typedef utils::bichart<node_type, std::allocator<node_type> > node_set_type;
  
  void clear()
  {
    nodes_.clear();
    
    agenda_.clear();
    heap_.clear();
    uniques_.clear();
  }
  
  // sort by greater item so that we can pop from less
  struct heap_compare
  {
    bool operator()(const score_span_pair_type& x, const score_span_pair_type& y) const
    {
      return x.first > y.first;
    }
  };

  void forward(const sentence_type& source,
	       const sentence_type& target,
	       model_type& theta,
	       const double beam)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    clear();
    
    // initialization
    nodes_.reserve(source_size + 1, target_size + 1);
    nodes_.resize(source_size + 1, target_size + 1);
    
    agenda.reserve(source_size + target_size + 1);
    agenda.resize(source_size + target_size + 1);
        
    // initialization
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = 0; trg <= target_size; ++ trg) 
	if (src < source_size || trg < target_size) {
	  
	  // epsilon at target
	  if (src < source_size)
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg)), theta);
	  
	  // epsilon at source
	  if (trg < target_size)
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg)), theta);
	  
	  // word-pair
	  if (src < source_size && trg < target_size)
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg + 1)), theta);
	}
    
    // iterate!
    const double infty = std::numeric_limits<double>::infinity();
    const size_type length_max = source_size + target_size;
    
    for (size_type length = 1; length != length_max; ++ length) 
      if (! agenda_[length].empty()) {
	span_pair_set_type& spans = agenda_[length];
	
	heap_.clear();
	heap_.reserve(spans.size());
	
	span_pair_set_type::const_iterator siter_end = spans.end();
	for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)  {
	  heap_.push_back(score_span_pair_type(nodes_(siter->source_.first_, siter->source_.last_,
						      siter->target_.first_, siter->target_.last_).score_,
					       *siter));
	  
	  std::push_heap(heap_.begin(), heap_.end(), heap_compare());
	}
	
	heap_type::iterator hiter_begin = heap.begin();
	heap_type::iterator hiter       = heap.end();
	heap_type::iterator hiter_end   = heap.end();
	
	// we will derive those with smaller reconstruction error...
	const double threshold = hiter_begin->first + beam;
	for (/**/; hiter_begin != hiter && hiter_begin->first > threshold; -- hiter)
	  std::pop_heap(hiter_begin, hiter, heap_compare());
	
	uniques_.clear();
	
	for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	  const span_pair_type& span_pair = iter->second;
	  
	  // we borrow the notation...
	  const difference_type l = length;
	  const difference_type s = span_pair_.source_.first_;
	  const difference_type t = span_pair_.source_.last_;
	  const difference_type u = span_pair_.target_.first_;
	  const difference_type v = span_pair_.target_.last_;
	  
	  const difference_type T = source_size;
	  const difference_type V = target_size;
	  
	  // remember, we have processed only upto length l. thus, do not try to combine with spans greather than l!
	  // also, keep used pair of spans in order to remove duplicates.
	  
	  for (difference_type S = utils::bithack::max(s - l, difference_type(0)); S <= s; ++ S) {
	    const difference_type L = l - (s - S);
	    
	    // straight
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == s); ++ U) {
	      // parent span: StUv
	      // span1: SsUu
	      // span2: stuv
	      
	      if (nodes_(S, s, U, u).score_ == infty) continue;
	      
	      const span_pair_type  span1(S, s, U, u);
	      const span_pair_type& span2(span_pair);
	      
	      if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
	      // compute score and add hyperedge
	      forward(span_pair_type(S, t, U, v), span1, span2, theta);
	    }

	    // inversion
	    for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: StuU
	      // span1: SsvU
	      // span2: stuv

	      if (nodes_(S, s, v, U).score_ == infty) continue;
	      
	      const span_pair_type  span1(S, s, v, U);
	      const span_pair_type& span2(span_pair);
	      
	      if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
	      // compute score and add hyperedge
	      forward(span_pair_type(S, t, u, U), span1, span2, theta);
	    }
	  }
	  
	  for (difference_type S = t; S <= utils::bithack::min(t + l, T); ++ S) {
	    const difference_type L = l - (S - t);
	    
	    // inversion
	    for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == t); ++ U) {
	      // parent span: sSUv
	      // span1: stuv
	      // span2: tSUu

	      if (nodes_(t, S, U, u).score_ == infty) continue;
	      
	      const span_pair_type& span1(span_pair);
	      const span_pair_type  span2(t, S, U, u);
	      
	      if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
	      // compute score and add hyperedge
	      forward(span_pair_type(s, S, U, v), span1, span2, theta);
	    }
	    
	    // straight
	    for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
	      // parent span: sSuU
	      // span1: stuv
	      // span2: tSvU
	      
	      if (nodes_(t, S, v, U).score_ == infty) continue;
	      
	      const span_pair_type& span1(span_pair);
	      const span_pair_type  span2(t, S, v, U);
	      
	      if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
	      forward(span_pair_type(s, S, u, U), span1, span2, theta);
	    }
	  }
	}
      }
  }

  // terminal!
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const span_pair_type& parent,
	       model_type& theta)
  {
    const size_type dimension = theta.dimension_;
    
    node_type& node = nodes_(parent.source_.first_, parent.source_.last_, parent.target_.first_, parent.target_.last_);
    
    const word_type& embedding_source = (! parent.source_.empty() ? source[parent.source_.first_] : vocab_type::EPSILON);
    const word_type& embedding_target = (! parent.target_.empty() ? target[parent.target_.first_] : vocab_type::EPSILON);
    
    node.score_ = 0.0;
    node.total_ = 0.0;
    node.output_norm_ = tensor_type(dimension * 2, 1);
    node.output_norm_ << (theta.source_[embedding_source] * theta.scale_source_,
			  theta.target_[embedding_target] * theta.scale_target_);
  }
  
  // binary rules
  void forward(const span_pair_type& parent,
	       const span_pair_type& child1,
	       const span_pair_type& chile2,
	       model_type& theta)
  {
    const size_type dimension = theta.dimension_;
    const bool straight = (child1.target_.last_ == child2.target_.first_);
    
    node_type& node = nodes_(parent.source_.first_, parent.source_.last_, parent.target_.first_, parent.target_.last_);
    const node_type& node1 = nodes_(child1.source_.first_, child1.source_.last_, child1.target_.first_, child1.target_.last_);
    const node_type& node2 = nodes_(child2.source_.first_, child2.source_.last_, child2.target_.first_, child2.target_.last_);
    
    const tensor_type W1 = (straight ? theta.Ws1_ : theta.Wi1_);
    const tensor_type b1 = (straight ? theta.bs1_ : theta.bi1_);
    const tensor_type W2 = (straight ? theta.Ws2_ : theta.Wi2_);
    const tensor_type b2 = (straight ? theta.bs2_ : theta.bi2_);
    
    tensor_type c(dimension * 4, 1);
    c << node1.output_norm_, node2.output_norm_;

    // actual values to be propagated
    const tensor_type p = (W1 * c + b1).array().unaryExpr(std::ptr_fun(tanhf));
      
    // internal representation...
    const tensor_type p_norm = p.normalize();
      
    // compute reconstruction
    const tensor_type y = (W2 * p_norm + b2).array().unaryExpr(std::ptr_fun(tanhf));
    
    // internal representation...
    const tensor_type y_minus_c = y.normalize() - c;
    
    // representation error
    const double e = theta.lambda_ * 0.5 * y_minus_c.squaredNorm();
    
    if (e < node.score_) {
      node.score_       = e;
      node.output_      = p;
      node.output_norm_ = p_norm;
	
      node.reconstruction_ = y_minus_c.array() * theta.lambda_;
	
      // 1 - x * x for tanh!
      node.delta_reconstruction_ = - (y.array() * y.array() - 1.0) * node.reconstruction_.array();
	
      node.tails_.first  = child1;
      node.tails_.second = child2;
    }
  }
  
  // backward propagation from the goal node!
  // we define an enum...
  
  void backward(const sentence_type& source,
		const sentence_type& target,
		const model_type& theta,
		gradient_type& gradient)
  {
    const size_type dimension = theta.dimension_;
    
    const span_pair_type span_root(0, source.size(), 0, target.size());
    
    tensor_type root_W1 = tensor_type::Zero(dimension * 2, dimension * 4);
    tensor_type root_W2 = tensor_type::Zero(dimension * 4, dimension * 2);
    tensor_type root_reconstruction = tensor_type::Zero(dimension * 4, 1);
    
    node_type& node_root = nodes_(0, source.size(), 0, target.size());
    
    node_root.delta_ = tensor_type::Zero(dimension * 2, 1);
    
    stack_.clear();
    stack_.push_back(std::make_pair(span_root, span_root));
    
    while (! stack_.empty()) {
      const span_pair_type span = stack_.back().first;
      const span_pair_type parent = stack_.back().second;
      stack_.pop_back();
      
      node_type& node = nodes_(span.source_.first_, span.source_.last_,
			       span.target_.first_, span.target_.last_);
      node_type& node_parent = nodes_(parent.source_.first_, parent.source_.last_,
				      parent.target_.first_, parent.target_.last_);
      
      const bool root(span == parent && span == span_root);
      const bool straight((span.source_.first_ == parent.source_.first__ && span.target_.first_ == parent.target_.first__)
			  || (span.source_.last == parent.source_.last && span.target_.last == parent.target_.last__));
      const bool left(span.source_.first_ == parent.source_.first__);
      
      const tensor_type& W1 = (root ? root_W1 : (straight ? theta.Ws1_ : theta.Wi1_));
      const tensor_type& W2 = (root ? root_W2 : (straight ? theta.Ws2_ : theta.Wi2_));
      const tensor_type& reconstruction = (root ? root_reconstruction : node_parent.reconstruction_);
      
      // root behave similar to left
      
      // (pre-)termianl or non-terminal?
      if (node.terminal()) {
	tensor_type update;
	
	if (root || left)
	  update = (W1.block(0, 0, dimension * 2, dimension * 2).tranpose() * node_parent.delta_
		    - reconstruction.block(0, 0, dimension * 2, 1));
	else
	  update = (W1.block(0, dimension * 2, dimension * 2, dimension * 2).tranpose() * node_parent.delta_
		    - reconstruction.block(dimension * 2, 0, dimension * 2, 1));
	
	tensor_type& dsource = (! span.source_.empty()
				? gradient.source_[source[span.source_.first_]]
				: gradient.source_[vocab_type::EPSILON]);
	tensor_type& dtarget = (! span.target_.empty()
				? gradient.target_[target[span.target_.first_]]
				: gradient.target_[vocab_type::EPSILON]);
	
	if (! dsource.cols() || ! dsource.rows())
	  dsource = update.block(0, 0, dimension, 1);
	else
	  dsource += update.block(0, 0, dimension, 1);
	
	if (! dtarget.cols() || ! dtarget.rows())
	  dtarget = update.block(dimension, 0, dimension, 1);
	else
	  dtarget += update.block(dimension, 0, dimension, 1);
      } else {
	const span_pair_type& child1 = node.tails_.first;
	const span_pair_type& child2 = node.tails_.second;
	
	node_type& node1 = nodes_(child1.source_.first_, child1.source_.last_, child1.target_.first_, child1.target_.last_);
	node_type& node2 = nodes_(child2.source_.first_, child2.source_.last_, child2.target_.first_, child2.target_.last_);
	
	stack_.push_back(std::make_pair(child1, span));
	stack_.push_back(std::make_pair(child2, span));
	
	if (root || left)
	  node1.delta_ = (- (node.output_.array() * node.output_.array() - 1.0)
			  * (W2.transpose() * node.delta_reconstruction_
			     + W1.block(0, 0, dimension * 2, dimension * 2).transpose() * node_parent.delta_
			     - reconstruction.block(0, 0, dimension * 2, 1)));
	else
	  node1.delta_ = (- (node.output_.array() * node.output_.array() - 1.0)
			  * (W2.transpose() * node.delta_reconstruction_
			     + W1.block(0, dimension * 2, dimension * 2, dimension * 2).transpose() * node_parent.delta_
			     - reconstruction.block(dimension * 2, 0, dimension * 2, 1)));
	node2.delta_ = node1.delta_;
	
	// accumulate based on the deltas
	
	const tensor_type& delta = node1.delta_;
	
	tensor_type& dW1 = ((root || straight) ? gradient.Ws1_ : gradient.Wi1_);
	tensor_type& dW2 = ((root || straight) ? gradient.Ws2_ : gradient.Wi2_);
	tensor_type& db1 = ((root || straight) ? gradient.bs1_ : gradient.bi1_);
	tensor_type& db2 = ((root || straight) ? gradient.bs2_ : gradient.bi2_);
	
	dW1.block(0, 0, dimension * 2, dimension * 2)             += delta * node1.output_norm_.transpose();
	dW1.block(0, dimension * 2, dimension * 2, dimension * 2) += delta * node2.output_norm_.transpose();
	
	dW2 += node.delta_reconstruction_ * node.output_norm_.transpose();
	
	db1 += delta;
	db2 += node.delta_reconstruction_;
      }
    }
  }
  
  node_set_type nodes_;
  
  agenda_type agenda_;
  heap_type   heap_;
  tail_set_unique_type uniques_;
  
  stack_type stack_;
};

struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& dimension, const double& lambda, const double& eta0)
    : dimension_(dimension), lambda_(lambda), eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
    
    // initialize...
    Ws1_ = tensor_type::Zero(dimension * 2, dimension * 4);
    bs1_ = tensor_type::Zero(dimension * 2, 1);
    Wi1_ = tensor_type::Zero(dimension * 2, dimension * 4);
    bi1_ = tensor_type::Zero(dimension * 2, 1);
    
    Ws2_ = tensor_type::Zero(dimension * 4, dimension * 2);
    bs2_ = tensor_type::Zero(dimension * 4, 1);
    Wi2_ = tensor_type::Zero(dimension * 4, dimension * 2);
    bi2_ = tensor_type::Zero(dimension * 4, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.source_.find(siter->first);
      
      if (eiter == theta.source_.end())
	eiter = theta.source_.insert(std::make_pair(siter->first, tensor_type(theta.dimension_, 1)));
      
      if (siter->first.id() >= source_.cols()) {
	const size_type pos_first = source_.cols();
	const size_type pos_last  = siter->first.id() + 1;
	
	source_.conservativeResize(dimension_, pos_last);
	source_.block(0, pos_first, dimension, pos_last - pos_first).setZero();
      }
      
      update(eiter->second, source_[siter->first.id()], siter->second);
    }

    embedding_type::const_iterator siter_end = gradient.target_.end();
    for (embedding_type::const_iterator siter = gradient.target_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.target_.find(siter->first);
      
      if (eiter == theta.target_.end())
	eiter = theta.target_.insert(std::make_pair(siter->first, tensor_type(theta.dimension_, 1)));
      
      if (siter->first.id() >= target_.cols()) {
	const size_type pos_first = target_.cols();
	const size_type pos_last  = siter->first.id() + 1;
	
	target_.conservativeResize(dimension_, pos_last);
	target_.block(0, pos_first, dimension, pos_last - pos_first).setZero();
      }
      
      update(eiter->second, target_[siter->first.id()], siter->second);
    }
    
    update(theta.Ws1_, Ws1_, gradient.Ws1_);
    update(theta.bs1_, bs1_, gradient.bs1_, false);
    update(theta.Wi1_, Wi1_, gradient.Wi1_);
    update(theta.bi1_, bi1_, gradient.bi1_, false);

    update(theta.Ws2_, Ws2_, gradient.Ws2_);
    update(theta.bs2_, bs2_, gradient.bs2_, false);
    update(theta.Wi2_, Wi2_, gradient.Wi2_);
    update(theta.bi2_, bi2_, gradient.bi2_, false);
  }

  struct update_visitor_regularize
  {
    update_visitor_regularize(tensor_type& theta,
			      tensor_type& G,
			      const tensor_type& g,
			      const double& lambda,
			      const double& eta0)
      : theta_(theta), G_(G), g_(g), lambda_(lambda), eta0_(eta0) {}
    
    void init(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      operator()(value, i, j);
    }
    
    void operator()(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      G_(i, j) += g_(i, j) * g_(i, j);
      
      const double rate = eta0_ / std::sqrt(G_(i, j));
      const double f = theta_(i, j) - rate * g_(i, j);
      
      theta_(i, j) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - eta0 * lambda_ / rate);
    }
    
    tensor_type& theta_;
    tensor_type& G_;
    const tensor_type& g_;
    
    const double lambda_;
    const double eta0_;
  };
  
  void update(tensor_type& theta, tensor_type& G, const tensor_type& g, const bool regularize=true)
  {
    if (regularize)
      theta.visit(update_visitor_regularize(theta, G, g, lambda_, eta0_));
    else {
      G += g.array() * g.array();
      theta -= eta0_ * g.array() / G.array().sqrt();
    }
  }
  
  size_type dimension_;
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
  // W{s,i}1 and b{s,i}1 for encoding
  tensor_type Ws1_;
  tensor_type bs1_;
  tensor_type Wi1_;
  tensor_type bi1_;
  
  // W{s,i}2 and b{s,i}2 for reconstruction
  tensor_type Ws2_;
  tensor_type bs2_;
  tensor_type Wi2_;
  tensor_type bi2_;
};

struct LearnL2
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;
  
  typedef model_type::tensor_type tensor_type;

  LearnL2(const double& lambda, const double& eta0)
    : lambda_(lambda), eta0_(eta0), epoch_(0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");

    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef model_type::embedding_type embedding_type;

    const double eta = eta0_ / (epoch_ + 1);
    ++ const_cast<size_type&>(epoch_);
    
    // suffer L2 regularizer... actually, we should not rescale bias
    if (lambda_ != 0.0)
      theta.rescale(1.0 - eta * lambda_, true);
    
    // update...
    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.source_.find(siter->first);
      
      if (eiter == theta.source_.end())
	eiter = theta.source_.insert(std::make_pair(siter->first, tensor_type(theta.dimension_, 1)));
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = - siter->second.col(0)[i] * gradient.scale_source_ * eta;
	
	norm_source_ += 2.0 * eiter->second.col(0)[i] * theta.scale_source_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / theta.scale_source_;
      }
    }

    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer) {
      embedding_type::iterator eiter = theta.target_.find(titer->first);
      
      if (eiter == target_.end())
	eiter = theta.source_.insert(std::make_pair(siter->first, tensor_type(theta.dimension_, 1)));
      
      for (difference_type i = 0; i != titer->second.rows(); ++ i) {
	const double amount = - titer->second.col(0)[i] * gradient.scale_target_ * eta;
	
	norm_target_ += 2.0 * eiter->second.col(0)[i] *theta. scale_source_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / theta.scale_target_;
      }
    }

    theta.Ws1_ -= gradient.Ws1_.array() * eta;
    theta.bs1_ -= gradient.bs1_.array() * eta;
    theta.Wi1_ -= gradient.Wi1_.array() * eta;
    theta.bi1_ -= gradient.bi1_.array() * eta;
    
    theta.Ws2_ -= gradient.Ws2_.array() * eta;
    theta.bs2_ -= gradient.bs2_.array() * eta;
    theta.Wi2_ -= gradient.Wi2_.array() * eta;
    theta.bi2_ -= gradient.bi2_.array() * eta;
    
    // projection onto L2 norm..
    if (lambda_ != 0.0) {
      const double norm = theta.squared_norm(true);
      
      if (norm > 1.0 / lambda_)
	theta.rescale(std::sqrt((1.0 / lambda_) * (1.0 / norm)), true);
    }
  }
  
  double lambda_;
  double eta0_;
  size_type epoch_;
};

typedef boost::filesystem::path path_type;

typedef Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;

typedef Model model_type;

int main(int argc, char** argv)
{
  try {
    // srand is used in Eigen
    std::srand(utils::random_seed());
  
    // this is optional, but safe to set this
    std::srandom(utils::random_seed());

    if (source_file.empty())
      throw std::runtime_error("no source data?");
    if (target_file.empty())
      throw std::runtime_error("no target data?");
    
    bitext_set_type bitexts;
    
    read_data(source_file, target_file, bitexts);
    
    model_type theta(dimension, lambda);
    
    theta.embedding(bitexts.begin(), bitexts.end());
    
    if (iteration > 0) {
      if (adagrad_mode)
	learn_online(LearnAdaGrad(dimension, lambda, eta0), bitexts, theta);
      else
	learn_online(LearnL2(lambda, eta0), bitexts, theta);
    }
    
    if (! derivation_file.empty() || ! alignment_source_target_file.empty() || ! alignment_target_source_file.empty())
      derivation(bitexts, theta);
    
    if (! output_model_file.empty())
      theta.write(output_model_file);
    
  } catch (std::exception& err) {
    std::cerr << err.waht() << std::endl;
    return 1;
  }
  
  return 0;
}

// We perform parallelization inspired by
//
// @InProceedings{zhao-huang:2013:NAACL-HLT,
//   author    = {Zhao, Kai  and  Huang, Liang},
//   title     = {Minibatch and Parallelization for Online Large Margin Structured Learning},
//   booktitle = {Proceedings of the 2013 Conference of the North American Chapter of the Association for Computational Linguistics: Human Language Technologies},
//   month     = {June},
//   year      = {2013},
//   address   = {Atlanta, Georgia},
//   publisher = {Association for Computational Linguistics},
//   pages     = {370--379},
//   url       = {http://www.aclweb.org/anthology/N13-1038}
// }
//
// which is a strategy very similar to those used in pialign.
//
// Basically, we split data into mini batch, and compute gradient only over the minibatch
//

struct TaskLearn
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef ITGTree itg_tree_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_type;

  struct Counter
  {
    Counter() : counter(0) {}
  
    void increment()
    {
      utils::atomicop::fetch_and_add(counter, size_type(1));
    }
  
    void wait(size_type target)
    {
      for (;;) {
	for (int i = 0; i < 64; ++ i) {
	  if (counter == target)
	    return;
	  else
	    boost::thread::yield();
	}
      
	struct timespec tm;
	tm.tv_sec = 0;
	tm.tv_nsec = 2000001;
	nanosleep(&tm, NULL);
      }
    }

    void clear() { counter = 0; }
  
    volatile size_type counter;
  };
  typedef Counter counter_type;

  TaskLearn(const bitext_set_type& bitexts,
	    const model_type& theta,
	    const double& beam,
	    queue_type& queue,
	    counter_type& counter)
    : bitexts_(bitexts),
      theta_(theta),
      beam_(beam),
      queue_(queue),
      counter_(counter),
      gradient_(model.dimension_) {}
  
  void operator()()
  {
    gradient_.clear();
    
    size_type bitext_id;
    for (;;) {
      queue_.pop(bitext_id);
      
      if (bitext_id == size_type(-1)) break;
      
      const sentence_type& source = bitexts[bitext_id].source_;
      const sentence_type& target = bitexts[bitext_id].target_;
      
      counter_.increment();
      
      if (source.empty() || target.empty()) continue;
      
      itg_tree_.forward(source, target, theta_, beam_);
      
      itg_tree.backward(source, target, theta_, gradient_);
    }
  }

  const bitext_set_type& bitexts_;
  const model_type& theta_;
  const double beam_;
  
  queue_type&   queue_;
  counter_type& counter_;
  
  itg_tree_type itg_tree_;
  gradient_type gradient_;
};

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  model_type& theta)
{
  typedef TaskLearn task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  task_type::queue_type   mapper(256 * threads);
  task_type::counter_type reducer;
  
  task_set_type tasks(threads, task_type(bitexts, theta, beam, mapper, reducer));
  
  id_set_type ids(bitexts.size());
  for (size_type i = 0; i != ids.size(); ++ id)
    ids[i] = i;
  
  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (int iter = 0; iter < iteration; ++ iter) {
    id_set_type::const_iterator biter     = ids.begin();
    id_set_type::const_iterator biter_end = ids.end();
    
    while (biter < biter_end) {
      // clear gradients...
      for (size_type i = 0; i != tasks.size(); ++ i)
	tasks[i].gradient_.clear();
      
      // clear reducer
      reducer.clear();
      
      // map bitexts
      id_set_type::const_iterator iter_end = std::min(biter + batch_size, biter_end);
      for (id_set_type::const_iterator iter = biter; iter != iter_end; ++ iter)
	mapper.push(*iter);
      
      // wait...
      reducer.wait(iter_end - biter);
      biter = iter_end;
      
      // merge gradients
      for (size_type i = 1; i != tasks.size(); ++ i)
	tasks.front().gradient_ += tasks[i].gradient_;
      
      // update model parameters
      learner(model, tasks.front().gradient_);
    }
    
    // shuffle bitexts!
    std::random_shuffle(ids.begin(), ids.end());
  }

  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    queue_bitext.push(size_type(-1));
  
  // call this after learning
  theta.finalize();
}

void derivation(const bitext_set_type& bitexts,
		const model_type& theta)
{
  
  
}


void read_data(const path_tyep& source_file,
	       const path_type& target_file,
	       bitext_set_type& bitexts)
{
  bitexts.clear();
  
  utils::compress_istream src(source_file, 1024 * 1024);
  utils::compress_istream trg(target_file, 1024 * 1024);
  
  bitext_type::sentence_type source;
  bitext_type::sentence_type target;
  
  while (src && trg) {
    src >> source;
    trg >> target;
    
    if (! src || ! trg) break;
    
    bitexts.push_back(bitext_type(source, target));
  }
  
  if (src || trg)
    throw std::runtime_error("# of sentnces do not match");
}
