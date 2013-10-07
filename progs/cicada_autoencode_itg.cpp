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

#include <Eigen/Core>

#include "cicada/symbol.hpp"

struct Model
{
  // model parameters
  
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;
  
  Model() {}
  Model(const size_type dimension) { initialize(dimension); }
  
  void initialize(const size_type dimension)
  {
    // intialize randomly...
    Ws1_ = tensor_type::Random(dimension, dimension * 2);
    bs1_ = tensor_type::Random(dimension, 1);
    Wi1_ = tensor_type::Random(dimension, dimension * 2);
    bi1_ = tensor_type::Random(dimension, 1);
    
    Ws2_ = tensor_type::Random(dimension * 2, dimension);
    bs2_ = tensor_type::Random(dimension * 2, 1);
    Wi2_ = tensor_type::Random(dimension * 2, dimension);
    bi2_ = tensor_type::Random(dimension * 2, 1);
  }
  
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
  
  typedef model_type::tensor_type tensor_type;
  
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
  
  struct Edge
  {
    tensor_type y1_;
    tensor_type y2_;
    tensor_type output1_;
    tensor_type output2_;
    span_pair_type child1_;
    span_pair_type child2_;
    
    bool straight() const { return child1_.target.last_ == child2_.target.first_; }
    bool inverted() const { return child2_.target.last_ == child1_.target.first_; }
  };
  
  typedef Edge edge_type;
  typedef utils::chunk_vector<edge_type, 4096 / sizeof(edge_type), std::allocator<edge_type> > edge_set_type;
  
  
};
