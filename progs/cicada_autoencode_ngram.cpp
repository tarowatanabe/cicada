//
//  Copyright(C) 2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

// 
// Lexical autoencoder
//
//  reconstruction
//    internal ------ sampled loss
//   embedding
//
// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by Pegasos
// SGD with L2 regularizer inspired by AdaGrad (default)
// SGD with L2/L2 regularizer from RDA (TODO)
//
// + batch algorithm using LBFGS (TODO)
//

#include <cstdlib>
#include <cmath>
#include <climits>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <set>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>

#include <Eigen/Core>

#include "cicada/symbol.hpp"
#include "cicada/sentence.hpp"
#include "cicada/vocab.hpp"
#include "cicada/alignment.hpp"

#include "utils/alloc_vector.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bichart.hpp"
#include "utils/bithack.hpp"
#include "utils/compact_map.hpp"
#include "utils/compact_set.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/mathop.hpp"
#include "utils/unordered_map.hpp"
#include "utils/repository.hpp"
#include "utils/program_options.hpp"
#include "utils/random_seed.hpp"
#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/vector2.hpp"
#include "utils/sampler.hpp"

#include <boost/random.hpp>
#include <boost/thread.hpp>
#include <boost/math/special_functions/expm1.hpp>
#include <boost/math/special_functions/log1p.hpp>


struct Model
{
  // model parameters
  
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  // we use float for the future compatibility with GPU :)
  typedef float parameter_type;
  typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic> tensor_type;

  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef utils::unordered_map<word_type, tensor_type,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, tensor_type> > >::type embedding_type;

  typedef boost::filesystem::path path_type;
  
  Model() : dimension_(0), order_(0), alpha_(0) {}
  Model(const size_type& dimension, const size_type& order) 
    : dimension_(dimension), order_(order), alpha_(0) { initialize(dimension, order); }
  Model(const size_type& dimension, const size_type& order, const double& alpha) 
    : dimension_(dimension), order_(order), alpha_(alpha) { initialize(dimension, order); }
  
  Model& operator-=(const Model& x)
  {
    embedding_type::const_iterator siter_end = x.embedding_.end();
    for (embedding_type::const_iterator siter = x.embedding_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = embedding_.find(siter->first);
      
      if (eiter == embedding_.end())
	eiter = embedding_.insert(std::make_pair(siter->first, tensor_type::Zero(dimension_, 1))).first;
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = - siter->second.col(0)[i] * x.scale_;
	
	norm_ += 2.0 * eiter->second.col(0)[i] * scale_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / scale_;
      }
    }
        
    Wl1_ -= x.Wl1_;
    bl1_ -= x.bl1_;
    Wl2_ -= x.Wl2_;
    bl2_ -= x.bl2_;
    
    Wc_ -= x.Wc_;
    bc_ -= x.bc_;
    
    return *this;
  }
  
  Model& operator+=(const Model& x)
  {
    embedding_type::const_iterator siter_end = x.embedding_.end();
    for (embedding_type::const_iterator siter = x.embedding_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = embedding_.find(siter->first);
      
      if (eiter == embedding_.end())
	eiter = embedding_.insert(std::make_pair(siter->first, tensor_type::Zero(dimension_, 1))).first;
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = siter->second.col(0)[i] * x.scale_;
	
	norm_ += 2.0 * eiter->second.col(0)[i] * scale_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / scale_;
      }
    }
    
    Wl1_ += x.Wl1_;
    bl1_ += x.bl1_;
    Wl2_ += x.Wl2_;
    bl2_ += x.bl2_;
    
    Wc_ += x.Wc_;
    bc_ += x.bc_;
    
    return *this;
  }
  
  void clear()
  {
    // embedding
    embedding_.clear();
    scale_ = 1.0;
    norm_ = 0.0;
    
    // matrices
    Wl1_.setZero();
    bl1_.setZero();
    Wl2_.setZero();
    bl2_.setZero();    

    Wc_.setZero();
    bc_.setZero();
  }

  void finalize()
  {
    if (scale_ != 1.0) {
      norm_ = 0.0;
      
      embedding_type::iterator siter_end = embedding_.end();
      for (embedding_type::iterator siter = embedding_.begin(); siter != siter_end; ++ siter) {
	siter->second *= scale_;
	norm_ += siter->second.squaredNorm();
      }
      
      scale_ = 1.0;
    }
  }
  
  void rescale(const double scaling, const bool ignore_bias)
  {
    scale_ *= scaling;

    if (scale_ == 0.0) {
      embedding_type::iterator siter_end = embedding_.end();
      for (embedding_type::iterator siter = embedding_.begin(); siter != siter_end; ++ siter)
	siter->second.setZero();
      
      scale_ = 1.0;
      norm_  = 0.0;
    } else
      norm_ *= scaling * scaling;

    
    Wl1_ *= scaling;
    Wl2_ *= scaling;

    Wc_ *= scaling;
    
    if (! ignore_bias) {
      bl1_ *= scaling;
      bl2_ *= scaling;

      bc_ *= scaling;
    }
  }

  double squared_norm(bool ignore_bias) const
  {
    double norm = norm_;
    
    norm += Wl1_.squaredNorm();
    norm += Wl2_.squaredNorm();

    norm += Wc_.squaredNorm();
    
    if (! ignore_bias) {
      norm += bl1_.squaredNorm();
      norm += bl2_.squaredNorm();

      norm += bc_.squaredNorm();
    }
    
    return norm;
  }
  
  void initialize(const size_type dimension, const size_type order)
  {
    if (dimension <= 0)
      throw std::runtime_error("invalid dimension");
    if (order <= 0)
      throw std::runtime_error("invalid order");
    
    // intialize randomly...
    dimension_ = dimension;
    order_     = order;
        
    // embedding
    embedding_.clear();
    scale_ = 1.0;
    norm_ = 0.0;
        
    // lexicon
    Wl1_ = tensor_type::Random(dimension, dimension * order);
    bl1_ = tensor_type::Random(dimension, 1);
    
    // lexicon reconstruction
    Wl2_ = tensor_type::Random(dimension * order, dimension);
    bl2_ = tensor_type::Random(dimension * order, 1);
    
    // classification
    Wc_ = tensor_type::Random(1, dimension);
    bc_ = tensor_type::Random(1, 1);
  }
  
  template <typename Iterator>
  void embedding(Iterator first, Iterator last)
  {
    tensor_type& bos     = embedding_[vocab_type::BOS];
    tensor_type& eos     = embedding_[vocab_type::EOS];
    
    if (! bos.rows() || ! bos.cols())
      bos = tensor_type::Random(dimension_, 1).normalized();
    if (! eos.rows() || ! eos.cols())
      eos = tensor_type::Random(dimension_, 1).normalized();
    
    for (/**/; first != last; ++ first) {
      tensor_type& we = embedding_[*first];

      if (! we.rows() || ! we.cols())
	we = tensor_type::Random(dimension_, 1).normalized();
    }
    
    // compute norm...
    scale_ = 1.0;
    norm_ = 0.0;
    
    embedding_type::const_iterator siter_end = embedding_.end();
    for (embedding_type::const_iterator siter = embedding_.begin(); siter != siter_end; ++ siter)
      norm_ += siter->second.squaredNorm();
  }
  
  struct real_policy : boost::spirit::karma::real_policies<parameter_type>
  {
    static unsigned int precision(parameter_type)
    {
      return 10;
    }
  };

  boost::spirit::karma::real_generator<parameter_type, real_policy> float10;


  void read_embedding(const path_type& path)
  {
    // we will overwrite embedding... thus we will not clear embedding_
    
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;
    typedef boost::fusion::tuple<std::string, parameter_set_type > embedding_parsed_type;
    typedef boost::spirit::istream_iterator iterator_type;
    
    qi::rule<iterator_type, std::string(), standard::blank_type>           word;
    qi::rule<iterator_type, embedding_parsed_type(), standard::blank_type> parser; 
    
    word   %= qi::lexeme[+(standard::char_ - standard::space)];
    parser %= word >> *qi::double_ >> (qi::eol | qi::eoi);
    
    utils::compress_istream is(path, 1024 * 1024);
    is.unsetf(std::ios::skipws);
    
    iterator_type iter(is);
    iterator_type iter_end;
    
    embedding_parsed_type parsed;
    
    while (iter != iter_end) {
      boost::fusion::get<0>(parsed).clear();
      boost::fusion::get<1>(parsed).clear();
      
      if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, parsed))
	if (iter != iter_end)
	  throw std::runtime_error("embedding parsing failed");
      
      if (boost::fusion::get<1>(parsed).size() != dimension_)
	throw std::runtime_error("invalid embedding size");
      
      embedding_[boost::fusion::get<0>(parsed)] = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), dimension_, 1);
    }
  }
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["dimension"] = utils::lexical_cast<std::string>(dimension_);
    rep["order"]     = utils::lexical_cast<std::string>(order_);
    rep["alpha"]     = utils::lexical_cast<std::string>(alpha_);
    
    const path_type embedding_file = rep.path("embedding.gz");
    
    utils::compress_ostream os(embedding_file, 1024 * 1024);
    std::ostream_iterator<char> iter(os);
    
    embedding_type::const_iterator siter_end = embedding_.end();
    for (embedding_type::const_iterator siter = embedding_.begin(); siter != siter_end; ++ siter) {
      karma::generate(iter, standard::string, siter->first);
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i)
	karma::generate(iter, karma::lit(' ') << float10, siter->second.row(i)[0]);
      
      karma::generate(iter, karma::lit('\n'));
    }
    
    // dump matrices...
    write(rep.path("Wl1.txt.gz"), rep.path("Wl1.bin"), Wl1_);
    write(rep.path("bl1.txt.gz"), rep.path("bl1.bin"), bl1_);

    write(rep.path("Wl2.txt.gz"), rep.path("Wl2.bin"), Wl2_);
    write(rep.path("bl2.txt.gz"), rep.path("bl2.bin"), bl2_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
  }

  void write(const path_type& path_text, const path_type& path_binary, const tensor_type& matrix) const
  {
    {
      utils::compress_ostream os(path_text, 1024 * 1024);
      os.precision(10);
      os << matrix;
    }
    
    {
      utils::compress_ostream os(path_binary, 1024 * 1024);
      
      const tensor_type::Index rows = matrix.rows();
      const tensor_type::Index cols = matrix.cols();
      
      os.write((char*) &rows, sizeof(tensor_type::Index));
      os.write((char*) &cols, sizeof(tensor_type::Index));
      os.write((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
    }
  }
  
  // dimension...
  size_type dimension_;
  size_type order_;
  
  // hyperparameter
  double alpha_;
  
  // Embedding
  embedding_type embedding_;
  double scale_;
  double norm_;
  
  // Wl1 and bl1 for encoding
  tensor_type Wl1_;
  tensor_type bl1_;
  
  // Wl2 and bl2 for reconstruction
  tensor_type Wl2_;
  tensor_type bl2_;  

  // Wc and bc for classification
  tensor_type Wc_;
  tensor_type bc_;
};


struct NGram
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type tensor_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;

  typedef std::vector<word_type, std::allocator<word_type> > word_set_type;
  
  template <typename Gen>
  std::pair<double, double> operator()(const sentence_type& sentence,
				       const word_set_type& words,
				       const model_type& theta,
				       gradient_type& gradient,
				       Gen& gen)
  {
    typedef model_type::embedding_type embedding_type;
    
    const size_type sentence_size = sentence.size();
    
    const size_type dimension = theta.dimension_;
    const size_type order     = theta.order_;
    
    double error = 0.0;
    double error_classification = 0.0;
    
    tensor_type input(dimension * order, 1);
    tensor_type input_prev(dimension * order, 1);
    tensor_type input_sampled(dimension * order, 1);
    
    // fill with BOS...
    embedding_type::const_iterator biter = theta.embedding_.find(vocab_type::BOS);
    if (biter == theta.embedding_.end())
      throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(vocab_type::BOS));
    
    if (biter->second.rows() != dimension)
      throw std::runtime_error("dimension does not match");
    
    for (size_type i = 0; i != order; ++ i)
      input_prev.block(dimension * i, 0, dimension, 1) = biter->second * theta.scale_;
    
    for (size_type pos = 0; pos <= sentence_size; ++ pos) {
      // copy previous input as a history
      input.block(0, 0, dimension * (order - 1), 1) = input_prev.block(dimension, 0, dimension * (order - 1), 1);
      input_sampled = input;
      
      const word_type& word = (pos == sentence_size ? vocab_type::EOS : sentence[pos]);
      
      const size_type pos_sampled = boost::random::uniform_int_distribution<size_t>(0, words.size() - 1)(gen);
      
      const word_type& word_sampled = words[pos_sampled];
      
      embedding_type::const_iterator piter = theta.embedding_.find(word);
      if (piter == theta.embedding_.end())
	throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(word));
      
      if (piter->second.rows() != dimension)
	throw std::runtime_error("wrong dimension");

      embedding_type::const_iterator miter = theta.embedding_.find(word_sampled);
      if (miter == theta.embedding_.end())
	throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(word_sampled));
      
      if (miter->second.rows() != dimension)
	throw std::runtime_error("wrong dimension");
      
      input.block(dimension * (order - 1), 0, dimension, 1)         = piter->second * theta.scale_;
      input_sampled.block(dimension * (order - 1), 0, dimension, 1) = miter->second * theta.scale_;
      
      const tensor_type p = (theta.Wl1_ * input + theta.bl1_).array().unaryExpr(std::ptr_fun(tanhf));
      const tensor_type p_norm = p.normalized();
      const tensor_type y = (theta.Wl2_ * p_norm + theta.bl2_).array().unaryExpr(std::ptr_fun(tanhf));
      const tensor_type y_minus_c = y.normalized() - input;
      
      const tensor_type p_sampled = (theta.Wl1_ * input_sampled + theta.bl1_).array().unaryExpr(std::ptr_fun(tanhf));
      const tensor_type p_sampled_norm = p_sampled.normalized();
      
      const double e = theta.alpha_ * 0.5 * y_minus_c.squaredNorm();
      
      const tensor_type reconstruction       = y_minus_c.array() * theta.alpha_;
      const tensor_type delta_reconstruction = - (y.array() * y.array() - 1.0) * reconstruction.array();
      
      const double y_p = std::tanh((theta.Wc_ * p_norm + theta.bc_)(0,0));
      const double y_m = std::tanh((theta.Wc_ * p_sampled_norm + theta.bc_)(0,0));
      
      const double e_classification = std::max(1.0 - (y_p - y_m), 0.0);
      
      const double delta_classification_p = - (1.0 - y_p * y_p) * (e_classification > 0.0);
      const double delta_classification_m =   (1.0 - y_m * y_m) * (e_classification > 0.0);
      
      // update error...
      error                += e;
      error_classification += e_classification;
	
      // backward...
      
      const tensor_type delta = (- (p.array() * p.array() - 1.0)
				 * (theta.Wl2_.transpose() * delta_reconstruction
				    + theta.Wc_.transpose() * delta_classification_p).array());
      
      const tensor_type delta_sampled = (- (p_sampled.array() * p_sampled.array() - 1.0)
					 * (theta.Wc_.transpose() * delta_classification_m).array());
      
      
      gradient.Wl1_ += delta * input.transpose();
      gradient.bl1_ += delta;
      
      gradient.Wl1_ += delta_sampled * input_sampled.transpose();
      gradient.bl1_ += delta_sampled;
      
      gradient.Wl2_ += delta_reconstruction * p_norm.transpose();
      gradient.bl2_ += delta_reconstruction;
      
      gradient.Wc_         += delta_classification_p * p_norm.transpose();
      gradient.bc_.array() += delta_classification_p;
      
      gradient.Wc_         += delta_classification_m * p_sampled_norm.transpose();
      gradient.bc_.array() += delta_classification_m;
      
      // update embedding
      const tensor_type delta_embedding_p = theta.Wl1_.transpose() * delta - reconstruction;
      const tensor_type delta_embedding_m = theta.Wl1_.transpose() * delta_sampled;
      
      for (size_type i = 0; i != order; ++ i) {
	const difference_type shift = difference_type(i) - order + 1;
	
	if (shift) {
	  const word_type& word = (difference_type(pos) + shift < 0 ? vocab_type::BOS : sentence[pos + shift]);
	  
	  tensor_type& dembedding = gradient.embedding_[word];
	  
	  if (! dembedding.cols() || ! dembedding.rows())
	    dembedding = tensor_type::Zero(dimension, 1);
	  
	  dembedding += delta_embedding_p.block(dimension * i, 0, dimension, 1);
	  dembedding += delta_embedding_m.block(dimension * i, 0, dimension, 1);
	} else {
	  tensor_type& dembedding_p = gradient.embedding_[word];
	  tensor_type& dembedding_m = gradient.embedding_[word_sampled];

	  if (! dembedding_p.cols() || ! dembedding_p.rows())
	    dembedding_p = tensor_type::Zero(dimension, 1);

	  if (! dembedding_m.cols() || ! dembedding_m.rows())
	    dembedding_m = tensor_type::Zero(dimension, 1);
	  
	  dembedding_p += delta_embedding_p.block(dimension * i, 0, dimension, 1);
	  dembedding_m += delta_embedding_m.block(dimension * i, 0, dimension, 1);
	}
      }
      
      input_prev.swap(input);
    }
    
    return std::make_pair(error, error_classification);
  }
};

struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& dimension, const size_type& order, const double& lambda, const double& eta0)
    : dimension_(dimension), order_(order), lambda_(lambda), eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
    
    // initialize...
    Wl1_ = tensor_type::Zero(dimension, dimension * order);
    bl1_ = tensor_type::Zero(dimension, 1);
    
    Wl2_ = tensor_type::Zero(dimension * order, dimension);
    bl2_ = tensor_type::Zero(dimension * order, 1);
    
    Wc_ = tensor_type::Zero(1, dimension);
    bc_ = tensor_type::Zero(1, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef model_type::embedding_type embedding_type;

    embedding_type::const_iterator siter_end = gradient.embedding_.end();
    for (embedding_type::const_iterator siter = gradient.embedding_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.embedding_.find(siter->first);
      
      if (eiter == theta.embedding_.end()) {
	std::cerr << "WARNING: this should not happen: embedding: " << siter->first << std::endl;
	eiter = theta.embedding_.insert(std::make_pair(siter->first, tensor_type::Zero(theta.dimension_, 1))).first;
      }
      
      if (siter->first.id() >= embedding_.cols()) {
	const size_type pos_first = embedding_.cols();
	const size_type pos_last  = siter->first.id() + 1;

	tensor_type& matrix = const_cast<tensor_type&>(embedding_);
	
	matrix.conservativeResize(dimension_, pos_last);
	matrix.block(0, pos_first, dimension_, pos_last - pos_first).setZero();
      }
      
      update(eiter->second, embedding_.col(siter->first.id()), siter->second, lambda_ != 0.0);
    }
    
    update(theta.Wl1_, Wl1_, gradient.Wl1_, lambda_ != 0.0);
    update(theta.bl1_, bl1_, gradient.bl1_, false);

    update(theta.Wl2_, Wl2_, gradient.Wl2_, lambda_ != 0.0);
    update(theta.bl2_, bl2_, gradient.bl2_, false);

    update(theta.Wc_, Wc_, gradient.Wc_, lambda_ != 0.0);
    update(theta.bc_, bc_, gradient.bc_, false);
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
      
      theta_(i, j) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
    }
    
    tensor_type& theta_;
    tensor_type& G_;
    const tensor_type& g_;
    
    const double lambda_;
    const double eta0_;
  };
  
  void update(tensor_type& theta, const tensor_type& __G, const tensor_type& g, const bool regularize=true) const
  {
    tensor_type& G = const_cast<tensor_type&>(__G);

    if (regularize) {
      update_visitor_regularize visitor(theta, G, g, lambda_, eta0_);
      
      theta.visit(visitor);
    } else {
      G.array() += g.array() * g.array();
      theta.array() -= eta0_ * g.array() / G.array().sqrt();
    }
  }
  
  size_type dimension_;
  size_type order_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type embedding_;
  
  // Wl1 and bl1 for encoding
  tensor_type Wl1_;
  tensor_type bl1_;
  
  // Wl2 and bl2 for reconstruction
  tensor_type Wl2_;
  tensor_type bl2_;

  // Wc and bc for classification
  tensor_type Wc_;
  tensor_type bc_;
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
    if (lambda_ <= 0.0)
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
    embedding_type::const_iterator siter_end = gradient.embedding_.end();
    for (embedding_type::const_iterator siter = gradient.embedding_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.embedding_.find(siter->first);
      
      if (eiter == theta.embedding_.end()) {
	std::cerr << "WARNING: this should not happen: " << siter->first << std::endl;
	eiter = theta.embedding_.insert(std::make_pair(siter->first, tensor_type::Zero(theta.dimension_, 1))).first;
      }
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = - siter->second.col(0)[i] * gradient.scale_ * eta;
	
	theta.norm_ += 2.0 * eiter->second.col(0)[i] * theta.scale_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / theta.scale_;
      }
    }


    theta.Wl1_.array() -= gradient.Wl1_.array() * eta;
    theta.bl1_.array() -= gradient.bl1_.array() * eta;

    theta.Wl2_.array() -= gradient.Wl2_.array() * eta;
    theta.bl2_.array() -= gradient.bl2_.array() * eta;
    
    theta.Wc_.array() -= gradient.Wc_.array() * eta;
    theta.bc_.array() -= gradient.bc_.array() * eta;
    
    // projection onto L2 norm..
    if (lambda_ != 0.0) {
      const double norm = theta.squared_norm(true);
      
      if (norm > 1.0 / lambda_)
	theta.rescale(std::sqrt(1.0 / lambda_) * std::sqrt(1.0 / norm), true);
    }
  }
  
  double lambda_;
  double eta0_;
  size_type epoch_;
};

typedef boost::filesystem::path path_type;

typedef cicada::Sentence sentence_type;
typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
typedef NGram::word_set_type word_set_type;

typedef Model model_type;

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type input_file;
path_type embedding_file;
path_type output_model_file;

double alpha = 0.001;
int dimension = 16;
int order = 3;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 1024;
double beam = 0.1;
double lambda = 1;
double eta0 = 1;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const sentence_set_type& sentences,
		  const word_set_type& words,
		  model_type& theta);
void read_data(const path_type& input_file,
	       sentence_set_type& sentences,
	       word_set_type& words);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (dimension <= 0)
      throw std::runtime_error("dimension must be positive");
    
    if (order <= 0)
      throw std::runtime_error("order size should be positive");
    
    if (alpha < 0.0 || alpha > 1.0)
      throw std::runtime_error("error factor should be: 0.0 <= alpha <= 1.0");

    
    if (int(optimize_sgd) + optimize_adagrad > 1)
      throw std::runtime_error("either one of optimize-{sgd,adagrad}");
    
    if (int(optimize_sgd) + optimize_adagrad == 0)
      optimize_adagrad = true;
    
    threads = utils::bithack::max(threads, 1);
    
    // srand is used in Eigen
    std::srand(utils::random_seed());
  
    // this is optional, but safe to set this
    ::srandom(utils::random_seed());
        
    if (input_file.empty())
      throw std::runtime_error("no data?");
    
    sentence_set_type sentences;
    word_set_type words;
    
    read_data(input_file, sentences, words);
    
    model_type theta(dimension, order, alpha);

    if (! embedding_file.empty()) {
      if (embedding_file != "-" && ! boost::filesystem::exists(embedding_file))
	throw std::runtime_error("no embedding: " + embedding_file.string());
      
      theta.read_embedding(embedding_file);
    }
    
    theta.embedding(words.begin(), words.end());
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension, order, lambda, eta0), sentences, words, theta);
      else
	learn_online(LearnL2(lambda, eta0), sentences, words, theta);
    }
    
    if (! output_model_file.empty())
      theta.write(output_model_file);
    
  } catch (std::exception& err) {
    std::cerr << err.what() << std::endl;
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


struct TaskAccumulate
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef NGram ngram_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::Symbol   word_type;
  typedef cicada::Vocab    vocab_type;
  
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

  TaskAccumulate(const sentence_set_type& sentences,
		 const word_set_type& words,
		 const model_type& theta,
		 const double& beam,
		 queue_type& queue,
		 counter_type& counter)
    : sentences_(sentences),
      words_(words),
      theta_(theta),
      beam_(beam),
      queue_(queue),
      counter_(counter),
      gradient_(theta.dimension_, theta.order_),
      error_(0),
      classification_(0),
      samples_(0) {}
  
  void operator()()
  {
    clear();

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    size_type sentence_id;
    for (;;) {
      queue_.pop(sentence_id);
      
      if (sentence_id == size_type(-1)) break;
      
      const sentence_type& sentence = sentences_[sentence_id];
      
      if (! sentence.empty()) {
	std::pair<double, double> errors = ngram_(sentence, words_, theta_, gradient_, generator);
	  
	error_          += errors.first;
	classification_ += errors.second;
	++ samples_;
      }
      
      counter_.increment();
    }
  }

  void clear()
  {
    gradient_.clear();
    error_ = 0.0;
    classification_ = 0.0;
    samples_ = 0;
  }

  const sentence_set_type& sentences_;
  const word_set_type& words_;
  const model_type& theta_;
  const double beam_;
  
  queue_type&            queue_;
  counter_type&          counter_;
  
  ngram_type ngram_;

  gradient_type gradient_;
  double        error_;
  double        classification_;
  size_type     samples_;
};

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

template <typename Learner>
void learn_online(const Learner& learner,
		  const sentence_set_type& sentences,
		  const word_set_type& words,
		  model_type& theta)
{
  typedef TaskAccumulate task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  task_type::queue_type   mapper(256 * threads);
  task_type::counter_type reducer;
  
  task_set_type tasks(threads, task_type(sentences,
					 words,
					 theta,
					 beam,
					 mapper,
					 reducer));
  
  id_set_type ids(sentences.size());
  for (size_type i = 0; i != ids.size(); ++ i)
    ids[i] = i;
  
  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;

    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);

    id_set_type::const_iterator biter     = ids.begin();
    id_set_type::const_iterator biter_end = ids.end();

    double error = 0.0;
    double classification = 0.0;
    size_type samples = 0;
    size_type num_bitext = 0;
    
    while (biter < biter_end) {
      // clear gradients...
      for (size_type i = 0; i != tasks.size(); ++ i)
	tasks[i].clear();
      
      // clear reducer
      reducer.clear();
      
      // map bitexts
      id_set_type::const_iterator iter_end = std::min(biter + batch_size, biter_end);
      for (id_set_type::const_iterator iter = biter; iter != iter_end; ++ iter) {
	mapper.push(*iter);
	
	++ num_bitext;
	if (debug) {
	  if (num_bitext % DEBUG_DOT == 0)
	    std::cerr << '.';
	  if (num_bitext % DEBUG_LINE == 0)
	    std::cerr << '\n';
	}
      }
      
      // wait...
      reducer.wait(iter_end - biter);
      biter = iter_end;
      
      // merge gradients
      error          += tasks.front().error_;
      classification += tasks.front().classification_;
      samples        += tasks.front().samples_;
      for (size_type i = 1; i != tasks.size(); ++ i) {
	tasks.front().gradient_ += tasks[i].gradient_;
	error          += tasks[i].error_;
	classification += tasks[i].classification_;
	samples        += tasks[i].samples_;
      }
      
      // update model parameters
      learner(theta, tasks.front().gradient_);
    }
    
    if (debug && ((num_bitext / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;
    if (debug)
      std::cerr << "# of bitexts: " << num_bitext << std::endl;

    if (debug)
      std::cerr << "reconstruction error: " << (error / samples) << std::endl
		<< "classification error: " << (classification / samples) << std::endl
		<< "parsed: " << samples << std::endl;
    
    // shuffle bitexts!
    std::random_shuffle(ids.begin(), ids.end());
  }

  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    mapper.push(size_type(-1));

  workers.join_all();
  
  // call this after learning
  theta.finalize();
}


void read_data(const path_type& input_file,
	       sentence_set_type& sentences,
	       word_set_type& words)
{
  typedef utils::compact_set<sentence_type::word_type,
			     utils::unassigned<sentence_type::word_type>, utils::unassigned<sentence_type::word_type>,
			     boost::hash<sentence_type::word_type>, std::equal_to<sentence_type::word_type>,
			     std::allocator<sentence_type::word_type> > unique_set_type;
  
  sentences.clear();

  unique_set_type uniques;
  
  utils::compress_istream is(input_file, 1024 * 1024);
  
  sentence_type sentence;
  
  while (is >> sentence) {
    sentences.push_back(sentence);
    
    uniques.insert(sentence.begin(), sentence.end());
  }
  
  words.clear();
  words.reserve(uniques.size());
  words.insert(words.end(), uniques.begin(), uniques.end());
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("input",    po::value<path_type>(&input_file),    "input file")
    
    ("embedding", po::value<path_type>(&embedding_file), "initial embedding")
    
    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("alpha",     po::value<double>(&alpha)->default_value(alpha),      "model parameter")
    ("dimension", po::value<int>(&dimension)->default_value(dimension), "dimension")
    ("order",     po::value<int>(&order)->default_value(order),         "context order size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("beam",              po::value<double>(&beam)->default_value(beam),          "beam width for parsing")
    ("lambda",            po::value<double>(&lambda)->default_value(lambda),      "regularization constant")
    ("eta0",              po::value<double>(&eta0)->default_value(eta0),          "\\eta_0 for decay")

    ("threads", po::value<int>(&threads), "# of threads")
    
    ("debug", po::value<int>(&debug)->implicit_value(1), "debug level")
    ("help", "help message");
  
  po::options_description desc_command;
  desc_command.add(opts_command);
  
  po::variables_map variables;
  po::store(po::parse_command_line(argc, argv, desc_command, po::command_line_style::unix_style & (~po::command_line_style::allow_guessing)), variables);
  
  po::notify(variables);

  if (variables.count("help")) {
    std::cout << argv[0] << " [options] [operations]\n"
	      << opts_command << std::endl;
    exit(0);
  }
}
