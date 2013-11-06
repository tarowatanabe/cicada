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
// < [ ((( word |||  word ))) ] <  ((( word ||| word ))) > >
//
// [ ] indicates straight and < > indicates inversion.

// we will try four learning algorithms:
//
// SGD with L2 regularizer inspired by Pegasos (default)
// SGD with L2 regularizer inspired by AdaGrad
// SGD with L2/L2 regularizer from RDA
//
// + batch algorithm using LBFGS
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
  void swap(Bitext& x, Bitext& y)
  {
    x.swap(y);
  }
};

namespace nn
{
  struct sigmoid
  {
    double operator()(const double& x) const
    {
      const double expx = std::exp(- x);
      return (expx == std::numeric_limits<double>::infinity() ? 0.0 : 1.0 / (expx + 1.0));
    }
  };

  struct dsigmoid
  {
    double operator()(const double& x) const
    {
      const double expx = std::exp(- x);

      if (expx == std::numeric_limits<double>::infinity())
	return 0.0;
      else {
	const double m = 1.0 / (expx + 1.0);
	return m * (1.0 - m);
      }
    }
  };
  
  struct tanh
  {
    double operator()(const double& x) const
    {
      return std::tanh(x);
    }
  };

  struct dtanh
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return Tp(1) - x * x;
    }
  };

  struct htanh
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return std::min(std::max(x, Tp(- 1)), Tp(1));
    }
  };
  
  struct dhtanh
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return Tp(- 1) <= x && x <= Tp(1);
    }
  };
  
  struct hinge
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return std::max(x, Tp(0));
    }
  };
  
  struct dhinge
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return x >= Tp(0);
    }
  };

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
  
  Model() : dimension_embedding_(0), dimension_hidden_(0), dimension_itg_(0),
	    window_(0), alpha_(0), beta_(0) {}
  Model(const size_type& dimension_embedding,
	const size_type& dimension_hidden,
	const size_type& dimension_itg,
	const size_type& window) 
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden), dimension_itg_(dimension_itg),
      window_(window), alpha_(0), beta_(0) { initialize(dimension_embedding, dimension_hidden, dimension_itg, window); }
  template <typename Gen>
  Model(const size_type& dimension_embedding,
	const size_type& dimension_hidden,
	const size_type& dimension_itg,
	const size_type& window, const double& alpha, const double& beta, Gen& gen) 
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden), dimension_itg_(dimension_itg),
      window_(window), alpha_(alpha), beta_(beta) { initialize(dimension_embedding, dimension_hidden, dimension_itg, window, gen); }
  
  Model& operator-=(const Model& x)
  {
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = source_.find(siter->first);
      
      if (eiter == source_.end())
	eiter = source_.insert(std::make_pair(siter->first, tensor_type::Zero(dimension_embedding_, 1))).first;
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = - siter->second.col(0)[i] * x.scale_source_;
	
	norm_source_ += 2.0 * eiter->second.col(0)[i] * scale_source_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / scale_source_;
      }
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      embedding_type::iterator eiter = target_.find(titer->first);
      
      if (eiter == target_.end())
	eiter = target_.insert(std::make_pair(titer->first, tensor_type::Zero(dimension_embedding_, 1))).first;
      
      for (difference_type i = 0; i != titer->second.rows(); ++ i) {
	const double amount = - titer->second.col(0)[i] * x.scale_target_;
	
	norm_target_ += 2.0 * eiter->second.col(0)[i] * scale_target_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / scale_target_;
      }
    }
    
    Ws1_ -= x.Ws1_;
    bs1_ -= x.bs1_;
    Wi1_ -= x.Wi1_;
    bi1_ -= x.bi1_;

    Ws2_ -= x.Ws2_;
    bs2_ -= x.bs2_;
    Wi2_ -= x.Wi2_;
    bi2_ -= x.bi2_;

    Wp1_ -= x.Wp1_;
    bp1_ -= x.bp1_;
    Wp2_ -= x.Wp2_;
    bp2_ -= x.bp2_;
    
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
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = source_.find(siter->first);
      
      if (eiter == source_.end())
	eiter = source_.insert(std::make_pair(siter->first, tensor_type::Zero(dimension_embedding_, 1))).first;
      
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
	eiter = target_.insert(std::make_pair(titer->first, tensor_type::Zero(dimension_embedding_, 1))).first;
      
      for (difference_type i = 0; i != titer->second.rows(); ++ i) {
	const double amount = titer->second.col(0)[i] * x.scale_target_;
	
	norm_target_ += 2.0 * eiter->second.col(0)[i] * scale_target_ * amount + amount * amount;
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
    
    Wp1_ += x.Wp1_;
    bp1_ += x.bp1_;
    Wp2_ += x.Wp2_;
    bp2_ += x.bp2_;

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

    Wp1_.setZero();
    bp1_.setZero();
    Wp2_.setZero();
    bp2_.setZero();    

    Wl1_.setZero();
    bl1_.setZero();
    Wl2_.setZero();
    bl2_.setZero();    
    
    Wc_.setZero();
    bc_.setZero();
  }

  void finalize()
  {
    if (scale_source_ != 1.0) {
      norm_source_ = 0.0;
      
      embedding_type::iterator siter_end = source_.end();
      for (embedding_type::iterator siter = source_.begin(); siter != siter_end; ++ siter) {
	siter->second *= scale_source_;
	norm_source_ += siter->second.squaredNorm();
      }
      
      scale_source_ = 1.0;
    }
    
    if  (scale_target_ != 1.0) {
      norm_target_ = 0.0;
      
      embedding_type::iterator titer_end = target_.end();
      for (embedding_type::iterator titer = target_.begin(); titer != titer_end; ++ titer) {
	titer->second *= scale_target_;
	norm_target_ += titer->second.squaredNorm();
      }
      
      scale_target_ = 1.0;
    }
  }

  void rescale(const double scaling, const bool ignore_bias)
  {
    scale_source_ *= scaling;
    scale_target_ *= scaling;

    if (scale_source_ == 0.0) {
      embedding_type::iterator siter_end = source_.end();
      for (embedding_type::iterator siter = source_.begin(); siter != siter_end; ++ siter)
	siter->second.setZero();
      
      scale_source_ = 1.0;
      norm_source_  = 0.0;
    } else
      norm_source_ *= scaling * scaling;

    if (scale_target_ == 0.0) {
      embedding_type::iterator titer_end = target_.end();
      for (embedding_type::iterator titer = target_.begin(); titer != titer_end; ++ titer)
	titer->second.setZero();
      
      scale_target_ = 1.0;
      norm_target_ = 0.0;
    } else
      norm_target_ *= scaling * scaling;
    
    Ws1_ *= scaling;
    Wi1_ *= scaling;
    
    Ws2_ *= scaling;
    Wi2_ *= scaling;

    Wp1_ *= scaling;
    Wp2_ *= scaling;

    Wl1_ *= scaling;
    Wl2_ *= scaling;

    Wc_ *= scaling;
    
    if (! ignore_bias) {
      bs1_ *= scaling;
      bi1_ *= scaling;
      
      bs2_ *= scaling;
      bi2_ *= scaling;

      bp1_ *= scaling;
      bp2_ *= scaling;

      bl1_ *= scaling;
      bl2_ *= scaling;

      bc_ *= scaling;
    }
  }

  double squared_norm(bool ignore_bias) const
  {
    double norm = norm_source_ + norm_target_;

    norm += Ws1_.squaredNorm();
    norm += Wi1_.squaredNorm();
    
    norm += Ws2_.squaredNorm();
    norm += Wi2_.squaredNorm();

    norm += Wp1_.squaredNorm();
    norm += Wp2_.squaredNorm();

    norm += Wl1_.squaredNorm();
    norm += Wl2_.squaredNorm();

    norm += Wc_.squaredNorm();
    
    if (! ignore_bias) {
      norm += bs1_.squaredNorm();
      norm += bi1_.squaredNorm();
      
      norm += bs2_.squaredNorm();
      norm += bi2_.squaredNorm();

      norm += bp1_.squaredNorm();
      norm += bp2_.squaredNorm();
      
      norm += bl1_.squaredNorm();
      norm += bl2_.squaredNorm();

      norm += bc_.squaredNorm();
    }
    
    return norm;
  }
  
  template <typename Gen>
  struct randomize
  {
    randomize(Gen& gen) : gen_(gen) {}

    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return boost::random::uniform_real_distribution<Tp>(-0.1, 0.1)(const_cast<Gen&>(gen_));
    }

    Gen& gen_;
  };

  template <typename Gen>
  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const size_type dimension_itg,
		  const size_type window,
		  Gen& gen)
  {
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_itg <= 0)
      throw std::runtime_error("invalid dimension");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    dimension_itg_       = dimension_itg;
    window_              = window;
    
    // embedding
    source_.clear();
    target_.clear();
    scale_source_ = 1.0;
    scale_target_ = 1.0;
    norm_source_ = 0.0;
    norm_target_ = 0.0;
    
    // score
    Ws1_ = tensor_type::Zero(dimension_itg, dimension_itg * 2).array().unaryExpr(randomize<Gen>(gen));
    bs1_ = tensor_type::Zero(dimension_itg, 1).array().unaryExpr(randomize<Gen>(gen));
    Wi1_ = tensor_type::Zero(dimension_itg, dimension_itg * 2).array().unaryExpr(randomize<Gen>(gen));
    bi1_ = tensor_type::Zero(dimension_itg, 1).array().unaryExpr(randomize<Gen>(gen));
    
    // reconstruction
    Ws2_ = tensor_type::Zero(dimension_itg * 2, dimension_itg).array().unaryExpr(randomize<Gen>(gen));
    bs2_ = tensor_type::Zero(dimension_itg * 2, 1).array().unaryExpr(randomize<Gen>(gen));
    Wi2_ = tensor_type::Zero(dimension_itg * 2, dimension_itg).array().unaryExpr(randomize<Gen>(gen));
    bi2_ = tensor_type::Zero(dimension_itg * 2, 1).array().unaryExpr(randomize<Gen>(gen));
    
    // phrases
    Wp1_ = tensor_type::Zero(dimension_itg, dimension_hidden).array().unaryExpr(randomize<Gen>(gen));
    bp1_ = tensor_type::Zero(dimension_itg, 1).array().unaryExpr(randomize<Gen>(gen));
    Wp2_ = tensor_type::Zero(dimension_hidden, dimension_itg).array().unaryExpr(randomize<Gen>(gen));
    bp2_ = tensor_type::Zero(dimension_hidden, 1).array().unaryExpr(randomize<Gen>(gen));
    
    // lexicon
    Wl1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1)).array().unaryExpr(randomize<Gen>(gen));
    bl1_ = tensor_type::Zero(dimension_hidden, 1).array().unaryExpr(randomize<Gen>(gen));
    Wl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden).array().unaryExpr(randomize<Gen>(gen));
    bl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1).array().unaryExpr(randomize<Gen>(gen));
    
    // classification
    Wc_ = tensor_type::Zero(1, dimension_itg).array().unaryExpr(randomize<Gen>(gen));
    bc_ = tensor_type::Zero(1, 1).array().unaryExpr(randomize<Gen>(gen));
  }

  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const size_type dimension_itg,
		  const size_type window)
  {
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_itg <= 0)
      throw std::runtime_error("invalid dimension");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    dimension_itg_       = dimension_itg;
    window_              = window;
    
    // embedding
    source_.clear();
    target_.clear();
    scale_source_ = 1.0;
    scale_target_ = 1.0;
    norm_source_ = 0.0;
    norm_target_ = 0.0;
    
    // score
    Ws1_ = tensor_type::Zero(dimension_itg, dimension_itg * 2);
    bs1_ = tensor_type::Zero(dimension_itg, 1);
    Wi1_ = tensor_type::Zero(dimension_itg, dimension_itg * 2);
    bi1_ = tensor_type::Zero(dimension_itg, 1);
    
    // reconstruction
    Ws2_ = tensor_type::Zero(dimension_itg * 2, dimension_itg);
    bs2_ = tensor_type::Zero(dimension_itg * 2, 1);
    Wi2_ = tensor_type::Zero(dimension_itg * 2, dimension_itg);
    bi2_ = tensor_type::Zero(dimension_itg * 2, 1);
    
    // phrases
    Wp1_ = tensor_type::Zero(dimension_itg, dimension_hidden);
    bp1_ = tensor_type::Zero(dimension_itg, 1);
    Wp2_ = tensor_type::Zero(dimension_hidden, dimension_itg);
    bp2_ = tensor_type::Zero(dimension_hidden, 1);
    
    // lexicon
    Wl1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1));
    bl1_ = tensor_type::Zero(dimension_hidden, 1);
    Wl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden);
    bl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);
    
    // classification
    Wc_ = tensor_type::Zero(1, dimension_itg);
    bc_ = tensor_type::Zero(1, 1);
  }
  
  template <typename IteratorSource, typename IteratorTarget, typename Gen>
  void embedding(IteratorSource sfirst, IteratorSource slast,
		 IteratorTarget tfirst, IteratorTarget tlast,
		 Gen& gen)
  {
    embedding(sfirst, slast, source_, scale_source_, norm_source_, gen);
    embedding(tfirst, tlast, target_, scale_target_, norm_target_, gen);
  }
  
  template <typename Iterator, typename Gen>
  void embedding(Iterator first, Iterator last, embedding_type& embedding, double& scale, double& norm, Gen& gen)
  {
    tensor_type& eps = embedding[vocab_type::EPSILON];
    tensor_type& bos = embedding[vocab_type::BOS];
    tensor_type& eos = embedding[vocab_type::EOS];
    
    if (! eps.rows() || ! eps.cols())
      eps = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen));
    if (! bos.rows() || ! bos.cols())
      bos = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen));
    if (! eos.rows() || ! eos.cols())
      eos = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen));
    
    for (/**/; first != last; ++ first) {
      tensor_type& we = embedding[*first];
      
      if (! we.rows() || ! we.cols())
	we = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen));
    }
    
    // compute norm...
    scale = 1.0;
    norm  = 0.0;
    
    embedding_type::const_iterator siter_end = embedding.end();
    for (embedding_type::const_iterator siter = embedding.begin(); siter != siter_end; ++ siter)
      norm += siter->second.squaredNorm();
  }
  
  struct real_policy : boost::spirit::karma::real_policies<parameter_type>
  {
    static unsigned int precision(parameter_type)
    {
      return 10;
    }
  };

  boost::spirit::karma::real_generator<parameter_type, real_policy> float10;

  void read_embedding(const path_type& source_file, const path_type& target_file)
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

    if (! source_file.empty()) {
      
      if (source_file != "-" && ! boost::filesystem::exists(source_file))
	throw std::runtime_error("no embedding: " + source_file.string());
      
      utils::compress_istream is(source_file, 1024 * 1024);
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
	
	if (boost::fusion::get<1>(parsed).size() != dimension_embedding_)
	  throw std::runtime_error("invalid embedding size");
	
	source_[boost::fusion::get<0>(parsed)] = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), dimension_embedding_, 1);
      }
    }

    if (! target_file.empty()) {
      
      if (target_file != "-" && ! boost::filesystem::exists(target_file))
	throw std::runtime_error("no embedding: " + target_file.string());
      
      utils::compress_istream is(target_file, 1024 * 1024);
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
	
	if (boost::fusion::get<1>(parsed).size() != dimension_embedding_)
	  throw std::runtime_error("invalid embedding size");
	
	target_[boost::fusion::get<0>(parsed)] = Eigen::Map<const tensor_type>(&(*boost::fusion::get<1>(parsed).begin()), dimension_embedding_, 1);
      }
    }
  }

  template <typename Rep, typename Tp>
  inline
  void key_value(const Rep& rep, const std::string& key, Tp& value)
  {
    typename Rep::const_iterator iter = rep.find(key);
    if (iter != rep.end())
      value = utils::lexical_cast<Tp>(iter->second);
  }

  void read(const path_type& path)
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    repository_type rep(path, repository_type::read);

    key_value(rep, "dimension-embedding", dimension_embedding_);
    key_value(rep, "dimension-hidden",    dimension_hidden_);
    key_value(rep, "dimension-itg",       dimension_itg_);
    key_value(rep, "window",    window_);
    key_value(rep, "alpha",     alpha_);
    key_value(rep, "beta",      beta_);
    
    read_embedding(rep.path("source.gz"), rep.path("target.gz"));
    
    // read matrix...
  }
  
  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["dimension-embedding"] = utils::lexical_cast<std::string>(dimension_embedding_);
    rep["dimension-hidden"]    = utils::lexical_cast<std::string>(dimension_hidden_);
    rep["dimension-itg"]       = utils::lexical_cast<std::string>(dimension_itg_);
    
    rep["window"]    = utils::lexical_cast<std::string>(window_);
    rep["alpha"]     = utils::lexical_cast<std::string>(alpha_);
    rep["beta"]      = utils::lexical_cast<std::string>(beta_);
    
    const path_type source_file = rep.path("source.gz");
    const path_type target_file = rep.path("target.gz");

    {
      utils::compress_ostream os(source_file, 1024 * 1024);
      std::ostream_iterator<char> iter(os);
      
      embedding_type::const_iterator siter_end = source_.end();
      for (embedding_type::const_iterator siter = source_.begin(); siter != siter_end; ++ siter) {
	karma::generate(iter, standard::string, siter->first);
	
	for (difference_type i = 0; i != siter->second.rows(); ++ i)
	  karma::generate(iter, karma::lit(' ') << float10, siter->second.row(i)[0]);
	
	karma::generate(iter, karma::lit('\n'));
      }
    }

    {
      utils::compress_ostream os(target_file, 1024 * 1024);
      std::ostream_iterator<char> iter(os);
      
      embedding_type::const_iterator titer_end = target_.end();
      for (embedding_type::const_iterator titer = target_.begin(); titer != titer_end; ++ titer) {
	karma::generate(iter, standard::string, titer->first);
	
	for (difference_type i = 0; i != titer->second.rows(); ++ i)
	  karma::generate(iter, karma::lit(' ') << float10, titer->second.row(i)[0]);
	
	karma::generate(iter, karma::lit('\n'));
      }
    }
    
    // dump matrices...
    // we will dump two: binary and text
    write(rep.path("Ws1.txt.gz"), rep.path("Ws1.bin"), Ws1_);
    write(rep.path("bs1.txt.gz"), rep.path("bs1.bin"), bs1_);
    write(rep.path("Wi1.txt.gz"), rep.path("Wi1.bin"), Wi1_);
    write(rep.path("bi1.txt.gz"), rep.path("bi1.bin"), bi1_);

    write(rep.path("Ws2.txt.gz"), rep.path("Ws2.bin"), Ws2_);
    write(rep.path("bs2.txt.gz"), rep.path("bs2.bin"), bs2_);
    write(rep.path("Wi2.txt.gz"), rep.path("Wi2.bin"), Wi2_);
    write(rep.path("bi2.txt.gz"), rep.path("bi2.bin"), bi2_);
    
    write(rep.path("Wp1.txt.gz"), rep.path("Wp1.bin"), Wp1_);
    write(rep.path("bp1.txt.gz"), rep.path("bp1.bin"), bp1_);
    write(rep.path("Wp2.txt.gz"), rep.path("Wp2.bin"), Wp2_);
    write(rep.path("bp2.txt.gz"), rep.path("bp2.bin"), bp2_);

    write(rep.path("Wl1.txt.gz"), rep.path("Wl1.bin"), Wl1_);
    write(rep.path("bl1.txt.gz"), rep.path("bl1.bin"), bl1_);
    write(rep.path("Wl2.txt.gz"), rep.path("Wl2.bin"), Wl2_);
    write(rep.path("bl2.txt.gz"), rep.path("bl2.bin"), bl2_);

    write(rep.path("Wc.txt.gz"), rep.path("Wc.bin"), Wc_);
    write(rep.path("bc.txt.gz"), rep.path("bc.bin"), bc_);
  }

  void read(const path_type& path,
	    tensor_type& matrix,
	    const tensor_type::Index rows_expected,
	    const tensor_type::Index cols_expected)
  {
    utils::compress_istream is(path, 1024 * 1024);
    
    tensor_type::Index rows;
    tensor_type::Index cols;
    
    is.read((char*) &rows, sizeof(tensor_type::Index));
    is.read((char*) &cols, sizeof(tensor_type::Index));

    if (rows != rows_expected)
      throw std::runtime_error("rows does noe match for " + path.string());
    if (cols != cols_expected)
      throw std::runtime_error("rows does noe match for " + path.string());
    
    matrix.resize(rows, cols);
    is.read((char*) matrix.data(), sizeof(tensor_type::Scalar) * rows * cols);
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
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type dimension_itg_;
  size_type window_;
  
  // hyperparameter
  double alpha_;
  double beta_;
  
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

  // Wp1 and bp1 for encoding
  tensor_type Wp1_;
  tensor_type bp1_;
  
  // Wp2 and bp2 for reconstruction
  tensor_type Wp2_;
  tensor_type bp2_;  
  
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


struct ITGTree
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type tensor_type;
  
  typedef Bitext bitext_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  typedef std::vector<word_type, std::allocator<word_type> > word_set_type;

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

    friend
    bool operator==(const Span& x, const Span& y)
    {
      return x.first_ == y.first_ && x.last_ == y.last_;
    }

    friend
    bool operator!=(const Span& x, const Span& y)
    {
      return x.first_ != y.first_ || x.last_ != y.last_;
    }

    friend
    std::ostream& operator<<(std::ostream& os, const Span& x)
    {
      os << x.first_ << ".." << x.last_;
      return os;
    }
    
    index_type first_;
    index_type last_;
  };
  
  typedef Span span_type;
  
  struct SpanPair
  {
    typedef Span span_type;
    typedef span_type value_type;
    
    typedef span_type::index_type index_type;
    
    SpanPair(const span_type& source, const span_type& target)
      : source_(source), target_(target) {}
    SpanPair(const index_type& w, const index_type& x, const index_type& y, const index_type& z)
      : source_(w, x), target_(y, z) {}
    template <typename _Span>
    SpanPair(const std::pair<_Span, _Span>& x)
      : source_(x.first), target_(x.second) {}
    SpanPair() : source_(), target_() {}

    bool empty() const { return source_.empty() && target_.empty(); }
    size_type size() const { return source_.size() + target_.size(); }

    friend
    bool operator==(const SpanPair& x, const SpanPair& y)
    {
      return x.source_ == y.source_ && x.target_ == y.target_;
    }

    friend
    bool operator!=(const SpanPair& x, const SpanPair& y)
    {
      return x.source_ != y.source_ || x.target_ != y.target_;
    }

    friend
    std::ostream& operator<<(std::ostream& os, const SpanPair& x)
    {
      os << x.source_ << "-" << x.target_;
      return os;
    }
    
    span_type source_;
    span_type target_;
  };
  
  typedef SpanPair span_pair_type;

  struct Hyperedge
  {
    span_pair_type span_;
    span_pair_type left_;
    span_pair_type right_;

    Hyperedge()
      : span_(), left_(), right_() {}
    Hyperedge(const span_pair_type& span)
      : span_(span), left_(), right_() {}
    Hyperedge(const span_pair_type& span, const span_pair_type& left, const span_pair_type& right)
      : span_(span), left_(left), right_(right) {}
    
    bool aligned() const { return ! span_.source_.empty() && ! span_.target_.empty(); }
    bool terminal() const { return left_.empty() && right_.empty(); }
    bool straight() const { return ! terminal() && left_.target_.last_  == right_.target_.first_; }
    bool inverted() const { return ! terminal() && left_.target_.first_ == right_.target_.last_; }
  };

  typedef Hyperedge hyperedge_type;

  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > span_pair_set_type;
  typedef std::vector<span_pair_set_type, std::allocator<span_pair_set_type> > agenda_type;

  typedef std::pair<span_pair_type, span_pair_type> span_pair_stack_type;

  typedef std::vector<span_pair_stack_type, std::allocator<span_pair_stack_type> > stack_type;

  typedef std::vector<hyperedge_type, std::allocator<hyperedge_type> > derivation_type;
  typedef std::vector<span_pair_type, std::allocator<span_pair_type> > stack_derivation_type;

  typedef std::vector<hyperedge_type, std::allocator<hyperedge_type> > hyperedge_set_type;
  typedef std::vector<hyperedge_set_type, std::allocator<hyperedge_set_type> > tree_type;

  typedef std::pair<double, span_pair_type> error_span_pair_type;
  typedef std::vector<error_span_pair_type, std::allocator<error_span_pair_type> > heap_type;

  typedef std::pair<span_pair_type, span_pair_type> tail_set_type;
  
  struct tail_set_unassigned
  {
    tail_set_type operator()() const
    {
      return tail_set_type(span_pair_type(span_type::index_type(-1), span_type::index_type(-1),
					  span_type::index_type(-1), span_type::index_type(-1)),
			   span_pair_type(span_type::index_type(-1), span_type::index_type(-1),
					  span_type::index_type(-1), span_type::index_type(-1)));
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
    
    Node()
      : error_(std::numeric_limits<double>::infinity()),
	total_(std::numeric_limits<double>::infinity()),
	error_classification_(0.0),
	total_classification_(0.0),
	cost_(std::numeric_limits<double>::infinity()) {}
    
    bool terminal() const { return tails_.first.empty() && tails_.second.empty(); }
    bool straight() const { return tails_.first.target_.last_ == tails_.second.target_.first_; }
    
    double      error_;
    double      total_;
    
    double      error_classification_;
    double      total_classification_;
    
    double      cost_;
    
    tensor_type output_;
    tensor_type output_norm_;
    tensor_type delta_;
    
    // sampled alternative
    tensor_type output_sampled_;
    tensor_type output_sampled_norm_;
    tensor_type delta_sampled_;
    
    // recosntruction
    tensor_type reconstruction_;
    tensor_type delta_reconstruction_;
    
    // classification (p for plus, m for minus)
    double delta_classification_p_;
    double delta_classification_m_;
    
    tail_set_type tails_;
  };
  
  typedef Node node_type;
  typedef utils::bichart<node_type, std::allocator<node_type> > node_set_type;

  struct RestCost
  {
    double cost_;
    double alpha_;
    double beta_;
    
    RestCost()
      : cost_(std::numeric_limits<double>::infinity()),
	alpha_(std::numeric_limits<double>::infinity()),
	beta_(std::numeric_limits<double>::infinity()) {}
  };
  
  typedef RestCost rest_cost_type;
  typedef std::vector<rest_cost_type, std::allocator<rest_cost_type> > rest_cost_set_type;
  
  struct Leaf
  {
    Leaf()
      : error_(std::numeric_limits<double>::infinity()),
	error_classification_(0.0),
	cost_(std::numeric_limits<double>::infinity()) {}
    
    double      error_;
    double      error_classification_;

    double      cost_;
    
    // input with context
    tensor_type input_;
    tensor_type input_sampled_;
    
    // hidden
    tensor_type hidden_;
    tensor_type hidden_norm_;
    
    // sampled hidden
    tensor_type hidden_sampled_;
    tensor_type hidden_sampled_norm_;
    
    // hidden reconstruction
    tensor_type hidden_reconstruction_;
    tensor_type hidden_delta_reconstruction_;
    
    // output + reconstruction
    tensor_type output_;
    tensor_type output_norm_;
    
    // sampled alternative
    tensor_type output_sampled_;
    tensor_type output_sampled_norm_;
    
    // reconstruction
    tensor_type reconstruction_;
    tensor_type delta_reconstruction_;
    
    // classification (p for plus, m for minus)
    double delta_classification_p_;
    double delta_classification_m_;

    // sampled classifications
    word_type   source_sampled_;
    word_type   target_sampled_;
  };

  typedef Leaf leaf_type;
  typedef utils::vector2<leaf_type, std::allocator<leaf_type> > leaf_set_type;
  
  void clear()
  {
    nodes_.clear();
    leaves_.clear();

    costs_source_.clear();
    costs_target_.clear();
    
    agenda_.clear();
    heap_.clear();
    uniques_.clear();

    tree_.clear();
    
    stack_.clear();
    stack_derivation_.clear();
  }
  
  // sort by greater item so that we can pop from less
  struct heap_compare
  {
    bool operator()(const error_span_pair_type& x, const error_span_pair_type& y) const
    {
      return x.first > y.first;
    }
  };

  void forward_backward(rest_cost_set_type& costs)
  {
    const size_type sentence_size = costs.size() - 1;
    
    // forward...
    costs[0].alpha_ = 0;
    for (size_type last = 1; last <= sentence_size; ++ last) {
      const size_type first = last - 1;
      
      costs[last].alpha_ = std::min(costs[last].alpha_, costs[first].alpha_ + costs[first].cost_);
    }
    
    // backward...
    costs[sentence_size].beta_ = 0;
    for (difference_type first = sentence_size - 1; first >= 0; -- first) {
      const size_type last = first + 1;
      
      costs[first].beta_ = std::min(costs[first].beta_, costs[first].cost_ + costs[last].beta_);
    }
  }

  template <typename Function, typename Derivative>
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const model_type& theta,
	       const double beam,
	       Function   func,
	       Derivative deriv)
  {
    const size_type source_size = source.size();
    const size_type target_size = target.size();

    // initialization
    clear();
    
    nodes_.reserve(source_size + 1, target_size + 1);
    nodes_.resize(source_size + 1, target_size + 1);

    leaves_.reserve(source_size + 1, target_size + 1);
    leaves_.resize(source_size + 1, target_size + 1);
    
    costs_source_.reserve(source_size + 1);
    costs_target_.reserve(target_size + 1);

    costs_source_.resize(source_size + 1);
    costs_target_.resize(target_size + 1);
    
    agenda_.reserve(source_size + target_size + 1);
    agenda_.resize(source_size + target_size + 1);
    
    // initialization
    //std::cerr << "initialize leaves" << std::endl;

    forward_leaves(source, target, theta, func, deriv);

    //std::cerr << "start initialization" << std::endl;

    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = 0; trg <= target_size; ++ trg) 
	if (src < source_size || trg < target_size) {

	  // epsilon at target
	  if (src < source_size) {
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg)), theta, func, deriv);
	    
	    costs_source_[src].cost_ = std::min(costs_source_[src].cost_, double(nodes_(src, src + 1, trg, trg).error_));
	  }
	  
	  // epsilon at source
	  if (trg < target_size) {
	    forward(source, target, span_pair_type(span_type(src, src), span_type(trg, trg + 1)), theta, func, deriv);
	    
	    costs_target_[trg].cost_ = std::min(costs_target_[trg].cost_, double(nodes_(src, src, trg, trg + 1).error_));
	  }
	  
	  // word-pair
	  if (src < source_size && trg < target_size) {
	    forward(source, target, span_pair_type(span_type(src, src + 1), span_type(trg, trg + 1)), theta, func, deriv);
	    
	    const double error = nodes_(src, src + 1, trg, trg + 1).error_;
	    
	    costs_source_[src].cost_ = std::min(costs_source_[src].cost_, error);
	    costs_target_[trg].cost_ = std::min(costs_target_[trg].cost_, error);
	  }
	}

    // estimate rest-costs
    forward_backward(costs_source_);
    forward_backward(costs_target_);

    //std::cerr << "start parsing" << std::endl;

    // iterate!
    const double infty = std::numeric_limits<double>::infinity();
    const size_type length_max = source_size + target_size;
    double beam_curr = beam;

    for (int iter = 0; iter != 10; ++ iter) {
      
      for (size_type length = 1; length != length_max; ++ length) 
	if (! agenda_[length].empty()) {
	  span_pair_set_type& spans = agenda_[length];
	  
	  heap_.clear();
	  heap_.reserve(spans.size());
	  
	  span_pair_set_type::const_iterator siter_end = spans.end();
	  for (span_pair_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)  {
	    const double error = (nodes_(siter->source_.first_, siter->source_.last_,
					 siter->target_.first_, siter->target_.last_).total_
				  + std::max(costs_source_[siter->source_.first_].alpha_ + costs_source_[siter->source_.last_].beta_,
					     costs_target_[siter->target_.first_].alpha_ + costs_target_[siter->target_.last_].beta_));
	    
	    heap_.push_back(error_span_pair_type(error, *siter));
	    
	    std::push_heap(heap_.begin(), heap_.end(), heap_compare());
	  }
	  
	  heap_type::iterator hiter_begin = heap_.begin();
	  heap_type::iterator hiter       = heap_.end();
	  heap_type::iterator hiter_end   = heap_.end();

	  if (length > 2) {
	    // we will derive those with smaller reconstruction error...
	    const double threshold = hiter_begin->first + beam_curr;
	    for (/**/; hiter_begin != hiter && hiter_begin->first < threshold; -- hiter)
	      std::pop_heap(hiter_begin, hiter, heap_compare());
	  } else
	    hiter = hiter_begin;

	  // clear uniques
	  uniques_.clear();
	  
	  //std::cerr << "length: " << length << " beam_curr: " << (hiter_end - hiter) << std::endl;
	  
	  for (heap_type::iterator iter = hiter ; iter != hiter_end; ++ iter) {
	    const span_pair_type& span_pair = iter->second;
	    
	    // we borrow the notation...
	    const difference_type l = length;
	    const difference_type s = span_pair.source_.first_;
	    const difference_type t = span_pair.source_.last_;
	    const difference_type u = span_pair.target_.first_;
	    const difference_type v = span_pair.target_.last_;
	  
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
	      
		if (nodes_(S, s, U, u).error_ == infty) continue;
	      
		const span_pair_type  span1(S, s, U, u);
		const span_pair_type& span2(span_pair);
	      
		if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
		// compute score and add hyperedge
		forward(span_pair_type(S, t, U, v), span1, span2, theta, func, deriv);
	      }

	      // inversion
	      for (difference_type U = v + (S == s); U <= utils::bithack::min(v + L, V); ++ U) {
		// parent span: StuU
		// span1: SsvU
		// span2: stuv

		if (nodes_(S, s, v, U).error_ == infty) continue;
	      
		const span_pair_type  span1(S, s, v, U);
		const span_pair_type& span2(span_pair);
	      
		if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
		// compute score and add hyperedge
		forward(span_pair_type(S, t, u, U), span1, span2, theta, func, deriv);
	      }
	    }
	  
	    for (difference_type S = t; S <= utils::bithack::min(t + l, T); ++ S) {
	      const difference_type L = l - (S - t);
	    
	      // inversion
	      for (difference_type U = utils::bithack::max(u - L, difference_type(0)); U <= u - (S == t); ++ U) {
		// parent span: sSUv
		// span1: stuv
		// span2: tSUu

		if (nodes_(t, S, U, u).error_ == infty) continue;
	      
		const span_pair_type& span1(span_pair);
		const span_pair_type  span2(t, S, U, u);
	      
		if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
		// compute score and add hyperedge
		forward(span_pair_type(s, S, U, v), span1, span2, theta, func, deriv);
	      }
	    
	      // straight
	      for (difference_type U = v + (S == t); U <= utils::bithack::min(v + L, V); ++ U) {
		// parent span: sSuU
		// span1: stuv
		// span2: tSvU
	      
		if (nodes_(t, S, v, U).error_ == infty) continue;
	      
		const span_pair_type& span1(span_pair);
		const span_pair_type  span2(t, S, v, U);
	      
		if (! uniques_.insert(std::make_pair(span1, span2)).second) continue;
	      
		forward(span_pair_type(s, S, u, U), span1, span2, theta, func, deriv);
	      }
	    }
	  }
	}
      
      if (nodes_(0, source.size(), 0, target.size()).error_ != infty) break;

      std::cerr << "parsing failed: " << beam_curr << std::endl;
      
      beam_curr *= 10;
    }
  }

  // leaves!
  template <typename Function, typename Derivative>
  void forward_leaves(const sentence_type& source,
		      const sentence_type& target,
		      const model_type& theta,
		      Function   func,
		      Derivative deriv)
  {
    typedef model_type::embedding_type embedding_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type dimension_embedding = theta.dimension_embedding_;
    const size_type dimension_hidden    = theta.dimension_hidden_;
    const size_type dimension_itg       = theta.dimension_itg_;
    const size_type window    = theta.window_;
    
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	leaf_type& leaf = leaves_(src, trg);
	
	leaf.input_ = tensor_type(dimension_embedding * 2 * (window * 2 + 1), 1);

	if (src == 0) {
	  embedding_type::const_iterator siter = theta.source_.find(vocab_type::EPSILON);
	  if (siter == theta.source_.end())
	    throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	  
	  if (siter->second.rows() != dimension_embedding)
	    throw std::runtime_error("dimensin does not for the source side");
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    leaf.input_.block(dimension_embedding * i, 0, dimension_embedding, 1) = siter->second * theta.scale_source_;
	} else {
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& embedding_source = (src + shift <= 0
						 ? vocab_type::BOS
						 : (src + shift > source_size
						    ? vocab_type::EOS
						    : source[src + shift - 1]));
	    
	    embedding_type::const_iterator siter = theta.source_.find(embedding_source);
	    if (siter == theta.source_.end())
	      throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	    
	    if (siter->second.rows() != dimension_embedding)
	      throw std::runtime_error("dimensin does not for the source side");
	    
	    leaf.input_.block(dimension_embedding * i, 0, dimension_embedding, 1) = siter->second * theta.scale_source_;
	  }
	}
	
	if (trg == 0) {
	  embedding_type::const_iterator titer = theta.target_.find(vocab_type::EPSILON);
	  if (titer == theta.target_.end())
	    throw std::runtime_error("no target embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	  
	  if (titer->second.rows() != dimension_embedding)
	    throw std::runtime_error("dimensin does not for the target side");
	  
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    leaf.input_.block(dimension_embedding * i + offset, 0, dimension_embedding, 1) = titer->second * theta.scale_target_;
	} else {
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& embedding_target = (trg + shift <= 0
						 ? vocab_type::BOS
						 : (trg + shift > target_size
						    ? vocab_type::EOS
						    : target[trg + shift - 1]));
	    
	    embedding_type::const_iterator titer = theta.target_.find(embedding_target);
	    if (titer == theta.target_.end())
	      throw std::runtime_error("no target embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	    
	    if (titer->second.rows() != dimension_embedding)
	      throw std::runtime_error("dimensin does not for the target side");
	    
	    leaf.input_.block(dimension_embedding * i + offset, 0, dimension_embedding, 1) = titer->second * theta.scale_target_;
	  }
	}
	
	{
	  const tensor_type& c = leaf.input_;
	  const tensor_type p = (theta.Wl1_ * c + theta.bl1_).array().unaryExpr(func);
	  const tensor_type p_norm = p.normalized();
	  const tensor_type y = (theta.Wl2_ * p_norm + theta.bl2_).array().unaryExpr(func);
	  
	  tensor_type y_normalized = y;
	  for (size_type i = 0; i != 2 * (window * 2 + 1); ++ i)
	    y_normalized.block(i * dimension_embedding, 0, dimension_embedding, 1).normalize();
	  
	  const tensor_type y_minus_c = y_normalized - c;
	  
	  const double e = theta.alpha_ * 0.5 * y_minus_c.squaredNorm();
	  
	  leaf.error_       = e;
	  leaf.cost_        = e;
	  
	  leaf.hidden_      = p;
	  leaf.hidden_norm_ = p_norm;
	  
	  leaf.hidden_reconstruction_       = y_minus_c.array() * theta.alpha_;
	  leaf.hidden_delta_reconstruction_ = y.array().unaryExpr(deriv) * leaf.hidden_reconstruction_.array();
	}
	
	{
	  const tensor_type& c = leaf.hidden_norm_;
	  const tensor_type p = (theta.Wp1_ * c + theta.bp1_).array().unaryExpr(func);
	  const tensor_type p_norm = p.normalized();
	  const tensor_type y = (theta.Wp2_ * p_norm + theta.bp2_).array().unaryExpr(func);
	  const tensor_type y_minus_c = y.normalized() - c;
	  
	  const double e = theta.alpha_ * 0.5 * y_minus_c.squaredNorm();
	  
	  leaf.error_ += e;
	  leaf.cost_  += e;
	  leaf.cost_  += - func((theta.Wc_ * p_norm + theta.bc_)(0,0)) * theta.beta_;
	  
	  leaf.output_      = p;
	  leaf.output_norm_ = p_norm;
	  
	  leaf.reconstruction_       = y_minus_c.array() * theta.alpha_;
	  leaf.delta_reconstruction_ = y.array().unaryExpr(deriv) * leaf.reconstruction_.array();
	}
      }
  }

  // terminal!
  template <typename Function, typename Derivative>
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const span_pair_type& parent,
	       const model_type& theta,
	       Function   func,
	       Derivative deriv)
  {
    node_type& node = nodes_(parent.source_.first_, parent.source_.last_, parent.target_.first_, parent.target_.last_);
    
    const leaf_type& leaf = leaves_(parent.source_.empty() ? 0 : parent.source_.first_ + 1,
				    parent.target_.empty() ? 0 : parent.target_.first_ + 1);
    
    node.error_ = leaf.error_;
    node.total_ = leaf.error_;
    node.cost_  = leaf.cost_;
    
    node.output_norm_ = leaf.output_norm_;
    
    agenda_[parent.size()].push_back(parent);
  }
  
  // binary rules
  template <typename Function, typename Derivative>
  void forward(const span_pair_type& parent,
	       const span_pair_type& child1,
	       const span_pair_type& child2,
	       const model_type& theta,
	       Function   func,
	       Derivative deriv)
  {
    const size_type dimension_embedding = theta.dimension_embedding_;
    const size_type dimension_hidden    = theta.dimension_hidden_;
    const size_type dimension_itg       = theta.dimension_itg_;
    const bool straight = (child1.target_.last_ == child2.target_.first_);
    
    node_type& node = nodes_(parent.source_.first_, parent.source_.last_, parent.target_.first_, parent.target_.last_);
    const node_type& node1 = nodes_(child1.source_.first_, child1.source_.last_, child1.target_.first_, child1.target_.last_);
    const node_type& node2 = nodes_(child2.source_.first_, child2.source_.last_, child2.target_.first_, child2.target_.last_);
    
    const tensor_type W1 = (straight ? theta.Ws1_ : theta.Wi1_);
    const tensor_type b1 = (straight ? theta.bs1_ : theta.bi1_);
    const tensor_type W2 = (straight ? theta.Ws2_ : theta.Wi2_);
    const tensor_type b2 = (straight ? theta.bs2_ : theta.bi2_);
    
    if (node1.output_norm_.rows() != dimension_itg)
      std::cerr << "dimension does not match for child1: "
		<< parent << " : " << child1 << " + " << child2 << std::endl
		<< "error1: " << node1.error_ << " error2: " << node2.error_ << std::endl;
    
    if (node2.output_norm_.rows() != dimension_itg)
      std::cerr << "dimension does not match for child2: "
		<< parent << " : " << child1 << " + " << child2 << std::endl
		<< "error1: " << node1.error_ << " error2: " << node2.error_ << std::endl;
    

    tensor_type c(dimension_itg * 2, 1);
    c << node1.output_norm_, node2.output_norm_;

    // actual values to be propagated
    const tensor_type p = (W1 * c + b1).array().unaryExpr(func);
      
    // internal representation...
    const tensor_type p_norm = p.normalized();
      
    // compute reconstruction
    const tensor_type y = (W2 * p_norm + b2).array().unaryExpr(func);

    tensor_type y_normalized = y;
    y_normalized.block(0, 0, dimension_itg, 1).normalize();
    y_normalized.block(dimension_itg, 0, dimension_itg, 1).normalize();
    
    // internal representation...
    const tensor_type y_minus_c = y_normalized - c;
    
    // representation error
    const double e     = theta.alpha_ * 0.5 * y_minus_c.squaredNorm();
    const double total = e + node1.total_ + node2.total_;
    const double cost  = node1.cost_ + node2.cost_ + e - func((theta.Wc_ * p_norm + theta.bc_)(0,0)) * theta.beta_;
    
    const double infty = std::numeric_limits<double>::infinity();
    
    if (cost < node.cost_) {

      if (node.error_ == infty)
	agenda_[parent.size()].push_back(parent);
      
      node.error_       = e;
      node.total_       = total;
      node.cost_        = cost;

      node.output_      = p;
      node.output_norm_ = p_norm;
      
      node.reconstruction_ = y_minus_c.array() * theta.alpha_;
      node.delta_reconstruction_ = y.array().unaryExpr(deriv) * node.reconstruction_.array();
      
      node.tails_.first  = child1;
      node.tails_.second = child2;
    }
    
  }
  
  // additional forward loop to derive "wrong translations"
  template <typename Function, typename Derivative, typename Gen>
  void forward(const sentence_type& source,
	       const sentence_type& target,
	       const word_set_type& sources,
	       const word_set_type& targets,
	       const model_type& theta,
	       Function   func,
	       Derivative deriv,
	       Gen& gen)
  {
    typedef model_type::embedding_type embedding_type;

    const size_type dimension_embedding = theta.dimension_embedding_;
    const size_type dimension_hidden    = theta.dimension_hidden_;
    const size_type dimension_itg       = theta.dimension_itg_;
    const size_type window    = theta.window_;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    tree_.reserve(source_size + target_size + 1);
    tree_.resize(source_size + target_size + 1);
    
    stack_derivation_.clear();
    stack_derivation_.push_back(span_pair_type(0, source_size, 0, target_size));
    
    while (! stack_derivation_.empty()) {
      const span_pair_type span = stack_derivation_.back();
      stack_derivation_.pop_back();
      
      const node_type& node = nodes_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_);
      
      if (node.terminal())
	tree_[span.size()].push_back(hyperedge_type(span));
      else {
	const span_pair_type& child1 = node.tails_.first;
	const span_pair_type& child2 = node.tails_.second;
	
	// this is the order we push in the backward operation
	stack_derivation_.push_back(child1);
	stack_derivation_.push_back(child2);
	
	tree_[span.size()].push_back(hyperedge_type(span, child1, child2));
      }
    }
    
    // now, we are ready to enumerate tree structure in a bottom up fashion
    const size_type length_max = source_size + target_size;
    for (size_type length = 1; length <= length_max; ++ length) 
      if (! tree_[length].empty()) {
	hyperedge_set_type::const_iterator eiter_end = tree_[length].end();
	for (hyperedge_set_type::const_iterator eiter = tree_[length].begin(); eiter != eiter_end; ++ eiter) {
	  const hyperedge_type& edge = *eiter;

	  const span_pair_type& parent = edge.span_;
	  const span_pair_type& child1 = edge.left_;
	  const span_pair_type& child2 = edge.right_;
	  
	  node_type& node = nodes_(parent.source_.first_, parent.source_.last_,
				   parent.target_.first_, parent.target_.last_);
	  
	  if (edge.terminal()) {
	    leaf_type& leaf = leaves_(edge.span_.source_.empty() ? 0 : edge.span_.source_.first_ + 1,
				      edge.span_.target_.empty() ? 0 : edge.span_.target_.first_ + 1);
	    
	    leaf.input_sampled_ = leaf.input_;
	    leaf.source_sampled_ = vocab_type::EPSILON;
	    leaf.target_sampled_ = vocab_type::EPSILON;
	    
	    if (! edge.span_.source_.empty()) {
	      const size_type pos = boost::random::uniform_int_distribution<size_t>(0, sources.size() - 1)(gen);

	      leaf.source_sampled_ = sources[pos];
	      
	      embedding_type::const_iterator siter = theta.source_.find(sources[pos]);
	      if (siter == theta.source_.end())
		throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	      
	      if (siter->second.rows() != dimension_embedding)
		throw std::runtime_error("dimensin does not for the source side");
	      
	      leaf.input_sampled_.block(dimension_embedding * window, 0, dimension_embedding, 1) = siter->second * theta.scale_source_;
	    }
	    
	    if (! edge.span_.target_.empty()) {
	      const size_type pos = boost::random::uniform_int_distribution<size_t>(0, targets.size() - 1)(gen);

	      leaf.target_sampled_ = targets[pos];
	      
	      embedding_type::const_iterator titer = theta.target_.find(targets[pos]);
	      if (titer == theta.target_.end())
		throw std::runtime_error("no target embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	      
	      if (titer->second.rows() != dimension_embedding)
		throw std::runtime_error("dimensin does not for the target side");
	      
	      const size_type offset = dimension_embedding * (window * 2 + 1);
	      
	      leaf.input_sampled_.block(offset + dimension_embedding * window, 0, dimension_embedding, 1) = titer->second * theta.scale_target_;
	    }
	    
	    leaf.hidden_sampled_      = (theta.Wl1_ * leaf.input_sampled_ + theta.bl1_).array().unaryExpr(func);
	    leaf.hidden_sampled_norm_ = leaf.hidden_sampled_.normalized();
	    
	    leaf.output_sampled_      = (theta.Wp1_ * leaf.hidden_sampled_ + theta.bp1_).array().unaryExpr(func);
	    leaf.output_sampled_norm_ = leaf.output_sampled_.normalized();
	    
	    node.output_sampled_norm_ = leaf.output_sampled_norm_;
	    
	    const double y_p = func((theta.Wc_ * node.output_norm_ + theta.bc_)(0,0));
	    const double y_m = func((theta.Wc_ * node.output_sampled_norm_ + theta.bc_)(0,0));
	    const double error = std::max(1.0 - (y_p - y_m), 0.0) * theta.beta_;
	    
	    leaf.error_classification_ = error;
	    node.error_classification_ = error;
	    node.total_classification_ = error;
	    
	    leaf.delta_classification_p_ = - deriv(y_p) * (error > 0.0) * theta.beta_;
	    leaf.delta_classification_m_ =   deriv(y_m) * (error > 0.0) * theta.beta_;
	  } else {
	    const bool straight = edge.straight();
	    
	    const node_type& node1 = nodes_(child1.source_.first_, child1.source_.last_,
					    child1.target_.first_, child1.target_.last_);
	    const node_type& node2 = nodes_(child2.source_.first_, child2.source_.last_,
					    child2.target_.first_, child2.target_.last_);
	    
	    const tensor_type W1 = (straight ? theta.Ws1_ : theta.Wi1_);
	    const tensor_type b1 = (straight ? theta.bs1_ : theta.bi1_);
	    
	    tensor_type c(dimension_itg * 2, 1);
	    c << node1.output_sampled_norm_, node2.output_sampled_norm_;
	    
	    node.output_sampled_      = (W1 * c + b1).array().unaryExpr(func);
	    node.output_sampled_norm_ = node.output_sampled_.normalized();
	    
	    const double y_p = func((theta.Wc_ * node.output_norm_ + theta.bc_)(0,0));
	    const double y_m = func((theta.Wc_ * node.output_sampled_norm_ + theta.bc_)(0,0));
	    const double error = std::max(1.0 - (y_p - y_m), 0.0) * theta.beta_;
	    
	    node.error_classification_ = error;
	    node.total_classification_ = error + node1.total_classification_ + node2.total_classification_;
	    
	    node.delta_classification_p_ = - deriv(y_p) * (error > 0.0) * theta.beta_;
	    node.delta_classification_m_ =   deriv(y_m) * (error > 0.0) * theta.beta_;
	  }
	}
      }
  }
  
  // backward propagation from the goal node!
  template <typename Function, typename Derivative>
  void backward(const sentence_type& source,
		const sentence_type& target,
		const model_type& theta,
		gradient_type& gradient,
		Function   func,
		Derivative deriv)
  {
    const size_type dimension_embedding = theta.dimension_embedding_;
    const size_type dimension_hidden    = theta.dimension_hidden_;
    const size_type dimension_itg       = theta.dimension_itg_;
    const size_type window    = theta.window_;

    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const tensor_type root_W1 = tensor_type::Ones(dimension_itg, dimension_itg * 2);
    const tensor_type root_reconstruction = tensor_type::Zero(dimension_itg * 2, 1);
    
    const span_pair_type span_root(0, source_size, 0, target_size);
    
    node_type& node_root = nodes_(0, source_size, 0, target_size);
    
    node_root.delta_         = tensor_type::Zero(dimension_itg, 1);
    node_root.delta_sampled_ = tensor_type::Zero(dimension_itg, 1);
    
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
      const bool straight((span.source_.first_ == parent.source_.first_ && span.target_.first_ == parent.target_.first_)
			  || (span.source_.last_ == parent.source_.last_ && span.target_.last_ == parent.target_.last_));
      const bool left(span.source_.first_ == parent.source_.first_);
      
      const bool straight_child = node.straight();
      
      const tensor_type& W1 = (root ? root_W1 : (straight ? theta.Ws1_ : theta.Wi1_));
      const tensor_type& reconstruction = (root ? root_reconstruction : node_parent.reconstruction_);
      
      // root behave similar to left
      
      // (pre-)termianl or non-terminal?
      if (node.terminal()) {
	leaf_type& leaf = leaves_(span.source_.empty() ? 0 : span.source_.first_ + 1,
				  span.target_.empty() ? 0 : span.target_.first_ + 1);
	
	tensor_type delta;
	
	if (root || left)
	  delta = (leaf.output_.array().unaryExpr(deriv)
		   * (theta.Wp2_.transpose() * leaf.delta_reconstruction_
		      + theta.Wc_.transpose() * leaf.delta_classification_p_
		      + W1.block(0, 0, dimension_itg, dimension_itg).transpose() * node_parent.delta_
		      - reconstruction.block(0, 0, dimension_itg, 1)).array());
	else
	  delta = (leaf.output_.array().unaryExpr(deriv)
		   * (theta.Wp2_.transpose() * leaf.delta_reconstruction_
		      + theta.Wc_.transpose() * leaf.delta_classification_p_
		      + W1.block(0, dimension_itg, dimension_itg, dimension_itg).transpose() * node_parent.delta_
		      - reconstruction.block(dimension_itg, 0, dimension_itg, 1)).array());
	
	gradient.Wp1_ += delta * leaf.hidden_.transpose();
	gradient.bp1_ += delta;
	
	gradient.Wp2_ += leaf.delta_reconstruction_ * leaf.output_norm_.transpose();
	gradient.bp2_ += leaf.delta_reconstruction_;
	
	gradient.Wc_         += leaf.delta_classification_p_ * leaf.output_norm_.transpose();
	gradient.bc_.array() += leaf.delta_classification_p_;
	
	const tensor_type delta_hidden = (leaf.hidden_.array().unaryExpr(deriv)
					  * (theta.Wl2_.transpose() * leaf.hidden_delta_reconstruction_
					     + theta.Wp1_.transpose() * delta
					     - leaf.reconstruction_).array());
	
	gradient.Wl1_ += delta_hidden * leaf.input_.transpose();
	gradient.bl1_ += delta_hidden;
	
	gradient.Wl2_ += leaf.hidden_delta_reconstruction_ * leaf.hidden_norm_.transpose();
	gradient.bl2_ += leaf.hidden_delta_reconstruction_;
	
	const tensor_type delta_embedding = theta.Wl1_.transpose() * delta_hidden - leaf.hidden_reconstruction_;
	
	if (span.source_.empty()) {
	  tensor_type& dsource = gradient.source_[vocab_type::EPSILON];
	  
	  if (! dsource.cols() || ! dsource.rows())
	    dsource = tensor_type::Zero(dimension_embedding, 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    dsource += delta_embedding.block(dimension_embedding * i, 0, dimension_embedding, 1);
	} else {
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& word = (span.source_.first_ + shift < 0
				     ? vocab_type::BOS
				     : (span.source_.first_ + shift >= source_size
					? vocab_type::EOS
					: source[span.source_.first_ + shift]));
	    
	    tensor_type& dsource = gradient.source_[word];
	    
	    if (! dsource.cols() || ! dsource.rows())
	      dsource = tensor_type::Zero(dimension_embedding, 1);
	    
	    dsource += delta_embedding.block(dimension_embedding * i, 0, dimension_embedding, 1);
	  }
	}
	
	if (span.target_.empty()) {
	  tensor_type& dtarget = gradient.target_[vocab_type::EPSILON];
	  
	  if (! dtarget.cols() || ! dtarget.rows())
	    dtarget = tensor_type::Zero(dimension_embedding, 1);

	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    dtarget += delta_embedding.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	} else {
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& word = (span.target_.first_ + shift < 0
				     ? vocab_type::BOS
				     : (span.target_.first_ + shift >= target_size
					? vocab_type::EOS
					: target[span.target_.first_ + shift]));
	    
	    tensor_type& dtarget = gradient.target_[word];
	    
	    if (! dtarget.cols() || ! dtarget.rows())
	      dtarget = tensor_type::Zero(dimension_embedding, 1);
	    
	    dtarget += delta_embedding.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	  }
	}
	
	{
	  tensor_type delta;
	  
	  if (root || left)
	    delta = (leaf.output_sampled_.array().unaryExpr(deriv)
		     * (theta.Wc_.transpose() * leaf.delta_classification_m_
			+ W1.block(0, 0, dimension_itg, dimension_itg).transpose()
			* node_parent.delta_sampled_).array());
	  else
	    delta = (leaf.output_sampled_.array().unaryExpr(deriv)
		     * (theta.Wc_.transpose() * leaf.delta_classification_m_
			+ W1.block(0, dimension_itg, dimension_itg, dimension_itg).transpose()
			* node_parent.delta_sampled_).array());
	  
	  gradient.Wp1_ += delta * leaf.hidden_sampled_.transpose();
	  gradient.bp1_ += delta;
	  
	  gradient.Wc_         += leaf.delta_classification_m_ * leaf.output_sampled_norm_.transpose();
	  gradient.bc_.array() += leaf.delta_classification_m_;
	  
	  const tensor_type delta_hidden = (leaf.hidden_.array().unaryExpr(deriv) * (theta.Wp1_.transpose() * delta).array());
	  
	  gradient.Wl1_ += delta_hidden * leaf.input_sampled_.transpose();
	  gradient.bl1_ += delta_hidden;
	  
	  const tensor_type delta_embedding = theta.Wl1_.transpose() * delta_hidden;
	  
	  if (span.source_.empty()) {
	  tensor_type& dsource = gradient.source_[vocab_type::EPSILON];
	  
	  if (! dsource.cols() || ! dsource.rows())
	    dsource = tensor_type::Zero(dimension_embedding, 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    dsource += delta_embedding.block(dimension_embedding * i, 0, dimension_embedding, 1);
	  } else {
	    for (size_type i = 0; i != window * 2 + 1; ++ i) {
	      const difference_type shift = difference_type(i) - window;
	      
	      const word_type& word = (span.source_.first_ + shift < 0
				       ? vocab_type::BOS
				       : (span.source_.first_ + shift >= source_size
					  ? vocab_type::EOS
					  : (shift == 0 
					     ? leaf.source_sampled_
					     : source[span.source_.first_ + shift])));
	      
	      tensor_type& dsource = gradient.source_[word];
	      
	      if (! dsource.cols() || ! dsource.rows())
		dsource = tensor_type::Zero(dimension_embedding, 1);
	      
	      dsource += delta_embedding.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    }
	  }
	  
	  
	  if (span.target_.empty()) {
	    tensor_type& dtarget = gradient.target_[vocab_type::EPSILON];
	    
	    if (! dtarget.cols() || ! dtarget.rows())
	      dtarget = tensor_type::Zero(dimension_embedding, 1);
	    
	    const size_type offset = dimension_embedding * (window * 2 + 1);
	    
	    for (size_type i = 0; i != window * 2 + 1; ++ i)
	      dtarget += delta_embedding.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	  } else {
	    const size_type offset = dimension_embedding * (window * 2 + 1);
	    
	    for (size_type i = 0; i != window * 2 + 1; ++ i) {
	      const difference_type shift = difference_type(i) - window;
	      
	      const word_type& word = (span.target_.first_ + shift < 0
				       ? vocab_type::BOS
				       : (span.target_.first_ + shift >= target_size
					  ? vocab_type::EOS
					  : (shift == 0
					     ? leaf.target_sampled_
					     : target[span.target_.first_ + shift])));
	    
	      tensor_type& dtarget = gradient.target_[word];
	    
	      if (! dtarget.cols() || ! dtarget.rows())
		dtarget = tensor_type::Zero(dimension_embedding, 1);
	    
	      dtarget += delta_embedding.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    }
	  }
	  
	}

      } else {
	const span_pair_type& child1 = node.tails_.first;
	const span_pair_type& child2 = node.tails_.second;
	
	node_type& node1 = nodes_(child1.source_.first_, child1.source_.last_, child1.target_.first_, child1.target_.last_);
	node_type& node2 = nodes_(child2.source_.first_, child2.source_.last_, child2.target_.first_, child2.target_.last_);
	
	stack_.push_back(std::make_pair(child1, span));
	stack_.push_back(std::make_pair(child2, span));
	
	const tensor_type& W2 = (straight_child ? theta.Ws2_ : theta.Wi2_);
	
	if (root || left)
	  node1.delta_ = (node.output_.array().unaryExpr(deriv)
			  * (W2.transpose() * node.delta_reconstruction_
			     + theta.Wc_.transpose() * node.delta_classification_p_
			     + W1.block(0, 0, dimension_itg, dimension_itg).transpose() * node_parent.delta_
			     - reconstruction.block(0, 0, dimension_itg, 1)).array());
	else
	  node1.delta_ = (node.output_.array().unaryExpr(deriv)
			  * (W2.transpose() * node.delta_reconstruction_
			     + theta.Wc_.transpose() * node.delta_classification_p_
			     + W1.block(0, dimension_itg, dimension_itg, dimension_itg).transpose() * node_parent.delta_
			     - reconstruction.block(dimension_itg, 0, dimension_itg, 1)).array());
	node2.delta_ = node1.delta_;
	
	// accumulate based on the deltas
	
	const tensor_type& delta = node1.delta_;
	
	tensor_type& dW1 = (straight_child ? gradient.Ws1_ : gradient.Wi1_);
	tensor_type& db1 = (straight_child ? gradient.bs1_ : gradient.bi1_);
	
	tensor_type& dW2 = (straight_child ? gradient.Ws2_ : gradient.Wi2_);
	tensor_type& db2 = (straight_child ? gradient.bs2_ : gradient.bi2_);
	
	dW1.block(0, 0, dimension_itg, dimension_itg)             += delta * node1.output_norm_.transpose();
	dW1.block(0, dimension_itg, dimension_itg, dimension_itg) += delta * node2.output_norm_.transpose();
	db1 += delta;
	
	dW2 += node.delta_reconstruction_ * node.output_norm_.transpose();
	db2 += node.delta_reconstruction_;
	
	gradient.Wc_         += node.delta_classification_p_ * node.output_norm_.transpose();
	gradient.bc_.array() += node.delta_classification_p_;
	
	{
	  if (root || left)
	    node1.delta_sampled_ = (node.output_sampled_.array().unaryExpr(deriv)
				    * (theta.Wc_.transpose() * node.delta_classification_m_
				       + W1.block(0, 0, dimension_itg, dimension_itg).transpose()
				       * node_parent.delta_sampled_).array());
	  else
	    node1.delta_sampled_ = (node.output_sampled_.array().unaryExpr(deriv)
				    * (theta.Wc_.transpose() * node.delta_classification_m_
				       + W1.block(0, dimension_itg, dimension_itg, dimension_itg).transpose()
				       * node_parent.delta_sampled_).array());
	  
	  node2.delta_sampled_ = node1.delta_sampled_;
	  
	  const tensor_type& delta = node1.delta_sampled_;
	  
	  dW1.block(0, 0, dimension_itg, dimension_itg)             += delta * node1.output_sampled_norm_.transpose();
	  dW1.block(0, dimension_itg, dimension_itg, dimension_itg) += delta * node2.output_sampled_norm_.transpose();
	  db1 += delta;

	  gradient.Wc_         += node.delta_classification_m_ * node.output_sampled_norm_.transpose();
	  gradient.bc_.array() += node.delta_classification_m_;
	}
      }
    }
  }

  void derivation(const sentence_type& source,
		  const sentence_type& target,
		  derivation_type& d)
  {
    d.clear();
    
    stack_derivation_.clear();
    stack_derivation_.push_back(span_pair_type(0, source.size(), 0, target.size()));
    
    while (! stack_derivation_.empty()) {
      const span_pair_type span = stack_derivation_.back();
      stack_derivation_.pop_back();
      
      const node_type& node = nodes_(span.source_.first_, span.source_.last_, span.target_.first_, span.target_.last_);
      
      if (node.terminal())
	d.push_back(hyperedge_type(span));
      else {
	const span_pair_type& child1 = node.tails_.first;
	const span_pair_type& child2 = node.tails_.second;
	
	d.push_back(hyperedge_type(span, child1, child2));
	
	// we will push in right-to-left order...!
	stack_derivation_.push_back(child2);
	stack_derivation_.push_back(child1);
      }
    }
  }
  
  node_set_type nodes_;
  leaf_set_type leaves_;

  tree_type tree_;

  rest_cost_set_type costs_source_;
  rest_cost_set_type costs_target_;
  
  agenda_type agenda_;
  heap_type   heap_;
  tail_set_unique_type uniques_;
  
  stack_type stack_;
  stack_derivation_type stack_derivation_;
};

struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& dimension_embedding,
	       const size_type& dimension_hidden,
	       const size_type& dimension_itg,
	        const size_type& window, const double& lambda, const double& eta0)
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden), dimension_itg_(dimension_itg),
      window_(window), lambda_(lambda), eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");
    
    // initialize...
    Ws1_ = tensor_type::Zero(dimension_itg, dimension_itg * 2);
    bs1_ = tensor_type::Zero(dimension_itg, 1);
    Wi1_ = tensor_type::Zero(dimension_itg, dimension_itg * 2);
    bi1_ = tensor_type::Zero(dimension_itg, 1);

    Ws2_ = tensor_type::Zero(dimension_itg * 2, dimension_itg);
    bs2_ = tensor_type::Zero(dimension_itg * 2, 1);
    Wi2_ = tensor_type::Zero(dimension_itg * 2, dimension_itg);
    bi2_ = tensor_type::Zero(dimension_itg * 2, 1);

    Wp1_ = tensor_type::Zero(dimension_itg, dimension_hidden);
    bp1_ = tensor_type::Zero(dimension_itg, 1);
    Wp2_ = tensor_type::Zero(dimension_hidden, dimension_itg);
    bp2_ = tensor_type::Zero(dimension_hidden, 1);
    
    Wl1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1));
    bl1_ = tensor_type::Zero(dimension_hidden, 1);
    Wl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden);
    bl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);

    Wc_ = tensor_type::Zero(1, dimension_itg);
    bc_ = tensor_type::Zero(1, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef model_type::embedding_type embedding_type;

    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.source_.find(siter->first);
      
      if (eiter == theta.source_.end()) {
	std::cerr << "WARNING: this should not happen" << std::endl;
	eiter = theta.source_.insert(std::make_pair(siter->first, tensor_type::Zero(theta.dimension_embedding_, 1))).first;
      }
      
      if (siter->first.id() >= source_.cols()) {
	const size_type pos_first = source_.cols();
	const size_type pos_last  = siter->first.id() + 1;

	tensor_type& matrix = const_cast<tensor_type&>(source_);
	
	matrix.conservativeResize(dimension_embedding_, pos_last);
	matrix.block(0, pos_first, dimension_embedding_, pos_last - pos_first).setZero();
      }
      
      update(eiter->second, source_.col(siter->first.id()), siter->second, lambda_ != 0.0);
    }

    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer) {
      embedding_type::iterator eiter = theta.target_.find(titer->first);
      
      if (eiter == theta.target_.end()) {
	std::cerr << "WARNING: this should not happen" << std::endl;
	eiter = theta.target_.insert(std::make_pair(titer->first, tensor_type::Zero(theta.dimension_embedding_, 1))).first;
      }
      
      if (titer->first.id() >= target_.cols()) {
	const size_type pos_first = target_.cols();
	const size_type pos_last  = titer->first.id() + 1;

	tensor_type& matrix = const_cast<tensor_type&>(target_);
	
	matrix.conservativeResize(dimension_embedding_, pos_last);
	matrix.block(0, pos_first, dimension_embedding_, pos_last - pos_first).setZero();
      }
      
      update(eiter->second, target_.col(titer->first.id()), titer->second, lambda_ != 0.0);
    }
    
    update(theta.Ws1_, Ws1_, gradient.Ws1_, lambda_ != 0.0);
    update(theta.bs1_, bs1_, gradient.bs1_, false);
    update(theta.Wi1_, Wi1_, gradient.Wi1_, lambda_ != 0.0);
    update(theta.bi1_, bi1_, gradient.bi1_, false);

    update(theta.Ws2_, Ws2_, gradient.Ws2_, lambda_ != 0.0);
    update(theta.bs2_, bs2_, gradient.bs2_, false);
    update(theta.Wi2_, Wi2_, gradient.Wi2_, lambda_ != 0.0);
    update(theta.bi2_, bi2_, gradient.bi2_, false);

    update(theta.Wp1_, Wp1_, gradient.Wp1_, lambda_ != 0.0);
    update(theta.bp1_, bp1_, gradient.bp1_, false);
    update(theta.Wp2_, Wp2_, gradient.Wp2_, lambda_ != 0.0);
    update(theta.bp2_, bp2_, gradient.bp2_, false);
    
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
  
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type dimension_itg_;
  size_type window_;
  
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

  // Wp1 and bp1 for encoding
  tensor_type Wp1_;
  tensor_type bp1_;
  
  // Wp2 and bp2 for reconstruction
  tensor_type Wp2_;
  tensor_type bp2_;

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
    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.source_.find(siter->first);
      
      if (eiter == theta.source_.end()) {
	std::cerr << "WARNING: this should not happen" << std::endl;
	eiter = theta.source_.insert(std::make_pair(siter->first, tensor_type::Zero(theta.dimension_embedding_, 1))).first;
      }
      
      for (difference_type i = 0; i != siter->second.rows(); ++ i) {
	const double amount = - siter->second.col(0)[i] * gradient.scale_source_ * eta;
	
	theta.norm_source_ += 2.0 * eiter->second.col(0)[i] * theta.scale_source_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / theta.scale_source_;
      }
    }

    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer) {
      embedding_type::iterator eiter = theta.target_.find(titer->first);
      
      if (eiter == theta.target_.end()) {
	std::cerr << "WARNING: this should not happen" << std::endl;
	eiter = theta.target_.insert(std::make_pair(titer->first, tensor_type::Zero(theta.dimension_embedding_, 1))).first;
      }
      
      for (difference_type i = 0; i != titer->second.rows(); ++ i) {
	const double amount = - titer->second.col(0)[i] * gradient.scale_target_ * eta;
	
	theta.norm_target_ += 2.0 * eiter->second.col(0)[i] * theta.scale_target_ * amount + amount * amount;
	eiter->second.col(0)[i] += amount / theta.scale_target_;
      }
    }

    theta.Ws1_.array() -= gradient.Ws1_.array() * eta;
    theta.bs1_.array() -= gradient.bs1_.array() * eta;
    theta.Wi1_.array() -= gradient.Wi1_.array() * eta;
    theta.bi1_.array() -= gradient.bi1_.array() * eta;
    
    theta.Ws2_.array() -= gradient.Ws2_.array() * eta;
    theta.bs2_.array() -= gradient.bs2_.array() * eta;
    theta.Wi2_.array() -= gradient.Wi2_.array() * eta;
    theta.bi2_.array() -= gradient.bi2_.array() * eta;

    theta.Wp1_.array() -= gradient.Wp1_.array() * eta;
    theta.bp1_.array() -= gradient.bp1_.array() * eta;
    theta.Wp2_.array() -= gradient.Wp2_.array() * eta;
    theta.bp2_.array() -= gradient.bp2_.array() * eta;

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

#if 0
      if (theta.scale_source_ < 0.001 || theta.scale_source_ > 1000 || theta.scale_target_ < 0.001 || theta.scale_target_ > 1000)
	theta.finalize();
#endif
    }
  }
  
  double lambda_;
  double eta0_;
  size_type epoch_;
};

typedef boost::filesystem::path path_type;

typedef Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;
typedef ITGTree::word_set_type word_set_type;

typedef Model model_type;

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type source_file;
path_type target_file;

path_type embedding_source_file;
path_type embedding_target_file;

path_type derivation_file;
path_type alignment_source_target_file;
path_type alignment_target_source_file;
path_type output_model_file;

double alpha = 0.99;
double beta = 0.01;
int dimension_embedding = 16;
int dimension_hidden = 128;
int dimension_itg = 64;
int window = 2;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 1024;
double beam = 0.1;
double lambda = 1e-5;
double eta0 = 1;

bool dump_mode = false;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const word_set_type& sources,
		  const word_set_type& targets,
		  model_type& theta);
void derivation(const bitext_set_type& bitexts,
		const model_type& theta);
void read_bitext(const path_type& source_file,
		 const path_type& target_file,
		 bitext_set_type& bitexts,
		 word_set_type& sources,
		 word_set_type& targets);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (dimension_embedding <= 0)
      throw std::runtime_error("dimension must be positive");
    if (dimension_hidden <= 0)
      throw std::runtime_error("dimension must be positive");
    if (dimension_itg <= 0)
      throw std::runtime_error("dimension must be positive");
    if (window < 0)
      throw std::runtime_error("window size should be >= 0");
    
    if (alpha < 0.0)
      throw std::runtime_error("alpha should be >= 0.0");
    if (beta < 0.0)
      throw std::runtime_error("beta should be >= 0.0");
    
    if (int(optimize_sgd) + optimize_adagrad > 1)
      throw std::runtime_error("either one of optimize-{sgd,adagrad}");
    
    if (int(optimize_sgd) + optimize_adagrad == 0)
      optimize_adagrad = true;
    
    threads = utils::bithack::max(threads, 1);
    
    // srand is used in Eigen
    std::srand(utils::random_seed());
  
    // this is optional, but safe to set this
    ::srandom(utils::random_seed());

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
        
    if (source_file.empty())
      throw std::runtime_error("no source data?");
    if (target_file.empty())
      throw std::runtime_error("no target data?");
    
    bitext_set_type bitexts;
    word_set_type sources;
    word_set_type targets;
    
    read_bitext(source_file, target_file, bitexts, sources, targets);
    
    model_type theta(dimension_embedding, dimension_hidden, dimension_itg, window, alpha, beta, generator);
    
    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());

      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta.read_embedding(embedding_source_file, embedding_target_file);
    }

    theta.embedding(sources.begin(), sources.end(), targets.begin(), targets.end(), generator);
    
    if (iteration > 0) {
      if (optimize_adagrad)
	learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, dimension_itg, window, lambda, eta0), bitexts, sources, targets, theta);
      else
	learn_online(LearnL2(lambda, eta0), bitexts, sources, targets, theta);
    }
    
    if (! derivation_file.empty() || ! alignment_source_target_file.empty() || ! alignment_target_source_file.empty())
      derivation(bitexts, theta);
    
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

struct OutputMapReduce
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef boost::filesystem::path path_type;
  
  typedef Bitext bitext_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  typedef ITGTree itg_tree_type;
  
  typedef itg_tree_type::span_type       span_type;
  typedef itg_tree_type::span_pair_type  span_pair_type;
  typedef itg_tree_type::hyperedge_type  hyperedge_type;
  typedef itg_tree_type::derivation_type derivation_type;

  struct bitext_derivation_type
  {
    size_type       id_;
    bitext_type     bitext_;
    derivation_type derivation_;
    
    bitext_derivation_type() : id_(size_type(-1)), bitext_(), derivation_() {}
    bitext_derivation_type(const size_type& id,
			   const bitext_type& bitext,
			   const derivation_type& derivation)
      : id_(id), bitext_(bitext), derivation_(derivation) {}
    
    void swap(bitext_derivation_type& x)
    {
      std::swap(id_, x.id_);
      bitext_.swap(x.bitext_);
      derivation_.swap(x.derivation_);
    }

    void clear()
    {
      id_ = size_type(-1);
      bitext_.clear();
      derivation_.clear();
    }
  };
  
  typedef bitext_derivation_type value_type;

  typedef utils::lockfree_list_queue<value_type, std::allocator<value_type> > queue_type;

  struct compare_value
  {
    bool operator()(const value_type& x, const value_type& y) const
    {
      return x.id_ < y.id_;
    }
  };
  typedef std::set<value_type, compare_value, std::allocator<value_type> > bitext_set_type;

};

namespace std
{
  inline
  void swap(OutputMapReduce::value_type& x,
	    OutputMapReduce::value_type& y)
  {
    x.swap(y);
  }
};

struct OutputDerivation : OutputMapReduce
{
  typedef std::vector<std::string, std::allocator<std::string> > stack_type;
	  
  OutputDerivation(const path_type& path,
		   queue_type& queue)
    : path_(path), queue_(queue) {}
  
  void operator()()
  {
    if (path_.empty()) {
      bitext_derivation_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_set_type bitexts;
      bitext_derivation_type bitext;
      size_type id = 0;
      
      utils::compress_ostream os(path_, 1024 * 1024);
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
	
	if (bitext.id_ == id) {
	  write(os, bitext);
	  ++ id;
	} else
	  bitexts.insert(bitext);
	
	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  write(os, *bitexts.begin());
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	write(os, *bitexts.begin());
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while writing derivation output?");
    }
  }
  
  void write(std::ostream& os, const value_type& bitext)
  {
    stack_.clear();
    
    derivation_type::const_iterator diter_end = bitext.derivation_.end();
    for (derivation_type::const_iterator diter = bitext.derivation_.begin(); diter != diter_end; ++ diter) {
      if (diter->terminal()) {
	const word_type& source = (! diter->span_.source_.empty()
				   ? bitext.bitext_.source_[diter->span_.source_.first_]
				   : vocab_type::EPSILON);
	const word_type& target = (! diter->span_.target_.empty()
				   ? bitext.bitext_.target_[diter->span_.target_.first_]
				   : vocab_type::EPSILON);
	
	os << "((( " << source << " ||| " << target << " )))";
	
	while (! stack_.empty() && stack_.back() != " ") {
	  os << stack_.back();
	  stack_.pop_back();
	}
	
	if (! stack_.empty() && stack_.back() == " ") {
	  os << stack_.back();
	  stack_.pop_back();
	}
      } else if (diter->straight()) {
	os << "[ ";
	stack_.push_back(" ]");
	stack_.push_back(" ");
      } else {
	os << "< ";
	stack_.push_back(" >");
	stack_.push_back(" ");
      }
    }
    
    os << '\n';
  }
  
  path_type   path_;
  queue_type& queue_;
  
  stack_type stack_;
};

struct OutputAlignment : OutputMapReduce
{
  typedef cicada::Alignment alignment_type;

  OutputAlignment(const path_type& path_source_target,
		  const path_type& path_target_source,
		  queue_type& queue)
    : path_source_target_(path_source_target),
      path_target_source_(path_target_source),
      queue_(queue) {}
  
  void operator()()
  {
    if (path_source_target_.empty() && path_target_source_.empty()) {
      bitext_derivation_type bitext;
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
      }
    } else {
      bitext_set_type bitexts;
      bitext_derivation_type bitext;
      size_type id = 0;
      
      std::auto_ptr<std::ostream> os_source_target(! path_source_target_.empty()
						   ? new utils::compress_ostream(path_source_target_, 1024 * 1024)
						   : 0);
      std::auto_ptr<std::ostream> os_target_source(! path_target_source_.empty()
						   ? new utils::compress_ostream(path_target_source_, 1024 * 1024)
						   : 0);
      
      for (;;) {
	queue_.pop_swap(bitext);
	
	if (bitext.id_ == size_type(-1)) break;
	
	if (bitext.id_ == id) {
	  write(os_source_target.get(), os_target_source.get(), bitext);
	  ++ id;
	} else
	  bitexts.insert(bitext);
	
	while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	  write(os_source_target.get(), os_target_source.get(), *bitexts.begin());
	  bitexts.erase(bitexts.begin());
	  ++ id;
	}
      }
      
      while (! bitexts.empty() && bitexts.begin()->id_ == id) {
	write(os_source_target.get(), os_target_source.get(), *bitexts.begin());
	bitexts.erase(bitexts.begin());
	++ id;
      }
      
      if (! bitexts.empty())
	throw std::runtime_error("error while writing derivation output?");
    }
  }


  
  void write(std::ostream* os_source_target, std::ostream* os_target_source, const value_type& bitext)
  {
    alignment_.clear();

    derivation_type::const_iterator diter_end = bitext.derivation_.end();
    for (derivation_type::const_iterator diter = bitext.derivation_.begin(); diter != diter_end; ++ diter)
      if (diter->terminal() && diter->aligned())
	alignment_.push_back(std::make_pair(diter->span_.source_.first_, diter->span_.target_.first_));
    
    if (os_source_target) {
      std::sort(alignment_.begin(), alignment_.end());
      *os_source_target << alignment_ << '\n';
    }
    
    if (os_target_source) {
      alignment_.inverse();
      
      std::sort(alignment_.begin(), alignment_.end());
      *os_target_source << alignment_ << '\n';
    }
  }

  path_type  path_source_target_;
  path_type  path_target_source_;
  queue_type& queue_;
  
  alignment_type alignment_;
};

struct TaskAccumulate
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef ITGTree itg_tree_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  typedef OutputMapReduce output_map_reduce_type;

  typedef output_map_reduce_type::bitext_derivation_type bitext_derivation_type;
  typedef output_map_reduce_type::queue_type queue_derivation_type;

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

  TaskAccumulate(const bitext_set_type& bitexts,
		 const word_set_type& sources,
		 const word_set_type& targets,
		 const model_type& theta,
		 const double& beam,
		 queue_type& queue,
		 counter_type& counter,
		 queue_derivation_type& queue_derivation,
		 queue_derivation_type& queue_alignment)
    : bitexts_(bitexts),
      sources_(sources),
      targets_(targets),
      theta_(theta),
      beam_(beam),
      queue_(queue),
      counter_(counter),
      queue_derivation_(queue_derivation),
      queue_alignment_(queue_alignment),
      gradient_(theta.dimension_embedding_, theta.dimension_hidden_, theta.dimension_itg_, theta.window_),
      error_(0),
      classification_(0),
      samples_(0) {}

  
  void operator()()
  {
    clear();

    bitext_derivation_type bitext_derivation;

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    size_type bitext_id;
    for (;;) {
      queue_.pop(bitext_id);
      
      if (bitext_id == size_type(-1)) break;
      
      const sentence_type& source = bitexts_[bitext_id].source_;
      const sentence_type& target = bitexts_[bitext_id].target_;

      bitext_derivation.id_     = bitext_id;
      bitext_derivation.bitext_ = bitexts_[bitext_id];
      bitext_derivation.derivation_.clear();
      
      if (! source.empty() && ! target.empty()) {

#if 0
	std::cerr << "source: " << source << std::endl
		  << "target: " << target << std::endl;
#endif
	
	itg_tree_.forward(source, target, theta_, beam_, nn::tanh(), nn::dtanh());

	const itg_tree_type::node_type& root = itg_tree_.nodes_(0, source.size(), 0, target.size());

	const bool parsed = (root.error_ != std::numeric_limits<double>::infinity());
	
	if (parsed) {
	  itg_tree_.forward(source, target, sources_, targets_, theta_, nn::tanh(), nn::dtanh(), generator);
	  
	  itg_tree_.backward(source, target, theta_, gradient_, nn::tanh(), nn::dtanh());
	  
	  error_          += root.total_;
	  classification_ += root.total_classification_;
	  ++ samples_;
	  
	  itg_tree_.derivation(source, target, bitext_derivation.derivation_);
	  
	} else
	  std::cerr << "failed parsing: " << std::endl
		    << "source: " << source << std::endl
		    << "target: " << target << std::endl;
      }
      
      counter_.increment();
      
      queue_derivation_.push(bitext_derivation);
      queue_alignment_.push(bitext_derivation);
    }
  }

  void clear()
  {
    gradient_.clear();
    error_ = 0.0;
    classification_ = 0.0;
    samples_ = 0;
  }

  const bitext_set_type& bitexts_;
  const word_set_type& sources_;
  const word_set_type& targets_;
  const model_type& theta_;
  const double beam_;
  
  queue_type&            queue_;
  counter_type&          counter_;
  queue_derivation_type& queue_derivation_;
  queue_derivation_type& queue_alignment_;
  
  itg_tree_type itg_tree_;

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
		  const bitext_set_type& bitexts,
		  const word_set_type& sources,
		  const word_set_type& targets,
		  model_type& theta)
{
  typedef TaskAccumulate task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputDerivation output_derivation_type;
  typedef OutputAlignment  output_alignment_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  task_type::queue_type   mapper(256 * threads);
  task_type::counter_type reducer;

  output_map_reduce_type::queue_type queue_derivation;
  output_map_reduce_type::queue_type queue_alignment;
  
  task_set_type tasks(threads, task_type(bitexts,
					 sources,
					 targets,
					 theta,
					 beam,
					 mapper,
					 reducer,
					 queue_derivation,
					 queue_alignment));
  
  id_set_type ids(bitexts.size());
  for (size_type i = 0; i != ids.size(); ++ i)
    ids[i] = i;
  
  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  for (int t = 0; t < iteration; ++ t) {
    if (debug)
      std::cerr << "iteration: " << (t + 1) << std::endl;

    const std::string iter_tag = '.' + utils::lexical_cast<std::string>(t + 1);

    boost::thread output_derivation(output_derivation_type(! derivation_file.empty() && dump_mode
							   ? add_suffix(derivation_file, iter_tag)
							   : path_type(),
							   queue_derivation));
    boost::thread output_alignment(output_alignment_type(! alignment_source_target_file.empty() && dump_mode
							 ? add_suffix(alignment_source_target_file, iter_tag)
							 : path_type(),
							 ! alignment_target_source_file.empty() && dump_mode
							 ? add_suffix(alignment_target_source_file, iter_tag)
							 : path_type(),
							 queue_alignment));

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

    queue_derivation.push(output_map_reduce_type::value_type());
    queue_alignment.push(output_map_reduce_type::value_type());

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
    
    output_derivation.join();
    output_alignment.join();
  }

  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    mapper.push(size_type(-1));

  workers.join_all();
  
  // call this after learning
  theta.finalize();
}

struct TaskDerivation
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef ITGTree itg_tree_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;

  typedef OutputMapReduce output_map_reduce_type;

  typedef output_map_reduce_type::bitext_derivation_type bitext_derivation_type;
  typedef output_map_reduce_type::queue_type queue_derivation_type;

  typedef utils::lockfree_list_queue<size_type, std::allocator<size_type> > queue_type;
  
  TaskDerivation(const bitext_set_type& bitexts,
		 const model_type& theta,
		 const double& beam,
		 queue_type& queue,
		 queue_derivation_type& queue_derivation,
		 queue_derivation_type& queue_alignment)
    : bitexts_(bitexts),
      theta_(theta),
      beam_(beam),
      queue_(queue),
      queue_derivation_(queue_derivation),
      queue_alignment_(queue_alignment) {}

  
  void operator()()
  {
    bitext_derivation_type bitext_derivation;
    
    size_type bitext_id;
    for (;;) {
      queue_.pop(bitext_id);
      
      if (bitext_id == size_type(-1)) break;
      
      const sentence_type& source = bitexts_[bitext_id].source_;
      const sentence_type& target = bitexts_[bitext_id].target_;
      
      bitext_derivation.id_     = bitext_id;
      bitext_derivation.bitext_ = bitexts_[bitext_id];
      bitext_derivation.derivation_.clear();
      
      if (! source.empty() && ! target.empty()) {

#if 0
	std::cerr << "source: " << source << std::endl
		  << "target: " << target << std::endl;
#endif
	
	itg_tree_.forward(source, target, theta_, beam_, nn::tanh(), nn::dtanh());
	
	const itg_tree_type::node_type& root = itg_tree_.nodes_(0, source.size(), 0, target.size());
	
	const bool parsed = (root.error_ != std::numeric_limits<double>::infinity());
	
	if (parsed)
	  itg_tree_.derivation(source, target, bitext_derivation.derivation_);
	else
	  std::cerr << "failed parsing: " << std::endl
		    << "source: " << source << std::endl
		    << "target: " << target << std::endl;
      }
      
      queue_derivation_.push(bitext_derivation);
      queue_alignment_.push(bitext_derivation);
    }
  }

  const bitext_set_type& bitexts_;
  const model_type& theta_;
  const double beam_;
  
  queue_type&            queue_;
  queue_derivation_type& queue_derivation_;
  queue_derivation_type& queue_alignment_;
  
  itg_tree_type itg_tree_;
};

void derivation(const bitext_set_type& bitexts,
		const model_type& theta)
{
  typedef TaskDerivation task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;
  
  typedef task_type::size_type size_type;
  
  typedef OutputMapReduce  output_map_reduce_type;
  typedef OutputDerivation output_derivation_type;
  typedef OutputAlignment  output_alignment_type;

  task_type::queue_type   mapper(256 * threads);
  output_map_reduce_type::queue_type queue_derivation;
  output_map_reduce_type::queue_type queue_alignment;
  
  task_set_type tasks(threads, task_type(bitexts,
					 theta,
					 beam,
					 mapper,
					 queue_derivation,
					 queue_alignment));

  boost::thread_group workers;
  for (size_type i = 0; i != tasks.size(); ++ i)
    workers.add_thread(new boost::thread(boost::ref(tasks[i])));
  
  boost::thread_group workers_dump;
  workers_dump.add_thread(new boost::thread(output_derivation_type(! derivation_file.empty()
								   ? derivation_file
								   : path_type(),
								   queue_derivation)));
  workers_dump.add_thread(new boost::thread(output_alignment_type(! alignment_source_target_file.empty()
								  ? alignment_source_target_file
								  : path_type(),
								  ! alignment_target_source_file.empty()
								  ? alignment_target_source_file
								  : path_type(),
								  queue_alignment)));
  
  for (size_type i = 0; i != bitexts.size(); ++ i)
    mapper.push(i);
  
  // termination
  for (size_type i = 0; i != tasks.size(); ++ i)
    mapper.push(size_type(-1));
  
  workers.join_all();

  queue_derivation.push(output_map_reduce_type::value_type());
  queue_alignment.push(output_map_reduce_type::value_type());
  
  workers_dump.join_all();
}

void read_bitext(const path_type& source_file,
		 const path_type& target_file,
		 bitext_set_type& bitexts,
		 word_set_type& sources,
		 word_set_type& targets)
{
  typedef utils::compact_set<bitext_type::word_type,
			     utils::unassigned<bitext_type::word_type>, utils::unassigned<bitext_type::word_type>,
			     boost::hash<bitext_type::word_type>, std::equal_to<bitext_type::word_type>,
			     std::allocator<bitext_type::word_type> > unique_set_type;
  
  bitexts.clear();

  unique_set_type uniques_source;
  unique_set_type uniques_target;
  
  utils::compress_istream src(source_file, 1024 * 1024);
  utils::compress_istream trg(target_file, 1024 * 1024);
  
  bitext_type::sentence_type source;
  bitext_type::sentence_type target;
  
  while (src && trg) {
    src >> source;
    trg >> target;
    
    if (! src || ! trg) break;
    
    bitexts.push_back(bitext_type(source, target));

    uniques_source.insert(source.begin(), source.end());
    uniques_target.insert(target.begin(), target.end());
  }
  
  if (src || trg)
    throw std::runtime_error("# of sentnces do not match");

  sources.clear();
  targets.clear();
  
  sources.reserve(uniques_source.size());
  targets.reserve(uniques_target.size());
  
  sources.insert(sources.end(), uniques_source.begin(), uniques_source.end());
  targets.insert(targets.end(), uniques_target.begin(), uniques_target.end());
}

void options(int argc, char** argv)
{
  namespace po = boost::program_options;
  
  po::options_description opts_command("command line options");
  opts_command.add_options()
    ("source",    po::value<path_type>(&source_file),    "source file")
    ("target",    po::value<path_type>(&target_file),    "target file")
    
    ("embedding-source", po::value<path_type>(&embedding_source_file), "initial source embedding")
    ("embedding-target", po::value<path_type>(&embedding_target_file), "initial target embedding")

    ("derivation",               po::value<path_type>(&derivation_file),               "output derivation")
    ("alignment-source-target",  po::value<path_type>(&alignment_source_target_file),  "output alignemnt for P(target | source)")
    ("alignment-target-source",  po::value<path_type>(&alignment_target_source_file),  "output alignemnt for P(source | target)")

    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("alpha",     po::value<double>(&alpha)->default_value(alpha),      "parameter for reconstruction error")
    ("beta",      po::value<double>(&beta)->default_value(beta),        "parameter for classificaiton error")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("dimension-itg",       po::value<int>(&dimension_itg)->default_value(dimension_itg),             "dimension for ITG layer")
    ("window",              po::value<int>(&window)->default_value(window),                           "context window size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("beam",              po::value<double>(&beam)->default_value(beam),          "beam width for parsing")
    ("lambda",            po::value<double>(&lambda)->default_value(lambda),      "regularization constant")
    ("eta0",              po::value<double>(&eta0)->default_value(eta0),          "\\eta_0 for decay")

    ("dump", po::bool_switch(&dump_mode), "dump intermediate derivations and alignments")
    
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
