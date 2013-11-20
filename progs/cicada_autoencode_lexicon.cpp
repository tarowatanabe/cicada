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
#include "utils/sampler.hpp"
#include "utils/resource.hpp"

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
  
  Model() : dimension_embedding_(0), dimension_hidden_(0), window_(0), alpha_(0), beta_(0) {}
  Model(const size_type& dimension_embedding, const size_type& dimension_hidden, const size_type& window) 
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden),
      window_(window), alpha_(0), beta_(0) { initialize(dimension_embedding, dimension_hidden, window); }
  template <typename Gen>
  Model(const size_type& dimension_embedding, const size_type& dimension_hidden, const size_type& window,
	const double& alpha, const double& beta, Gen& gen) 
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden),
      window_(window), alpha_(alpha), beta_(beta) { initialize(dimension_embedding, dimension_hidden, window, gen); }

  Model& operator-=(const Model& x)
  {
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      tensor_type& embedding = source_[siter->first];

      if (! embedding.rows())
	embedding = - siter->second;
      else
	embedding -= siter->second;
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      tensor_type& embedding = target_[titer->first];
      
      if (! embedding.rows())
	embedding = - titer->second;
      else
	embedding -= titer->second;
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
    embedding_type::const_iterator siter_end = x.source_.end();
    for (embedding_type::const_iterator siter = x.source_.begin(); siter != siter_end; ++ siter) {
      tensor_type& embedding = source_[siter->first];

      if (! embedding.rows())
	embedding = siter->second;
      else
	embedding += siter->second;
    }

    embedding_type::const_iterator titer_end = x.target_.end();
    for (embedding_type::const_iterator titer = x.target_.begin(); titer != titer_end; ++ titer) {
      tensor_type& embedding = target_[titer->first];
      
      if (! embedding.rows())
	embedding = titer->second;
      else
	embedding += titer->second;
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
    source_.clear();
    target_.clear();
    
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
  }

  template <typename Gen>
  struct randomize
  {
    randomize(Gen& gen, const double range=0.01) : gen_(gen), range_(range) {}
    
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return boost::random::uniform_real_distribution<Tp>(-range_, range_)(const_cast<Gen&>(gen_));
    }
    
    Gen& gen_;
    double range_;
  };

  template <typename Gen>
  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const size_type window,
		  Gen& gen)
  {
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    window_              = window;
    
    // embedding
    source_.clear();
    target_.clear();
    
    // intialize randomly...

    const double range_l = std::sqrt(6.0 / (dimension_hidden_ + dimension_embedding * 2 * (window * 2 + 1)));
    const double range_c = std::sqrt(6.0 / (dimension_hidden_ + 1));
    
    // lexicon
    Wl1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1)).array().unaryExpr(randomize<Gen>(gen, range_l));
    bl1_ = tensor_type::Zero(dimension_hidden, 1);
    
    // lexicon reconstruction
    Wl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden).array().unaryExpr(randomize<Gen>(gen, range_l));
    bl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);

    // classification
    Wc_ = tensor_type::Zero(1, dimension_hidden).array().unaryExpr(randomize<Gen>(gen, range_c));
    bc_ = tensor_type::Ones(1, 1).array();
  }

  void initialize(const size_type dimension_embedding,
		  const size_type dimension_hidden,
		  const size_type window)
  {
    if (dimension_embedding <= 0)
      throw std::runtime_error("invalid dimension");
    if (dimension_hidden <= 0)
      throw std::runtime_error("invalid dimension");
    
    dimension_embedding_ = dimension_embedding;
    dimension_hidden_    = dimension_hidden;
    window_              = window;
    
    // embedding
    source_.clear();
    target_.clear();
        
    // lexicon
    Wl1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1));
    bl1_ = tensor_type::Zero(dimension_hidden, 1);
    
    // lexicon reconstruction
    Wl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden);
    bl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);

    // classification
    Wc_ = tensor_type::Zero(1, dimension_hidden);
    bc_ = tensor_type::Zero(1, 1);
  }
  
  
  
  template <typename IteratorSource, typename IteratorTarget, typename Gen>
  void embedding(IteratorSource sfirst, IteratorSource slast,
		 IteratorTarget tfirst, IteratorTarget tlast,
		 Gen& gen)
  {
    embedding(sfirst, slast, source_, gen);
    embedding(tfirst, tlast, target_, gen);
  }
  
  template <typename Iterator, typename Gen>
  void embedding(Iterator first, Iterator last, embedding_type& embedding, Gen& gen)
  {
    tensor_type& eps = embedding[vocab_type::EPSILON];
    tensor_type& bos = embedding[vocab_type::BOS];
    tensor_type& eos = embedding[vocab_type::EOS];

    const double range_e = std::sqrt(6.0 / (dimension_embedding_ + 1));
    
    if (! eps.rows() || ! eps.cols())
      eps = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen, range_e));
    if (! bos.rows() || ! bos.cols())
      bos = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen, range_e));
    if (! eos.rows() || ! eos.cols())
      eos = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen, range_e));
    
    for (/**/; first != last; ++ first) {
      tensor_type& we = embedding[*first];
      
      if (! we.rows() || ! we.cols())
	we = tensor_type::Zero(dimension_embedding_, 1).array().unaryExpr(randomize<Gen>(gen, range_e));
    }
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

  void write(const path_type& path) const
  {
    // we use a repository structure...
    typedef utils::repository repository_type;
    
    namespace karma = boost::spirit::karma;
    namespace standard = boost::spirit::standard;
    
    repository_type rep(path, repository_type::write);
    
    rep["dimension-embedding"] = utils::lexical_cast<std::string>(dimension_embedding_);
    rep["dimension-hidden"]    = utils::lexical_cast<std::string>(dimension_hidden_);
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
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type window_;
  
  // hyperparameter
  double alpha_;
  double beta_;
  
  // Embedding
  embedding_type source_;
  embedding_type target_;
  
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

struct Dictionary
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model    model_type;
  
  typedef uint64_t count_type;
  
  typedef cicada::Sentence  sentence_type;
  typedef cicada::Alignment alignment_type;
  typedef cicada::Symbol    word_type;
  typedef cicada::Vocab     vocab_type;
  
  struct Dict
  {
    typedef utils::compact_map<word_type, double,
			       utils::unassigned<word_type>, utils::unassigned<word_type>,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, double> > > logprob_set_type;
    typedef utils::compact_map<word_type, count_type,
			       utils::unassigned<word_type>, utils::unassigned<word_type>,
			       boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<std::pair<const word_type, count_type> > > count_set_type;
    
    typedef std::vector<word_type, std::allocator<word_type> >   word_set_type;
    typedef boost::random::discrete_distribution<>               distribution_type;
    
    Dict() {}

    template <typename Tp>
    struct compare_pair
    {
      bool operator()(const Tp& x, const Tp& y) const
      {
	return (x.second > y.second
		|| (x.second == y.second
		    && static_cast<const std::string&>(x.first) < static_cast<const std::string&>(y.first)));
      }
    };
    
    void initialize()
    {
      typedef std::pair<word_type, double> word_prob_type;
      typedef std::vector<word_prob_type, std::allocator<word_prob_type> > word_prob_set_type;
      typedef std::vector<double, std::allocator<double> > prob_set_type;
      
      logprobs_.clear();
      words_.clear();
      
      word_prob_set_type word_probs(counts_.begin(), counts_.end());
      std::sort(word_probs.begin(), word_probs.end(), compare_pair<word_prob_type>());
      
      prob_set_type probs;
      words_.reserve(word_probs.size());
      probs.reserve(word_probs.size());
      
      word_prob_set_type::const_iterator witer_end = word_probs.end();
      for (word_prob_set_type::const_iterator witer = word_probs.begin(); witer != witer_end; ++ witer) {
	words_.push_back(witer->first);
	probs.push_back(witer->second);
	logprobs_[witer->first] = std::log(witer->second);
      }
      
      // initialize distribution
      distribution_ = distribution_type(probs.begin(), probs.end());
    }
    
    double logprob(const word_type& word) const
    {
      logprob_set_type::const_iterator liter = logprobs_.find(word);
      if (liter != logprobs_.end())
	return liter->second;
      else
	return - std::numeric_limits<double>::infinity();
    }
    
    template <typename Gen>
    word_type draw(Gen& gen) const
    {
      return words_[distribution_(gen)];
    }

    count_type& operator[](const word_type& word) { return counts_[word]; }
    
    count_set_type    counts_;
    logprob_set_type  logprobs_;
    word_set_type     words_;
    distribution_type distribution_;
  };

  typedef Dict dict_type;
  typedef utils::alloc_vector<dict_type, std::allocator<dict_type> > dict_set_type;
  
  Dictionary() {}

  dict_type& operator[](const word_type& word) { return dicts_[word.id()]; }

  void swap(Dictionary& x) { dicts_.swap(x.dicts_); }

  void initialize()
  {
    for (size_type i = 0; i != dicts_.size(); ++ i)
      if (dicts_.exists(i))
	dicts_[i].initialize();
  }
  
  void clear()
  {
    dicts_.clear();
  }
  
  double logprob(const word_type& source, const word_type& target) const
  {
    if (dicts_.exists(source.id()))
      return dicts_[source.id()].logprob(target);
    else
      return dicts_[vocab_type::UNK.id()].logprob(target);
  }
  
  template <typename Gen>
  word_type draw(const word_type& source, Gen& gen) const
  {
    if (dicts_.exists(source.id()))
      return dicts_[source.id()].draw(gen);
    else
      return dicts_[vocab_type::UNK.id()].draw(gen);
  }

  dict_set_type dicts_;
};

struct Lexicon
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;
  typedef Dictionary dictionary_type;

  typedef model_type::parameter_type parameter_type;
  typedef model_type::tensor_type tensor_type;
  
  typedef Bitext bitext_type;
  
  typedef bitext_type::word_type word_type;
  typedef bitext_type::sentence_type sentence_type;
  typedef bitext_type::vocab_type vocab_type;
  
  Lexicon(const dictionary_type& dict_source_target,
	  const dictionary_type& dict_target_source)
    : dict_source_target_(dict_source_target),
      dict_target_source_(dict_target_source) {}
  
  template <typename Function, typename Derivative, typename Gen>
  std::pair<double, double> operator()(const sentence_type& source,
				       const sentence_type& target,
				       const model_type& theta,
				       gradient_type& gradient,
				       Function   func,
				       Derivative deriv,
				       Gen& gen)
  {
    typedef model_type::embedding_type embedding_type;
    
    const size_type source_size = source.size();
    const size_type target_size = target.size();
    
    const size_type dimension_embedding = theta.dimension_embedding_;
    const size_type dimension_hidden    = theta.dimension_hidden_;
    const size_type window    = theta.window_;
    
    double error = 0.0;
    double error_classification = 0.0;

    tensor_type input(dimension_embedding * 2 * (window * 2 + 1), 1);
    tensor_type input_sampled(dimension_embedding * 2 * (window * 2 + 1), 1);
    
    for (size_type src = 0; src <= source_size; ++ src)
      for (size_type trg = (src == 0); trg <= target_size; ++ trg) {
	
	// forward...

	if (! src) {
	  embedding_type::const_iterator siter = theta.source_.find(vocab_type::EPSILON);
	  if (siter == theta.source_.end())
	    throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	  
	  if (siter->second.rows() != dimension_embedding)
	    throw std::runtime_error("dimensin does not for the source side");
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    input.block(dimension_embedding * i, 0, dimension_embedding, 1) = siter->second;
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
	      throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(embedding_source));
	    
	    if (siter->second.rows() != dimension_embedding)
	      throw std::runtime_error("dimensin does not for the source side");
	    
	    input.block(dimension_embedding * i, 0, dimension_embedding, 1) = siter->second;
	  }
	}
	
	if (! trg) {
	  embedding_type::const_iterator titer = theta.target_.find(vocab_type::EPSILON);
	  if (titer == theta.target_.end())
	    throw std::runtime_error("no target embedding for " + static_cast<const std::string&>(vocab_type::EPSILON));
	  
	  if (titer->second.rows() != dimension_embedding)
	    throw std::runtime_error("dimensin does not for the target side");
	  
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i)
	    input.block(dimension_embedding * i + offset, 0, dimension_embedding, 1) = titer->second;
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
	      throw std::runtime_error("no target embedding for " + static_cast<const std::string&>(embedding_target));
	    
	    if (titer->second.rows() != dimension_embedding)
	      throw std::runtime_error("dimensin does not for the target side");
	    
	    input.block(dimension_embedding * i + offset, 0, dimension_embedding, 1) = titer->second;
	  }
	}
	
	// compute sampled input...
	
	input_sampled = input;
	
	word_type source_sampled = vocab_type::EPSILON;
	word_type target_sampled = vocab_type::EPSILON;
	
	if (src) {
	  source_sampled = dict_target_source_.draw(trg ? target[trg - 1] : vocab_type::EPSILON, gen);
	  
	  embedding_type::const_iterator siter = theta.source_.find(source_sampled);
	  if (siter == theta.source_.end())
	    throw std::runtime_error("no source embedding for " + static_cast<const std::string&>(source_sampled));
	  
	  if (siter->second.rows() != dimension_embedding)
	    throw std::runtime_error("dimensin does not for the source side");
	  
	  input_sampled.block(dimension_embedding * window, 0, dimension_embedding, 1) = siter->second;
	}
	
	if (trg) {
	  target_sampled = dict_source_target_.draw(src ? source[src - 1] : vocab_type::EPSILON, gen);
	  
	  embedding_type::const_iterator titer = theta.target_.find(target_sampled);
	  if (titer == theta.target_.end())
	    throw std::runtime_error("no target embedding for " + static_cast<const std::string&>(target_sampled));
	  
	  if (titer->second.rows() != dimension_embedding)
	    throw std::runtime_error("dimensin does not for the target side");
	  
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  input_sampled.block(offset + dimension_embedding * window, 0, dimension_embedding, 1) = titer->second;
	}

#if 0
	std::cerr << "input:   rows: " << input.rows() << " cols: " << input.cols() << std::endl
		  << "sampled: rows: " << input_sampled.rows() << " cols: " << input_sampled.cols() << std::endl;
#endif
	
	const tensor_type p = (theta.Wl1_ * input + theta.bl1_).array().unaryExpr(func);
	const tensor_type p_norm = p.normalized();
	const tensor_type y = (theta.Wl2_ * p_norm + theta.bl2_).array().unaryExpr(func);
	
	tensor_type y_normalized = y;
	for (size_type i = 0; i != 2 * (window * 2 + 1); ++ i)
	  y_normalized.block(i * dimension_embedding, 0, dimension_embedding, 1).normalize();
	
	const tensor_type y_minus_c = y_normalized - input;
	
	const tensor_type p_sampled = (theta.Wl1_ * input_sampled + theta.bl1_).array().unaryExpr(func);
	const tensor_type p_sampled_norm = p_sampled.normalized();
	
	const double e = theta.alpha_ * 0.5 * y_minus_c.squaredNorm();

	//std::cerr << "error: " << e << std::endl;
	
	const tensor_type reconstruction       = y_minus_c.array() * theta.alpha_;
	const tensor_type delta_reconstruction = y.array().unaryExpr(deriv) * reconstruction.array();
	
	const double y_p = (theta.Wc_ * p_norm + theta.bc_)(0,0);
	const double y_m = (theta.Wc_ * p_sampled_norm + theta.bc_)(0,0);
	
	const double e_classification = std::max(1.0 - (y_p - y_m), 0.0) * theta.beta_;
	
	//std::cerr << "classification: " << e_classification << std::endl;
	
	const double delta_classification_p = - (e_classification > 0.0) * theta.beta_;
	const double delta_classification_m =   (e_classification > 0.0) * theta.beta_;

	// update error...
	error                += e;
	error_classification += e_classification;

	// backward...
	
	const tensor_type delta = (p.array().unaryExpr(deriv)
				   * (theta.Wl2_.transpose() * delta_reconstruction
				      + theta.Wc_.transpose() * delta_classification_p).array());
	
	const tensor_type delta_sampled = (p_sampled.array().unaryExpr(deriv)
					   * (theta.Wc_.transpose() * delta_classification_m).array());
	
#if 0
	std::cerr << "delta:   rows: " << delta.rows() << " cols: " << delta.cols() << std::endl
		  << "sampled: rows: " << delta_sampled.rows() << " cols: " << delta_sampled.cols() << std::endl;
#endif
	
	gradient.Wl1_ += delta * input.transpose();
	gradient.bl1_ += delta;
	
	gradient.Wl1_ += delta_sampled * input_sampled.transpose();
	gradient.bl1_ += delta_sampled;
	
	gradient.Wl2_ += delta_reconstruction * p_norm.transpose();
	gradient.bl2_ += delta_reconstruction;

	//std::cerr << "update classification" << std::endl;
	
	gradient.Wc_         += delta_classification_p * p_norm.transpose();
	gradient.bc_.array() += delta_classification_p;
	
	gradient.Wc_         += delta_classification_m * p_sampled_norm.transpose();
	gradient.bc_.array() += delta_classification_m;
	
	//std::cerr << "update embedding" << std::endl;
	
	// update embedding
	const tensor_type delta_embedding_p = theta.Wl1_.transpose() * delta - reconstruction;
	const tensor_type delta_embedding_m = theta.Wl1_.transpose() * delta_sampled;
	
	if (! src) {
	  tensor_type& dsource = gradient.source_[vocab_type::EPSILON];
	  
	  if (! dsource.cols() || ! dsource.rows())
	    dsource = tensor_type::Zero(dimension_embedding, 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    dsource += delta_embedding_p.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    dsource += delta_embedding_m.block(dimension_embedding * i, 0, dimension_embedding, 1);
	  }
	} else {
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& word = (src + shift <= 0
				     ? vocab_type::BOS
				     : (src + shift > source_size
					? vocab_type::EOS
					: source[src + shift - 1]));
	    
	    tensor_type& dsource = gradient.source_[word];
	    
	    if (! dsource.cols() || ! dsource.rows())
	      dsource = tensor_type::Zero(dimension_embedding, 1);
	    
	    dsource += delta_embedding_p.block(dimension_embedding * i, 0, dimension_embedding, 1);

	    if (shift != 0)
	      dsource += delta_embedding_m.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    else {
	      tensor_type& dsource = gradient.source_[source_sampled];
	      
	      if (! dsource.cols() || ! dsource.rows())
		dsource = tensor_type::Zero(dimension_embedding, 1);
	      
	      dsource += delta_embedding_m.block(dimension_embedding * i, 0, dimension_embedding, 1);
	    }
	  }
	}

	if (! trg) {
	  tensor_type& dtarget = gradient.target_[vocab_type::EPSILON];
	  
	  if (! dtarget.cols() || ! dtarget.rows())
	    dtarget = tensor_type::Zero(dimension_embedding, 1);
	  
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    dtarget += delta_embedding_p.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    dtarget += delta_embedding_m.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	  }
	} else {
	  const size_type offset = dimension_embedding * (window * 2 + 1);
	  
	  for (size_type i = 0; i != window * 2 + 1; ++ i) {
	    const difference_type shift = difference_type(i) - window;
	    
	    const word_type& word = (trg + shift <= 0
				     ? vocab_type::BOS
				     : (trg + shift > target_size
					? vocab_type::EOS
					: target[trg + shift - 1]));
	    
	    tensor_type& dtarget = gradient.target_[word];
	    
	    if (! dtarget.cols() || ! dtarget.rows())
	      dtarget = tensor_type::Zero(dimension_embedding, 1);
	    
	    dtarget += delta_embedding_p.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);

	    if (shift != 0)
	      dtarget += delta_embedding_m.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    else {
	      tensor_type& dtarget = gradient.target_[target_sampled];
	      
	      if (! dtarget.cols() || ! dtarget.rows())
		dtarget = tensor_type::Zero(dimension_embedding, 1);
	      
	      dtarget += delta_embedding_m.block(dimension_embedding * i + offset, 0, dimension_embedding, 1);
	    }
	  }
	}
      }
    
    return std::make_pair(error, error_classification);
  }

  const dictionary_type& dict_source_target_;
  const dictionary_type& dict_target_source_;
};

struct LearnAdaGrad
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef Model model_type;
  typedef Model gradient_type;

  typedef cicada::Symbol   word_type;
  
  typedef model_type::tensor_type tensor_type;
  
  LearnAdaGrad(const size_type& dimension_embedding,
	       const size_type& dimension_hidden,
	       const size_type& window,
	       const double& lambda,
	       const double& eta0)
    : dimension_embedding_(dimension_embedding), dimension_hidden_(dimension_hidden),
      window_(window), lambda_(lambda), eta0_(eta0)
  {
    if (lambda_ < 0.0)
      throw std::runtime_error("invalid regularization");
    
    if (eta0_ <= 0.0)
      throw std::runtime_error("invalid learning rate");

    const size_type vocabulary_size = word_type::allocated();
    source_ = tensor_type::Zero(dimension_embedding, vocabulary_size);
    target_ = tensor_type::Zero(dimension_embedding, vocabulary_size);
    
    // initialize...
    Wl1_ = tensor_type::Zero(dimension_hidden, dimension_embedding * 2 * (window * 2 + 1));
    bl1_ = tensor_type::Zero(dimension_hidden, 1);
    
    Wl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), dimension_hidden);
    bl2_ = tensor_type::Zero(dimension_embedding * 2 * (window * 2 + 1), 1);    
    
    Wc_ = tensor_type::Zero(1, dimension_hidden);
    bc_ = tensor_type::Zero(1, 1);
  }
  
  void operator()(model_type& theta, const gradient_type& gradient) const
  {
    typedef model_type::embedding_type embedding_type;

    embedding_type::const_iterator siter_end = gradient.source_.end();
    for (embedding_type::const_iterator siter = gradient.source_.begin(); siter != siter_end; ++ siter) {
      embedding_type::iterator eiter = theta.source_.find(siter->first);
      
      if (eiter == theta.source_.end()) {
	std::cerr << "WARNING: this should not happen: source: " << siter->first << std::endl;
	eiter = theta.source_.insert(std::make_pair(siter->first, tensor_type::Zero(theta.dimension_embedding_, 1))).first;
      }
      
      update(siter->first, eiter->second, const_cast<tensor_type&>(source_), siter->second, lambda_ != 0.0);
    }

    embedding_type::const_iterator titer_end = gradient.target_.end();
    for (embedding_type::const_iterator titer = gradient.target_.begin(); titer != titer_end; ++ titer) {
      embedding_type::iterator eiter = theta.target_.find(titer->first);
      
      if (eiter == theta.target_.end()) {
	std::cerr << "WARNING: this should not happen: target: "  << titer->first << std::endl;
	eiter = theta.target_.insert(std::make_pair(titer->first, tensor_type::Zero(theta.dimension_embedding_, 1))).first;
      }
      
      update(titer->first, eiter->second, const_cast<tensor_type&>(target_), titer->second, lambda_ != 0.0);
    }
    
    update(theta.Wl1_, const_cast<tensor_type&>(Wl1_), gradient.Wl1_, lambda_ != 0.0);
    update(theta.bl1_, const_cast<tensor_type&>(bl1_), gradient.bl1_, false);

    update(theta.Wl2_, const_cast<tensor_type&>(Wl2_), gradient.Wl2_, lambda_ != 0.0);
    update(theta.bl2_, const_cast<tensor_type&>(bl2_), gradient.bl2_, false);

    update(theta.Wc_, const_cast<tensor_type&>(Wc_), gradient.Wc_, lambda_ != 0.0);
    update(theta.bc_, const_cast<tensor_type&>(bc_), gradient.bc_, false);
  }

  template <typename Theta, typename GradVar, typename Grad>
  struct update_visitor_regularize
  {
    update_visitor_regularize(Eigen::MatrixBase<Theta>& theta,
			      Eigen::MatrixBase<GradVar>& G,
			      const Eigen::MatrixBase<Grad>& g,
			      const double& lambda,
			      const double& eta0)
      : theta_(theta), G_(G), g_(g), lambda_(lambda), eta0_(eta0) {}
    
    void init(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      operator()(value, i, j);
    }
    
    void operator()(const tensor_type::Scalar& value, tensor_type::Index i, tensor_type::Index j)
    {
      if (g_(i, j) == 0) return;
      
      G_(i, j) += g_(i, j) * g_(i, j);
      
      const double rate = eta0_ / std::sqrt(G_(i, j));
      const double f = theta_(i, j) - rate * g_(i, j);
      
      theta_(i, j) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
    }
    
    Eigen::MatrixBase<Theta>&      theta_;
    Eigen::MatrixBase<GradVar>&    G_;
    const Eigen::MatrixBase<Grad>& g_;
    
    const double lambda_;
    const double eta0_;
  };
  
  struct learning_rate
  {
    template <typename Tp>
    Tp operator()(const Tp& x) const
    {
      return (x == 0.0 ? 0.0 : 1.0 / std::sqrt(x));
    }
  };
  
  template <typename Theta, typename GradVar, typename Grad>
  void update(Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const bool regularize=true) const
  {
    if (regularize) {
      update_visitor_regularize<Theta, GradVar, Grad> visitor(theta, G, g, lambda_, eta0_);
      theta.visit(visitor);
    } else {
      G.array() += g.array().square();
      theta.array() -= eta0_ * g.array() * G.array().unaryExpr(learning_rate());
    }
  }

  template <typename Theta, typename GradVar, typename Grad>
  void update(const word_type& word,
	      Eigen::MatrixBase<Theta>& theta,
	      Eigen::MatrixBase<GradVar>& G,
	      const Eigen::MatrixBase<Grad>& g,
	      const bool regularize=true) const
  {
    if (regularize) {
      for (int row = 0; row != g.rows(); ++ row) 
	if (g(row, 0) != 0) {
	  G(row, word.id()) += g(row, 0) * g(row, 0);
	  
	  const double rate = eta0_ / std::sqrt(G(row, word.id()));
	  const double f = theta(row, 0) - rate * g(row, 0);
	  
	  theta(row, 0) = utils::mathop::sgn(f) * std::max(0.0, std::fabs(f) - rate * lambda_);
	}
    } else {
      G.col(word.id()).array() += g.array().square();
      theta.array() -= eta0_ * g.array() * G.col(word.id()).array().unaryExpr(learning_rate());
    }
  }
  
  size_type dimension_embedding_;
  size_type dimension_hidden_;
  size_type window_;
  
  double lambda_;
  double eta0_;
  
  // embedding
  tensor_type source_;
  tensor_type target_;
  
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


typedef boost::filesystem::path path_type;

typedef Bitext bitext_type;
typedef std::vector<bitext_type, std::allocator<bitext_type> > bitext_set_type;

typedef Model model_type;
typedef Dictionary dictionary_type;

static const size_t DEBUG_DOT  = 10000;
static const size_t DEBUG_WRAP = 100;
static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

path_type source_file;
path_type target_file;

path_type embedding_source_file;
path_type embedding_target_file;

path_type output_model_file;

double alpha = 0.99;
double beta = 0.01;
int dimension_embedding = 32;
int dimension_hidden = 128;
int window = 2;

bool optimize_sgd = false;
bool optimize_adagrad = false;

int iteration = 10;
int batch_size = 1024;
double lambda = 1e-3;
double eta0 = 1;
int cutoff = 3;

int threads = 2;

int debug = 0;

template <typename Learner>
void learn_online(const Learner& learner,
		  const bitext_set_type& bitexts,
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta);
void read_bitext(const path_type& source_file,
		 const path_type& target_file,
		 bitext_set_type& bitexts,
		 dictionary_type& dict_source_target,
		 dictionary_type& dict_target_source);

void options(int argc, char** argv);

int main(int argc, char** argv)
{
  try {
    options(argc, argv);

    if (dimension_embedding <= 0)
      throw std::runtime_error("dimension must be positive");
    if (dimension_hidden <= 0)
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
    
    dictionary_type dict_source_target;
    dictionary_type dict_target_source;
    
    read_bitext(source_file, target_file, bitexts, dict_source_target, dict_target_source);
    
    model_type theta(dimension_embedding, dimension_hidden, window, alpha, beta, generator);

    if (! embedding_source_file.empty() || ! embedding_target_file.empty()) {
      if (embedding_source_file != "-" && ! boost::filesystem::exists(embedding_source_file))
	throw std::runtime_error("no embedding: " + embedding_source_file.string());

      if (embedding_target_file != "-" && ! boost::filesystem::exists(embedding_target_file))
	throw std::runtime_error("no embedding: " + embedding_target_file.string());
      
      theta.read_embedding(embedding_source_file, embedding_target_file);
    }

    const dictionary_type::dict_type::word_set_type& sources = dict_target_source[cicada::Vocab::EPSILON].words_;
    const dictionary_type::dict_type::word_set_type& targets = dict_source_target[cicada::Vocab::EPSILON].words_;
    
    theta.embedding(sources.begin(), sources.end(), targets.begin(), targets.end(), generator);
    
    if (iteration > 0)
      learn_online(LearnAdaGrad(dimension_embedding, dimension_hidden, window, lambda, eta0),
		   bitexts,
		   dict_source_target,
		   dict_target_source,
		   theta);
    
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

  typedef Lexicon lexicon_type;
  
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

  TaskAccumulate(const bitext_set_type& bitexts,
		 const dictionary_type& dict_source_target,
		 const dictionary_type& dict_target_source,
		 const model_type& theta,
		 queue_type& queue,
		 counter_type& counter)
    : bitexts_(bitexts),
      theta_(theta),
      queue_(queue),
      counter_(counter),
      lexicon_(dict_source_target, dict_target_source),
      gradient_(theta.dimension_embedding_, theta.dimension_hidden_, theta.window_),
      error_(0),
      classification_(0),
      samples_(0) {}
  
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
      return Tp(- 1) < x && x < Tp(1);
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
      return x > Tp(0);
    }
  };

  void operator()()
  {
    clear();

    boost::mt19937 generator;
    generator.seed(utils::random_seed());
    
    size_type bitext_id;
    for (;;) {
      queue_.pop(bitext_id);
      
      if (bitext_id == size_type(-1)) break;
      
      const sentence_type& source = bitexts_[bitext_id].source_;
      const sentence_type& target = bitexts_[bitext_id].target_;
      
      if (! source.empty() && ! target.empty()) {

#if 0
	std::cerr << "source: " << source << std::endl
		  << "target: " << target << std::endl;
#endif
	
	std::pair<double, double> errors = lexicon_(source, target, theta_, gradient_, htanh(), dhtanh(),
						    generator);
	  
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

  const bitext_set_type& bitexts_;
  const model_type& theta_;
  
  queue_type&            queue_;
  counter_type&          counter_;
  
  lexicon_type lexicon_;

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
		  const dictionary_type& dict_source_target,
		  const dictionary_type& dict_target_source,
		  model_type& theta)
{
  typedef TaskAccumulate task_type;
  typedef std::vector<task_type, std::allocator<task_type> > task_set_type;

  typedef task_type::size_type size_type;

  typedef std::vector<size_type, std::allocator<size_type> > id_set_type;
  
  task_type::queue_type   mapper(256 * threads);
  task_type::counter_type reducer;
  
  task_set_type tasks(threads, task_type(bitexts,
					 dict_source_target,
					 dict_target_source,
					 theta,
					 mapper,
					 reducer));
  
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

    id_set_type::const_iterator biter     = ids.begin();
    id_set_type::const_iterator biter_end = ids.end();

    double error = 0.0;
    double classification = 0.0;
    size_type samples = 0;
    size_type num_bitext = 0;

    utils::resource start;
    
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

    utils::resource end;
    
    if (debug && ((num_bitext / DEBUG_DOT) % DEBUG_WRAP))
      std::cerr << std::endl;
    if (debug)
      std::cerr << "# of bitexts: " << num_bitext << std::endl;

    if (debug)
      std::cerr << "reconstruction error: " << (error / samples) << std::endl
		<< "classification error: " << (classification / samples) << std::endl;
    
    if (debug)
      std::cerr << "cpu time:    " << end.cpu_time() - start.cpu_time() << std::endl
		<< "user time:   " << end.user_time() - start.user_time() << std::endl;

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


void read_bitext(const path_type& source_file,
		 const path_type& target_file,
		 bitext_set_type& bitexts,
		 dictionary_type& dict_source_target,
		 dictionary_type& dict_target_source)
{
  typedef cicada::Vocab vocab_type;
  typedef cicada::Symbol word_type;
  typedef bitext_type::sentence_type sentence_type;

  bitexts.clear();
  dict_source_target.clear();
  dict_target_source.clear();

  utils::compress_istream src(source_file, 1024 * 1024);
  utils::compress_istream trg(target_file, 1024 * 1024);
  
  sentence_type source;
  sentence_type target;
  
  while (src && trg) {
    src >> source;
    trg >> target;
    
    if (! src || ! trg) break;
    
    bitexts.push_back(bitext_type(source, target));
    
    sentence_type::const_iterator siter_begin = source.begin();
    sentence_type::const_iterator siter_end   = source.end();
    sentence_type::const_iterator titer_begin = target.begin();
    sentence_type::const_iterator titer_end   = target.end();
    
    {
      dictionary_type::dict_type& dict = dict_source_target[vocab_type::EPSILON];
      
      for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	++ dict[*titer];
      
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter) {
	dictionary_type::dict_type& dict = dict_source_target[*siter];
	
	for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
	  ++ dict[*titer];
      }
    }

    {
      dictionary_type::dict_type& dict = dict_target_source[vocab_type::EPSILON];
      
      for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	++ dict[*siter];
      
      for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer) {
	dictionary_type::dict_type& dict = dict_target_source[*titer];
	
	for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
	  ++ dict[*siter];
      }
    }
  }
  
  if (src || trg)
    throw std::runtime_error("# of sentnces do not match");

  if (cutoff > 1) {
    typedef dictionary_type::dict_type::count_set_type word_set_type;
    
    word_set_type words_source;
    word_set_type words_target;
    
    const word_set_type& counts_source = dict_target_source[vocab_type::EPSILON].counts_;
    const word_set_type& counts_target = dict_source_target[vocab_type::EPSILON].counts_;
    
    word_set_type::const_iterator siter_end = counts_source.end();
    for (word_set_type::const_iterator siter = counts_source.begin(); siter != siter_end; ++ siter)
      if (siter->second >= cutoff)
	words_source.insert(*siter);
    
    word_set_type::const_iterator titer_end = counts_target.end();
    for (word_set_type::const_iterator titer = counts_target.begin(); titer != titer_end; ++ titer)
      if (titer->second >= cutoff)
	words_target.insert(*titer);
    
    dictionary_type dict_source_target_new;
    dictionary_type dict_target_source_new;
    
    for (word_type::id_type i = 0; i != dict_source_target.dicts_.size(); ++ i)
      if (dict_source_target.dicts_.exists(i)) {
	word_type source(i);
	
	if (source != vocab_type::EPSILON && words_source.find(source) == words_source.end())
	  source = vocab_type::UNK;
	
	dictionary_type::dict_type& dict = dict_source_target_new[source];
	
	word_set_type::const_iterator titer_end = dict_source_target[i].counts_.end();
	for (word_set_type::const_iterator titer = dict_source_target[i].counts_.begin(); titer != titer_end; ++ titer)
	  if (words_target.find(titer->first) == words_target.end())
	    dict[vocab_type::UNK] += titer->second;
	  else
	    dict[titer->first] += titer->second;
      }
    
    for (word_type::id_type i = 0; i != dict_target_source.dicts_.size(); ++ i)
      if (dict_target_source.dicts_.exists(i)) {
	word_type target(i);
	
	if (target != vocab_type::EPSILON && words_target.find(target) == words_target.end())
	  target = vocab_type::UNK;
	
	dictionary_type::dict_type& dict = dict_target_source_new[target];
	
	word_set_type::const_iterator siter_end = dict_target_source[i].counts_.end();
	for (word_set_type::const_iterator siter = dict_target_source[i].counts_.begin(); siter != siter_end; ++ siter)
	  if (words_source.find(siter->first) == words_source.end())
	    dict[vocab_type::UNK] += siter->second;
	  else
	    dict[siter->first] += siter->second;
      }

    dict_source_target.swap(dict_source_target_new);
    dict_target_source.swap(dict_target_source_new);
    
    bitext_set_type::iterator biter_end = bitexts.end();
    for (bitext_set_type::iterator biter = bitexts.begin(); biter != biter_end; ++ biter) {

      sentence_type::iterator siter_end = biter->source_.end();
      for (sentence_type::iterator siter = biter->source_.begin(); siter != siter_end; ++ siter)
	if (words_source.find(*siter) == words_source.end())
	  *siter = vocab_type::UNK;

      sentence_type::iterator titer_end = biter->target_.end();
      for (sentence_type::iterator titer = biter->target_.begin(); titer != titer_end; ++ titer)
	if (words_target.find(*titer) == words_target.end())
	  *titer = vocab_type::UNK;	
    }
    
  }

  dict_source_target[vocab_type::BOS][vocab_type::BOS] = 1;
  dict_source_target[vocab_type::EOS][vocab_type::EOS] = 1;
  
  dict_source_target.initialize();
  dict_target_source.initialize();
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
    
    ("output-model", po::value<path_type>(&output_model_file), "output model parameter")
    
    ("alpha",     po::value<double>(&alpha)->default_value(alpha),      "parameter for reconstruction error")
    ("beta",      po::value<double>(&beta)->default_value(beta),        "parameter for classificaiton error")
    
    ("dimension-embedding", po::value<int>(&dimension_embedding)->default_value(dimension_embedding), "dimension for embedding")
    ("dimension-hidden",    po::value<int>(&dimension_hidden)->default_value(dimension_hidden),       "dimension for hidden layer")
    ("window",              po::value<int>(&window)->default_value(window),                           "context window size")
    
    ("optimize-sgd",     po::bool_switch(&optimize_sgd),     "SGD (Pegasos) optimizer")
    ("optimize-adagrad", po::bool_switch(&optimize_adagrad), "AdaGrad optimizer")
    
    ("iteration",         po::value<int>(&iteration)->default_value(iteration),   "max # of iterations")
    ("batch",             po::value<int>(&batch_size)->default_value(batch_size), "mini-batch size")
    ("cutoff",            po::value<int>(&cutoff)->default_value(cutoff),         "cutoff count for vocabulary (<= 1 to keep all)")
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
