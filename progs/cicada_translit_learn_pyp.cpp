//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

//
// pyp-transliteration
//

// this is a simple monotone transliteration model...
// we have only one restaurant + two priors for source/target, resampled from the model.
//
// the lambda sampling is taken from:
// @InProceedings{mochihashi-yamada-ueda:2009:ACLIJCNLP,
//   author    = {Mochihashi, Daichi  and  Yamada, Takeshi  and  Ueda, Naonori},
//   title     = {Bayesian Unsupervised Word Segmentation with Nested Pitman-Yor Language Modeling},
//   booktitle = {Proceedings of the Joint Conference of the 47th Annual Meeting of the ACL and the 4th International Joint Conference on Natural Language Processing of the AFNLP},
//   month     = {August},
//   year      = {2009},
//   address   = {Suntec, Singapore},
//   publisher = {Association for Computational Linguistics},
//   pages     = {100--108},
//   url       = {http://www.aclweb.org/anthology/P/P09/P09-1012}
// }
//


#include <map>
#include <iterator>

#include <cicada/sentence.hpp>
#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include "utils/vector2.hpp"
#include "utils/piece.hpp"
#include "utils/resource.hpp"
#include "utils/program_options.hpp"
#include "utils/compress_stream.hpp"
#include "utils/mathop.hpp"
#include "utils/bithack.hpp"
#include "utils/lockfree_list_queue.hpp"
#include "utils/chinese_restaurant_process.hpp"
#include "utils/unordered_map.hpp"
#include "utils/unordered_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/sampler.hpp"
#include "utils/repository.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/simple_vector.hpp"

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <boost/thread.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/fusion/tuple.hpp>

struct PYPTranslit
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef std::string  word_type;
  typedef std::string  phrase_type;
  typedef utils::piece piece_type;
  
  typedef utils::piece segment_type;
  struct segment_pair_type
  {
    segment_type source;
    segment_type target;
    
    segment_pair_type() : source(), target() {}
    segment_pair_type(const segment_type& __source,
		      const segment_type& __target)
      : source(__source), target(__target) {}
  };
  
  typedef std::vector<segment_pair_type, std::allocator<segment_pair_type> > derivation_type;
  
  struct phrase_pair_type
  {
    phrase_type source;
    phrase_type target;
    
    phrase_pair_type() : source(), target() {}
    phrase_pair_type(const phrase_type& __source,
		     const phrase_type& __target)
      : source(__source), target(__target) {}
    
    friend
    bool operator==(const phrase_pair_type& x, const phrase_pair_type& y)
    {
      return (x.source < y.source || (!(y.source < x.source) && x.target < y.target));
    }
    
    friend
    size_t hash_value(phrase_pair_type const& x)
    {
      typedef utils::hashmurmur<size_t> hasher_type;
      
      return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
    }
  };
  
  typedef utils::chinese_restaurant_process<phrase_pair_type, boost::hash<phrase_pair_type>, std::equal_to<phrase_pair_type>, std::allocator<phrase_pair_type > > table_type;
  
  struct length_base_type
  {
    typedef utils::vector2<double, std::allocator<double> > matrix_type;
    
    length_base_type(const double& __lambda,
		     const double& __strength_shape,
		     const double& __strength_rate)
      : lambda(__lambda),
	strength_shape(__strength_shape),
	strength_rate(__strength_rate)
    {
      initialize(__lambda);
    }
    
    double prob(const size_type source_size, const size_type target_size) const
    {
      if (source_size == 0 || target_size == 0) return 0.0;
      
      if (source_size < 32 && target_size < 32)
	return matrix(source_size, target_size);
      else
	return std::exp(utils::mathop::log_poisson(target_size, lambda * source_size))
    }
    
    void initialize(const double& __lambda)
    {
      lambda = __lambda;
      
      matrix.clear();
      matrix.resize(32, 32);
      
      for (size_type source_size = 1; source_size != 32; ++ source_size)
	for (size_type target_size = 1; target_size != 32; ++ target_size)
	  matrix(source_size, target_size) = std::exp(utils::mathop::log_poisson(target_size, lambda * source_size));
    }
    
    matrix_type matrix;
    
    double lambda;
    double strength_shape;
    double strength_rate;
  };
  
  
  PYPTranslit(const double __discount,
	      const double __strength,
	      const double __discount_alpha,
	      const double __discount_beta,
	      const double __strength_shape,
	      const double __strength_rate,
	      const double __lambda,
	      const double __lambda_strength_shape,
	      const double __lambda_strength_rate)
    : table(__discount, __strength, __discount_alpha, __discount_beta, __strength_shape, __strength_rate),
      base(__lambda, __lambda_strength_shape, __lambda_strength_rate)
  { }
  
  
  template <typename Iterator, typename Sampler>
  bool increment(Iterator first, Iterator last, Sampler& sampler, const double temperature=1.0)
  {
    for (/**/; first != last; ++ first)
      increment(first->source, first->target, sampler, temperature);
  }

  template <typename Iterator, typename Sampler>
  bool decrement(Iterator first, Iterator last, Sampler& sampler)
  {
    for (/**/; first != last; ++ first)
      decrement(first->source, first->target, sampler);
  }
  
  template <typename Sampler>
  bool increment(const piece_type& source, const piece_type& target, Sampler& sampler, const double temperature=1.0)
  {
    return table.increment(phrase_pair_type(source, target), base.prob(source.size(), target.size()), sampler, temperature);
  }
  
  template <typename Sampler>
  bool decrement(const piece_type& source, const piece_type& target, Sampler& sampler)
  {
    return table.increment(phrase_pair_type(source, target), sampler);
  }
  
  double prob(const piece_type& source, const piece_type& target) const
  {
    return table.prob(phrase_pair_type(source, target), base.prob(source.size(), target.size()));
  }
  
  double log_likelihood() const
  {
    return table.log_likelihood();
  }
  
  double log_likelihood(const double& discount, const double& srength) const
  {
    if (strength <= - discount) return - std::numeric_limits<double>::infinity();
    
    return table.log_likelihood(discount, strength);
  }
  
  template <typename Sampler>
  void sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    // we will resample lambda in base...
    double ratio = 0.0;
    typename table_type::const_iterator titer_end = table.end();
    for (typename table_type::const_iterator titer = table.begin(); titer != titer_end; ++ titer)
      ratio += (double(titer->first.target.size()) / titer->first.source.size()) * titer->second.size_table();
    
    base.lambda = sampler.gamma(base.strength_shape + ratio, base.strength_rate + table.size_table());
    
    table.sample_parameters(sampler, num_loop, num_iterations);
  }
  
  template <typename Sampler>
  void slice_sample_parameters(Sampler& sampler, const int num_loop = 2, const int num_iterations = 8)
  {
    // we will resample lambda in base...
    double ratio = 0.0;
    typename table_type::const_iterator titer_end = table.end();
    for (typename table_type::const_iterator titer = table.begin(); titer != titer_end; ++ titer)
      ratio += (double(titer->first.target.size()) / titer->first.source.size()) * titer->second.size_table();
    
    base.lambda = sampler.gamma(base.strength_shape + ratio, base.strength_rate + table.size_table());
    
    table.slice_sample_parameters(sampler, num_loop, num_iterations);
  }
  
  
public: 
  table_type        table;
  length_base_type  base;
};
