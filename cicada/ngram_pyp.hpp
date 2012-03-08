//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_PYP__HPP__
#define __CICADA__NGRAM_PYP__HPP__ 1

// PYPLM!

#include <stdint.h>

#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/packed_vector.hpp>
#include <utils/succinct_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/array_power2.hpp>
#include <utils/spinlock.hpp>
#include <utils/bithack.hpp>
#include <utils/mathop.hpp>

namespace cicada
{
  class NGramPYP : public utils::hashmurmur<size_t>
  {
  public:
    typedef Symbol                  word_type;
    typedef Vocab                   vocab_type;

    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    typedef word_type::id_type id_type;
    typedef uint64_t           count_type;
    typedef float              logprob_type;
    typedef double             prob_type;
    
    typedef boost::filesystem::path   path_type;
    
  private:
    typedef utils::hashmurmur<size_t> hasher_type;
    
    typedef utils::packed_vector_mapped<id_type, std::allocator<id_type> >       id_set_type;
    typedef utils::packed_vector_mapped<count_type, std::allocator<count_type> > count_set_type;
    typedef utils::succinct_vector_mapped<std::allocator<int32_t> >              position_set_type;
    typedef std::vector<size_type, std::allocator<size_type> >                   offset_set_type;
    typedef std::vector<double, std::allocator<double> >                         parameter_set_type;
    
    struct spinlock_type
    {
      typedef utils::spinlock             mutex_type;
      typedef mutex_type::scoped_lock     lock_type;
      typedef mutex_type::scoped_try_lock trylock_type;
      
      spinlock_type() : mutex() {}
      spinlock_type(const spinlock_type&) {}
      spinlock_type& operator=(const spinlock_type&) { return *this; }
      
      mutex_type mutex;
    };
    
    struct cache_pos_type
    {
      size_type pos;
      size_type pos_next;
      id_type id;
      
      cache_pos_type() : pos(size_type(-1)), pos_next(size_type(-1)), id(id_type(-1)) {}
    };
    
    typedef utils::array_power2<cache_pos_type, 1024 * 128, std::allocator<cache_pos_type> > cache_pos_set_type;

  public:
    NGramPYP() {}
    NGramPYP(const path_type& path)  { open(path); }

    path_type path() const { return index_.path().parent_path(); }
    size_type size() const { return index_.size(); }
    bool empty() const { return index_.empty(); }

    bool is_open() const { return index_.is_open(); }
    
    void open(const path_type& path);
    void close() { clear(); }
    void clear()
    {
      index_.clear();
      count_.clear();
      total_.clear();
      positions_.clear();
      offsets_.clear();
      vocab_.clear();
      caches_pos_.clear();
    }

  public:
    static NGramPYP& create(const path_type& path);
    
  public:
    const vocab_type& vocab() const { return vocab_; }
    const int& order() const { return order_; }
    
    size_type parent(size_type pos) const
    {
      return (pos < offsets_[1] ? size_type(-1) : positions_.select(pos + 1 - offsets_[1], true) + (offsets_[1] + 1) - pos - 1);
    }
    
    bool has_child(size_type pos) const
    {
      return children_first(pos) != children_last(pos);
    }
    
    size_type children_first(size_type pos) const
    {
      if (pos == size_type(-1) || pos == 0)
	return (~pos) & offsets_[1];
      else
	return children_last(pos - 1);
    }
    
    size_type children_last(size_type pos) const
    {
      if (pos == size_type(-1) || pos >= offsets_[offsets_.size() - 2]) {
	const size_type is_root_mask = size_type(pos == size_type(-1)) - 1;
	return ((~is_root_mask) & offsets_[1]) | (is_root_mask & offsets_.back());
      }
      
      position_set_type::size_type last = positions_.select(pos + 2 - 1, false);
      
      const size_type last_mask = size_type(last == position_set_type::size_type(-1)) - 1;
      
      return ((~last_mask) & offsets_.back()) | (last_mask & ((last + 1 + offsets_[1] + 1) - (pos + 2)));
    }
    
    size_type find(size_type pos, const id_type& id) const
    {
      if (id == id_type(-1))
	return size_type(-1);
      else {
	// lock here!
	spinlock_type::trylock_type lock(const_cast<spinlock_type::mutex_type&>(spinlock_pos_.mutex));
	
	if (lock) {
	  const size_type cache_pos = hasher_type::operator()(id, pos) & (caches_pos_.size() - 1);
	  cache_pos_type& cache = const_cast<cache_pos_type&>(caches_pos_[cache_pos]);
	  if (cache.pos != pos || cache.id != id) {
	    cache.pos      = pos;
	    cache.id       = id;
	    cache.pos_next = __find(pos, id);
	  }
	  return cache.pos_next;
	} else
	  return __find(pos, id);
      }
    }

    template <typename Iterator>
    logprob_type operator()(Iterator first, Iterator last) const
    {
      return logprob(first, last);
    }

    template <typename Iterator>
    logprob_type logprob(Iterator first, Iterator last) const
    {
      return utils::mathop::log(prob(first, last));
    }

    template <typename Iterator>
    prob_type prob(Iterator first, Iterator last) const
    {
      if (first == last) return 1.0;
      
      typedef typename std::iterator_traits<Iterator>::value_type value_type;
      
      return __prob_dispatch(std::max(first, last - order_), last, value_type());
    }
    
    template <typename Iterator>
    bool exists(Iterator first, Iterator last) const
    {
      typedef std::reverse_iterator<Iterator> reverse_iterator;

      if (first == last) return false;
      
      // we will traverse from the back!
      reverse_iterator begin(last);
      reverse_iterator end(first);
      
      size_type node = size_type(-1);
      for (reverse_iterator iter = begin; iter != end; ++ iter) {
	node = find(node, *iter);
	
	if (node == size_type(-1))
	  return false;
      }
      return true;
    }
    
    template <typename Iterator>
    prob_type __prob_dispatch(Iterator first, Iterator last, id_type) const
    {
      typedef std::reverse_iterator<Iterator> reverse_iterator;
      
      const id_type word = *(last - 1);
      
      double p = __prob_dispatch(0, size_type(-1), word, p0_);
      
      // we will traverse from the back!
      reverse_iterator begin(last - 1);
      reverse_iterator end(first);
      
      size_type node = size_type(-1);
      int order = 1;
      for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
	node = find(node, *iter);
	
	if (node == size_type(-1))
	  return p;
	else
	  p = __prob_dispatch(order, node, word, p);
      }
      
      return p;
    }

    template <typename Iterator, typename _Word>
    prob_type __prob_dispatch(Iterator first, Iterator last, _Word) const
    {
      typedef std::reverse_iterator<Iterator> reverse_iterator;
      
      const id_type word = vocab_[*(last - 1)];
      
      double p = __prob_dispatch(0, size_type(-1), word, p0_);
      
      // we will traverse from the back!
      reverse_iterator begin(last - 1);
      reverse_iterator end(first);
      
      size_type node = size_type(-1);
      int order = 1;
      for (reverse_iterator iter = begin; iter != end; ++ iter, ++ order) {
	node = find(node, vocab_[*iter]);
	
	if (node == size_type(-1))
	  return p;
	else
	  p = __prob_dispatch(order, node, word, p);
      }
      
      return p;
    }
    
    prob_type __prob_dispatch(int order, size_type pos, const id_type& word, const double p0) const
    {
      const size_type pos_total = utils::bithack::branch(pos == size_type(-1), size_type(0), pos + 1);
      
      const count_type customer_size = total_[(pos_total << 1)];
      const count_type table_size    = total_[(pos_total << 1) + 1];
      
      const double r = table_size * discount_[order] + strength_[order];
      
      if (! customer_size)
	return r * p0 / (double(customer_size) + strength_[order]);
      
      const size_type pos_dish = find(pos, word);
      
      if (pos_dish == size_type(-1))
	return r * p0 / (double(customer_size) + strength_[order]);
      else {
	const count_type customer_dish_size = count_[(pos_dish << 1)];
	const count_type table_dish_size    = count_[(pos_dish << 1) + 1];
	
	return (double(customer_dish_size) - discount_[order] * table_dish_size + r * p0) / (double(customer_size) + strength_[order]);
      }
    }

  private:
    size_type __find(size_type pos, const id_type& id) const
    {
      const size_type pos_first = children_first(pos);
      const size_type pos_last  = children_last(pos);
      
      if (pos_first == pos_last) return size_type(-1);
      
      const size_type child = lower_bound(pos_first, pos_last, id);
      
      return utils::bithack::branch(child != pos_last && !(id < index_[child]), child, size_type(-1));
    }
    
    size_type lower_bound(size_type first, size_type last, const id_type& id) const
    {
      // otherwise...
      size_type length = last - first;
      
      if (length <= 128) {
	for (/**/; first != last && index_[first] < id; ++ first);
	return first;
      } else {
	while (length > 0) {
	  const size_t half  = length >> 1;
	  const size_t middle = first + half;
	  
	  const bool is_less = index_[middle] < id;
	  
	  first  = utils::bithack::branch(is_less, middle + 1, first);
	  length = utils::bithack::branch(is_less, length - half - 1, half);
	}
	return first;
      }
    }
    
  private:
    id_set_type       index_;
    count_set_type    count_;
    count_set_type    total_;
    position_set_type positions_;
    offset_set_type   offsets_;
    
    double     p0_;
    count_type counts0_;
    
    parameter_set_type discount_;
    parameter_set_type strength_;
    
    double discount_alpha_;
    double discount_beta_;
    double strength_shape_;
    double strength_rate_;
    
    vocab_type        vocab_;
    
    int               order_;
    
    spinlock_type      spinlock_pos_;
    cache_pos_set_type caches_pos_;
  };
  
};

#endif
