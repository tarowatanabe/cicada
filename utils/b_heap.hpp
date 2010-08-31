#ifndef __UTILS__B_HEAP__HPP__
#define __UTILS__B_HEAP__HPP__ 1

#include <vector>
#include <functional>

#include <utils/bithack.hpp>

#include <boost/integer/static_log2.hpp>

namespace utils
{
  
  // priority_queue like interface...
  template <typename _Tp, typename _Sequence=std::vector<_Tp>, typename _Compare=std::less<typename _Sequence::value_type>, size_t _PageSize=8 >
  class b_heap : public _Compare
  {
  public:
    typedef typename _Sequence::value_type                value_type;
    typedef typename _Sequence::reference                 reference;
    typedef typename _Sequence::const_reference           const_reference;
    typedef typename _Sequence::size_type                 size_type;
    typedef          _Sequence                            container_type;

    typedef typename _Sequence::const_iterator const_iterator;
    typedef typename _Sequence::iterator       iterator;
    
  protected:
    _Sequence c;
    
  public:
    b_heap(size_type x) : c(x + 1) { clear(); }
    b_heap() { c.push_back(value_type()); }
    b_heap(size_type x, const _Compare& __comp) : _Compare(__comp), c(x + 1) { clear(); }
    b_heap(const _Compare& __comp) : _Compare(__comp) { c.push_back(value_type()); }
    
    bool empty() const { return c.size() == 1; }
    size_type size() const { return c.size() - 1; }
    
    void clear()
    {
      c.clear();
      c.push_back(value_type());
    }
    
    void reserve(size_type x) { c.reserve(x + 1); }
    

    // a conventional interface...
    const_reference top()
    {
      return c[1];
    }
    
    void push(const value_type& x)
    {
      c.push_back(x);
      siftup(c.size() - 1);
    }
    
    void pop()
    {
      std::swap(c[1], c.back());
      c.pop_back();
      
      siftdown(1);
    }
    
    const_iterator begin() const { return c.begin() + 1; }
    const_iterator end() const { return c.end(); }
    
  private:
    void siftdown(size_type index)
    {
      while (not_leaf(index)) {
	size_type max_child_index = max_child(index);
	if (! _Compare::operator()(c[max_child_index], c[index])) {
	  std::swap(c[max_child_index], c[index]);
	  index = max_child_index;
	} else
	  return;
      }
    }

    void siftup(size_type index)
    {
      while (index != 1) {
	size_type parent = parent_index(index);
	
	if (_Compare::operator()(c[parent], c[index])) {
	  std::swap(c[parent], c[index]);
	  index = parent;
	} else
	  return;
      }
    }
    
    bool not_leaf(size_type index) const
    {
      const size_type first_child = first_child_index(index);
      return first_child < c.size();
    }
    
    size_type max_child(size_type index) const
    {
      const size_t first_index = first_child_index(index);
      const size_t end_index = first_index + 2;
      
      const_iterator min_element = std::max_element(c.begin() + first_index,
						    std::min(c.begin() + end_index,
							     c.end()),
						    static_cast<const _Compare&>(*this));
      return min_element - c.begin();
    }
    
    // some definitions...
    static const bool      bh_is_power2   = utils::bithack::static_is_power2<_PageSize>::result;
    static const size_type bh_next_power2 = utils::bithack::static_next_largest_power2<_PageSize>::result;
    static const size_type bh_power2      = (bh_is_power2 ? _PageSize : bh_next_power2);
    
    static const size_type bh_psize = (bh_power2 >= 8 ? bh_power2 : size_type(8));
    static const size_type bh_shift = utils::bithack::static_floor_log2<bh_psize>::result;
    static const size_type bh_mask  = bh_psize - 1;
    
    static const size_type bh_half   = bh_psize >> 1;
    static const size_type bh_hshift = bh_shift - 2;
    static const size_type bh_hmask  = bh_mask >> 1;

    
    static inline size_type bh_pg(size_type index)
    {
      return index >> bh_shift;
    }
    
    static inline size_type bh_po(size_type index)
    {
      return index & bh_mask;
    }
    
    static inline size_type parent_index(size_type idx)
    {
      size_type po = bh_po(idx);
      if (idx < bh_psize || po > 3)
	return (idx & ~bh_mask) | (po >> 1);
      else if (po < 2) {
	size_type ip = (idx - bh_psize) >> bh_shift;
	ip += (ip & ~bh_hmask);
	ip |= bh_psize >> 1;
	return ip;
      } else
	return idx - 2;
    }

    static inline size_type first_child_index(size_type idx)
    {
      if (idx > bh_mask && !(idx & (bh_mask - 1)))
	return idx + 2;
      else if (idx & (bh_psize >> 1)) {
	size_type i1 = (idx & ~bh_mask) >> 1;
	i1 |= idx & (bh_mask >> 1);
	i1 += 1;
	i1 <<= bh_shift;
	return i1;
      } else
	return idx + (idx & bh_mask);
    }
  };
};

#endif
