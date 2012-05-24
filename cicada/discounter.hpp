// -*- mode: c++ -*-

//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__DISCOUNTER__HPP__
#define __CICADA__DISCOUNTER__HPP__ 1

namespace cicada
{
  struct Discounter
  {
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    typedef uint64_t           count_type;
    
    Discounter()
      : discount1(-1.0), discount2(-1.0), discount3plus(-1.0), modified(false) {}
    template <typename Iterator>
    Discounter(Iterator first, Iterator last)
      : discount1(-1.0), discount2(-1.0), discount3plus(-1.0), modified(false) { estimate(first, last); }
    
    template <typename Iterator>
    void estimate(Iterator first, Iterator last)
    {
      clear();
      
      if (first == last) return;
      Iterator iter1 = first;
      ++ first;
      if (first == last) return;
      Iterator iter2 = first;
      estimate_kn_smooth(iter1, iter2);
      ++ first;
      if (first == last) return;
      Iterator iter3 = first;
      ++ first;
      if (first == last) return;
      Iterator iter4 = first;
      estimate_modified_kn_smooth(iter1, iter2, iter3, iter4);
    }
    
    template <typename Iterator>
    void estimate_kn_smooth(Iterator iter1, Iterator iter2)
    {
      if (iter1->second == 0 || iter2->second == 0)
	return;
      
      discount1 = double(iter1->second) / (double(iter1->second) + 2.0 * iter2->second);
      mincount1 = iter1->first;
      mincount2 = iter2->first;
      modified = false;
    }
    
    template <typename Iterator>
    void estimate_modified_kn_smooth(Iterator iter1, Iterator iter2, Iterator iter3, Iterator iter4)
    {
      if (iter1->second == 0 || iter2->second == 0 || iter3->second == 0 || iter4->second == 0)
	return;
      
      const double y = double(iter1->second) / (double(iter1->second) + 2.0 * iter2->second);
      
      discount1     = 1.0 - 2.0 * y * iter2->second / iter1->second;
      discount2     = 2.0 - 3.0 * y * iter3->second / iter2->second;
      discount3plus = 3.0 - 4.0 * y * iter4->second / iter3->second;
      
      mincount1 = iter1->first;
      mincount2 = iter2->first;
      mincount3 = iter3->first;
      mincount4 = iter4->first;
      modified = true;
    }
    
    double lower_order_weight(const count_type total, const count_type observed, const count_type min2, const count_type min3) const
    { 
      if (modified)
	return (discount1 >= 0.0 && discount2 >= 0.0 && discount3plus > 0.0
		? (discount1 * (observed - min2) + discount2 * (min2 - min3) + discount3plus * min3) / total // modified-KN smoothing
		: double(observed) / (total + observed)); // witten-bell smoothing
      else
	return (discount1 >= 0.0
		? (discount1 * observed) / total    // KN smoothing
		: double(observed) / (total + observed)); // witten-bell smoothing
    }
    
    double discount(const count_type count, const count_type total, const count_type observed) const
    {
      if (count <= 0) return 1.0;
      
      if (modified) {
        if (discount1 >= 0.0 && discount2 >= 0.0 && discount3plus >= 0.0) { // modified-KN smoothing
          if (count == mincount1)
            return (double(count) - discount1) / count;
          else if (count == mincount2)
            return (double(count) - discount2) / count;
          else
            return (double(count) - discount3plus) / count;
        } else
          return (count < mincount1 ? 0.0 : double(total) / (total + observed)); // witten-bell smoothing
      } else
        return (discount1 >= 0.0
                ? (count < mincount1 ? 0.0 : (double(count) - discount1) / count)  // KN smoothing
                : (count < mincount1 ? 0.0 : double(total) / (total + observed))); // witten-bell smoothing
    }
    
    void clear()
    {
      // fall back to witten bell
      discount1 = -1.0;
      discount2 = -1.0;
      discount3plus = -1.0;
      modified = false;
    }
    
    double discount1;
    double discount2;
    double discount3plus;
    
    count_type mincount1;
    count_type mincount2;
    count_type mincount3;
    count_type mincount4;
    
    bool modified;
  };

};

#endif
