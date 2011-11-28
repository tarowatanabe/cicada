//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__MERT_KBEST_IMPL__HPP__
#define __CICADA__MERT_KBEST_IMPL__HPP__ 1

#include <vector>
#include <climits>

#include "cicada/weight_vector.hpp"
#include "cicada/dot_product.hpp"

#include "cicada_kbest_impl.hpp"

struct EnvelopeKBest
{
  typedef cicada::WeightVector<double> weight_set_type;
  
  EnvelopeKBest(const weight_set_type& __origin,
		const weight_set_type& __direction)
    
    : origin(__origin),
      direction(__direction) {}

  struct line_type
  {
    line_type() : x(- std::numeric_limits<double>::infinity()), m(0), y(0), hypothesis(0) {}
    line_type(const double& __m, const double& __y, const hypothesis_type& __hypothesis)
      : x(- std::numeric_limits<double>::infinity()), m(__m), y(__y), hypothesis(&__hypothesis) {}
    
    double x;
    double m;
    double y;
    
    const hypothesis_type* hypothesis;
  };
  
  typedef std::vector<line_type, std::allocator<line_type> > line_set_type;

  struct compare_slope
  {
    bool operator()(const line_type& x, const line_type& y) const
    {
      return x.m < y.m;
    }
  };

  void operator()(line_set_type& lines)
  {
    std::sort(lines.begin(), lines.end(), compare_slope());
  
    int j = 0;
    int K = lines.size();
      
    for (int i = 0; i < K; ++ i) {
      line_type line = lines[i];
      line.x = - std::numeric_limits<double>::infinity();
	
      if (0 < j) {
	if (lines[j - 1].m == line.m) { // parallel line...
	  if (line.y <= lines[j - 1].y) continue;
	  -- j;
	}
	while (0 < j) {
	  line.x = (line.y - lines[j - 1].y) / (lines[j - 1].m - line.m);
	  if (lines[j - 1].x < line.x) break;
	  -- j;
	}
	
	if (0 == j)
	  line.x = - std::numeric_limits<double>::infinity();
      }
      
      lines[j++] = line;
    }
    
    lines.resize(j);
  }

  template <typename Iterator>
  void operator()(Iterator first, Iterator last, line_set_type& lines)
  {
    lines.clear();
    
    for (/**/; first != last; ++ first) {
      const hypothesis_type& hyp = *first;
      
      const double m = cicada::dot_product(direction, hyp.features.begin(), hyp.features.end(), 0.0);
      const double y = cicada::dot_product(origin,    hyp.features.begin(), hyp.features.end(), 0.0);
      
      lines.push_back(line_type(m, y, hyp));
    }
    
    operator()(lines);
  }

  template <typename Iterator>
  void operator()(Iterator first1, Iterator last1, Iterator first2, Iterator last2, line_set_type& lines)
  {
    lines.clear();
    
    for (/**/; first1 != last1; ++ first1) {
      const hypothesis_type& hyp = *first1;
      
      for (Iterator iter2 = first2; iter2 != last2; ++ iter2) {
	const hypothesis_type& hyp2 = *iter2;
	
	const double m = (cicada::dot_product(direction, hyp.features.begin(), hyp.features.end(), 0.0)
			  - cicada::dot_product(direction, hyp2.features.begin(), hyp2.features.end(), 0.0));
	const double y = (cicada::dot_product(origin,    hyp.features.begin(), hyp.features.end(), 0.0)
			  - cicada::dot_product(origin,    hyp2.features.begin(), hyp2.features.end(), 0.0));
	
	lines.push_back(line_type(m, y, hyp));
      }
    }
    
    operator()(lines);
  }

  template <typename Iterator, typename OutputIterator>
  void operator()(Iterator first, Iterator last, OutputIterator result)
  {
    for (/**/; first != last; ++ first) {
      const hypothesis_type& hyp = *first;
      
      const double m = cicada::dot_product(direction, hyp.features.begin(), hyp.features.end(), 0.0);
      const double y = cicada::dot_product(origin,    hyp.features.begin(), hyp.features.end(), 0.0);
      
      *result = line_type(m, y, hyp);
      ++ result;
    }
  }

  template <typename Iterator, typename OutputIterator>
  void operator()(Iterator first1, Iterator last1, Iterator first2, Iterator last2, OutputIterator result)
  {
    for (/**/; first1 != last1; ++ first1) {
      const hypothesis_type& hyp = *first1;
      
      for (Iterator iter2 = first2; iter2 != last2; ++ iter2) {
	const hypothesis_type& hyp2 = *iter2;
	
	const double m = (cicada::dot_product(direction, hyp.features.begin(), hyp.features.end(), 0.0)
			  - cicada::dot_product(direction, hyp2.features.begin(), hyp2.features.end(), 0.0));
	const double y = (cicada::dot_product(origin,    hyp.features.begin(), hyp.features.end(), 0.0)
			  - cicada::dot_product(origin,    hyp2.features.begin(), hyp2.features.end(), 0.0));
	
	*result = line_type(m, y, hyp);
	++ result;
      }
    }
  }
  
  inline
  void operator()(const hypothesis_set_type& kbests, line_set_type& lines)
  {
    operator()(kbests.begin(), kbests.end(), lines);
  }

  inline
  void operator()(const hypothesis_set_type& kbests, const hypothesis_set_type& oracles, line_set_type& lines)
  {
    operator()(kbests.begin(), kbests.end(), oracles.begin(), oracles.end(), lines);
  }
  
  const weight_set_type& origin;
  const weight_set_type& direction;
};


#endif
