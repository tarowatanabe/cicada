// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__QUANTIZER__HPP__
#define __CICADA__QUANTIZER__HPP__ 1

#include <map>
#include <vector>
#include <algorithm>
#include <numeric>

#include <utils/mathop.hpp>

namespace cicada
{
  
  struct Quantizer
  {
    typedef size_t size_type;
    
    template <typename Counts, typename CodeBook, typename CodeMap>
    static inline
    void quantize(Counts& counts, CodeBook& codebook, CodeMap& codemap, const int debug=0)
    {
      typedef typename Counts::key_type    logprob_type;
      typedef typename Counts::mapped_type count_type;
      
      typedef std::pair<logprob_type, count_type> logprob_count_type;
      typedef std::vector<logprob_count_type, std::allocator<logprob_count_type> > logprob_set_type;
      typedef std::multimap<size_type, logprob_set_type, std::greater<size_type>, std::allocator<std::pair<const size_type, logprob_set_type> > > quantized_type;
      
      quantized_type quantized;
      codemap.clear();
      size_type num_center = 256;
      
      
      // zeros...
      {
	typename Counts::iterator citer = counts.begin();
	if (citer != counts.end() && citer->first <= boost::numeric::bounds<logprob_type>::lowest()) {
	  quantized.insert(std::make_pair(citer->second, logprob_set_type(1, *citer)));
	  counts.erase(citer);
	  -- num_center;
	}
      }
      
      // ones...?
      {
	typename Counts::iterator citer = counts.find(0.0);
	if (citer != counts.end()) {
	  quantized.insert(std::make_pair(citer->second, logprob_set_type(1, *citer)));
	  counts.erase(citer);
	  -- num_center;
	}
      }

      // high..
      if (! counts.empty()) {
	typename Counts::iterator citer = counts.end();
	-- citer;
	if (citer->first >= boost::numeric::bounds<logprob_type>::highest()) {
	  quantized.insert(std::make_pair(citer->second, logprob_set_type(1, *citer)));
	  counts.erase(citer);
	  -- num_center;
	}
      }
      
      // split...
      typename Counts::const_iterator citer_zero = counts.upper_bound(0.0);
      
      size_type num_backward = 0;
      for (typename Counts::const_iterator citer = citer_zero; citer != counts.end(); ++ citer, ++ num_backward);
      const size_type num_forward = counts.size() - num_backward;
      
      size_type num_center_forward = (num_backward == 0 ? num_center : size_t(double(num_center) * (double(num_forward) / counts.size())));
      if (num_backward > 0 && num_center_forward == num_center)
	-- num_center_forward;
      size_type num_center_backward = num_center - num_center_forward;
      
      logprob_set_type logprobs;
      
      if (num_forward > 0) {// forward...
	const size_type interval = num_forward / num_center_forward;
	std::vector<size_type, std::allocator<size_type> > bins(num_center_forward, interval);
	
	const size_type total_bins = std::accumulate(bins.begin(), bins.end(), size_type(0));
	if (num_forward > total_bins) {
	  for (size_type i = 0; i < num_forward - total_bins; ++ i)
	    ++ bins[i];
	}
	
	typename Counts::const_iterator citer_end = citer_zero;
	typename Counts::const_iterator citer = counts.begin();
	
	for (size_t i = 0; i < num_center_forward && citer != citer_end; ++ i) {
	  size_type num_data = 0;
	  logprobs.clear();
	  
	  for (size_type bin = 0; bin < bins[i] && citer != citer_end; ++ bin, ++ citer) {
	    logprobs.push_back(*citer);
	    num_data += citer->second;
	  }
	  quantized.insert(std::make_pair(num_data, logprobs));	  
	}
      }
    
      if (num_backward > 0) {// backward...
	const size_type interval = num_backward / num_center_backward;
	std::vector<size_type, std::allocator<size_type> > bins(num_center_backward, interval);
	
	const size_type total_bins = std::accumulate(bins.begin(), bins.end(), size_type(0));
	if (num_backward > total_bins) {
	  for (size_type i = 0; i < num_backward - total_bins; ++ i)
	    ++ bins[i];
	}
	
	typename Counts::const_iterator citer_end = counts.end();
	typename Counts::const_iterator citer = citer_zero;
	
	for (size_type i = 0; i < num_center_backward && citer != citer_end; ++ i) {
	  size_type num_data = 0;
	  logprobs.clear();
	  
	  for (size_type bin = 0; bin < bins[i] && citer != citer_end; ++ bin, ++ citer) {
	    logprobs.push_back(*citer);
	    num_data += citer->second;
	  }
	  quantized.insert(std::make_pair(num_data, logprobs));
	}
      }
      
      // we will now create codebook...
      // id assignment is performed in sorted order...
      size_type code_id = 0;
      typename quantized_type::const_iterator qiter_end = quantized.end();
      for (typename quantized_type::const_iterator qiter = quantized.begin(); qiter != qiter_end; ++ qiter, ++ code_id) {
	
	if (code_id > 255)
	  throw std::runtime_error("invalid code-id?");

	const size_type& count = qiter->first;
	const logprob_set_type& logprobs = qiter->second;
	
	double logsum = boost::numeric::bounds<double>::lowest();
	typename logprob_set_type::const_iterator liter_end = logprobs.end();
	for (typename logprob_set_type::const_iterator liter = logprobs.begin(); liter != liter_end; ++ liter) {
	  logsum = utils::mathop::logsum(logsum, double(liter->first) + utils::mathop::log(double(liter->second)));
	  codemap[liter->first] = code_id;
	}

	if (debug)
	  std::cerr << "code: " << code_id << " count: " << count << std::endl;
	
	codebook[code_id] = logsum - utils::mathop::log(double(count));
      }
    }
  };
  
};

#endif
