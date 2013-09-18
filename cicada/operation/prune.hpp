// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PRUNE__HPP__
#define __CICADA__OPERATION__PRUNE__HPP__ 1

#include <cicada/operation.hpp>

#include <utils/sampler.hpp>

namespace cicada
{
  namespace operation
  {
    class Prune : public Operation
    {
    private:
      typedef utils::sampler<boost::mt19937> sampler_type;
      
    public:
      Prune(const std::string& parameter, const int __debug);

      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);

      const weights_path_type* weights;
      const weight_set_type*   weights_assigned;
  
      size_t kbest;
      size_t edge;
      double beam;
      double density;
      double scale;
      bool   sample;
      bool   uniform;

      sampler_type sampler;
  
      bool weights_one;
      bool weights_fixed;

      feature_set_type weights_extra;

      bool semiring_tropical;
      bool semiring_logprob;
      bool semiring_log;
  
      int debug;
    };

  };
};


#endif
