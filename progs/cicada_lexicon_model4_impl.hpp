//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA_LEXICON_MODEL4_IMPL__HPP__
#define __CICADA_LEXICON_MODEL4_IMPL__HPP__ 1

#include <numeric>
#include <set>

#include "cicada_lexicon_impl.hpp"

#include "utils/vector2_aligned.hpp"
#include "utils/vector3_aligned.hpp"
#include "utils/mathop.hpp"
#include "utils/aligned_allocator.hpp"

#include "kuhn_munkres.hpp"
#include "itg_alignment.hpp"
#include "dependency_hybrid.hpp"
#include "dependency_degree2.hpp"
#include "dependency_mst.hpp"

struct LearnModel4 : public LearnBase
{
  LearnModel4(const LearnBase& __base)
    : LearnBase(__base) {}
  
  
  void operator()(const sentence_type& source, const sentence_type& target, const alignment_type& alignment)
  {
    
  }
  
  void operator()(const sentence_type& source, const sentence_type& target)
  {
    
  }
};

struct LearnModel4Posterior : public LearnBase
{
  LearnModel4Posterior(const LearnBase& __base)
    : LearnBase(__base) {}
};

struct LearnModel4Symmetric : public LearnBase
{
  LearnModel4Symmetric(const LearnBase& __base)
    : LearnBase(__base) {}
};

struct LearnModel4SymmetricPosterior : public LearnBase
{
  LearnModel4SymmetricPosterior(const LearnBase& __base)
    : LearnBase(__base) {}

};

struct ViterbiModel4 : public ViterbiBase
{
  ViterbiModel4(const ttable_type& __ttable_source_target,
		const ttable_type& __ttable_target_source,
		const atable_type& __atable_source_target,
		const atable_type& __atable_target_source,
		const classes_type& __classes_source,
		const classes_type& __classes_target)
    : ViterbiBase(__ttable_source_target, __ttable_target_source,
		  __atable_source_target, __atable_target_source,
		  __classes_source, __classes_target) {}
  
};

struct PosteriorModel4 : public ViterbiBase
{
  
};

struct ITGModel4 : public ViterbiBase
{

};

struct MaxMatchModel4 : public ViterbiBase
{
  
};

struct DependencyModel4 : public ViterbiBase
{

  
};

template <typename Analyzer>
struct __DependencyModel4Base : public DependencyModel4
{
  
};

typedef __DependencyModel4Base<DependencyHybrid>            DependencyHybridModel4;
typedef __DependencyModel4Base<DependencyHybridSingleRoot>  DependencyHybridSingleRootModel4;
typedef __DependencyModel4Base<DependencyDegree2>           DependencyDegree2Model4;
typedef __DependencyModel4Base<DependencyDegree2SingleRoot> DependencyDegree2SingleRootModel4;
typedef __DependencyModel4Base<DependencyMST>               DependencyMSTModel4;
typedef __DependencyModel4Base<DependencyMSTSingleRoot>     DependencyMSTSingleRootModel4;


struct PermutationModel4 : public ViterbiBase
{
  
};

#endif
