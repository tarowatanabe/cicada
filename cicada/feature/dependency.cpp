
#include "feature/dependency.hpp"

namespace cicada
{
  namespace feature
  {
    
    Dependency::Dependency()
      : attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent") {}
    
    Dependency::Dependency(const Dependency& x)
      : base_type(static_cast<const base_type&>(x)),
	attr_dependency_pos("dependency-pos"),
	attr_dependency_head("dependency-head"),
	attr_dependency_dependent("dependency-dependent") {}
    
    Dependency& Dependency::operator=(const Dependency& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      return *this;
    }
    
  };
};
