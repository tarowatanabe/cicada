
#include "feature/dependency.hpp"

#include "parameter.hpp"

#include "utils/piece.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/simple_vector.hpp"

#include <boost/fusion/tuple.hpp>
#include <boost/array.hpp>

namespace cicada
{
  namespace feature
  {
    class DependencyImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef Dependency::feature_set_type feature_set_type;
      typedef Dependency::attribute_set_type attribute_set_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      struct __attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
      {
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
      };

      DependencyImpl()
	: attr_dependency_pos("dependency-pos"),
	  attr_dependency_head("dependency-head"),
	  attr_dependency_dependent("dependency-dependent") {}
      
      
      attribute_type attr_dependency_pos;
      attribute_type attr_dependency_head;
      attribute_type attr_dependency_dependent;
    };

    Dependency::Dependency(const std::string& parameter)
      : pimpl(new impl_type()),
	
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "dependency")
	throw std::runtime_error("is this really dependency feature function? " + parameter);
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	
      }
      
    }
    
    Dependency::Dependency(const Dependency& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type()) {}

    Dependency::~Dependency() { if (pimpl) delete pimpl; }
    
    Dependency& Dependency::operator=(const Dependency& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
  };
};
