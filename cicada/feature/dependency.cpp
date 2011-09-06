
#include <vector>

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

      typedef Dependency::feature_function_type feature_function_type;

      typedef feature_function_type::symbol_type symbol_type;
      typedef feature_function_type::vocab_type  vocab_type;

      typedef feature_function_type::hypergraph_type hypergraph_type;
      typedef feature_function_type::lattice_type    lattice_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef std::pair<int, int> lattice_edge_type;
      typedef std::vector<lattice_edge_type, std::allocator<lattice_edge_type> > lattice_edge_set_type;
      
      typedef std::pair<symbol_type, symbol_type> terminal_pos_type;
      typedef std::vector<terminal_pos_type, std::allocator<terminal_pos_type> > terminal_pos_set_type;
      
      struct __attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
      {
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
      };

      DependencyImpl(const int __order)
	: order(__order),
	  lattice(0),
	  forced_feature(false),
	  attr_dependency_pos("dependency-pos"),
	  attr_dependency_head("dependency-head"),
	  attr_dependency_dependent("dependency-dependent") {}

      void dependency_score(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features)
      {
	
	
      }	

      void clear()
      {
	
      }
      
      void assign(const lattice_type& __lattice)
      {
	lattice = &__lattice;

	edges.clear();
	terminals.clear();
	
	edges.push_back(std::make_pair(-1, 0));
	terminals.push_back(std::make_pair(vocab_type::EPSILON, vocab_type::X));
	
	for (size_type pos = 0; pos != lattice->size(); ++ pos) {
	  lattice_type::arc_set_type::const_iterator aiter_end = lattice->operator[](pos).end();
	  for (lattice_type::arc_set_type::const_iterator aiter = lattice->operator[](pos).begin(); aiter != aiter_end; ++ aiter) {
	    edges.push_back(std::make_pair(pos, pos + aiter->distance));
	    terminals.push_back(std::make_pair(aiter->label.terminal(), aiter->label.pos()));
	    
	    if (terminals.back().second.empty())
	      terminals.back().second = vocab_type::X;
	  }
	}
      }

      int order;

      const lattice_type*   lattice;
      lattice_edge_set_type edges;
      terminal_pos_set_type terminals;
      
      bool forced_feature;
      
      attribute_type attr_dependency_pos;
      attribute_type attr_dependency_head;
      attribute_type attr_dependency_dependent;
    };

    Dependency::Dependency(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "dependency")
	throw std::runtime_error("is this really dependency feature function? " + parameter);
      
      int order = 2;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for dependency: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (order <= 0)
	throw std::runtime_error("we do not support zero or negative orders");
      
      pimpl = new impl_type(order);
      
      base_type::__state_size = sizeof(int) * (1 << (order) - 1);
      base_type::__feature_name = "dependency";
      base_type::__sparse_feature = true;
    }
    
    Dependency::Dependency(const Dependency& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl)) {}
    
    Dependency::~Dependency() { if (pimpl) delete pimpl; }
    
    Dependency& Dependency::operator=(const Dependency& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    
    void Dependency::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->dependency_score(state, states, edge, features);
    }
    
    void Dependency::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {

    }
    
    void Dependency::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {

    }
    
    void Dependency::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {

    }
    
    void Dependency::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    void Dependency::initialize()
    {
      pimpl->clear();
    }
    
    void Dependency::assign(const size_type& id,
			    const hypergraph_type& hypergraph,
			    const lattice_type& lattice,
			    const span_set_type& spans,
			    const sentence_set_type& targets,
			    const ngram_count_set_type& ngram_counts)
    {
      pimpl->assign(lattice);
    }
    
  };
};
