//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/span.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"

#include "utils/indexed_set.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/chart.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class SpanImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;

      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Lattice    lattice_type;
      typedef cicada::SpanVector span_set_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;

      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<int, int> span_type;

      typedef utils::chart<std::string, std::allocator<std::string> > label_chart_type;
      
#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<span_type, std::string, utils::hashmurmur<size_t>, std::equal_to<span_type>,
				      std::allocator<std::pair<const span_type, std::string> > > label_map_type;
#else
      typedef sgi::hash_map<span_type, std::string, utils::hashmurmur<size_t>, std::equal_to<span_type>,
			    std::allocator<std::pair<const span_type, std::string> > > label_map_type;
      
#endif

      SpanImpl()
	: forced_feature(false), attr_span_first("span-first"), attr_span_last("span-last") {}
      
      
      struct __rule_span : public boost::static_visitor<int>
      {
	int operator()(const attribute_set_type::int_type& x) const { return x; }
	template <typename Tp>
	int operator()(const Tp& x) const { throw std::runtime_error("no phrasal span with integer?"); }
      };
      
      int rule_span(const attribute_set_type& attrs, const attribute_type& attr) const
      {
	attribute_set_type::const_iterator iter = attrs.find(attr);
	if (iter == attrs.end())
	  throw std::runtime_error("no rule span attribute?");
	
	return boost::apply_visitor(__rule_span(), iter->second);
      }

      
      void span_score(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      feature_set_type& features) const
      {
	if (label_map.empty()) return;

	span_type span(rule_span(edge.attributes, attr_span_first),
		       rule_span(edge.attributes, attr_span_last));
	
	int* context = reinterpret_cast<int*>(state);
	context[0] = span.first;
	context[1] = span.second;
	
	std::string rule_string = "span:" + span_label(span) + '(';
	
	int pos_non_terminal = 0;
	phrase_type::const_iterator piter_end = edge.rule->rhs.end();
	for (phrase_type::const_iterator piter = edge.rule->rhs.begin(); piter != piter_end; ++ piter)
	  if (piter->is_non_terminal()) {
	    int antecedent_index = piter->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = pos_non_terminal;
	    
	    const int* antecedent_context = reinterpret_cast<const int*>(states[antecedent_index]);
	    
	    rule_string += span_label(std::make_pair(antecedent_context[0], antecedent_context[1]));
	    
	    span.first = antecedent_context[1];
	    
	    ++ pos_non_terminal;
	  } else if (span.first >= span.second) {
	    std::cerr << "WARNING: invalid span?" << std::endl;
	  } else {
	    rule_string += '<' + span_label(std::make_pair(span.first, span.first + 1)) + '>';
	    ++ span.first;
	  }
	
	rule_string += ')';
	
	if (forced_feature || feature_set_type::feature_type::exists(rule_string))
	  features[rule_string] += 1.0;
      }
      
      std::string strip_label(const std::string& label) const
      {
	const size_t size = label.size();
	if (size > 2 && label[0] == '[' && label[size - 1] == ']')
	  return label.substr(1, size - 2);
	else
	  return label;
      }

      const std::string& span_label(const span_type& span) const
      {
	std::string& label = const_cast<std::string&>(label_chart(span.first, span.second));
	
	if (! label.empty())
	  return label;
	    
	// exact match...
	label_map_type::const_iterator niter = label_map.find(std::make_pair(span.first, span.second));
	if (niter != label_map.end()) {
	  label = niter->second;
	  return label;
	}
	
	// try binary combination...
	for (int middle = span.first + 1; middle < span.second; ++ middle) {
	  label_map_type::const_iterator piter = label_map.find(std::make_pair(span.first, middle));
	  label_map_type::const_iterator niter = label_map.find(std::make_pair(middle, span.second));
	  
	  if (piter != label_map.end() && niter != label_map.end()) {
	    label = '[' + strip_label(piter->second) + '+' + strip_label(niter->second) + ']';
	    return label;
	  }
	}
	    
	// try right-substitution...
	for (int last_super = span.second + 1; last_super < static_cast<int>(label_chart.size()); ++ last_super) {
	  label_map_type::const_iterator siter = label_map.find(std::make_pair(span.first, last_super));
	  label_map_type::const_iterator riter = label_map.find(std::make_pair(span.second, last_super));
	  
	  if (siter != label_map.end() && riter != label_map.end()) {
	    label = '[' + strip_label(siter->second) + '/' + strip_label(riter->second) + ']';
	    return label;
	  }
	}
	
	// try left-subtraction...
	for (int first_super = span.first - 1; first_super >= 0; -- first_super) {
	  label_map_type::const_iterator siter = label_map.find(std::make_pair(first_super, span.second));
	  label_map_type::const_iterator liter = label_map.find(std::make_pair(first_super, span.first));
	  
	  if (siter != label_map.end() && liter != label_map.end()) { 
	    label = '[' + strip_label(liter->second) + '\\' + strip_label(siter->second) + ']';
	    return label;
	  }
	}
	
	// try tripple combination...
	for (int middle1 = span.first + 1; middle1 < span.second; ++ middle1)
	  for (int middle2 = middle1 + 1; middle2 < span.second; ++ middle2) {
	    label_map_type::const_iterator iter1 = label_map.find(std::make_pair(span.first, middle1));
	    label_map_type::const_iterator iter2 = label_map.find(std::make_pair(middle1, middle2));
	    label_map_type::const_iterator iter3 = label_map.find(std::make_pair(middle2, span.second));
	    
	    if (iter1 != label_map.end() && iter2 != label_map.end() && iter3 != label_map.end()) {
	      label = '[' + strip_label(iter1->second) + '+' + strip_label(iter2->second) + '+' + strip_label(iter3->second) + ']';
	      return label;
	    }
	  }
	
	// try longest left and longest right
	{
	  label_map_type::const_iterator liter = label_map.end();
	  for (int last_left = span.second - 1; span.first < last_left && liter == label_map.end(); -- last_left)
	    liter = label_map.find(std::make_pair(span.first, last_left));
	  
	  label_map_type::const_iterator riter = label_map.end();
	  for (int first_right = span.first + 1; first_right < span.second && riter == label_map.end(); ++ first_right)
	    riter = label_map.find(std::make_pair(first_right, span.second));
	  
	  if (liter != label_map.end() && riter != label_map.end()) {
	    label = '[' + strip_label(liter->second) + ".." + strip_label(riter->second) + ']';
	    return label;
	  }
	}
	
	static const std::string __default = "[x]";
	label = __default;
	return label;
      }
      
      
      void assign(const hypergraph_type& hypergraph,
		  const lattice_type& lattice,
		  const span_set_type& spans)
      {
	span_type span_goal(std::numeric_limits<int>::max(), std::numeric_limits<int>::min());
	
	hypergraph_type::node_type::edge_set_type::const_iterator eiter_end = hypergraph.nodes[hypergraph.goal].edges.end();
	for (hypergraph_type::node_type::edge_set_type::const_iterator eiter = hypergraph.nodes[hypergraph.goal].edges.begin(); eiter != eiter_end; ++ eiter) {
	  const hypergraph_type::edge_type& edge = hypergraph.edges[*eiter];
	  
	  span_goal.first  = utils::bithack::min(span_goal.first,  rule_span(edge.attributes, attr_span_first));
	  span_goal.second = utils::bithack::max(span_goal.second, rule_span(edge.attributes, attr_span_last));
	}
	
	label_chart.clear();
	label_chart.reserve(span_goal.second + 1);
	label_chart.resize(span_goal.second + 1);
	
	label_map.clear();
	span_set_type::const_iterator siter_end = spans.end();
	for (span_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)
	  if (! siter->label.empty())
	    label_map[std::make_pair(siter->first, siter->last)] = siter->label;
      }

      label_chart_type label_chart;
      label_map_type   label_map;

      bool forced_feature;

      attribute_type attr_span_first;
      attribute_type attr_span_last;
    };
    

    
    Span::Span(const std::string& parameter)
      : pimpl(new impl_type())
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "span")
	throw std::runtime_error("is this really span feature function? " + parameter);
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for span: " << piter->first << "=" << piter->second << std::endl;
      
      
      base_type::__state_size = sizeof(int) * 2;
      base_type::__feature_name = "span";
      base_type::__sparse_feature = true;
    }
    
    Span::~Span() { std::auto_ptr<impl_type> tmp(pimpl); }

    Span::Span(const Span& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    Span& Span::operator=(const Span& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Span::apply(state_ptr_type& state,
		     const state_ptr_set_type& states,
		     const edge_type& edge,
		     feature_set_type& features,
		     feature_set_type& estimates,
		     const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->span_score(state, states, edge, features);
    }
    
    void Span::apply_coarse(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
    {}
    void Span::apply_predict(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Span::apply_scan(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  const int dot,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {}
    void Span::apply_complete(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }


    void Span::assign(const size_type& id,
		      const hypergraph_type& hypergraph,
		      const lattice_type& lattice,
		      const span_set_type& spans,
		      const sentence_set_type& targets,
		      const ngram_count_set_type& ngram_counts)
    {
      pimpl->assign(hypergraph, lattice, spans);
    }

    void Span::initialize()
    {
      
    }
  };
};
