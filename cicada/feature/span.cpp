
#include <utility>
#include <memory>

#include "cicada/feature/span.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/span_node.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie.hpp"
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
      typedef cicada::SpanVector span_set_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<int, int> span_node_type;
      typedef std::vector<span_node_type, std::allocator<span_node_type> > span_node_set_type;

      typedef utils::chart<std::string, std::allocator<std::string> > label_chart_type;
      
#ifdef HAVE_TR1_UNORDERED_MAP
      typedef std::tr1::unordered_map<span_node_type, std::string, utils::hashmurmur<size_t>, std::equal_to<span_node_type>,
				      std::allocator<std::pair<const span_node_type, std::string> > > label_map_type;
#else
      typedef sgi::hash_map<span_node_type, std::string, utils::hashmurmur<size_t>, std::equal_to<span_node_type>,
			    std::allocator<std::pair<const span_node_type, std::string> > > label_map_type;
      
#endif

      SpanImpl()
	: forced_feature(false) {}
      
      
      void span_score(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      feature_set_type& features) const
      {
	if (label_map.empty() || spans_node.empty()) return;
	
	if (edge.head >= spans_node.size()) {
	  std::cerr << "WARNING: node id exceeds our limit!" << std::endl;
	  return;
	}

	
	span_node_type span = spans_node[edge.head];
	int* context = reinterpret_cast<int*>(state);
	context[0] = span.first;
	context[1] = span.second;
	
	std::string rule_string = "span:" + span_label(span) + '(';
	
	int pos_non_terminal = 0;
	phrase_type::const_iterator piter_end = edge.rule->source.end();
	for (phrase_type::const_iterator piter = edge.rule->source.begin(); piter != piter_end; ++ piter)
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

      const std::string& span_label(const span_node_type& span) const
      {
	if (! label_chart.empty())
	  return label_chart(span.first, span.second);
	
	const_cast<label_chart_type&>(label_chart).reserve(spans_node.back().second + 1);
	const_cast<label_chart_type&>(label_chart).resize(spans_node.back().second + 1);
	
	for (int length = 1; length != label_chart.size(); ++ length)
	  for (int first = 0; first + length != label_chart.size(); ++ first) {
	    const int last = first + length;
	    
	    std::string& label = const_cast<std::string&>(label_chart(first, last));
	    
	    // exact match...
	    label_map_type::const_iterator niter = label_map.find(std::make_pair(first, last));
	    if (niter != label_map.end()) {
	      label = niter->second;
	      continue;
	    }
	    
	    // try binary combination...
	    for (int middle = first + 1; middle != last; ++ middle) {
	      label_map_type::const_iterator piter = label_map.find(std::make_pair(first, middle));
	      label_map_type::const_iterator niter = label_map.find(std::make_pair(middle, last));
	      
	      if (piter != label_map.end() && niter != label_map.end()) {
		label = '[' + strip_label(piter->second) + '+' + strip_label(niter->second) + ']';
		break;
	      }
	    }
	    
	    if (! label.empty()) continue;
	    
	    // try right-substitution...
	    for (int last_super = last + 1; last_super < label_chart.size(); ++ last_super) {
	      label_map_type::const_iterator siter = label_map.find(std::make_pair(first, last_super));
	      label_map_type::const_iterator riter = label_map.find(std::make_pair(last, last_super));
	      
	      if (siter != label_map.end() && riter != label_map.end()) {
		label = '[' + strip_label(siter->second) + '/' + strip_label(riter->second) + ']';
		break;
	      }
	    }

	    if (! label.empty()) continue;

	    // try left-subtraction...
	    for (int first_super = first - 1; first_super >= 0; -- first_super) {
	      label_map_type::const_iterator siter = label_map.find(std::make_pair(first_super, last));
	      label_map_type::const_iterator liter = label_map.find(std::make_pair(first_super, first));
	  
	      if (siter != label_map.end() && liter != label_map.end()) { 
		label = '[' + strip_label(liter->second) + '\\' + strip_label(siter->second) + ']';
		break;
	      }
	    }
	
	    if (! label.empty()) continue;

	    // try again with entries in label_chart!
	    // try binary combination...
	    for (int middle = first + 1; middle != last; ++ middle) {
	      const std::string& prev = label_chart(first, middle);
	      const std::string& next = label_chart(middle, last);
	      
	      if (! prev.empty() && ! next.empty()) {
		label = '[' + strip_label(prev) + "++" + strip_label(next) + ']';
		break;
	      }
	    }

	    if (! label.empty()) continue;
	
	    // try right-substitution...
	    for (int last_super = last + 1; last_super < label_chart.size(); ++ last_super) {
	      const std::string& super = label_chart(first, last_super);
	      const std::string& right = label_chart(last, last_super);

	      if (! super.empty() && ! right.empty()) {
		label = '[' + strip_label(super) + "//" + strip_label(right) + ']';
		break;
	      }
	    }

	    if (! label.empty()) continue;

	    // try left-subtraction...
	    for (int first_super = first - 1; first_super >= 0; -- first_super) {
	      const std::string& super = label_chart(first_super, last);
	      const std::string& left  = label_chart(first_super, first);
	  
	      if (! super.empty() && ! left.empty()) { 
		label = '[' + strip_label(left) + "\\\\" + strip_label(super) + ']';
		break;
	      }
	    }

	    if (! label.empty()) continue;
	
	    // try longest left and longest right
	    {
	      label_map_type::const_iterator liter = label_map.end();
	      for (int last_left = last - 1; first < last_left && liter == label_map.end(); -- last_left)
		liter = label_map.find(std::make_pair(first, last_left));
	  
	      label_map_type::const_iterator riter = label_map.end();
	      for (int first_right = first + 1; first_right < last && riter == label_map.end(); ++ first_right)
		riter = label_map.find(std::make_pair(first_right, last));
	      
	      if (liter != label_map.end() && riter != label_map.end()) {
		label = '[' + strip_label(liter->second) + ".." + strip_label(riter->second) + ']';
		break;
	      }
	    }
	    
	    if (! label.empty()) continue;
	    
	    static const std::string __default = "[x]";
	    
	    label = __default;
	  }
	
	
	return label_chart(span.first, span.second);
      }
      
      
      void assign(const hypergraph_type& hypergraph)
      {
	spans_node.clear();
	spans_node.reserve(hypergraph.nodes.size());
	spans_node.resize(hypergraph.nodes.size());
	
	cicada::span_node(hypergraph, spans_node);
	
	const span_node_type& span_goal = spans_node[hypergraph.goal];
	
	label_chart.clear();
      }

      void assign(const span_set_type& spans)
      {
	label_map.clear();
	
	span_set_type::const_iterator siter_end = spans.end();
	for (span_set_type::const_iterator siter = spans.begin(); siter != siter_end; ++ siter)
	  if (! siter->label.empty())
	    label_map[std::make_pair(siter->first, siter->last)] = siter->label;
      }

      label_chart_type label_chart;
      label_map_type   label_map;
      
      span_node_set_type spans_node;

      bool forced_feature;
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
    
    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }
    
    void Span::operator()(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates) const
    {
      const std::string& __feature_prefix = base_type::feature_name();
      for (feature_set_type::iterator fiter = features.begin(); fiter != features.end(); /**/)
	if (equal_prefix(__feature_prefix, fiter->first))
	  features.erase(fiter ++);
	else
	  ++ fiter;

      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->span_score(state, states, edge, features);
    }
    
    void Span::operator()(const state_ptr_type& state,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates) const
    {
      
    }

    void Span::assign(const hypergraph_type& hypergraph)
    {
      pimpl->assign(hypergraph);
    }

    void Span::assign(const span_set_type& spans)
    {
      pimpl->assign(spans);
    }

    void Span::initialize()
    {
      
    }
  };
};
