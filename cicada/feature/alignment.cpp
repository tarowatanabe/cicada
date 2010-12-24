//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "feature/alignment.hpp"
#include "parameter.hpp"

namespace cicada
{
  namespace feature
  {
    namespace align
    {

      struct __attribute_integer : public boost::static_visitor<cicada::AttributeVector::int_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
      
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -2; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -2; }
      };
    
      static const cicada::Attribute __attr_source_size("source-size");
      static const cicada::Attribute __attr_target_size("target-size");
      static const cicada::Attribute __attr_source_position("source-position");
      static const cicada::Attribute __attr_target_position("target-position");

      // relative position...
      
      RelativePosition::RelativePosition(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: normalizers(), sentence(0), source_size(0), target_size(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (param.name() != "relative-position")
	  throw std::runtime_error("is this really relative position feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "cluster") == 0) {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (strcasecmp(piter->first.c_str(), "stemmer") == 0)
	    normalizers.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else
	    std::cerr << "WARNING: unsupported parameter for relative-position: " << piter->first << "=" << piter->second << std::endl;
	}
	
	__state_size = 0;
	__feature_name = "relative-position";
      }
      
      void RelativePosition::operator()(const feature_function_type& feature_function,
					state_ptr_type& state,
					const state_ptr_set_type& states,
					const edge_type& edge,
					feature_set_type& features,
					feature_set_type& estimates,
					const bool final) const
      {
	attribute_set_type::const_iterator siter = edge.attributes.find(__attr_source_position);
	attribute_set_type::const_iterator titer = edge.attributes.find(__attr_target_position);
	
	if (siter == edge.attributes.end() || titer == edge.attributes.end()) return;
	
	const int source_pos = boost::apply_visitor(__attribute_integer(), siter->second);
	const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	
	if (source_pos < 0 || target_pos < 0) return;
	
	const double value = std::fabs(static_cast<double>(source_pos) / source_size
				       - static_cast<double>(target_pos) / target_size);
	
	if (value != 0.0)
	  features[feature_function.feature_name()] = value;
	else
	  features.erase(feature_function.feature_name());
	
	if (! normalizers.empty() && sentence) {
	  normalizer_set_type::const_iterator niter_end = normalizers.end();
	  for (normalizer_set_type::const_iterator niter = normalizers.begin(); niter != niter_end; ++ niter) {
	    const symbol_type& wc = niter->operator()(sentence->operator[](target_pos));
	    
	    if (wc == sentence->operator[](target_pos)) continue;
	    
	    const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					  + ':'
					  + static_cast<const std::string&>(wc));
	    if (value != 0.0)
	      features[feature] = value;
	    else
	      features.erase(value);
	  }
	}
      }
      
      void RelativePosition::operator()(const feature_function_type& feature_function,
					const size_type& id,
					const hypergraph_type& hypergraph,
					const lattice_type& lattice,
					const span_set_type& spans,
					const sentence_set_type& targets,
					const ngram_count_set_type& ngram_counts)
      {
	source_size = 0;
	target_size = 0;
	
	sentence = 0;
	if (! targets.empty())
	  sentence = &targets.front();
      
	if (! hypergraph.is_valid()) return;

	const hypergraph_type::node_type& node = hypergraph.nodes[hypergraph.goal];
	const hypergraph_type::edge_type& edge = hypergraph.edges[node.edges.front()];
      
	attribute_set_type::const_iterator siter = edge.attributes.find(__attr_source_size);
	attribute_set_type::const_iterator titer = edge.attributes.find(__attr_target_size);
      
	if (siter != edge.attributes.end())
	  source_size = boost::apply_visitor(__attribute_integer(), siter->second);
	if (titer != edge.attributes.end())
	  target_size = boost::apply_visitor(__attribute_integer(), titer->second);
	
	source_size = std::max(0, source_size);
	target_size = std::max(0, target_size);
      }

      
      NullJump::NullJump(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (param.name() != "null-jump")
	  throw std::runtime_error("is this really relative position feature function? " + parameter);

	__state_size   = sizeof(int);
	__feature_name = "null-jump";

	feature_none_none = "null-jump:none-none";
	feature_none_word = "null-jump:none-word";
	feature_word_none = "null-jump:word-none";
	feature_word_word = "null-jump:word-word";
      }
      
      // define state-full features...
      void NullJump::operator()(const feature_function_type& feature_function,
				state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
      {
	// how do we define distortion covered by forest...
	
	if (states.empty()) {
	  
	  attribute_set_type::const_iterator titer = edge.attributes.find(__attr_target_position);
	  if (titer == edge.attributes.end())
	    throw std::runtime_error("we do not support non alignment forest");
	  
	  *reinterpret_cast<int*>(state) = (boost::apply_visitor(__attribute_integer(), titer->second) < 0);
	} else if (states.size() == 1)
	  *reinterpret_cast<int*>(state) = *reinterpret_cast<const int*>(states[0]);
	else {
	  // we assume penn-treebank style grammar....
	  
	  int prev = *reinterpret_cast<const int*>(states.front());
	  state_ptr_set_type::const_iterator siter_end = states.end();
	  for (state_ptr_set_type::const_iterator siter = states.begin() + 1; siter != siter_end; ++ siter) {
	    const int next = *reinterpret_cast<const int*>(*siter);
	    
	    if (prev)
	      features[next ? feature_none_none : feature_none_word] += 1.0;
	    else 
	      features[next ? feature_word_none : feature_word_word] += 1.0;
	    
	    prev = next;
	  }
	  *reinterpret_cast<int*>(state) = prev;
	} 
      }

      TargetBigram::TargetBigram(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: normalizers(), sentence(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (param.name() != "target-bigram")
	  throw std::runtime_error("is this really relative position feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "cluster") == 0) {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (strcasecmp(piter->first.c_str(), "stemmer") == 0)
	    normalizers.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else
	    std::cerr << "WARNING: unsupported parameter for target-bigram: " << piter->first << "=" << piter->second << std::endl;
	}
	
	__state_size = sizeof(symbol_type);
	__feature_name = "target-bigram";
      }
      
      
      void TargetBigram::operator()(const feature_function_type& feature_function,
				    state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
      {
	attribute_set_type::const_iterator titer = edge.attributes.find(__attr_target_position);
	
	if (titer == edge.attributes.end()) return;
	
	const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	
	
	
      }
      
      void TargetBigram::operator()(const feature_function_type& feature_function,
				    const size_type& id,
				    const hypergraph_type& hypergraph,
				    const lattice_type& lattice,
				    const span_set_type& spans,
				    const sentence_set_type& targets,
				    const ngram_count_set_type& ngram_counts)
      {
	sentence = 0;
	if (! targets.empty())
	  sentence = &targets.front();
      }

      
    };
  };
};
