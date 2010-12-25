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

            
      Path::Path(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (param.name() != "path")
	  throw std::runtime_error("is this really path feature function? " + parameter);
	
	__state_size   = sizeof(int);
	__feature_name = "path";
      }
      
      // define state-full features...
      void Path::operator()(const feature_function_type& feature_function,
			    state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
      {
	
	
      }
      
      NullPath::NullPath(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (param.name() != "null-path")
	  throw std::runtime_error("is this really null-path feature function? " + parameter);

	__state_size   = sizeof(int);
	__feature_name = "null-path";

	feature_none_none = "null-path:none-none";
	feature_none_word = "null-path:none-word";
	feature_word_none = "null-path:word-none";
	feature_word_word = "null-path:word-word";
      }
      
      // define state-full features...
      void NullPath::operator()(const feature_function_type& feature_function,
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
	  throw std::runtime_error("is this really target bigram feature function? " + parameter);
	
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
	
	__state_size = sizeof(symbol_type) * 2;
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
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	
	if (! sentence) {
	  context[0] = vocab_type::EPSILON;
	  context[1] = vocab_type::EPSILON;
	  return;
	}
	
	if (states.empty()) {
	  attribute_set_type::const_iterator titer = edge.attributes.find(__attr_target_position);
	  if (titer == edge.attributes.end())
	    throw std::runtime_error("we do not support non alignment forest");
	  
	  const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	  if (target_pos >= 0) {
	    context[0] = sentence->operator[](target_pos);
	    context[1] = sentence->operator[](target_pos);
	  } else {
	    context[0] = vocab_type::EPSILON;
	    context[1] = vocab_type::EPSILON;
	  }
	} else if (states.size() == 1) {
	  const symbol_type* context_antecedent = reinterpret_cast<const symbol_type*>(states[0]);
	  std::copy(context_antecedent, context_antecedent + 2, context);
	} else {
	  const symbol_type* context_antecedent = reinterpret_cast<const symbol_type*>(states.front());
	  
	  context[0] = context_antecedent[0];
	  symbol_type prev = context_antecedent[1];
	  state_ptr_set_type::const_iterator siter_end = states.end();
	  for (state_ptr_set_type::const_iterator siter = states.begin() + 1; siter != siter_end; ++ siter) {
	    const symbol_type* context_antecedent = reinterpret_cast<const symbol_type*>(states.front());
	    
	    // bigram between prev and context_antecedent[0]
	    
	    if (prev == vocab_type::EPSILON) {
	      if (context_antecedent[0] == vocab_type::EPSILON) {
		// do nothing...
	      } else {
		context[0] = context_antecedent[0];
		prev       = context_antecedent[1];
	      }
	    } else {
	      if (context_antecedent[0] == vocab_type::EPSILON) {
		// do nothing...
	      } else {
		// bigram between prev and context_antecedent[0]

		const symbol_type& next = context_antecedent[1];

		const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					      + ':' + static_cast<const std::string&>(prev)
					      + ':' + static_cast<const std::string&>(next));

		features[feature] += 1.0;
		
		normalizer_set_type::const_iterator niter_end = normalizers.end();
		for (normalizer_set_type::const_iterator niter = normalizers.begin(); niter != niter_end; ++ niter) {
		  const symbol_type prev_norm = niter->operator()(prev);
		  const symbol_type next_norm = niter->operator()(next);

		  if (prev_norm != prev || next_norm != next) {
		    const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
						  + ':' + static_cast<const std::string&>(prev_norm)
						  + ':' + static_cast<const std::string&>(next_norm));
		    
		    features[feature] += 1.0;
		  }
		}
		
		prev = next;
	      }
	    }
	  }
	  
	  context[1] = prev;
	}
	
	if (final) {
	  if (context[0] != vocab_type::EPSILON) {
	    const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					  + ':' + static_cast<const std::string&>(vocab_type::BOS)
					  + ':' + static_cast<const std::string&>(context[0]));
	    
	    features[feature] += 1.0;
	    
	    normalizer_set_type::const_iterator niter_end = normalizers.end();
	    for (normalizer_set_type::const_iterator niter = normalizers.begin(); niter != niter_end; ++ niter) {
	      const symbol_type norm = niter->operator()(context[0]);
	      
	      if (norm == context[0]) continue;
	      
	      const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					    + ':' + static_cast<const std::string&>(vocab_type::BOS)
					    + ':' + static_cast<const std::string&>(norm));
	      
	      features[feature] += 1.0;
	    }
	  }
	  
	  if (context[1] != vocab_type::EPSILON) {
	    const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					  + ':' + static_cast<const std::string&>(context[1])
					  + ':' + static_cast<const std::string&>(vocab_type::EOS));
	    
	    features[feature] += 1.0;
	    
	    normalizer_set_type::const_iterator niter_end = normalizers.end();
	    for (normalizer_set_type::const_iterator niter = normalizers.begin(); niter != niter_end; ++ niter) {
	      const symbol_type norm = niter->operator()(context[1]);
	      
	      if (norm == context[1]) continue;
	      
	      const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					    + ':' + static_cast<const std::string&>(norm)
					    + ':' + static_cast<const std::string&>(vocab_type::EOS));
	      
	      features[feature] += 1.0;
	    }
	  }
	}
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

      WordPair::WordPair(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: normalizers_source(), normalizers_target(), sentence(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (param.name() != "word-pair")
	  throw std::runtime_error("is this really word pair feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "cluster-source") == 0) {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers_source.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (strcasecmp(piter->first.c_str(), "cluster-target") == 0) {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers_target.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (strcasecmp(piter->first.c_str(), "stemmer-source") == 0)
	    normalizers_source.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else if (strcasecmp(piter->first.c_str(), "stemmer-target") == 0)
	    normalizers_target.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else
	    std::cerr << "WARNING: unsupported parameter for word-pair: " << piter->first << "=" << piter->second << std::endl;
	}
	
	__state_size = 0;
	__feature_name = "word-pair";
      }
      
      
      void WordPair::operator()(const feature_function_type& feature_function,
				state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
      {
	if (! states.empty()) return;
	
	attribute_set_type::const_iterator titer = edge.attributes.find(__attr_target_position);
	
	if (titer == edge.attributes.end())
	  throw std::runtime_error("we do not support non alignment forest");
	
	const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	
	if (target_pos < -1)
	  throw std::runtime_error("we do not support non alignment forest");
	
	const symbol_type& target = (target_pos >= 0 && sentence ? sentence->operator[](target_pos) : vocab_type::EPSILON);
	
	rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
	  
	  const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					+ ':' + static_cast<const std::string&>(*riter)
					+ ':' + static_cast<const std::string&>(target));

	  features[feature] += 1.0;
	  
	  normalizer_set_type::const_iterator siter_end = normalizers_source.end();
	  for (normalizer_set_type::const_iterator siter = normalizers_source.begin(); siter != siter_end; ++ siter) {
	    const symbol_type source_norm = siter->operator()(*riter);
	    
	    normalizer_set_type::const_iterator titer_end = normalizers_target.end();
	    for (normalizer_set_type::const_iterator titer = normalizers_target.begin(); titer != titer_end; ++ titer) {
	      const symbol_type target_norm = titer->operator()(target);
	      
	      if (*riter != source_norm || target != target_norm) {
		const feature_type feature = (static_cast<const std::string&>(feature_function.feature_name())
					      + ':' + static_cast<const std::string&>(source_norm)
					      + ':' + static_cast<const std::string&>(target_norm));
		
		features[feature] += 1.0;
	      }
	    }
	  }
	}
      }
      
      void WordPair::operator()(const feature_function_type& feature_function,
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
