//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "feature/alignment.hpp"
#include "parameter.hpp"

#include "utils/piece.hpp"
#include "utils/bithack.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/unordered_map.hpp"
#include "utils/simple_vector.hpp"

#include <boost/fusion/tuple.hpp>
#include <boost/array.hpp>

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
    
      Base::Base()
	: attr_source_size("source-size"),
	  attr_target_size("target-size"),
	  attr_source_position("source-position"),
	  attr_target_position("target-position") {}

      // relative position...
      
      RelativePosition::RelativePosition(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: normalizers(), sentence(0), source_size(0), target_size(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (utils::ipiece(param.name()) != "relative-position")
	  throw std::runtime_error("is this really relative position feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "cluster") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "stemmer")
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
					const bool final) const
      {
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_source_position);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	
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
      
	attribute_set_type::const_iterator siter = edge.attributes.find(attr_source_size);
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_size);
      
	if (siter != edge.attributes.end())
	  source_size = boost::apply_visitor(__attribute_integer(), siter->second);
	if (titer != edge.attributes.end())
	  target_size = boost::apply_visitor(__attribute_integer(), titer->second);
	
	source_size = std::max(0, source_size);
	target_size = std::max(0, target_size);
      }

      class PathImpl : public Base
      {
      public:
	typedef boost::fusion::tuple<symbol_type, symbol_type, symbol_type, symbol_type, int> item_type;
	typedef boost::array<feature_type, 5> features_type;
	
	typedef utils::unordered_map<item_type, features_type, utils::hashmurmur<size_t>, std::equal_to<item_type>,
				     std::allocator<std::pair<const item_type, features_type> > >::type feature_map_type;
	
      public:
	void operator()(const symbol_type& source_prev,
			const symbol_type& target_prev,
			const symbol_type& target_curr,
			const symbol_type& target_next,
			const int distance,
			feature_set_type& features)
	{
	  const item_type item(source_prev, target_prev, target_curr, target_next, distance);
	  
	  feature_map_type::iterator iter = maps.find(item);
	  if (iter == maps.end()) {
	    iter = maps.insert(std::make_pair(item, features_type())).first;

	    const std::string distance_str  = utils::lexical_cast<std::string>(distance);
	    const std::string distance_mark = (distance > 0 ? "P" : (distance < 0 ? "N" : "Z"));
	    
	    iter->second[0] = "path:" + distance_str;

	    iter->second[1] = ("path:"
			       + static_cast<const std::string&>(source_prev)
			       + '+' + static_cast<const std::string&>(target_curr)
			       + ':' + distance_mark);
	    
	    iter->second[2] = ("path:"
			       + static_cast<const std::string&>(source_prev)
			       + '+' + static_cast<const std::string&>(target_prev)
			       + '+' + static_cast<const std::string&>(target_curr)
			       + '+' + static_cast<const std::string&>(target_next)
			       + ':' + distance_mark);
	    
	    iter->second[3] = ("path:"
			       + static_cast<const std::string&>(source_prev)
			       + '+' + static_cast<const std::string&>(target_curr)
			       + ':' + distance_str);
	    
	    iter->second[4] = ("path:"
			       + static_cast<const std::string&>(source_prev)
			       + '+' + static_cast<const std::string&>(target_prev)
			       + '+' + static_cast<const std::string&>(target_curr)
			       + '+' + static_cast<const std::string&>(target_next)
			       + ':' + distance_str);
	    
	    std::sort(iter->second.begin(), iter->second.end());
	  }
	  
	  features_type::const_iterator fiter_end = iter->second.end();
	  for (features_type::const_iterator fiter = iter->second.begin(); fiter != fiter_end; ++ fiter)
	    features[*fiter] += 1.0;
	}
	
	void clear()
	{
	  maps.clear();
	}
	
      public:
        feature_map_type maps;
      };

            
      Path::Path(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: pimpl(new PathImpl()), normalizers_source(), normalizers_target(), sentence(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (utils::ipiece(param.name()) != "path")
	  throw std::runtime_error("is this really path feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "cluster-source") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers_source.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "cluster-target") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers_target.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "stemmer-source")
	    normalizers_source.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else if (utils::ipiece(piter->first) == "stemmer-target")
	    normalizers_target.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else
	    std::cerr << "WARNING: unsupported parameter for path: " << piter->first << "=" << piter->second << std::endl;
	}

	if (normalizers_source.size() > 1)
	  throw std::runtime_error("we do no support multiple normalizers (yet)");
	if (normalizers_target.size() > 1)
	  throw std::runtime_error("we do no support multiple normalizers (yet)");
	
	__state_size   = sizeof(int) + sizeof(symbol_type);
	__feature_name = "path";
      }

      Path::Path(const Path& x)
	: Base(static_cast<const Base&>(x)),
	  pimpl(new PathImpl()),
	  normalizers_source(x.normalizers_source),
	  normalizers_target(x.normalizers_target),
	  sentence(0) {}
      
      Path::~Path() { if (pimpl) delete pimpl; }
      
      Path& Path::operator=(const Path& x)
      {
	static_cast<Base&>(*this) = static_cast<const Base&>(x);
	pimpl->clear();
	normalizers_source = x.normalizers_source;
	normalizers_target = x.normalizers_target;
	return *this;
      }
      
      void Path::operator()(const feature_function_type& feature_function,
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
	
	if (pimpl->maps.size() > 1024 * 1024)
	  pimpl->maps.clear();
      }
	
      // define state-full features...
      void Path::operator()(const feature_function_type& feature_function,
			    state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    const bool final) const
      {
	if (states.empty()) {
	  attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	  attribute_set_type::const_iterator siter = edge.attributes.find(attr_source_position);
	  
	  if (titer == edge.attributes.end())
	    throw std::runtime_error("we do not support non alignment forest");
	  if (siter == edge.attributes.end())
	    throw std::runtime_error("we do not support non alignment forest");

	  const int pos_target = boost::apply_visitor(__attribute_integer(), titer->second);
	  const int pos_source = boost::apply_visitor(__attribute_integer(), siter->second);
	  
	  int*         state_target = reinterpret_cast<int*>(state);
	  symbol_type* state_source = reinterpret_cast<symbol_type*>(state_target + 1);
	  
	  *state_target = pos_target;
	  *state_source = (normalizers_source.empty() ? vocab_type::EPSILON : normalizers_source.front()(edge.rule->rhs.front()));
	  
	  if (pos_source == 0 && pos_target >= 0) {
	    // fire path from BOS!
	    const int distance = pos_target + 1;
	    
	    const symbol_type target_prev = (normalizers_target.empty() || ! sentence || pos_target == 0
					     ? vocab_type::BOS
					     : normalizers_target.front()(sentence->operator[](pos_target - 1)));
	    const symbol_type target_curr = (normalizers_target.empty() || ! sentence
					     ? vocab_type::EPSILON
					     : normalizers_target.front()(sentence->operator[](pos_target)));
	    const symbol_type target_next = (normalizers_target.empty() || ! sentence || pos_target + 1 >= sentence->size()
					     ? vocab_type::EOS
					     : normalizers_target.front()(sentence->operator[](pos_target + 1)));
	    
	    pimpl->operator()(vocab_type::BOS, target_prev, target_curr, target_next, distance, features);
	  }
	} else if (states.size() == 1) {
	  int*         state_target = reinterpret_cast<int*>(state);
	  symbol_type* state_source = reinterpret_cast<symbol_type*>(state_target + 1);
	  
	  const int*         ant_target = reinterpret_cast<const int*>(states[0]);
	  const symbol_type* ant_source = reinterpret_cast<const symbol_type*>(ant_target + 1);
	  
	  *state_target = *ant_target;
	  *state_source = *ant_source;
	} else {
	  int         prev_target = *reinterpret_cast<const int*>(states.front());
	  symbol_type prev_source = *reinterpret_cast<const symbol_type*>(reinterpret_cast<const int*>(states.front()) + 1);
	  state_ptr_set_type::const_iterator siter_end = states.end();
	  for (state_ptr_set_type::const_iterator siter = states.begin() + 1; siter != siter_end; ++ siter) {
	    const int         next_target = *reinterpret_cast<const int*>(*siter);
	    const symbol_type next_source = *reinterpret_cast<const symbol_type*>(reinterpret_cast<const int*>(*siter) + 1);
	    
	    if (next_target >= 0) {
	      // fire feature!
	      
	      const symbol_type target_prev = (normalizers_target.empty() || ! sentence || next_target == 0
					       ? vocab_type::BOS
					       : normalizers_target.front()(sentence->operator[](next_target - 1)));
	      const symbol_type target_curr = (normalizers_target.empty() || ! sentence
					       ? vocab_type::EPSILON
					       : normalizers_target.front()(sentence->operator[](next_target)));
	      const symbol_type target_next = (normalizers_target.empty() || ! sentence || next_target + 1 >= sentence->size()
					       ? vocab_type::EOS
					       : normalizers_target.front()(sentence->operator[](next_target + 1)));

	      
	      if (prev_target < 0) {
		const int distance = next_target + 1; // distance from BOS
		
		pimpl->operator()(vocab_type::BOS, target_prev, target_curr, target_next, distance, features);
	      } else {
		const int distance = next_target - prev_target;
		
		pimpl->operator()(prev_source, target_prev, target_curr, target_next, distance, features);
	      }
	      
	      prev_target = next_target;
	      prev_source = next_source;
	    }
	  }
	  
	  *reinterpret_cast<int*>(state)                                     = prev_target;
	  *reinterpret_cast<symbol_type*>(reinterpret_cast<int*>(state) + 1) = prev_source;
	}
	
	if (final && sentence) {
	  const int         prev_target = *reinterpret_cast<const int*>(states.front());
	  const symbol_type prev_source = *reinterpret_cast<const symbol_type*>(reinterpret_cast<const int*>(states.front()) + 1);

	  // fire feature for EOS...
	  const int distance = sentence->size() - prev_target;
	  
	  const symbol_type target_prev = (prev_target < 0
					   ? vocab_type::BOS
					   : (normalizers_target.empty() || ! sentence
					      ? vocab_type::EPSILON
					      : normalizers_target.front()(sentence->operator[](prev_target))));
	  
	  pimpl->operator()(prev_source, target_prev, vocab_type::EOS, vocab_type::EOS, distance, features);
	}
      }
      
      NullPath::NullPath(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (utils::ipiece(param.name()) != "null-path")
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
				const bool final) const
      {
	// how do we define distortion covered by forest...
	
	if (states.empty()) {
	  attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	  attribute_set_type::const_iterator siter = edge.attributes.find(attr_source_position);
	  
	  if (titer == edge.attributes.end())
	    throw std::runtime_error("we do not support non alignment forest");
	  if (siter == edge.attributes.end())
	    throw std::runtime_error("we do not support non alignment forest");
	  
	  const int pos_target = boost::apply_visitor(__attribute_integer(), titer->second);
	  const int pos_source = boost::apply_visitor(__attribute_integer(), siter->second);
	  
	  *reinterpret_cast<int*>(state) = (pos_target < 0);
	  
	  // fire feature for BOS
	  if (pos_source == 0)
	    features[pos_target < 0 ? feature_word_none : feature_word_word] += 1.0;
	      
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
	
	// fire featuer for EOS
	if (final)
	  features[*reinterpret_cast<int*>(state) ? feature_none_word : feature_word_word] += 1.0;
      }
      
      class FertilityLocalImpl
      {
      public:
	void clear() {}
      };
      
      FertilityLocal::FertilityLocal(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: pimpl(new impl_type()), normalizers(), sentence(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (utils::ipiece(param.name()) != "fertility-local")
	  throw std::runtime_error("is this really local fertility feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "cluster") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "stemmer")
	    normalizers.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else
	    std::cerr << "WARNING: unsupported parameter for fertility-local: " << piter->first << "=" << piter->second << std::endl;
	}
	
	__state_size = sizeof(int) * 2 + sizeof(symbol_type);
	__feature_name = "fertility-local";
      }

      FertilityLocal::FertilityLocal(const FertilityLocal& x)
	: pimpl(new impl_type()),
	  normalizers(x.normalizers),
	  sentence(0) {}
      
      FertilityLocal::~FertilityLocal() { if (pimpl) delete pimpl; }
      
      FertilityLocal& FertilityLocal::operator=(const FertilityLocal& x)
      {
	static_cast<Base&>(*this) = static_cast<const Base&>(x);
	pimpl->clear();
	normalizers = x.normalizers;
	return *this;
      }
      
      void FertilityLocal::operator()(const feature_function_type& feature_function,
				      state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      const bool final) const
      {
	
	
      }
      
      void FertilityLocal::operator()(const feature_function_type& feature_function,
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


      TargetBigram::TargetBigram(const std::string& parameter, size_type& __state_size, feature_type& __feature_name)
	: normalizers(), sentence(0)
      {
	typedef cicada::Parameter parameter_type;
	
	const parameter_type param(parameter);
	
	if (utils::ipiece(param.name()) != "target-bigram")
	  throw std::runtime_error("is this really target bigram feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "cluster") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "stemmer")
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
				    const bool final) const
      {
	symbol_type* context = reinterpret_cast<symbol_type*>(state);
	
	if (! sentence) {
	  context[0] = vocab_type::EPSILON;
	  context[1] = vocab_type::EPSILON;
	  return;
	}
	
	if (states.empty()) {
	  attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
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
	
	if (utils::ipiece(param.name()) != "word-pair")
	  throw std::runtime_error("is this really word pair feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "cluster-source") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers_source.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "cluster-target") {
	    if (! boost::filesystem::exists(piter->second))
	      throw std::runtime_error("no cluster file: " + piter->second);
	    
	    normalizers_target.push_back(normalizer_type(&cicada::Cluster::create(piter->second)));
	  } else if (utils::ipiece(piter->first) == "stemmer-source")
	    normalizers_source.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else if (utils::ipiece(piter->first) == "stemmer-target")
	    normalizers_target.push_back(normalizer_type(&cicada::Stemmer::create(piter->second)));
	  else
	    std::cerr << "WARNING: unsupported parameter for word-pair: " << piter->first << "=" << piter->second << std::endl;
	}

	if (normalizers_source.size() != normalizers_target.size())
	  throw std::runtime_error("# of normalizers do not match");
	
	__state_size = 0;
	__feature_name = "word-pair";
      }
      
      
      void WordPair::operator()(const feature_function_type& feature_function,
				state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				const bool final) const
      {
	if (! states.empty()) return;
	
	attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	
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
	  normalizer_set_type::const_iterator titer_end = normalizers_target.end();
	  normalizer_set_type::const_iterator siter = normalizers_source.begin();
	  normalizer_set_type::const_iterator titer = normalizers_target.begin();
	  
	  for (/**/; siter != siter_end; ++ siter, ++ titer) {
	    const symbol_type source_norm = siter->operator()(*riter);
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
