//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "bleu_multi.hpp"

#include "cicada/parameter.hpp"

namespace cicada
{
  namespace feature
  {
    BleuMulti::BleuMulti(const std::string& parameter)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "bleu-multi" && param.name() != "bleu-multiple")
	throw std::runtime_error("this is not bleu-multi feature: " + parameter);
      
      std::string order;
      std::string exact;
      std::string tokenizer;
      
      int size = 0;
      
      std::vector<std::string, std::allocator<std::string> > names;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "order") == 0)
	  order = piter->second;
	else if (strcasecmp(piter->first.c_str(), "exact") == 0)
	  exact = piter->second;
	else if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	  tokenizer = piter->second;
	else if (strcasecmp(piter->first.c_str(), "name") == 0)
	  names.push_back(piter->second);
	else if (strcasecmp(piter->first.c_str(), "size") == 0)
	  size = boost::lexical_cast<int>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for bleu-multi: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (size <= 0)
	throw std::runtime_error("invalid bleu-multi size: " + boost::lexical_cast<std::string>(size));
      
      if (static_cast<int>(names.size()) > size)
	throw std::runtime_error("feature-name size and bleu size do not match");
      
      if (static_cast<int>(names.size()) < size)
	names.resize(size);
      
      base_type::__state_size = 0;
      for (size_t i = 0; i != names.size(); ++ i) {
	if (names[i].empty())
	  names[i] = "bleu-multi:" + boost::lexical_cast<std::string>(i);
	
	parameter_type param_bleu(param);
	param_bleu.name() = "bleu";
	param_bleu.erase("name");
	param_bleu.erase("size");
	
	param_bleu.push_back(std::make_pair("name", names[i]));
	
	bleus.push_back(feature_function_type::create(boost::lexical_cast<std::string>(param_bleu)));
	
	base_type::__state_size += bleus.back()->state_size();
      }
      
      
      base_type::__feature_name = "bleu-multi";
      
    }
    
    BleuMulti::BleuMulti(const BleuMulti& x)
      : base_type(static_cast<const base_type&>(x)),
	bleus(x.bleus)
    {
      for (size_t i = 0; i != bleus.size(); ++ i)
	bleus[i] = bleus[i]->clone();
    }
    
    BleuMulti& BleuMulti::operator=(const BleuMulti& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);

      bleus = x.bleus;
      
      for (size_t i = 0; i != bleus.size(); ++ i)
	bleus[i] = bleus[i]->clone();
      
      return *this;
    }
    
    void BleuMulti::apply(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {
      state_ptr_type     state_bleu;
      state_ptr_set_type states_bleu(states.size());
      
      size_t size = 0;
      for (size_t i = 0; i != bleus.size(); ++ i) {
	state_bleu = static_cast<char*>(state) + size;
	for (size_t k = 0; k != states.size(); ++ k)
	  states_bleu[k] = static_cast<char*>(states[k]) + size;
	
	bleus[i]->apply(state_bleu, states_bleu, edge, features, estimates, final);
	size += bleus[i]->state_size();
      }
    }
    void BleuMulti::apply_coarse(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {

    }
    void BleuMulti::apply_predict(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
    }
    void BleuMulti::apply_scan(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       const int dot,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {
    }
    void BleuMulti::apply_complete(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {
      state_ptr_type     state_bleu;
      state_ptr_set_type states_bleu(states.size());
      
      size_t size = 0;
      for (size_t i = 0; i != bleus.size(); ++ i) {
	state_bleu = static_cast<char*>(state) + size;
	for (size_t k = 0; k != states.size(); ++ k)
	  states_bleu[k] = static_cast<char*>(states[k]) + size;
	
	bleus[i]->apply_complete(state_bleu, states_bleu, edge, features, estimates, final);
	size += bleus[i]->state_size();
      }
    }
    
    void BleuMulti::assign(const size_type& id,
			   const hypergraph_type& hypergraph,
			   const lattice_type& lattice,
			   const span_set_type& spans,
			   const sentence_set_type& targets,
			   const ngram_count_set_type& ngram_counts)
    {
      if (targets.size() != bleus.size())
	throw std::runtime_error("target size do not match with bleu-multi size");

      for (size_t i = 0; i != bleus.size(); ++ i)
	bleus[i]->assign(id, hypergraph, lattice, spans, sentence_set_type(1, targets[i]), ngram_counts);
    }
    
    void BleuMulti::initialize()
    {
      for (size_t i = 0; i != bleus.size(); ++ i)
	bleus[i]->initialize();
    }
  };
};
