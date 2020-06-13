//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <stdexcept>
#include <set>
#include <memory>

#include "model.hpp"

#include "utils/simple_vector.hpp"
#include "utils/bithack.hpp"

namespace cicada
{
  class StateAllocator : public std::allocator<char>
  {
  public:
    typedef Model::state_type   state_type;

    typedef state_type::pointer pointer;
    typedef size_t              size_type;
    typedef ptrdiff_t           difference_type;
  
    typedef std::allocator<char> allocator_type;
    typedef utils::simple_vector<pointer, std::allocator<pointer> > state_set_type;
  
  public:
    static const size_type chunk_size  = 1024 * 16;
    static const size_type chunk_mask  = chunk_size - 1;
  
  public:

    StateAllocator()
      : states(), state_iterator(0),
	cache(0),
	state_size(0),
	state_alloc_size(0),
	state_chunk_size(0) {}
  
    StateAllocator(size_type __state_size)
      : states(), state_iterator(0),
	cache(0),
	state_size(__state_size),
	state_alloc_size(0),
	state_chunk_size(0)
    {
      if (state_size != 0) {
	// sizeof(char*) aligned size..
	const size_type alignment_size = utils::bithack::max(sizeof(pointer), size_type(16));
	const size_type alignment_mask = ~(alignment_size - 1);
	
	state_alloc_size = (state_size + alignment_size - 1) & alignment_mask;
	state_chunk_size = state_alloc_size * chunk_size;
      }
    }
  
    StateAllocator(const StateAllocator& x) 
      : allocator_type(static_cast<const allocator_type&>(x)),
	states(), state_iterator(0),
	cache(0),
	state_size(x.state_size),
	state_alloc_size(x.state_alloc_size),
	state_chunk_size(x.state_chunk_size) {}
  
    StateAllocator& operator=(const StateAllocator& x)
    { 
      clear();
      
      static_cast<allocator_type&>(*this) = static_cast<const allocator_type&>(x);
      
      state_size = x.state_size;
      state_alloc_size = x.state_alloc_size;
      state_chunk_size = x.state_chunk_size;

      return *this;
    }
  
    ~StateAllocator() { clear(); }
  
  public:
  
    state_type allocate()
    {
      if (state_size == 0) return 0;
    
      if (cache) {
	pointer state = cache;
	cache = *reinterpret_cast<pointer*>(state);
	
	// clear buffer...
	std::fill(state, state + state_alloc_size, 0);
	
	return state_type(state);
      }
      
      const size_type chunk_pos = state_iterator & chunk_mask;
    
      if (chunk_pos == 0) {
	states.push_back(allocator_type::allocate(state_chunk_size));
	std::uninitialized_fill(states.back(), states.back() + state_chunk_size, 0);

	//std::cerr << "allocated: " << ((void*) states.back()) << " size: " << states.size() << std::endl;
      }
    
      ++ state_iterator;
    
      return state_type(states.back() + chunk_pos * state_alloc_size);
    }
  
    void deallocate(const state_type& state)
    {
      if (state.empty() || state_size == 0 || states.empty()) return;
      
      *reinterpret_cast<pointer*>(const_cast<state_type&>(state).base) = cache;
      cache = state.base;
    }

    state_type clone(const state_type& state)
    {
      if (state_size == 0)
	return state_type();
      
      state_type state_new = allocate();
      
      std::copy(state.base, state.base + state_alloc_size, state_new.base);
      
      return state_new;
    }
  
    void clear()
    {
      //std::cerr << "allocated states size: " << states.size() << std::endl;

      state_set_type::iterator siter_end = states.end();
      for (state_set_type::iterator siter = states.begin(); siter != siter_end; ++ siter) {
	//std::cerr << "de-allocated: " << ((void*) *siter) << std::endl;
	allocator_type::deallocate(*siter, state_chunk_size);
      }
    
      states.clear();
      state_iterator = 0;
      cache = 0;
    }
  
  private:  
    state_set_type states;
    size_type state_iterator;
    pointer   cache;

    size_type state_size;
    size_type state_alloc_size;
    size_type state_chunk_size;
  };

  
  Model::Model()
    : models(),
      allocator(new state_allocator_type()),
      offsets(),
      states_size(0) {}

  Model::Model(const feature_function_ptr_type& x)
    : models(1, x),
      allocator(new state_allocator_type()),
      offsets(),
      states_size(0) {}
  
  Model::Model(const Model& x)
    : models(x.models),
      allocator(new state_allocator_type(*x.allocator)),
      offsets(x.offsets),
      states_size(x.states_size) {}
  
  Model::~Model() { std::unique_ptr<state_allocator_type> tmp(allocator); }
  
  Model& Model::operator=(const Model& x)
  {
    models      = x.models;
    *allocator  = *x.allocator;
    offsets     = x.offsets;
    states_size = x.states_size;
    
    return *this;
  }

  
  Model::state_type Model::apply(const state_set_type& node_states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const
  {
    features.rehash(features.size() + models.size());

    state_type state = allocator->allocate();
    
    feature_function_type::state_ptr_set_type states(edge.tails.size());

    //std::cerr << "apply features for: " << *(edge.rule) << std::endl;

    for (size_t i = 0; i != models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];
      
      if (feature_function.state_size())
	for (size_t k = 0; k != states.size(); ++ k)
	  states[k] = node_states[edge.tails[k]].base + offsets[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply(state_feature, states, edge, features, final);
    }
    
    //std::cerr << "apply features end" << std::endl;
    
    return state;
  }

  Model::state_type Model::apply_coarse(const state_set_type& node_states,
					const edge_type& edge,
					feature_set_type& features,
					const bool final) const
  {
    features.rehash(features.size() + models.size());
    
    state_type state = allocator->allocate();
    
    feature_function_type::state_ptr_set_type states(edge.tails.size());

    //std::cerr << "apply features for: " << *(edge.rule) << std::endl;
    
    for (size_t i = 0; i != models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];
      
      if (feature_function.state_size())
	for (size_t k = 0; k != states.size(); ++ k)
	  states[k] = node_states[edge.tails[k]].base + offsets[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply_coarse(state_feature, states, edge, features, final);
    }
    
    //std::cerr << "apply features end" << std::endl;
    
    return state;
  }

  void Model::apply_predict(state_type& state,
			    const state_set_type& node_states,
			    const edge_type& edge,
			    feature_set_type& features,
			    const bool final) const
  {
    features.rehash(features.size() + models.size());
    
    if (state.empty())
      state = allocator->allocate();

    feature_function_type::state_ptr_set_type states(edge.tails.size());
    
    for (size_t i = 0; i != models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply_predict(state_feature, states, edge, features, final);
    }
  }
  
  void Model::apply_scan(state_type& state,
			 const state_set_type& node_states,
			 const edge_type& edge,
			 const int dot,
			 feature_set_type& features,
			 const bool final) const
  {
    features.rehash(features.size() + models.size());

    feature_function_type::state_ptr_set_type states(edge.tails.size());
    
    for (size_t i = 0; i != models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply_scan(state_feature, states, edge, dot, features, final);
    }
  }

  void Model::apply_complete(state_type& state,
			     const state_set_type& node_states,
			     const edge_type& edge,
			     feature_set_type& features,
			     const bool final) const
  {
    features.rehash(features.size() + models.size());

    feature_function_type::state_ptr_set_type states(edge.tails.size());
    
    for (size_t i = 0; i != models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];

      if (feature_function.state_size())
	for (size_t k = 0; k != states.size(); ++ k)
	  states[k] = node_states[edge.tails[k]].base + offsets[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply_complete(state_feature, states, edge, features, final);
    }
  }
  
  void Model::deallocate(const state_type& state) const
  {
    const_cast<state_allocator_type&>(*allocator).deallocate(state);
  }

  Model::state_type Model::clone(const state_type& state) const
  {
    return const_cast<state_allocator_type&>(*allocator).clone(state);
  }
  
  
  void Model::initialize()
  {
    typedef std::set<feature_type, std::less<feature_type>, std::allocator<feature_type> > feature_unique_type;
    
    feature_unique_type feature_names;
    
    offsets.clear();
    offsets.reserve(models.size());
    offsets.resize(models.size());
    states_size = 0;

    const size_t alignment_size = utils::bithack::max(sizeof(void*), size_type(16));
    const size_t alignment_mask = ~(alignment_size - 1);
    
    for (size_t i = 0; i != models.size(); ++ i) {
      offsets[i] = states_size;
      
      //std::cerr << "offset: " << i << " = " << offsets[i] << std::endl;
      
      // multiple of pointer size...
      if (models[i]->state_size())
	states_size += (models[i]->state_size() + alignment_size - 1) & alignment_mask;
      
      if (feature_names.find(models[i]->feature_name()) != feature_names.end())
	throw std::runtime_error("you have already registered feature: " + static_cast<const std::string&>(models[i]->feature_name()));
      feature_names.insert(models[i]->feature_name());
      
      models[i]->initialize();
    }
    
    *allocator = state_allocator_type(states_size);
  }
  
};
