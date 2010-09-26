#include <stdexcept>
#include <set>
#include <memory>

#include "model.hpp"

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
    typedef std::vector<pointer, std::allocator<pointer> > state_set_type;
  
  public:
    static const size_type chunk_size  = 1024 * 16;
    static const size_type chunk_mask  = chunk_size - 1;
  
  public:
  
    StateAllocator(size_type __state_size = 0)
      : states(), state_iterator(0),
	cache(0),
	state_size(__state_size),
	state_alloc_size(0),
	state_chunk_size(0)
    {
      if (state_size != 0) {
	// sizeof(char*) aligned size..
	const size_type pointer_size = sizeof(pointer);
	const size_type pointer_mask = ~(pointer_size - 1);
	
	state_alloc_size = (state_size + pointer_size - 1) & pointer_mask;
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
      static_cast<allocator_type&>(*this) = static_cast<const allocator_type&>(x);
    
      clear();
    
      state_size = x.state_size;
      state_alloc_size = x.state_alloc_size;
      state_chunk_size = x.state_chunk_size;
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
      }
    
      ++ state_iterator;
    
      return state_type(states.back() + chunk_pos * state_alloc_size);
    }
  
    void deallocate(const state_type& state)
    {
      if (state_size == 0 || states.empty()) return;
      
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
      state_set_type::iterator siter_end = states.end();
      for (state_set_type::iterator siter = states.begin(); siter != siter_end; ++ siter)
	allocator_type::deallocate(*siter, state_chunk_size);
    
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
  
  Model::Model(const Model& x)
    : models(x.models),
      allocator(new state_allocator_type(*x.allocator)),
      offsets(x.offsets),
      states_size(x.states_size) {}
  
  Model::~Model() { std::auto_ptr<state_allocator_type> tmp(allocator); }
  
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
				 feature_set_type& estimates,
				 const bool final) const
  {
    state_type state = allocator->allocate();
    
    feature_function_type::state_ptr_set_type states(edge.tails.size());

    //std::cerr << "apply features for: " << *(edge.rule) << std::endl;
    
    for (int i = 0; i < models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];
      
      if (feature_function.state_size())
	for (int k = 0; k < states.size(); ++ k)
	  states[k] = node_states[edge.tails[k]].base + offsets[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply(state_feature, states, edge, features, estimates, final);
    }
    
    //std::cerr << "apply features end" << std::endl;
    
    return state;
  }

  Model::state_type Model::apply_coarse(const state_set_type& node_states,
					const edge_type& edge,
					feature_set_type& features,
					feature_set_type& estimates,
					const bool final) const
  {
    state_type state = allocator->allocate();
    
    feature_function_type::state_ptr_set_type states(edge.tails.size());

    //std::cerr << "apply features for: " << *(edge.rule) << std::endl;
    
    for (int i = 0; i < models.size(); ++ i) {
      const feature_function_type& feature_function = *models[i];
      
      if (feature_function.state_size())
	for (int k = 0; k < states.size(); ++ k)
	  states[k] = node_states[edge.tails[k]].base + offsets[i];
      
      feature_function_type::state_ptr_type state_feature = state.base + offsets[i];
      
      feature_function.apply_coarse(state_feature, states, edge, features, estimates, final);
    }
    
    //std::cerr << "apply features end" << std::endl;
    
    return state;
  }

  void Model::apply_predict(state_type& state,
			    const state_set_type& node_states,
			    const edge_type& edge,
			    const int dot,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final=false) const
  {
    
    
  }
  
  void Model::apply_scan(state_type& state,
			 const state_set_type& node_states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final=false) const
  {
    
    
  }

  void Model::apply_complete(state_type& state,
			     const state_set_type& node_states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final=false) const
  {
    
    
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
    
    for (int i = 0; i < models.size(); ++ i) {
      offsets[i] = states_size;
      
      //std::cerr << "offset: " << i << " = " << offsets[i] << std::endl;
      
      states_size += models[i]->state_size();
      
      if (feature_names.find(models[i]->feature_name()) != feature_names.end())
	throw std::runtime_error("you have already registered feature: " + static_cast<const std::string&>(models[i]->feature_name()));
      feature_names.insert(models[i]->feature_name());
      
      models[i]->initialize();
    }
    
    *allocator = state_allocator_type(states_size);
  }
  
};
