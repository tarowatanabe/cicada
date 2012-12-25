
#include <iostream>

#include "neuron/neuron.hpp"

#include "utils/random_seed.hpp"

int main(int argc, char** argv)
{
  std::srand(utils::random_seed());

  cicada::neuron::Sequential m1;
  
  // 10 input, 3 output
  m1.push_back(new cicada::neuron::Linear(10, 3));
  m1.push_back(new cicada::neuron::Tanh());

  std::cout << "weight" << std::endl
	    << dynamic_cast<cicada::neuron::Linear&>(*m1.front()).weight << std::endl
	    << "bias" << std::endl
	    << dynamic_cast<cicada::neuron::Linear&>(*m1.front()).bias << std::endl;

  

  const cicada::neuron::Layer::tensor_type input = cicada::neuron::Layer::tensor_type::Random(10, 1);
  
  std::cout << "input" << std::endl
	    << input << std::endl;
  
  m1.forward(input);

  std::cout << "output" << std::endl
	    << m1.data_output << std::endl
	    << "linear" << std::endl
	    << m1.front()->data_output << std::endl;
}
