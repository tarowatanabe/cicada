
#include <iostream>
#include <sstream>

#include "neuron/neuron.hpp"

#include "utils/random_seed.hpp"

template <typename L>
void test_dimenstion(const char* name)
{
  L sum1(false);
  L sum2(true);
  
  const cicada::neuron::Layer::tensor_type input = cicada::neuron::Layer::tensor_type::Random(3, 7);
  
  sum1.forward(input);
  sum2.forward(input);
  
  std::cout << name << " input:" << std::endl
	    << input << std::endl
	    << name << " forward1:" << std::endl
	    << sum1.data_output << std::endl
	    << name << " forward2:" << std::endl
	    << sum2.data_output << std::endl;
  
  const cicada::neuron::Layer::tensor_type grad1 = cicada::neuron::Layer::tensor_type::Random(7, 1);
  const cicada::neuron::Layer::tensor_type grad2 = cicada::neuron::Layer::tensor_type::Random(3, 1);
  
  sum1.backward(input, grad1);
  sum2.backward(input, grad2);
  
  std::cout << name << " grad1:"<< std::endl
	    << grad1 << std::endl
	    << name << " backward1:" << std::endl
	    << sum1.gradient_input << std::endl
	    << name << " grad2:"<< std::endl
	    << grad2 << std::endl
	    << name << " backward2:" << std::endl
	    << sum2.gradient_input << std::endl;
}

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

  std::cout << m1 << std::endl
	    << m1.front() << std::endl
	    << m1.back() << std::endl;
  
  std::ostringstream sstr;
  sstr << m1;

  cicada::neuron::Layer::layer_ptr_type layers(cicada::neuron::Layer::construct(sstr.str()));

  std::cout << "re-generated: " << layers << std::endl;

  const cicada::neuron::Layer::tensor_type input = cicada::neuron::Layer::tensor_type::Random(10, 1);
  
  std::cout << "input" << std::endl
	    << input << std::endl;
  
  m1.forward(input);

  std::cout << "output" << std::endl
	    << m1.data_output << std::endl
	    << "linear" << std::endl
	    << m1.front()->data_output << std::endl;
  
  // summation

  test_dimenstion<cicada::neuron::Sum>("sum");
  test_dimenstion<cicada::neuron::Mean>("mean");
  test_dimenstion<cicada::neuron::Max>("max");
  test_dimenstion<cicada::neuron::Min>("min");
}
