
#include <cstdlib>

#include <iostream>

#include <utils/logprob.hpp>
#include <utils/mathop.hpp>

#include <boost/array.hpp>

int main(int argc, char** argv)
{
  
  srandom(time(0) * getpid());
 
  for (int i = 0; i < 4; ++ i) {

    double mul = 1.0;
    utils::logprob logmul = 1.0;
    utils::logprob logmul_from = 1.0;
    for (int samples = 0; samples < 1024; ++ samples) {
      const double value = double(random()) / random();
      utils::logprob logprob = value;
      
      logmul *= logprob;
      mul *= value;
      logmul_from *= utils::logprob::from_log(std::log(value));
    }

    
    
    std::cerr << "mul: " << mul << " logmul: " << logmul << " logmul-from: " << logmul_from << " log: " << utils::mathop::log(mul) << std::endl;

    boost::array<utils::logprob, 16> logprobs;
    for (int i = 0; i < 16; ++ i)
      if (logprobs[i] != utils::logprob())
	std::cerr << "differ..." << std::endl;

    boost::array<utils::logprob, 16>* logprobs_ptr = new boost::array<utils::logprob, 16>();
    for (int i = 0; i < 16; ++ i)
      if (logprobs_ptr->operator[](i) != utils::logprob())
	std::cerr << "differ..." << std::endl;
    
    delete logprobs_ptr;
  }
}
