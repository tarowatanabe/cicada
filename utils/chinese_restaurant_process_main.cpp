#include <iostream>
#include <vector>
#include <string>

#include "chinese_restaurant_process.hpp"
#include "sampler.hpp"

typedef utils::chinese_restaurant_process<std::string> crp_type;
typedef utils::sampler<boost::mt19937> sampler_type;

std::ostream& operator<<(std::ostream& os, const crp_type& crp)
{
  os << "PYP(discount=" << crp.discount() << ",strength=" << crp.strength() << ")" << '\n'
     << "customers = " << crp.size_customer() << '\n';
  
  crp_type::const_iterator citer_end = crp.end();
  for (crp_type::const_iterator citer = crp.begin(); citer != citer_end; ++ citer) {
    os << '\t' << citer->first << " customer: " << citer->second.size_customer() << " tables: " << citer->second.size_table() << '\n';
    
    crp_type::mapped_type::const_iterator iter_end = citer->second.end();
    for (crp_type::mapped_type::const_iterator iter = citer->second.begin(); iter != iter_end; ++ iter) {
      
    }
  }
  
  return os;
}


int main(int argc, char** argv)
{
  sampler_type sampler;
  
  {
    crp_type crp(0.1, 5);
    
    double base = 0.25;
    int total = 0;
    total += crp.increment("hoge", base, sampler);
    total += crp.increment("foo", base, sampler);
    total += crp.increment("foo", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
    total += crp.increment("bar", base, sampler);
     
    std::cout << "total: " << total << std::endl;
    std::cout << "  P(hoge)=" << crp.prob("hoge", base) << std::endl;
    std::cout << "  P(foo)="  << crp.prob("foo", base) << std::endl;
    std::cout << "  P(bar)="  << crp.prob("bar", base) << std::endl;
    std::cout << "  P(boom)=" << crp.prob("boom", base) << std::endl;
    std::cout << "total=" << (crp.prob("hoge", base) + crp.prob("foo", base) + crp.prob("bar", base) + crp.prob("boom", base)) << std::endl;
    std::cout << "log-likelihood=" << crp.log_likelihood() << std::endl;
    std::cout << crp << std::flush;
    
    total -= crp.decrement("hoge", sampler);
    total -= crp.decrement("bar", sampler);
    std::cout << crp << std::flush;
    total -= crp.decrement("bar", sampler);
    std::cout << crp << std::flush;
  }
  
  {
    typedef utils::chinese_restaurant_process<int> crp_type;
    
    crp_type crp(0.5, 1);
    
    double tot = 0;
    double xt = 0;
    int cust = 10;
    std::vector<int> hist(cust + 1, 0);
    
    for (int i = 0; i < cust; ++ i)
      crp.increment(1, 1.0, sampler);
    
    const int samples = 100000;
    const bool simulate = true;
    for (int k = 0; k < samples; ++ k) {
      if (! simulate) {
        crp.clear();
        for (int i = 0; i < cust; ++ i)
	  crp.increment(1, 1.0, sampler);
	
      } else {
        const int da = sampler() * cust;
        const bool a = sampler() < 0.5;
	
        if (a) {
          for (int i = 0; i < da; ++ i)	
	    crp.increment(1, 1.0, sampler);
          for (int i = 0; i < da; ++ i)
	    crp.decrement(1, sampler);
          xt += 1.0;
        } else {
          for (int i = 0; i < da; ++ i)
	    crp.decrement(1, sampler);
	  
          for (int i = 0; i < da; ++ i)
	    crp.increment(1, 1.0, sampler);
        }
      }
      int c = crp.size_table(1);
      ++ hist[c];
      tot += c;
    }
    std::cout << "P(a) = " << (xt / samples) << std::endl;
    std::cout << "E[num tables] = " << (tot / samples) << std::endl;
    double error = std::fabs((tot / samples) - 5.4);
    std::cout << "error = " << error << std::endl;
    for (int i = 1; i <= cust; ++ i)
      std::cout << i << ' ' << hist[i] << std::endl;
  }
  
  
  {
    crp_type crp(0.1,50.0,1,1,1,1);
    crp.increment("foo", 1.0, sampler);
    std::cout << crp.log_likelihood() << std::endl;
  
  }
}
