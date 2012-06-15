#include <iostream>
#include <string>
#include <set>

#include <utils/stree_set.hpp>

#include <boost/random.hpp>

typedef std::set<int> std_set_type;
typedef utils::stree_set<int>::type vec_set_type;

void verify(const std_set_type& stdset, const vec_set_type& vecset)
{
   if (stdset.size() != vecset.size() || ! std::equal(stdset.begin(), stdset.end(), vecset.begin()))
     std::cerr << "different stdset and vecset" << std::endl;

   for (std_set_type::const_iterator iter = stdset.begin(); iter != stdset.end(); ++ iter) {
     vec_set_type::const_iterator viter = vecset.find(*iter);
     if (viter == vecset.end())
       std::cerr << "no entry: " << *iter << std::endl;
     else if (*viter != *iter)
       std::cerr << "different value: " << *iter << " " << *viter << std::endl;
   }
     
     
}


int main(int argc, char** argv)
{
  boost::mt19937 generator;
  generator.seed(time(0) * getpid());
  boost::random_number_generator<boost::mt19937> gen(generator);
    
  std_set_type stdset;
  vec_set_type vecset;
  
  
  for (int i = 0; i != 1024 * 16; ++ i) {
    const int value = gen(1024 * 1024 * 1024);
    stdset.insert(value);
    vecset.insert(value);
  }
  
  verify(stdset, vecset);
  
  for (int i = 0; i != 1024 * 16; ++ i) {
    const int value = gen(1024 * 1024 * 1024);
    stdset.insert(value);
    vecset.insert(value);
  }

  verify(stdset, vecset);

  vecset.clear();
  vecset.insert(stdset.begin(), stdset.end());
  
  verify(stdset, vecset);
}
