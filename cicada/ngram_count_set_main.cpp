//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "ngram_count_set.hpp"

#include <cicada/msgpack/ngram_count_set.hpp>

#include "msgpack_main_impl.hpp"

int main(int argc, char** argv)
{
  const char tmp[] = "{[\"morning\"]:5.0, [\"bad\", \"morning\"]:3.0, [\"good\", \"morning\"]:2.0}";

  cicada::NGramCountSet counts;
  counts.assign(tmp);
  std::cout << "count: " << counts << std::endl;
  
  msgpack_test(counts);

  cicada::NGramCountSet input;
  while (std::cin >> input) {
    std::cout << "input: " << input << std::endl;
        
    msgpack_test(input);
  }
}

