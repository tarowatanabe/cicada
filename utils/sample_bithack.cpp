//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "utils/bithack.hpp"

#include <boost/numeric/conversion/bounds.hpp>

int main(int argc, char** argv)
{
  std::cout << utils::bithack::is_power2(0) << ' ' << utils::bithack::next_largest_power2(0) << std::endl;
  std::cout << utils::bithack::is_power2(1) << ' ' << utils::bithack::next_largest_power2(1) << std::endl;
  std::cout << utils::bithack::is_power2(2) << ' ' << utils::bithack::next_largest_power2(2) << std::endl;
  std::cout << utils::bithack::is_power2(3) << ' ' << utils::bithack::next_largest_power2(3) << std::endl;
  std::cout << utils::bithack::is_power2(4) << ' ' << utils::bithack::next_largest_power2(4) << std::endl;

  std::cout << utils::bithack::min(-56, 72) << std::endl;
  std::cout << utils::bithack::max(-56, 72) << std::endl;
  
  std::cout << utils::bithack::min(size_t(-56), size_t(72)) << std::endl;
  std::cout << utils::bithack::max(size_t(-56), size_t(72)) << std::endl;

  std::cout << utils::bithack::min(56, 72) << std::endl;
  std::cout << utils::bithack::max(56, 72) << std::endl;

  std::cout << utils::bithack::min(56, boost::numeric::bounds<int>::highest()) << std::endl;
  std::cout << utils::bithack::max(56, boost::numeric::bounds<int>::highest()) << std::endl;

  std::cout << utils::bithack::min(56, boost::numeric::bounds<int>::lowest()) << std::endl;
  std::cout << utils::bithack::max(56, boost::numeric::bounds<int>::lowest()) << std::endl;

  std::cout << int(utils::bithack::abs((int8_t) 56)) << std::endl;
  std::cout << int(utils::bithack::abs((int8_t) -56)) << std::endl;

  std::cout << utils::bithack::abs((int16_t) 56) << std::endl;
  std::cout << utils::bithack::abs((int16_t) - 56) << std::endl;

  std::cout << utils::bithack::abs((int32_t) 56) << std::endl;
  std::cout << utils::bithack::abs((int32_t) - 56) << std::endl;

  std::cout << utils::bithack::abs((int64_t) 56) << std::endl;
  std::cout << utils::bithack::abs((int64_t) - 56) << std::endl;
  

  std::cout << utils::bithack::bit_count((int8_t) 86) << std::endl;
  std::cout << utils::bithack::bit_count((int16_t) 86) << std::endl;
  std::cout << utils::bithack::bit_count((int32_t) 86) << std::endl;
  std::cout << utils::bithack::bit_count((int64_t) 86) << std::endl;
  
  std::cout << utils::bithack::static_most_significant_bit<86>::result << std::endl;
  std::cout << utils::bithack::most_significant_bit((int8_t) 86) << std::endl;
  
  std::cout << "log_2 8  " << utils::bithack::floor_log2(8) << " " << utils::bithack::static_floor_log2<8>::result << std::endl;
  std::cout << "log_2 16 " << utils::bithack::floor_log2(16) << std::endl;
  std::cout << "log_2 32 " << utils::bithack::floor_log2(32) << std::endl;
  std::cout << "log_2 64 " << utils::bithack::floor_log2(64) << std::endl;

  std::cerr << "5 < 64 ? 0 : 5 = " << utils::bithack::branch(5 < 64, 0, 5) << std::endl;
  std::cerr << "5 < 3 ? 0 : 5 = " << utils::bithack::branch(5 < 3, 0, 5) << std::endl;
}
