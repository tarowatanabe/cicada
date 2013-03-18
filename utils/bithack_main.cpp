//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <iostream>

#include "utils/bithack.hpp"

#include <boost/numeric/conversion/bounds.hpp>

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE bithack_test

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(power2)
{
  BOOST_CHECK(utils::bithack::is_power2(0));
  BOOST_CHECK(utils::bithack::is_power2(1));
  BOOST_CHECK(utils::bithack::is_power2(2));
  BOOST_CHECK(! utils::bithack::is_power2(3));
  BOOST_CHECK(utils::bithack::is_power2(4));
  BOOST_CHECK(! utils::bithack::is_power2(5));
  BOOST_CHECK(! utils::bithack::is_power2(6));
  BOOST_CHECK(! utils::bithack::is_power2(7));
  BOOST_CHECK(utils::bithack::is_power2(8));

  BOOST_CHECK_EQUAL( 1, utils::bithack::next_largest_power2(0));
  BOOST_CHECK_EQUAL( 2, utils::bithack::next_largest_power2(1));
  BOOST_CHECK_EQUAL( 4, utils::bithack::next_largest_power2(2));
  BOOST_CHECK_EQUAL( 4, utils::bithack::next_largest_power2(3));
  BOOST_CHECK_EQUAL( 8, utils::bithack::next_largest_power2(4));
  BOOST_CHECK_EQUAL( 8, utils::bithack::next_largest_power2(5));
  BOOST_CHECK_EQUAL( 8, utils::bithack::next_largest_power2(6));
  BOOST_CHECK_EQUAL( 8, utils::bithack::next_largest_power2(7));
  BOOST_CHECK_EQUAL(16, utils::bithack::next_largest_power2(8));
}

BOOST_AUTO_TEST_CASE(min_max)
{
  //
  // revise this implementation... so that we can make a meaningful comparison
  // 
  BOOST_CHECK_EQUAL(-56, utils::bithack::min(-56, 72));
  BOOST_CHECK_EQUAL( 72, utils::bithack::max(-56, 72));
  BOOST_CHECK_EQUAL(size_t(72), utils::bithack::min(size_t(-56), size_t(72)));
  BOOST_CHECK_EQUAL(size_t(-56), utils::bithack::max(size_t(-56), size_t(72)));
  BOOST_CHECK_EQUAL( 56, utils::bithack::min(56, 72));
  BOOST_CHECK_EQUAL( 72, utils::bithack::max(56, 72));

  BOOST_CHECK_EQUAL(56, utils::bithack::min(56, boost::numeric::bounds<int>::highest()));
  BOOST_CHECK_EQUAL(boost::numeric::bounds<int>::highest(), utils::bithack::max(56, boost::numeric::bounds<int>::highest()));
}


#if 0
int main(int argc, char** argv)
{


  std::cout << "min/max 56 and lowest" << std::endl;

  std::cout << utils::bithack::min(56, boost::numeric::bounds<int>::lowest()) << std::endl;
  std::cout << utils::bithack::max(56, boost::numeric::bounds<int>::lowest()) << std::endl;

  std::cout << "min/max 56 and smallest" << std::endl;

  std::cout << utils::bithack::min(56, boost::numeric::bounds<int>::smallest()) << std::endl;
  std::cout << utils::bithack::max(56, boost::numeric::bounds<int>::smallest()) << std::endl;

  std::cout << "int8_t abs 56 and -56" << std::endl;
 
  std::cout << int(utils::bithack::abs((int8_t) 56)) << std::endl;
  std::cout << int(utils::bithack::abs((int8_t) -56)) << std::endl;

  std::cout << "int16_t abs 56 and -56" << std::endl;
  
  std::cout << utils::bithack::abs((int16_t) 56) << std::endl;
  std::cout << utils::bithack::abs((int16_t) - 56) << std::endl;

  std::cout << "int32_t abs 56 and -56" << std::endl;

  std::cout << utils::bithack::abs((int32_t) 56) << std::endl;
  std::cout << utils::bithack::abs((int32_t) - 56) << std::endl;

  std::cout << "int64_t abs 56 and -56" << std::endl;
  
  std::cout << utils::bithack::abs((int64_t) 56) << std::endl;
  std::cout << utils::bithack::abs((int64_t) - 56) << std::endl;
  
  std::cout << "bits int8_t/int16_t/int32_t/int64_t 86" << std::endl;

  std::cout << utils::bithack::bit_count((int8_t) 86) << std::endl;
  std::cout << utils::bithack::bit_count((int16_t) 86) << std::endl;
  std::cout << utils::bithack::bit_count((int32_t) 86) << std::endl;
  std::cout << utils::bithack::bit_count((int64_t) 86) << std::endl;
  
  std::cout << "significant bit 86" << std::endl;

  std::cout << utils::bithack::static_most_significant_bit<86>::result << std::endl;
  std::cout << utils::bithack::most_significant_bit((int8_t) 86) << std::endl;
  
  std::cout << "log_2 8  " << utils::bithack::floor_log2(8) << " " << utils::bithack::static_floor_log2<8>::result << std::endl;
  std::cout << "log_2 16 " << utils::bithack::floor_log2(16) << std::endl;
  std::cout << "log_2 32 " << utils::bithack::floor_log2(32) << std::endl;
  std::cout << "log_2 64 " << utils::bithack::floor_log2(64) << std::endl;

  std::cerr << "5 < 64 ? 0 : 5 = " << utils::bithack::branch(5 < 64, 0, 5) << std::endl;
  std::cerr << "5 < 3 ? 0 : 5 = " << utils::bithack::branch(5 < 3, 0, 5) << std::endl;
}
#endif
