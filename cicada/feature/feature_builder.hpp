// -*- mode: c++ -*-
//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__FEATURE_BUILDER__HPP__
#define __CICADA__FEATURE__FEATURE_BUILDER__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE

#include <boost/spirit/include/karma.hpp>

#include <cstring>

#include <string>
#include <vector>
#include <iterator>

#include <cicada/Attribute.hpp>
#include <cicada/Feature.hpp>
#include <cicada/Symbol.hpp>

#include <utils/piece.hpp>

namespace cicada
{
  namespace feature
  {
    struct FeatureBuilder
    {
      typedef Attribute attribute_type;
      typedef Symbol    symbol_type;
      typedef Feature   feature_type;

      typedef std::vector<char, std::allocator<char> > buffer_type;
	
      typedef buffer_type::iterator       iterator;
      typedef buffer_type::const_iterator const_iterator;

      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const symbol_type& x)
      {
	builder.buffer.insert(builder.buffer.end(), x.begin(), x.end());
	return builder;
      }

      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const attribute_type& x)
      {
	builder.buffer.insert(builder.buffer.end(), x.begin(), x.end());
	return builder;
      }

      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const feature_type& x)
      {
	builder.buffer.insert(builder.buffer.end(), x.begin(), x.end());
	return builder;
      }

	
      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const std::string& x)
      {
	builder.buffer.insert(builder.buffer.end(), x.begin(), x.end());
	return builder;
      }

      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const char* x)
      {
	builder.buffer.insert(builder.buffer.end(), x, x + std::strlen(x));
	return builder;
      }
      
      template <size_t N>
      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const char (&x)[N])
      {
	builder.buffer.insert(builder.buffer.end(), x, x + N);
	return builder;
      }
      
      friend
      FeatureBuilder& operator<<(FeatureBuilder& builder, const int& x)
      {
	namespace karma = boost::spirit::karma;
	namespace standard = boost::spirit::standard;
	
	std::back_insert_iterator<buffer_type> iter(builder.buffer);
	karma::generate(iter, karma::int_, x);
	return builder;
      }

      void clear() { buffer.clear(); }
	
      bool exists() const
      {
	return feature_type::exists(utils::piece(buffer.begin(), buffer.end()));
      }
      
      operator feature_type() const { return feature_type(utils::piece(buffer.begin(), buffer.end())); }
      
      operator std::string() const { return std::string(buffer.begin(), buffer.end()); }
	
    private:
      buffer_type buffer;
    };
  };
};

#endif
