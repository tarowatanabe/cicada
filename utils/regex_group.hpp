// -*- mode: c++ -*-

#ifndef __UTILS__REGEX_GROUP__HPP__
#define __UTILS__REGEX_GROUP__HPP__ 1

#include <stdint.h>

#include <stdexcept>
#include <string>
#include <vector>

#include <unicode/unistr.h>
#include <unicode/regex.h>

namespace utils
{
  
  // regex substituter
  class __regex_substitution
  {
  public:
    __regex_substitution() : pattern(), matcher() { clear(); }
    __regex_substitution(const std::string& __pattern,
			 const std::string& __substitution,
			 const bool __case_insensitive=false,
			 const bool __replace_all=false)
      : pattern(), matcher()
    {
      __initialize(UnicodeString::fromUTF8(__pattern), UnicodeString::fromUTF8(__substitution), __case_insensitive, __replace_all);
    }
    __regex_substitution(const UnicodeString& __pattern,
			 const UnicodeString& __substitution,
			 const bool __case_insensitive=false,
			 const bool __replace_all=false)
      : pattern(), matcher()
    {
      __initialize(__pattern, __substitution, __case_insensitive, __replace_all);
    }
    __regex_substitution(const __regex_substitution& x)
      : pattern(), matcher() { assign(x); }
    ~__regex_substitution() { clear(); }
    
    __regex_substitution& operator=(const __regex_substitution& x)
    {
      assign(x);
      return *this;
    }
    
    UnicodeString& operator()(UnicodeString& data) const
    {
      if (replace_all) {
	while (1) {
	  const_cast<RegexMatcher*>(matcher)->reset(data);
	  if (! const_cast<RegexMatcher*>(matcher)->find()) break;
	  
	  UErrorCode status = U_ZERO_ERROR;
	  data = const_cast<RegexMatcher*>(matcher)->replaceAll(substitution, status);
	  if (U_FAILURE(status))
	    throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
	}
      } else {
	UErrorCode status = U_ZERO_ERROR;
	const_cast<RegexMatcher*>(matcher)->reset(data);
	data = const_cast<RegexMatcher*>(matcher)->replaceAll(substitution, status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("RegexMatcher::replaceAll(): ") + u_errorName(status));
      }
      return data;
    }
    
    void assign(const __regex_substitution& x)
    {
      clear();
      
      if (x.pattern)
	pattern = x.pattern->clone();
      
      if (pattern) {
	UErrorCode status = U_ZERO_ERROR;
	matcher = pattern->matcher(status);
	if (U_FAILURE(status))
	  throw std::runtime_error(std::string("Pattern::Mathcer() ") + u_errorName(status));
      }
      substitution = x.substitution;
      replace_all = x.replace_all;
    }
    
    void clear()
    {
      if (matcher)
	delete matcher;
      if (pattern)
	delete pattern;
      
      matcher = 0;
      pattern = 0;
      substitution = UnicodeString();
      replace_all = false;
    }

  private:
    
    void __initialize(const UnicodeString& __pattern,
		      const UnicodeString& __substitution,
		      const bool __case_insensitive=false,
		      const bool __replace_all=false)
    {
      clear();
      
      UErrorCode status = U_ZERO_ERROR;
      pattern = RegexPattern::compile(__pattern, __case_insensitive ? UREGEX_CASE_INSENSITIVE : 0, status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RegexPattern::compile() ") + u_errorName(status));
      
      status = U_ZERO_ERROR;
      matcher = pattern->matcher(status);
      if (U_FAILURE(status))
	throw std::runtime_error(std::string("RegexPattern::matcher() ") + u_errorName(status));
      
      substitution = __substitution;
      replace_all = __replace_all;
    }
    
  private:
    RegexPattern* pattern;
    RegexMatcher* matcher;
    UnicodeString substitution;
    bool replace_all;
  };
  

  template <typename __RegexGroup, typename __Regex>
  struct __regex_substitution_assign
  {
    typedef __RegexGroup regex_group_type;
    typedef __Regex      regex_type;
    typedef __regex_substitution_assign<regex_group_type, regex_type> self_type;
    
    __regex_substitution_assign(regex_group_type& __owner) : owner(__owner) {}
    
    self_type& operator()(const char* pattern, const char* substitution, const bool insensitive=false, const bool subst_all=false)
    {
      return operator()(std::string(pattern), std::string(substitution), insensitive, subst_all);
    }
    
    self_type& operator()(const std::string& pattern, const std::string& substitution, const bool insensitive=false, const bool subst_all=false)
    {
      owner.push_back(regex_type(pattern, substitution, insensitive, subst_all));
      return *this;
    }
    
    regex_group_type& owner;
  };

  
  class regex_group
  {
  public:
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef __regex_substitution regex_type;
    typedef __regex_substitution_assign<regex_group, regex_type> assign_type;
    
    typedef std::vector<regex_type, std::allocator<regex_type> > regex_set_type;
    
  public:
    regex_group() : group() {}
    
  public:
    
    UnicodeString& operator()(UnicodeString& data) const
    {
      regex_set_type::const_iterator riter_end = group.end();
      for (regex_set_type::const_iterator riter = group.begin(); riter != riter_end; ++ riter)
	riter->operator()(data);
      return data;
    }
    
    assign_type insert()
    {
      return assign_type(*this);
    }
    
    size_t size() const { return group.size(); }
    bool empty() const { return group.empty(); }
    void clear() { group.clear(); }
    
    void push_back(const regex_type& x)
    {
      group.push_back(x);
    }
    
  private:
    regex_set_type group;
  };
  
};

#endif
