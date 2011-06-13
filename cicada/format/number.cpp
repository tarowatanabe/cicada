
#include <vector>

#include <unicode/format.h>
#include <unicode/rbnf.h>
#include <unicode/unistr.h>
#include <unicode/locid.h>
#include <unicode/bytestream.h>

#include "number.hpp"

namespace cicada
{
  namespace foramt
  {
    class NumberImpl
    {
    public:
      typedef icu::NumberFormat parser_type;
      typedef icu::NumberFormat generator_type;
      
      typedef std::vector<parser_type*, std::allocator<parser_type*> >       parser_set_type;
      typedef std::vector<generator_type*, std::allocator<generator_type*> > generator_set_type;
      
      typedef cicada::Format::phrase_type     phrase_type;
      typedef cicada::Format::phrase_set_type phrase_set_type;

    public:
      NumberImpl() {}
      NumberImpl(const NumberImpl& x) { assign(x); }
      NumberImpl& operator=(const NumberImpl& x)
      {
	assign(x);
	return *this;
      }
      
      ~NumberImpl() { clear(); }
      
    public:
      void operator()(const phrase_type& input, phrase_set_type& output) const
      {
	typedef std::set<std::string, std::less<std::string>, std::allocator<std::string> > phrase_unique_type;

	UnicodeString uphrase = UnicodeString::fromUTF8(input);
	UnicodeString ugenerated;
	std::string   generated;

	phrase_unique_type uniques;
	
	parser_set_type::const_iterator piter_end = parsers.end();
	for (parser_set_type::const_iterator piter = parsers.begin(); piter != piter_end; ++ piter) {
	  const parser_type* parser = *piter;
	  icu::Formattable   formattable;
	  icu::UErrorCode    status(U_ZERO_ERROR);
	  
	  parser->parse(uphrase, formattable, status);
	  
	  if (U_FAILURE(status)) continue;
	  
	  generator_set_type::const_iterator giter_end = generators.end();
	  for (generator_set_type::const_iterator giter = generators.begin(); giter != giter_end; ++ giter) {
	    const generator_type* generator = *giter;
	    icu::UErrorCode    status(U_ZERO_ERROR);
	    
	    ugenerated.remove();
	    generator->format(formattable, ugenerated, status);
	    
	    if (U_FAILURE(status)) continue;
	    
	    generated.clear();
	    ugenerated.toUTF8String(generated);
	    uniques.insert(generated);
	  }
	}
	
	output.clear();
	output.insert(output.end(), uniques.begin(), uniques.end());
      }
      
    public:
      void assign(const NumberImpl& x)
      {
	clear();
	
	parser_set_type::const_iterator piter_end = x.parsers.end();
	for (parser_set_type::const_iterator piter = x.parsers.begin(); piter != piter_end; ++ piter)
	  parsers.push_back(piter->clone());
	
	generator_set_type::const_iterator giter_end = x.generators.end();
	for (generator_set_type::const_iterator giter = x.generators.begin(); giter != giter_end; ++ giter)
	  generators.push_back(giter->clone());
      }
      
      void clear()
      {
	parser_set_type::iterator piter_end = parsers.end();
	for (parser_set_type::iterator piter = parsers.begin(); piter != piter_end; ++ piter)
	  delete *piter;
	
	generator_set_type::iterator giter_end = generators.end();
	for (generator_set_type::iterator giter = generators.begin(); giter != giter_end; ++ giter)
	  delete *giter;
	
	parsers.clear();
	generators.clear();
      }
      
    private:
      parser_set_type    parsers;
      generator_set_tyep generators;
    };
    
    
  };
};
