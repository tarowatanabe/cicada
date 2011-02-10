//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/thread.hpp>

#include <iostream>
#include <sstream>

#include "global_lexicon.hpp"

#include "utils/compress_stream.hpp"
#include "utils/repository.hpp"
#include "utils/tempfile.hpp"
#include "utils/lexical_cast.hpp"

namespace cicada
{
  void GlobalLexicon::write(const path_type& file) const
  {
    if (path() == file) return;
    
    typedef utils::repository repository_type;
    
    repository_type rep(file, repository_type::write);
    
    lexicon.write(rep.path("lexicon"));
    vocab.write(rep.path("vocab"));
  }

  template <typename Lexicon, typename Word, typename Iterator>
  inline
  void dump(Lexicon& lexicon, const Word& target, Iterator first, Iterator last)
  {
    typedef typename Lexicon::key_type key_type;

    key_type codes[2];

    codes[0] = target.id();
    for (/**/; first != last; ++ first) {
      codes[1] = first->first.id();
      lexicon.insert(codes, 2, first->second);
    }
  }
  
  void GlobalLexicon::open(const path_type& path)
  {
    typedef utils::repository repository_type;

    clear();

    repository_type rep(path, repository_type::read);
    
    if (boost::filesystem::exists(rep.path("lexicon")) && boost::filesystem::exists(rep.path("vocab"))) {
      lexicon.open(rep.path("lexicon"));
      vocab.open(rep.path("vocab"));
    } else {
      typedef boost::fusion::tuple<std::string, std::string, double > lexicon_parsed_type;
      typedef boost::spirit::istream_iterator iterator_type;
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      qi::rule<iterator_type, std::string(), standard::blank_type>         word;
      qi::rule<iterator_type, lexicon_parsed_type(), standard::blank_type> parser; 
      
      word   %= qi::lexeme[+(standard::char_ - standard::space)];
      parser %= word >> word >> qi::double_ >> (qi::eol | qi::eoi);
      
      repository_type::const_iterator iter = rep.find("size");
      if (iter == rep.end())
	throw std::runtime_error("no size?");

      const size_type size = boost::lexical_cast<size_type>(iter->second);
      
      const path_type tmp_dir = utils::tempfile::tmp_dir();
      
      const path_type lexicon_path = utils::tempfile::directory_name(tmp_dir / "cicada.lexicon.XXXXXX");
      const path_type vocab_path   = utils::tempfile::directory_name(tmp_dir / "cicada.vocab.XXXXXX");
      
      utils::tempfile::insert(lexicon_path);
      utils::tempfile::insert(vocab_path);

      lexicon.open(lexicon_path, lexicon_type::WRITE);
      
      for (size_t node = 0; node < size; ++ node) {
	typedef std::pair<word_type, weight_type> word_weight_type;
	typedef std::vector<word_weight_type, std::allocator<word_weight_type> > word_weight_set_type;

	const path_type path = rep.path(std::string("lexicon.") + boost::lexical_cast<std::string>(node) + ".gz");
	
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error(std::string("no lexicon file: ") + path.file_string());
	
	utils::compress_istream is(path, 1024 * 1024);
	is.unsetf(std::ios::skipws);

	iterator_type iter(is);
	iterator_type iter_end;
  
	lexicon_parsed_type lexicon_parsed;
	
	word_weight_set_type weights;
	word_type target_prev;

	while (iter != iter_end) {
	  boost::fusion::get<0>(lexicon_parsed).clear();
	  boost::fusion::get<1>(lexicon_parsed).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, lexicon_parsed))
	    if (iter != iter_end)
	      throw std::runtime_error("global lexicon parsing failed");
	  
	  const word_type target(boost::fusion::get<0>(lexicon_parsed));
	  const word_type source(boost::fusion::get<1>(lexicon_parsed));
	  const weight_type weight(boost::fusion::get<2>(lexicon_parsed));
	  
	  if (target != target_prev) {
	    if (! weights.empty())
	      dump(lexicon, target_prev, weights.begin(), weights.end());
	    
	    weights.clear();
	    target_prev = target;
	  }
	  
	  weights.push_back(std::make_pair(source, weight));
	}
	
	if (! weights.empty())
	  dump(lexicon, target_prev, weights.begin(), weights.end());
	weights.clear();
      }
      
      lexicon.close();
      word_type::write(vocab_path);
      
      ::sync();
      
      while (! lexicon_type::exists(lexicon_path))
	boost::thread::yield();
      while (! vocab_type::exists(vocab_path))
	boost::thread::yield();
      
      lexicon.open(lexicon_path);
      vocab.open(vocab_path);
    }
  }

};
