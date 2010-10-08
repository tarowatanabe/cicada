#include <iostream>
#include <sstream>

#include "global_lexicon.hpp"

#include "utils/compress_stream.hpp"
#include "utils/repository.hpp"
#include "utils/tempfile.hpp"
#include "utils/space_separator.hpp"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


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
      typedef boost::tokenizer<utils::space_separator> tokenizer_type;
      
      repository_type::const_iterator iter = rep.find("size");
      if (iter == rep.end())
	throw std::runtime_error("no size?");

      const size_type size = atoi(iter->second.c_str());
      
      const path_type tmp_dir = utils::tempfile::tmp_dir();
      
      const path_type lexicon_path = utils::tempfile::directory_name(tmp_dir / "cicada.lexicon.XXXXXX");
      const path_type vocab_path   = utils::tempfile::directory_name(tmp_dir / "cicada.vocab.XXXXXX");
      
      utils::tempfile::insert(lexicon_path);
      utils::tempfile::insert(vocab_path);

      lexicon.open(lexicon_path, lexicon_type::WRITE);
      
      for (int node = 0; node < size; ++ node) {
	typedef std::pair<word_type, weight_type> word_weight_type;
	typedef std::vector<word_weight_type, std::allocator<word_weight_type> > word_weight_set_type;

	const path_type path = rep.path(std::string("lexicon.") + boost::lexical_cast<std::string>(node) + ".gz");
	
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error(std::string("no lexicon file: ") + path.file_string());
	
	utils::compress_istream is(path, 1024 * 1024);
	
	std::string line;
	word_weight_set_type weights;
	word_type target_prev;
	
	while (std::getline(is, line)) {
	  tokenizer_type tokenizer(line);
	  
	  tokenizer_type::iterator iter = tokenizer.begin();
	  if (iter == tokenizer.end()) continue;
	  
	  const std::string word1 = *iter;
	  ++ iter;
	  if (iter == tokenizer.end()) continue;
	  
	  const std::string word2 = *iter;
	  ++ iter;
	  if (iter == tokenizer.end()) continue;
	  
	  const std::string word3 = *iter;
	  
	  const word_type target(word1);
	  const word_type source(word2);
	  
	  const weight_type weight(boost::lexical_cast<weight_type>(word3));
	  
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
      
      lexicon.open(lexicon_path);
      vocab.open(vocab_path);
    }
  }

};
