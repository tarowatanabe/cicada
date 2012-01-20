//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include "lexicon.hpp"
#include "parameter.hpp"

#include "utils/compress_stream.hpp"
#include "utils/repository.hpp"
#include "utils/tempfile.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/spinlock.hpp"
#include "utils/sgi_hash_map.hpp"
#include "utils/thread_specific_ptr.hpp"

#include <boost/thread.hpp>

namespace cicada
{
  void Lexicon::open(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const path_type path = param.name();
    
    typedef utils::repository repository_type;
    
    clear();

    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error("no file? " + path.string());
    
    if (boost::filesystem::is_directory(path)) {
      repository_type rep(path, repository_type::read);
      
      lexicon.open(rep.path("lexicon"));
      vocab.open(rep.path("vocab"));
      
      repository_type::const_iterator siter = rep.find("smooth");
      if (siter == rep.end())
	throw std::runtime_error("no smoothing parameter...?");
      
      smooth = boost::lexical_cast<weight_type>(siter->second);
    } else {
      bool feature_mode = false;
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "feature")
	  feature_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
      
      // feature-mode will index by the series of tokens + the last token

      const path_type tmp_dir = utils::tempfile::tmp_dir();
      const path_type path_tmp = utils::tempfile::directory_name(tmp_dir / "cicada.lexicon.XXXXXX");
      
      utils::tempfile::insert(path_tmp);
      repository_type rep(path_tmp, repository_type::write);
      
      const path_type lexicon_path = rep.path("lexicon");
      const path_type vocab_path   = rep.path("vocab");
      
      lexicon.open(lexicon_path, lexicon_type::WRITE);
      
      if (feature_mode) {
	typedef std::vector<word_id_type, std::allocator<word_id_type> > code_set_type;
	typedef std::vector<std::string, std::allocator<std::string> > lexicon_parsed_type;
	typedef boost::spirit::istream_iterator iterator_type;
	
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	qi::rule<iterator_type, std::string(), standard::blank_type>         word;
	qi::rule<iterator_type, lexicon_parsed_type(), standard::blank_type> parser; 
	
	word   %= qi::lexeme[+(standard::char_ - standard::space)];
	parser %= *word >> (qi::eol | qi::eoi);
	
	utils::compress_istream is(path, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iterator_type iter(is);
	iterator_type iter_end;
	
	lexicon_parsed_type lexicon_parsed;
	code_set_type codes;

	smooth = std::numeric_limits<weight_type>::infinity();
	
	while (iter != iter_end) {
	  lexicon_parsed.clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, lexicon_parsed))
	    if (iter != iter_end)
	      throw std::runtime_error("global lexicon parsing failed");
	  
	  if (lexicon_parsed.size() < 2) continue;
	  
	  const weight_type weight(utils::lexical_cast<weight_type>(lexicon_parsed.back()));
	  
	  codes.clear();
	  lexicon_parsed_type::const_iterator liter_end = lexicon_parsed.end() - 1;
	  for (lexicon_parsed_type::const_iterator liter = lexicon_parsed.begin(); liter != liter_end; ++ liter)
	    codes.push_back(word_type(*liter).id());
	  
	  lexicon.insert(&(*codes.begin()), codes.size(), weight);
	  
	  smooth = std::min(smooth, weight);
	}

      } else {
	typedef boost::fusion::tuple<std::string, std::string, weight_type > lexicon_parsed_type;
	typedef boost::spirit::istream_iterator iterator_type;
	
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	qi::rule<iterator_type, std::string(), standard::blank_type>         word;
	qi::rule<iterator_type, lexicon_parsed_type(), standard::blank_type> parser; 
	
	word   %= qi::lexeme[+(standard::char_ - standard::space)];
	parser %= word >> word >> qi::float_ >> (qi::eol | qi::eoi); // weight type!
  
	utils::compress_istream is(path, 1024 * 1024);
	is.unsetf(std::ios::skipws);
	
	iterator_type iter(is);
	iterator_type iter_end;
	
	lexicon_parsed_type lexicon_parsed;
	
	word_id_type codes[2];
	smooth = std::numeric_limits<weight_type>::infinity();
	
	while (iter != iter_end) {
	  boost::fusion::get<0>(lexicon_parsed).clear();
	  boost::fusion::get<1>(lexicon_parsed).clear();
	  
	  if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, lexicon_parsed))
	    if (iter != iter_end)
	      throw std::runtime_error("global lexicon parsing failed");
	  
	  codes[0] = word_type(boost::fusion::get<1>(lexicon_parsed)).id(); // source
	  codes[1] = word_type(boost::fusion::get<0>(lexicon_parsed)).id(); // target
	  
	  lexicon.insert(codes, 2, boost::fusion::get<2>(lexicon_parsed));

	  smooth = std::min(smooth, boost::fusion::get<2>(lexicon_parsed));
	}
      }

      if (smooth == std::numeric_limits<weight_type>::infinity())
	smooth = 1e-40;
      
      rep["smooth"] = boost::lexical_cast<std::string>(smooth);

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
  
  void Lexicon::write(const path_type& file) const
  {
    if (path() == file) return;
    
    typedef utils::repository repository_type;
    
    repository_type rep(file, repository_type::write);
    
    lexicon.write(rep.path("lexicon"));
    vocab.write(rep.path("vocab"));
    
    rep["smooth"] = boost::lexical_cast<std::string>(smooth);
  }
  
  template <typename Tp>
  struct hash_string : public utils::hashmurmur<size_t>
  {
    size_t operator()(const Tp& x) const
    {
      return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
    }
  };

#ifdef HAVE_TR1_UNORDERED_MAP
  typedef std::tr1::unordered_map<std::string, Lexicon, hash_string<std::string>, std::equal_to<std::string>,
				  std::allocator<std::pair<const std::string, Lexicon> > > lexicon_map_type;
#else
  typedef sgi::hash_map<std::string, Lexicon, hash_string<std::string>, std::equal_to<std::string>,
			std::allocator<std::pair<const std::string, Lexicon> > > lexicon_map_type;
#endif

  namespace impl
  {
    typedef utils::spinlock             mutex_type;
    typedef mutex_type::scoped_lock     lock_type;
    
    static mutex_type       __lexicon_mutex;
    static lexicon_map_type __lexicon_map;
  };

#ifdef HAVE_TLS
  static __thread lexicon_map_type* __lexicons_tls = 0;
  static boost::thread_specific_ptr<lexicon_map_type> __lexicons;
#else
  static utils::thread_specific_ptr<lexicon_map_type> __lexicons;
#endif

  Lexicon& Lexicon::create(const std::string& parameter)
  {
    
    
#ifdef HAVE_TLS
    if (! __lexicons_tls) {
      __lexicons.reset(new lexicon_map_type());
      __lexicons_tls = __lexicons.get();
    }
    lexicon_map_type& lexicons_map = *__lexicons_tls;    
#else
    if (! __lexicons.get())
      __lexicons.reset(new lexicon_map_type());
    
    lexicon_map_type& lexicons_map = *__lexicons;
#endif
    
    lexicon_map_type::iterator iter = lexicons_map.find(parameter);
    if (iter == lexicons_map.end()) {
      impl::lock_type lock(impl::__lexicon_mutex);
      
      lexicon_map_type::iterator iter_global = impl::__lexicon_map.find(parameter);
      if (iter_global == impl::__lexicon_map.end())
	iter_global = impl::__lexicon_map.insert(std::make_pair(parameter, Lexicon(parameter))).first;
      
      iter = lexicons_map.insert(*iter_global).first;
    }
    
    return iter->second;
  }


};
