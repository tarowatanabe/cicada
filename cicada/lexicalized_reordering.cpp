//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <iostream>
#include <sstream>

#include "lexicalized_reordering.hpp"
#include "quantizer.hpp"
#include "parameter.hpp"

#include "utils/compress_stream.hpp"
#include "utils/repository.hpp"
#include "utils/tempfile.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>


namespace cicada
{
  
  
  void LexicalizedReordering::ScoreSet::read(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    clear();

    repository_type rep(path, repository_type::read);
    
    if (boost::filesystem::exists(rep.path("quantized"))) {
      quantized.open(rep.path("quantized"));
      
      const path_type score_map_file = rep.path("score-map");
      
      if (! boost::filesystem::exists(score_map_file))
	throw std::runtime_error(std::string("no map file? ") + score_map_file.file_string());
      
      std::ifstream is(score_map_file.file_string().c_str());
      is.read((char*) &(*maps.begin()), sizeof(score_type) * maps.size());
    } else
      score.open(rep.path("score"));
  }
  
  void LexicalizedReordering::ScoreSet::write(const path_type& file) const
  {
    typedef utils::repository repository_type;
    
    if (path() == file) return;
    
    repository_type rep(file, repository_type::write);
    
    if (quantized.is_open()) {
      quantized.write(rep.path("quantized"));
      
      std::ofstream os(rep.path("score-map").file_string().c_str());
      os.write((char*) &(*maps.begin()), sizeof(score_type) * maps.size());
    } else
      score.write(rep.path("score"));
  }

  void LexicalizedReordering::quantize()
  {
    typedef score_type base_type;

    typedef std::map<base_type, size_type, std::less<base_type>, std::allocator<std::pair<const base_type, size_type> > > counts_type;
    typedef std::map<base_type, quantized_type, std::less<base_type>, std::allocator<std::pair<const base_type, quantized_type> > > codemap_type;
    typedef boost::array<base_type, 256> codebook_type;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    counts_type      counts;
    codebook_type    codebook;
    codemap_type     codemap;
    
    for (size_t feature = 0; feature < scores.size(); ++ feature)
      if (scores[feature].score.is_open()) {
	
	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.score.quantized.XXXXXX");
	utils::tempfile::insert(path);
	
	boost::iostreams::filtering_ostream os;
	os.push(utils::packed_sink<quantized_type, std::allocator<quantized_type> >(path));
	
	counts.clear();
	codemap.clear();
	std::fill(codebook.begin(), codebook.end(), 0.0);
	
	score_set_type::score_set_type::const_iterator liter_end = scores[feature].score.end();
	for (score_set_type::score_set_type::const_iterator liter = scores[feature].score.begin(); liter != liter_end; ++ liter)
	  ++ counts[*liter];
	
	Quantizer::quantize(counts, codebook, codemap);
	
	for (score_set_type::score_set_type::const_iterator liter = scores[feature].score.begin(); liter != liter_end; ++ liter) {
	  codemap_type::const_iterator citer = codemap.find(*liter);
	  if (citer == codemap.end())
	    throw std::runtime_error("no codemap?");
	  
	  os.write((char*) &(citer->second), sizeof(quantized_type));
	}
	
	for (int i = 0; i < 256; ++ i)
	  scores[feature].maps[i] = codebook[i];
	
	os.pop();
	utils::tempfile::permission(path);
	
	scores[feature].quantized.open(path);
	scores[feature].score.clear();
      }
  }
  
  void LexicalizedReordering::write(const path_type& file) const
  {
    if (path() == file) return;
    
    typedef utils::repository repository_type;
    
    repository_type rep(file, repository_type::write);

    phrases.write(rep.path("phrases"));
    
    const size_type feature_size = scores.size();
    for (size_t feature = 0; feature < feature_size; ++ feature) {
      std::ostringstream stream_score;
      stream_score << "score-" << std::setfill('0') << std::setw(6) << feature;
      
      scores[feature].write(rep.path(stream_score.str()));
    }
    
    vocab.write(rep.path("vocab"));
    
    rep["feature-size"] = boost::lexical_cast<std::string>(feature_size);
    
    rep["fe"]            = utils::lexical_cast<std::string>(fe);
    rep["bidirectional"] = utils::lexical_cast<std::string>(bidirectional);
    rep["monotonicity"]  = utils::lexical_cast<std::string>(monotonicity);
  }
  
  template <typename Rep>
  bool assign_parameter(const Rep& rep, const std::string& name)
  {
    typename Rep::const_iterator iter = rep.find(name);
    if (iter == rep.end())
      throw std::runtime_error("no parameter: " + name);
    return utils::lexical_cast<bool>(iter->second);
  }

  typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
  typedef std::vector<float, std::allocator<float> > scores_parsed_type;

  typedef boost::fusion::tuple<phrase_parsed_type, scores_parsed_type>                     phrase_scores_parsed_type;
  typedef boost::fusion::tuple<phrase_parsed_type, phrase_parsed_type, scores_parsed_type> phrase_pair_scores_parsed_type;
  
  template <typename Iterator>
  struct lexicalized_reordering_parser : boost::spirit::qi::grammar<Iterator, phrase_scores_parsed_type(), boost::spirit::standard::space_type>
  {
    lexicalized_reordering_parser() : lexicalized_reordering_parser::base_type(phrase_scores)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::lexeme;
      using qi::hold;
      using standard::char_;
      using standard::space;
      using qi::float_;
      
      phrase %= +lexeme[+(char_ - space) - "|||"];
      scores %= +float_;
      phrase_scores %= phrase >> "|||" >> scores;
    }
    
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type>        phrase;
    boost::spirit::qi::rule<Iterator, scores_parsed_type(), boost::spirit::standard::space_type>        scores;
    boost::spirit::qi::rule<Iterator, phrase_scores_parsed_type(), boost::spirit::standard::space_type> phrase_scores;
  };
  
  template <typename Iterator>
  struct lexicalized_reordering_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_scores_parsed_type(), boost::spirit::standard::space_type>
  {
    lexicalized_reordering_pair_parser() : lexicalized_reordering_pair_parser::base_type(phrase_scores)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::lexeme;
      using qi::hold;
      using standard::char_;
      using standard::space;
      using qi::float_;
      
      phrase %= +lexeme[+(char_ - space) - "|||"];
      scores %= +float_;
      phrase_scores %= phrase >> "|||" >> phrase >> "|||" >> scores;
    }
    
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type>        phrase;
    boost::spirit::qi::rule<Iterator, scores_parsed_type(), boost::spirit::standard::space_type>        scores;
    boost::spirit::qi::rule<Iterator, phrase_pair_scores_parsed_type(), boost::spirit::standard::space_type> phrase_scores;
  };
  
  
  void LexicalizedReordering::open(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    typedef utils::repository repository_type;

    clear();
    
    const parameter_type param(parameter);
    
    const path_type path = param.name();
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error("no lexicalized reordering file" + param.name());
    
    if (boost::filesystem::is_directory(path)) {
      repository_type rep(path, repository_type::read);
      
      phrases.open(rep.path("phrases"));
      vocab.open(rep.path("vocab"));
      
      repository_type::const_iterator iter = rep.find("feature-size");
      if (iter == rep.end())
	throw std::runtime_error("no feature size?");
      const size_type feature_size = atoi(iter->second.c_str());
      
      scores.reserve(feature_size);
      scores.resize(feature_size);
      
      for (size_t feature = 0; feature < feature_size; ++ feature) {
	std::ostringstream stream_score;
	stream_score << "score-" << std::setfill('0') << std::setw(6) << feature;
	
	scores[feature].read(rep.path(stream_score.str()));
      }

      fe            = assign_parameter(rep, "fe");
      bidirectional = assign_parameter(rep, "bidirectional");
      monotonicity  = assign_parameter(rep, "monotonicity");
      
      return;
    }
    
    // we will read from path....
    
    for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
      if (strcasecmp(piter->first.c_str(), "fe") == 0)
	fe = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "bidirectional") == 0)
	bidirectional = utils::lexical_cast<bool>(piter->second);
      else if (strcasecmp(piter->first.c_str(), "monotonicity") == 0)
	monotonicity = utils::lexical_cast<bool>(piter->second);
      else
	std::cerr << "WARNING: unsupported parameter for lexicalized-reordering: " << piter->first << "=" << piter->second << std::endl;
    }
    
    typedef utils::compress_ostream ostream_type;
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    typedef std::pair<ostream_ptr_type, path_type> stream_path_type;
    typedef std::vector<stream_path_type, std::allocator<stream_path_type> > stream_path_set_type;
    typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > sequence_type;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    const path_type path_phrases = utils::tempfile::directory_name(tmp_dir / "cicada.phrase.XXXXXX");
    const path_type path_vocab   = utils::tempfile::directory_name(tmp_dir / "cicada.vocab.XXXXXX");

    stream_path_set_type score_streams;
    phrases.open(path_phrases, phrase_db_type::WRITE);
    
    id_type id_phrase = 0;
    int     feature_size = -1;
    
    if (! fe) {
      utils::compress_istream is(path, 1024 * 1024);
      std::string line;
      
      lexicalized_reordering_parser<std::string::const_iterator> parser;
      phrase_scores_parsed_type parsed;
      sequence_type phrase;
      
      while (std::getline(is, line)) {
	boost::fusion::get<0>(parsed).clear();
	boost::fusion::get<1>(parsed).clear();
	
	std::string::const_iterator iter_end = line.end();
	std::string::const_iterator iter = line.begin();

	const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::space, parsed);
	
	if (! result || iter != iter_end) continue;
	
	phrase.clear();
	phrase_parsed_type::const_iterator piter_end = boost::fusion::get<0>(parsed).end();
	for (phrase_parsed_type::const_iterator piter = boost::fusion::get<0>(parsed).begin(); piter != piter_end; ++ piter)
	  phrase.push_back(word_type(*piter).id());
	
	phrases.insert(&(*phrase.begin()), phrase.size(), id_phrase);
	++ id_phrase;
	
	if (feature_size < 0) {
	  feature_size = boost::fusion::get<1>(parsed).size();
	  
	  if (bidirectional) {
	    if (monotonicity) {
	      if (feature_size != 4)
		throw std::runtime_error("we are bidirectional-monotonicity and expecting 4 features");
	    } else {
	      if (feature_size != 6)
		throw std::runtime_error("we are bidirectional and expecting 6 features");
	    }
	  } else {
	    if (monotonicity) {
	      if (feature_size != 4)
		throw std::runtime_error("we are monotonicity and expecting 2 features");
	    } else {
	      if (feature_size != 6)
		throw std::runtime_error("we are expecting 3 features");
	    }
	  }
	  
	  score_streams.reserve(feature_size);
	  score_streams.resize(feature_size);
	  
	  for (int feature = 0; feature < feature_size; ++ feature) {
	    score_streams[feature].second = utils::tempfile::file_name(tmp_dir / "cicada.feature.XXXXXX");
	    utils::tempfile::insert(score_streams[feature].second);
	    
	    score_streams[feature].first.reset(new utils::compress_ostream(score_streams[feature].second, 1024 * 1024));
	  }
	} else if (feature_size != static_cast<int>(boost::fusion::get<1>(parsed).size()))
	  throw std::runtime_error("invalid # of features...");
	
	for (int feature = 0; feature < feature_size; ++ feature)
	  score_streams[feature].first->write((char*) &boost::fusion::get<1>(parsed)[feature], sizeof(score_type));
      }
      
    } else {
      utils::compress_istream is(path, 1024 * 1024);
      std::string line;

      lexicalized_reordering_pair_parser<std::string::const_iterator> parser;
      phrase_pair_scores_parsed_type parsed;
      sequence_type phrase;
      
      while (std::getline(is, line)) {
	boost::fusion::get<0>(parsed).clear();
	boost::fusion::get<1>(parsed).clear();
	boost::fusion::get<2>(parsed).clear();
	
	std::string::const_iterator iter_end = line.end();
	std::string::const_iterator iter = line.begin();

	const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, parser, boost::spirit::standard::space, parsed);
	
	if (! result || iter != iter_end) continue;
	
	phrase.clear();
	{
	  phrase_parsed_type::const_iterator piter_end = boost::fusion::get<0>(parsed).end();
	  for (phrase_parsed_type::const_iterator piter = boost::fusion::get<0>(parsed).begin(); piter != piter_end; ++ piter)
	    phrase.push_back(word_type(*piter).id());
	}
	phrase.push_back(vocab_type::EMPTY.id());
	{
	  phrase_parsed_type::const_iterator piter_end = boost::fusion::get<1>(parsed).end();
	  for (phrase_parsed_type::const_iterator piter = boost::fusion::get<1>(parsed).begin(); piter != piter_end; ++ piter)
	    phrase.push_back(word_type(*piter).id());
	}
	
	phrases.insert(&(*phrase.begin()), phrase.size(), id_phrase);
	++ id_phrase;
	
	if (feature_size < 0) {
	  feature_size = boost::fusion::get<2>(parsed).size();
	  
	  if (bidirectional) {
	    if (monotonicity) {
	      if (feature_size != 4)
		throw std::runtime_error("we are bidirectional-monotonicity and expecting 4 features");
	    } else {
	      if (feature_size != 6)
		throw std::runtime_error("we are bidirectional and expecting 6 features");
	    }
	  } else {
	    if (monotonicity) {
	      if (feature_size != 4)
		throw std::runtime_error("we are monotonicity and expecting 2 features");
	    } else {
	      if (feature_size != 6)
		throw std::runtime_error("we are expecting 3 features");
	    }
	  }
	  
	  score_streams.reserve(feature_size);
	  score_streams.resize(feature_size);
	  
	  for (int feature = 0; feature < feature_size; ++ feature) {
	    score_streams[feature].second = utils::tempfile::file_name(tmp_dir / "cicada.feature.XXXXXX");
	    utils::tempfile::insert(score_streams[feature].second);
	    
	    score_streams[feature].first.reset(new utils::compress_ostream(score_streams[feature].second, 1024 * 1024));
	  }
	} else if (feature_size != static_cast<int>(boost::fusion::get<2>(parsed).size()))
	  throw std::runtime_error("invalid # of features...");
	
	for (int feature = 0; feature < feature_size; ++ feature)
	  score_streams[feature].first->write((char*) &boost::fusion::get<2>(parsed)[feature], sizeof(score_type));
      }
    }
    
    phrases.close();
    phrases.open(path_phrases);
    
    word_type::write(path_vocab);
    vocab.open(path_vocab);
    
    scores.reserve(feature_size);
    scores.resize(feature_size);
    for (int feature = 0; feature < feature_size; ++ feature) {
      score_streams[feature].first->reset();
      utils::tempfile::permission(score_streams[feature].second);
      scores[feature].score.open(score_streams[feature].second);
    }
  }
  
};
