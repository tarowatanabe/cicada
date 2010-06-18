#include <stdexcept>

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/filesystem.hpp>
#include <utils/compress_stream.hpp>

template <typename Path, typename Weights>
void read_weights(const Path& path, Weights& weights)
{
  if (path.empty()) return;
  
  if (path != "-" && ! boost::filesystem::exists(path))
    throw std::runtime_error("no feture weights?" + path.file_string());
      
  utils::compress_istream is(path);
  is >> weights;
}
		  

template <typename Iterator>
inline
bool parse_id(size_t& id, Iterator& iter, Iterator end)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  namespace phoenix = boost::phoenix;
  
  using qi::phrase_parse;
  using qi::_1;
  using qi::ulong_;
  using standard::space;
  
  using phoenix::ref;
  
  return phrase_parse(iter, end, ulong_ [ref(id) = _1] >> "|||", space);
}

template <typename Iterator>
inline
bool parse_separator(Iterator& iter, Iterator end)
{
  namespace qi = boost::spirit::qi;
  namespace standard = boost::spirit::standard;
  namespace phoenix = boost::phoenix;
  
  using qi::phrase_parse;
  using standard::space;
  
  return phrase_parse(iter, end, "|||", space);
}

template <typename HyperGraph, typename Lattice, typename SentenceSet, typename Sentence>
inline
bool parse_line(const std::string& line,
		size_t& id,
		HyperGraph& hypergraph,
		Lattice& lattice,
		Lattice& target,
		SentenceSet& target_sentences,
		Sentence& sentence,
		const bool input_id,
		const bool input_lattice,
		const bool input_forest, 
		const bool input_bitext)
{
  std::string::const_iterator iter = line.begin();
  std::string::const_iterator end = line.end();
  
  if (input_id)
    if (! parse_id(id, iter, end))
      throw std::runtime_error("invalid id-prefixed format");
  
  if (input_lattice) {
    if (! lattice.assign(iter, end))
      throw std::runtime_error("invalid lattive format");
  } else if (input_forest) {
    if (! hypergraph.assign(iter, end))
	  throw std::runtime_error("invalid hypergraph format");
  } else {
    if (! sentence.assign(iter, end))
      throw std::runtime_error("invalid sentence format");
    
    lattice = Lattice(sentence);
  }
  
  if (input_bitext) {
    target_sentences.clear();
    
    while (parse_separator(iter, end)) {
      target_sentences.push_back(Sentence());
      
      if (! target_sentences.back().assign(iter, end))
	throw std::runtime_error("invalid sentence format");
    }
    
    if (target_sentences.empty())
      throw std::runtime_error("no bitext?");
    
    target = Lattice(target_sentences.front());
  }
  
  return iter == end;
}
