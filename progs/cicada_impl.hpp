
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
  if (! path.empty()) return;
  
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
