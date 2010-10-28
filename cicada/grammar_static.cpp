
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "grammar_static.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "succinct_db/succinct_trie_database.hpp"
#include "succinct_db/succinct_hash.hpp"

#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/tempfile.hpp"
#include "utils/group_aligned_code.hpp"
#include "utils/byte_aligned_code.hpp"
#include "utils/simple_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/arc_list.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/filesystem.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

namespace cicada
{
  
  // TODO: we need to index by index-stripped non-terminals!
  
  struct GrammarStaticImpl : public utils::hashmurmur<uint64_t>
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol  symbol_type;
    typedef cicada::Symbol  word_type;
    typedef cicada::Feature feature_type;
    typedef cicada::Vocab   vocab_type;
    
    typedef Transducer::rule_type     rule_type;
    typedef Transducer::rule_ptr_type rule_ptr_type;
    typedef Transducer::rule_set_type rule_set_type;
    
    typedef rule_type::feature_set_type feature_set_type;
    
    typedef rule_type::symbol_set_type symbol_set_type;
    typedef rule_type::symbol_set_type phrase_type;

    
    
    typedef float          score_type;
    typedef uint8_t        quantized_type;
    
    typedef uint32_t       id_type;

    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;
    
    typedef boost::filesystem::path path_type;

    typedef char byte_type;
    typedef char mapped_type;
    
    typedef std::allocator<std::pair<word_type::id_type, mapped_type> > rule_alloc_type;
    
    typedef succinctdb::succinct_trie_database<word_type::id_type, mapped_type, rule_alloc_type > rule_db_type;
    typedef succinctdb::succinct_hash_mapped<byte_type, std::allocator<byte_type> > phrase_db_type;

    typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;

    class ScoreSet
    {
    public:
      typedef utils::map_file<score_type, std::allocator<score_type> >                     score_set_type;
      typedef utils::packed_vector_mapped<quantized_type, std::allocator<quantized_type> > quantized_set_type;
      typedef boost::array<score_type, 256>                                                score_map_type;
      
      ScoreSet() {}
      ScoreSet(const path_type& path) { read(path); }
      
      void read(const path_type& path);
      void write(const path_type& file) const;
      
      void clear()
      {
	score.clear();
	quantized.clear();
      }
      void close() { clear(); }
      
      score_type operator[](size_type pos) const
      {
	return (quantized.is_open()
		? maps[quantized[pos]]
		: score[pos]);
      }
      
      path_type path() const
      {
	return (quantized.is_open()
		? quantized.path().parent_path()
		: score.path().parent_path());
      }
      bool empty() const { return quantized.empty() && score.empty(); }
      size_type size() const
      {
	return (quantized.is_open()
		? quantized.size()
		: score.size());
      }
      
      score_set_type     score;
      quantized_set_type quantized;
      score_map_type     maps;
    };
    
    typedef ScoreSet score_set_type;
    typedef std::vector<score_set_type, std::allocator<score_set_type> > score_db_type;
    
    // caching...
    template <typename Tp>
    struct __cache_pos
    {
      Tp        value;
      size_type pos;
      
      __cache_pos() : value(), pos(size_type(-1)) {}
    };
    
    typedef utils::arc_list<size_type, rule_set_type, 16,
			    std::equal_to<size_type>,
			    std::allocator<std::pair<size_type, rule_set_type> > > cache_rule_set_type;
    typedef __cache_pos<phrase_type>            cache_phrase_type;
    
    typedef utils::array_power2<cache_rule_set_type, 1024 * 16, std::allocator<cache_rule_set_type> > cache_rule_map_type;
    typedef utils::array_power2<cache_phrase_type,   1024 * 16, std::allocator<cache_phrase_type> >   cache_phrase_set_type;
            

  public:
    GrammarStaticImpl(const std::string& parameter) : max_span(15) { read(parameter); }

    GrammarStaticImpl(const GrammarStaticImpl& x)
      : rule_db(x.rule_db),
	source_db(x.source_db),
	target_db(x.target_db),
	score_db(x.score_db),
	vocab(x.vocab),
	feature_names(x.feature_names),
	max_span(x.max_span) {}

    GrammarStaticImpl& operator=(const GrammarStaticImpl& x)
    {
      clear();
      
      rule_db       = x.rule_db;
      source_db     = x.source_db;
      target_db     = x.target_db;
      score_db      = x.score_db;
      vocab         = x.vocab;
      feature_names = x.feature_names;
      max_span      = x.max_span;
      
      return *this;
    }
    
  public:

    void clear()
    {
      rule_db.clear();
      source_db.clear();
      target_db.clear();
      score_db.clear();
      vocab.clear();
      feature_names.clear();

      cache_rule_sets.clear();
      cache_sources.clear();
      cache_targets.clear();

      max_span = 15;
    }
    
    size_type find(const word_type& word) const
    {
      size_type node = 0;
      return find(word, node);
    }
    
    size_type find(const word_type& word, size_type node) const
    {
      const word_type::id_type id = vocab[word];
      return rule_db.find(&id, 1, node);
    }
    
    template <typename Iterator>
    size_type find(Iterator first, Iterator last) const
    {
      size_type node = 0;
      return find(first, last, node);
    }
    
    template <typename Iterator>
    size_type find(Iterator first, Iterator last, size_type node) const
    {
      for (/**/; first != last && is_valid(node); ++ first)
	node = find(*first, node);
      return node;
    }
    
    // valid implies that you can continue searching from node...
    bool is_valid(size_type node) const { return rule_db.is_valid(node); }
    bool has_children(size_type node) const { return rule_db.has_children(node); }
    
    // exists implies data associated with the node exists...
    bool exists(size_type node) const { return rule_db.exists(node); }
    
    const rule_set_type& read_rule_set(size_type node) const
    {
      const size_type cache_pos = hasher_type::operator()(node) & (cache_rule_sets.size() - 1);
      
      cache_rule_set_type& cache = const_cast<cache_rule_set_type&>(cache_rule_sets[cache_pos]);
      
      std::pair<cache_rule_set_type::iterator, bool> result = cache.find(node);
      if (! result.second) {
	typedef std::vector<byte_type, std::allocator<byte_type> >  code_set_type;
	
	rule_set_type& options = result.first->second;
	options.clear();
	
	rule_db_type::cursor cursor_end = rule_db.cend(node);
	for (rule_db_type::cursor cursor = rule_db.cbegin(node); cursor != cursor_end; ++ cursor) {
	  const size_type pos = cursor.node();
	  
	  code_set_type codes(rule_db[pos].begin(), rule_db[pos].end());
	  
	  code_set_type::const_iterator hiter = codes.begin();
	  code_set_type::const_iterator citer = codes.begin();
	  code_set_type::const_iterator citer_end = codes.end();
	  
	  id_type pos_feature;
	  id_type pos_source;
	  id_type pos_target;
	  
	  size_type code_pos = 0;
	  
	  const size_type offset_feature = utils::group_aligned_decode(pos_feature, &(*hiter), code_pos);
	  citer = hiter + offset_feature;
	  hiter += offset_feature & (- size_type((code_pos & 0x03) == 0x03));
	  ++ code_pos;
	  
	  const size_type offset_source = utils::group_aligned_decode(pos_source, &(*hiter), code_pos);
	  citer = hiter + offset_source;
	  hiter += offset_source & (- size_type((code_pos & 0x03) == 0x03));
	  ++ code_pos;
	  
	  const phrase_type rule_source = read_phrase(pos_source, cache_sources, source_db);
	  const size_type   rule_arity = rule_source.arity();
	  
	  while (citer != citer_end) {
	    
	    id_type id_lhs;
	    const size_type offset_lhs = utils::group_aligned_decode(id_lhs, &(*hiter), code_pos);
	    citer = hiter + offset_lhs;
	    hiter += offset_lhs & (- size_type((code_pos & 0x03) == 0x03));
	    ++ code_pos;
	    
	    const size_type offset_target = utils::group_aligned_decode(pos_target, &(*hiter), code_pos);
	    citer = hiter + offset_target;
	    hiter += offset_target & (- size_type((code_pos & 0x03) == 0x03));
	    ++ code_pos;
	    
	    const symbol_type lhs = vocab[id_lhs];
	    const phrase_type rule_target = read_phrase(pos_target, cache_targets, target_db);
	    
	    rule_ptr_type rule(new rule_type(lhs, rule_source, rule_target, rule_arity));
	    rule->sort_source_index();
	    
	    for (size_t feature = 0; feature < score_db.size(); ++ feature) {
	      const score_type score = score_db[feature][pos_feature];
	      
	      // ignore zero score...
	      if (score == 0.0) continue;
	      
	      // when zero, we will use inifinity...
	      rule->features[feature_names[feature]] = (score <= boost::numeric::bounds<score_type>::lowest()
							? - std::numeric_limits<feature_set_type::mapped_type>::infinity()
							: (score >= boost::numeric::bounds<score_type>::highest()
							   ? std::numeric_limits<feature_set_type::mapped_type>::infinity()
							   : feature_set_type::mapped_type(score)));
	    }
	    
	    ++ pos_feature;
	    
	    options.push_back(rule);
	  }
	}
      }
      
      return result.first->second;
    }

  private:    

    const phrase_type& read_phrase(size_type pos,
				   const cache_phrase_set_type& cache_phrases,
				   const phrase_db_type& phrase_db) const
    {
      const size_type cache_pos = hasher_type::operator()(pos) & (cache_phrases.size() - 1);
      
      cache_phrase_type& cache = const_cast<cache_phrase_type&>(cache_phrases[cache_pos]);
      if (cache.pos != pos) {
	typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
	typedef std::vector<byte_type, std::allocator<byte_type> >  code_set_type;
	typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > id_set_type;

	code_set_type codes(phrase_db[pos].begin(), phrase_db[pos].end());
	
	code_set_type::const_iterator hiter = codes.begin();
	code_set_type::const_iterator citer = codes.begin();
	code_set_type::const_iterator citer_end = codes.end();
	
	id_set_type phrase_id;
	word_type::id_type id;
	
	for (size_type code_pos = 0; citer != citer_end; ++ code_pos) {
	  const size_type offset = utils::group_aligned_decode(id, &(*hiter), code_pos);
	  
	  citer = hiter + offset;
	  hiter += offset & (- size_type((code_pos & 0x03) == 0x03));
	  
	  phrase_id.push_back(id);
	}
	
	sequence_type phrase(phrase_id.size());
	sequence_type::iterator piter = phrase.begin();
	id_set_type::const_iterator iiter_end = phrase_id.end();
	for (id_set_type::const_iterator iiter = phrase_id.begin(); iiter != iiter_end; ++ iiter, ++ piter)
	  *piter = vocab[*iiter];
	
	cache.pos = pos;
	cache.value = phrase_type(phrase.begin(), phrase.end());
      }
      return cache.value;
    }

  public:
    size_type feature_size() const { return score_db.size(); }
    bool empty() const { return score_db.empty(); }
    path_type path() const { return rule_db.path().parent_path(); }
    bool is_open() const { return ! score_db.empty(); }


    void quantize();
    void read(const std::string& parameter);
    void write(const path_type& path) const;

    void read_text(const std::string& path);
    void read_binary(const path_type& path);

  private:
    
    rule_db_type    rule_db;
    
    phrase_db_type  source_db;
    phrase_db_type  target_db;
    
    score_db_type   score_db;
    
    vocab_type      vocab;

    feature_name_set_type feature_names;
    
    // caching..
    cache_rule_map_type   cache_rule_sets;
    
    cache_phrase_set_type cache_sources;
    cache_phrase_set_type cache_targets;

  public:
    int max_span;
  };


  void GrammarStaticImpl::ScoreSet::read(const path_type& path)
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
  
  void GrammarStaticImpl::ScoreSet::write(const path_type& file) const
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
  

  void GrammarStaticImpl::quantize()
  {
    typedef score_type base_type;

    typedef std::map<base_type, size_type, std::less<base_type>, std::allocator<std::pair<const base_type, size_type> > > counts_type;
    typedef std::map<base_type, quantized_type, std::less<base_type>, std::allocator<std::pair<const base_type, quantized_type> > > codemap_type;
    typedef boost::array<base_type, 256> codebook_type;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    counts_type      counts;
    codebook_type    codebook;
    codemap_type     codemap;
    
    for (size_t feature = 0; feature < score_db.size(); ++ feature)
      if (score_db[feature].score.is_open()) {
	
	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.score.quantized.XXXXXX");
	utils::tempfile::insert(path);
	
	boost::iostreams::filtering_ostream os;
	os.push(utils::packed_sink<quantized_type, std::allocator<quantized_type> >(path));
	
	counts.clear();
	codemap.clear();
	std::fill(codebook.begin(), codebook.end(), 0.0);
	
	score_set_type::score_set_type::const_iterator liter_end = score_db[feature].score.end();
	for (score_set_type::score_set_type::const_iterator liter = score_db[feature].score.begin(); liter != liter_end; ++ liter)
	  ++ counts[*liter];
	
	Quantizer::quantize(counts, codebook, codemap);
	
	for (score_set_type::score_set_type::const_iterator liter = score_db[feature].score.begin(); liter != liter_end; ++ liter) {
	  codemap_type::const_iterator citer = codemap.find(*liter);
	  if (citer == codemap.end())
	    throw std::runtime_error("no codemap?");
	  
	  os.write((char*) &(citer->second), sizeof(quantized_type));
	}
	
	for (int i = 0; i < 256; ++ i)
	  score_db[feature].maps[i] = codebook[i];
	
	os.pop();
	utils::tempfile::permission(path);
	
	score_db[feature].quantized.open(path);
	score_db[feature].score.clear();
      }
  }
  
  void GrammarStaticImpl::read(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
	  
    const parameter_type param(parameter);
    
    const path_type path = param.name();
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no grammar file") + param.name());
    
    if (boost::filesystem::is_directory(path))
      read_binary(path);
    else
      read_text(parameter);
    
    parameter_type::const_iterator siter = param.find("max-span");
    if (siter != param.end())
      max_span = boost::lexical_cast<int>(siter->second);
  }
  
  void GrammarStaticImpl::write(const path_type& file) const
  {
        typedef utils::repository repository_type;
    
    if (file == path()) return;
    
    repository_type rep(file, repository_type::write);
    
    rule_db.write(rep.path("rule"));
    source_db.write(rep.path("source"));
    target_db.write(rep.path("target"));
    
    vocab.write(rep.path("vocab"));
    
    const size_type feature_size = score_db.size();
    for (size_t feature = 0; feature < feature_size; ++ feature) {
      std::ostringstream stream_score;
      stream_score << "score-" << std::setfill('0') << std::setw(6) << feature;
      
      score_db[feature].write(rep.path(stream_score.str()));

      const std::string name(std::string("feature") + boost::lexical_cast<std::string>(feature));
      
      rep[name] = feature_names[feature];
    }
    
    rep["feature-size"] = boost::lexical_cast<std::string>(feature_size);
  }
  
  void GrammarStaticImpl::read_binary(const path_type& path)
  {
    typedef utils::repository repository_type;
	
    repository_type rep(path, repository_type::read);
    
    rule_db.open(rep.path("rule"));
    source_db.open(rep.path("source"));
    target_db.open(rep.path("target"));
    
    vocab.open(rep.path("vocab"));
    
    repository_type::const_iterator iter = rep.find("feature-size");
    if (iter == rep.end())
      throw std::runtime_error("no feature size?");

    const size_type feature_size = atoi(iter->second.c_str());
    
    feature_names.reserve(feature_size);
    feature_names.resize(feature_size);
    score_db.reserve(feature_size);
    score_db.resize(feature_size);
    
    for (size_t feature = 0; feature < feature_size; ++ feature) {
      std::ostringstream stream_score;
      stream_score << "score-" << std::setfill('0') << std::setw(6) << feature;
      
      score_db[feature].read(rep.path(stream_score.str()));
      
      const std::string name(std::string("feature") + boost::lexical_cast<std::string>(feature));
      repository_type::const_iterator iter = rep.find(name);
      if (iter == rep.end())
	throw std::runtime_error(std::string("no feature name?: ") + name);
      
      feature_names[feature] = iter->second;
    }
    
#if 0
    std::cerr << "phrase-table: " << path
		<< " feature size: " << feature_size
		<< std::endl;
#endif
  }
  
  
  typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
  typedef std::vector<float, std::allocator<float> > scores_parsed_type;
  typedef boost::fusion::tuple<std::string, phrase_parsed_type, phrase_parsed_type, scores_parsed_type > rule_parsed_type;

  template <typename Iterator>
  struct rule_grammar_parser_static : boost::spirit::qi::grammar<Iterator, rule_parsed_type(), boost::spirit::standard::space_type>
  {
    
    rule_grammar_parser_static() : rule_grammar_parser_static::base_type(rule_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
    
      using qi::phrase_parse;
      using qi::lexeme;
      using qi::attr;
      using qi::hold;
      using standard::char_;
      using qi::float_; // FLOAT!
      using qi::_1;
      using standard::space;
      
      lhs %= (lexeme[char_('[') >> +(char_ - space - ']') >> char_(']')]);
      phrase %= *(lexeme[+(char_ - space) - "|||"]);
      rule_grammar %= (hold[lhs >> "|||"] | attr("")) >> phrase >> "|||" >> phrase >> "|||" >> (+float_);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> lhs;
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), boost::spirit::standard::space_type> rule_grammar;
  };

  struct ScoreSetStream
  {
    typedef boost::filesystem::path path_type;
    typedef utils::compress_ostream ostream_type;
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    
    ScoreSetStream() : ostream() {}
    ScoreSetStream(const path_type& _path)
      : ostream(new ostream_type(_path, 1024 * 1024)), path(_path) {}
    
    ostream_ptr_type ostream;
    path_type        path;
  };


  template <typename Phrase, typename Codes>
  inline
  void encode_phrase(const Phrase& phrase, Codes& codes)
  {
    typedef size_t size_type;
    
    codes.clear();
    codes.resize(phrase.size() * 8, 0);
    
    typename Codes::iterator hiter = codes.begin();
    typename Codes::iterator citer = codes.begin();
    
    typename Phrase::const_iterator piter_begin = phrase.begin();
    typename Phrase::const_iterator piter_end = phrase.end();
    
    for (typename Phrase::const_iterator piter = piter_begin; piter != piter_end; ++ piter) {
      const size_type offset = utils::group_aligned_encode(piter->id(), &(*hiter), piter - piter_begin);
      
      citer = hiter + offset;
      hiter += offset & (- size_type(((piter - piter_begin) & 0x03) == 0x03));
    }
    codes.resize(citer - codes.begin());
  }

  template <typename Iterator, typename Buffer>
  inline
  void encode_index(Iterator first, Iterator last, Buffer& buffer)
  {
    buffer.clear();
    for (/**/; first != last; ++ first)
      buffer.push_back(first->non_terminal().id());
  }
  
  void GrammarStaticImpl::read_text(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const path_type path = param.name();

    typedef succinctdb::succinct_hash<byte_type, std::allocator<byte_type> > phrase_db_type;
    
    typedef ScoreSetStream score_stream_type;
    typedef std::vector<score_stream_type, std::allocator<score_stream_type> > score_stream_set_type;
    
    // feature-id, lhs-id, target-id
    typedef boost::tuple<id_type, symbol_type::id_type, id_type> rule_option_type;
    typedef std::vector<rule_option_type, std::allocator<rule_option_type> > rule_option_set_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
    
    typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > id_set_type;
    
    typedef std::vector<byte_type, std::allocator<byte_type> >  code_set_type;
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no file? ") + path.file_string());
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    const path_type path_rule   = utils::tempfile::directory_name(tmp_dir / "cicada.rule.XXXXXX");
    const path_type path_source = utils::tempfile::directory_name(tmp_dir / "cicada.source.XXXXXX");
    const path_type path_target = utils::tempfile::directory_name(tmp_dir / "cicada.target.XXXXXX");
    const path_type path_vocab  = utils::tempfile::directory_name(tmp_dir / "cicada.vocab.XXXXXX");
    
    utils::tempfile::insert(path_rule);
    utils::tempfile::insert(path_source);
    utils::tempfile::insert(path_target);
    utils::tempfile::insert(path_vocab);
    
    rule_db.open(path_rule, rule_db_type::WRITE);
    phrase_db_type        sources_db(1024 * 1024 * 4);
    phrase_db_type        targets_db(1024 * 1024 * 4);
    
    score_stream_set_type score_streams;
    
    id_type id_rule = 0;
    
    int feature_size = -1;
    
    sequence_type source_prev;
    sequence_type source;
    sequence_type target;
    rule_parsed_type rule;
    
    id_set_type source_index;
    
    code_set_type codes_source;
    code_set_type codes_target;
    code_set_type codes_option;
    
    rule_option_set_type rule_options;
    
    const std::string sep("|||");
    
    utils::compress_istream is(path, 1024 * 1024);
    
    std::string line;
    
    // we will construct this parser everytimt...
    rule_grammar_parser_static<std::string::const_iterator> rule_parser;

    size_type arity_source = 0;

    while (std::getline(is, line)) {
      boost::fusion::get<0>(rule).clear();
      boost::fusion::get<1>(rule).clear();
      boost::fusion::get<2>(rule).clear();
      boost::fusion::get<3>(rule).clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);
      
      if (! result || iter != iter_end) continue;
      
      source.clear();
      source.insert(source.end(), boost::fusion::get<1>(rule).begin(), boost::fusion::get<1>(rule).end());
      
      target.clear();
      target.insert(target.end(), boost::fusion::get<2>(rule).begin(), boost::fusion::get<2>(rule).end());
      
      if (source != source_prev) {
	
	if (! rule_options.empty()) {
	  // encode options...
	   
	  {
	    codes_option.clear();
	    codes_option.resize(rule_options.size() * 16 + 16, 0);
	     
	    code_set_type::iterator hiter = codes_option.begin();
	    code_set_type::iterator citer = codes_option.begin();
	    size_type pos = 0;

	    const id_type id_feature = boost::get<0>(rule_options.front());

	    encode_phrase(source_prev, codes_source);
	     
	    const id_type id_source = sources_db.insert(&(*codes_source.begin()), codes_source.size(),
							hasher_type::operator()(codes_source.begin(), codes_source.end(), 0));
	     
	    const size_type offset_feature = utils::group_aligned_encode(id_feature, &(*hiter), pos);
	    citer = hiter + offset_feature;
	    hiter += offset_feature & (- size_type((pos & 0x03) == 0x03));
	    ++ pos;
	     
	    const size_type offset_source = utils::group_aligned_encode(id_source, &(*hiter), pos);
	    citer = hiter + offset_source;
	    hiter += offset_source & (- size_type((pos & 0x03) == 0x03));
	    ++ pos;

	    rule_option_set_type::const_iterator piter_end = rule_options.end();
	    for (rule_option_set_type::const_iterator piter = rule_options.begin(); piter != piter_end; ++ piter) {

	      const symbol_type::id_type id_lhs = boost::get<1>(*piter);
	      const id_type              id_target = boost::get<2>(*piter);

	      const size_type offset_lhs = utils::group_aligned_encode(id_lhs, &(*hiter), pos);
	      citer = hiter + offset_lhs;
	      hiter += offset_lhs & (- size_type((pos & 0x03) == 0x03));
	      ++ pos;
	       
	      const size_type offset_target = utils::group_aligned_encode(id_target, &(*hiter), pos);
	      citer = hiter + offset_target;
	      hiter += offset_target & (- size_type((pos & 0x03) == 0x03));
	      ++ pos;
	    }
	     
	    codes_option.resize(citer - codes_option.begin());
	  }
	   
	  // encode source.....
	  encode_index(source_prev.begin(), source_prev.en(), source_index);

	  // insert...
	  rule_db.insert(&(*source_index.begin()), source_index.size(), &(*codes_option.begin()), codes_option.size());
	}

	rule_options.clear();
	source_prev = source;

	arity_source = 0;
	for (sequence_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter)
	  arity_source += siter->is_non_terminal();
      }
       
      size_type arity_target = 0;
      for (sequence_type::const_iterator titer = target.begin(); titer != target.end(); ++ titer)
	arity_target += titer->is_non_terminal();
       
      if (arity_source != arity_target)
	throw std::runtime_error("# of non-terminals do not match...");
       
      // lhs...
      const std::string& lhs = boost::fusion::get<0>(rule);
      const word_type::id_type id_lhs = word_type(lhs.empty() ? vocab_type::X : word_type(lhs)).id();

#if 0
      std::cerr << "lhs: " << lhs;
      std::cerr << " source: ";
      std::copy(source_prev.begin(), source_prev.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << "target: ";
      std::copy(target.begin(), target.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << "features: ";
      std::copy(boost::fusion::get<3>(rule).begin(), boost::fusion::get<3>(rule).end(), std::ostream_iterator<score_type>(std::cerr, " "));
      std::cerr << std::endl;
#endif
      
      // scores...
      if (feature_size < 0) {
	feature_size = boost::fusion::get<3>(rule).size();
	
	score_streams.reserve(feature_size);
	score_streams.resize(feature_size);
	
	for (int feature = 0; feature < feature_size; ++ feature) {
	  score_streams[feature].path = utils::tempfile::file_name(tmp_dir / "cicada.feature.XXXXXX");
	  utils::tempfile::insert(score_streams[feature].path);
	  
	  score_streams[feature].ostream.reset(new utils::compress_ostream(score_streams[feature].path, 1024 * 1024));
	}
      } else if (feature_size != static_cast<int>(boost::fusion::get<3>(rule).size()))
	throw std::runtime_error("invalid # of features...");
      
      for (int feature = 0; feature < feature_size; ++ feature)
	score_streams[feature].ostream->write((char*) &boost::fusion::get<3>(rule)[feature], sizeof(score_type));
       
      // encode target...
      encode_phrase(target, codes_target);
       
      const id_type id_target = targets_db.insert(&(*codes_target.begin()), codes_target.size(),
						  hasher_type::operator()(codes_target.begin(), codes_target.end(), 0));
       
      // put into rule_options...
      rule_options.push_back(boost::make_tuple(id_rule ++, id_lhs, id_target));
    }
     
    if (! rule_options.empty()) {
      // encode options...
	   
      {
	codes_option.clear();
	codes_option.resize(rule_options.size() * 16 + 16, 0);
	     
	code_set_type::iterator hiter = codes_option.begin();
	code_set_type::iterator citer = codes_option.begin();
	size_type pos = 0;

	const id_type id_feature = boost::get<0>(rule_options.front());

	encode_phrase(source_prev, codes_source);
	     
	const id_type id_source = sources_db.insert(&(*codes_source.begin()), codes_source.size(),
						    hasher_type::operator()(codes_source.begin(), codes_source.end(), 0));
	     
	const size_type offset_feature = utils::group_aligned_encode(id_feature, &(*hiter), pos);
	citer = hiter + offset_feature;
	hiter += offset_feature & (- size_type((pos & 0x03) == 0x03));
	++ pos;
	     
	const size_type offset_source = utils::group_aligned_encode(id_source, &(*hiter), pos);
	citer = hiter + offset_source;
	hiter += offset_source & (- size_type((pos & 0x03) == 0x03));
	++ pos;

	rule_option_set_type::const_iterator piter_end = rule_options.end();
	for (rule_option_set_type::const_iterator piter = rule_options.begin(); piter != piter_end; ++ piter) {

	  const symbol_type::id_type id_lhs = boost::get<1>(*piter);
	  const id_type              id_target = boost::get<2>(*piter);

	  const size_type offset_lhs = utils::group_aligned_encode(id_lhs, &(*hiter), pos);
	  citer = hiter + offset_lhs;
	  hiter += offset_lhs & (- size_type((pos & 0x03) == 0x03));
	  ++ pos;
	       
	  const size_type offset_target = utils::group_aligned_encode(id_target, &(*hiter), pos);
	  citer = hiter + offset_target;
	  hiter += offset_target & (- size_type((pos & 0x03) == 0x03));
	  ++ pos;
	}
	     
	codes_option.resize(citer - codes_option.begin());
      }
	   
      // encode source.. we will use index-stripped indexing!
      encode_index(source_prev.begin(), source_prev.en(), source_index);

      // insert...
      rule_db.insert(&(*source_index.begin()), source_index.size(), &(*codes_option.begin()), codes_option.size());
    }

    // source phrases...
    sources_db.write(path_source);
    sources_db.clear();
    source_db.open(path_source);
    
    // target phrases...
    targets_db.write(path_target);
    targets_db.clear();
    target_db.open(path_target);

    // rules....
    rule_db.close();
    rule_db.open(path_rule);
     
    // vocabulary...
    word_type::write(path_vocab);
    vocab.open(path_vocab);
        
    // scores...
    score_db.reserve(feature_size);
    score_db.resize(feature_size);
    
    feature_names.clear();
    feature_names.reserve(feature_size);
    feature_names.resize(feature_size, feature_type());
     
    for (int feature = 0; feature < feature_size; ++ feature) {
      score_streams[feature].ostream->reset();
      utils::tempfile::permission(score_streams[feature].path);
      score_db[feature].score.open(score_streams[feature].path);

      const std::string name(std::string("feature") + boost::lexical_cast<std::string>(feature));

      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	feature_names[feature] = feature_type(piter->second);
      
      // default name...!
      if (feature_names[feature] == feature_type())
	feature_names[feature] = std::string("rule-table-") + boost::lexical_cast<std::string>(feature);
    }
  }
  
  
  GrammarStatic::GrammarStatic(const std::string& parameter)
    : pimpl(new impl_type(parameter)) {}

  GrammarStatic::~GrammarStatic() { std::auto_ptr<impl_type> tmp(pimpl); }

  GrammarStatic::GrammarStatic(const GrammarStatic& x)
    : pimpl(new impl_type(*x.pimpl)) {}

  GrammarStatic& GrammarStatic::operator=(const GrammarStatic& x)
  {
    *pimpl = *x.pimpl;
    return *this;
  }

  GrammarStatic::transducer_ptr_type GrammarStatic::clone() const
  {
    return transducer_ptr_type(new GrammarStatic(*this));
  }

  bool GrammarStatic::valid_span(int first, int last, int distance) const
  {
    return distance <= pimpl->max_span;
  }
  
  GrammarStatic::id_type GrammarStatic::root() const
  {
    return 0;
  }
  
  GrammarStatic::id_type GrammarStatic::next(const id_type& node, const symbol_type& symbol) const
  {
    const impl_type::size_type pos = pimpl->find(symbol.non_terminal(), node);
    
    return (pimpl->is_valid(pos) ? id_type(pos) : id_type(0));
  }
  
  bool GrammarStatic::has_next(const id_type& node) const
  {
    return (pimpl->is_valid(node) && pimpl->has_children(node));
  }
  
  const GrammarStatic::rule_set_type& GrammarStatic::rules(const id_type& node) const
  {
    static const rule_set_type __empty;
    
    return (pimpl->is_valid(node) && pimpl->exists(node) ? pimpl->read_rule_set(node) : __empty);
  }
  
  void GrammarStatic::quantize()
  {
    pimpl->quantize();
  }
  
  void GrammarStatic::write(const path_type& path)
  {
    pimpl->write(path);
  }
};
