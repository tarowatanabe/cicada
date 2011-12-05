//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "grammar_static.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "feature_vector_codec.hpp"
#include "attribute_vector_codec.hpp"

#include "succinct_db/succinct_trie_database.hpp"
#include "succinct_db/succinct_hash.hpp"

#include "utils/repository.hpp"
#include "utils/compress_stream.hpp"
#include "utils/tempfile.hpp"
#include "utils/group_aligned_code.hpp"
#include "utils/simple_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/arc_list.hpp"
#include "utils/packed_device.hpp"
#include "utils/packed_vector.hpp"
#include "utils/vertical_coded_device.hpp"
#include "utils/vertical_coded_vector.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/bithack.hpp"
#include "utils/byte_aligned_code.hpp"
#include "utils/json_string_parser.hpp"

#include <boost/lexical_cast.hpp>

#include <boost/filesystem.hpp>

#include <boost/spirit/include/qi.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/thread.hpp>

#include <google/dense_hash_map>

namespace cicada
{
  
  // TODO: we need to index by index-stripped non-terminals!
  
  struct GrammarStaticImpl : public utils::hashmurmur<uint64_t>
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;

    typedef uint64_t off_type;
    
    typedef cicada::Symbol    symbol_type;
    typedef cicada::Symbol    word_type;
    typedef cicada::Feature   feature_type;
    typedef cicada::Attribute attribute_type;
    typedef cicada::Vocab     vocab_type;
    
    typedef Transducer::rule_type          rule_type;
    typedef Transducer::rule_ptr_type      rule_ptr_type;
    typedef Transducer::rule_pair_type     rule_pair_type;
    typedef Transducer::rule_pair_set_type rule_pair_set_type;
    
    typedef Transducer::feature_set_type feature_set_type;
    
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

    template <typename Key>
    struct KeyVocab
    {
      typedef typename Key::id_type id_type;
      typedef succinctdb::succinct_hash_mapped<byte_type, std::allocator<byte_type> > data_type;
      
      KeyVocab() : data() {}
      KeyVocab(const path_type& path) : data(path) {}
      
      void write(const path_type& path) const { data.write(path); }
      void read(const path_type& path) { data.open(path); }
      void open(const path_type& path) { data.open(path); }
      void clear() { data.clear(); }
      
      bool empty() const { return data.empty(); }
      path_type path() const { return data.path(); }
      
      Key operator[](const id_type& id) const
      {
	return Key(data[id].begin(), data[id].end());
      }
      
      data_type data;
    };
    
    struct PackedData
    {
      typedef uint64_t off_type;
      typedef utils::map_file<byte_type, std::allocator<byte_type> >           data_type;
      typedef utils::vertical_coded_vector_mapped<off_type, std::allocator<off_type> > offset_type;
      
      PackedData() : data(), offset() {}
      PackedData(const path_type& path) : data(), offset() { open(path); }
      
      bool empty() const { return data.empty(); }
      path_type path() const { return data.path().parent_path(); }
      void clear()
      {
	data.clear();
	offset.clear();
      }
      
      utils::piece operator[](size_t i) const
      {
	if (i == 0)
	  return utils::piece(data.begin(), data.begin() + offset[i]);
	else
	  return utils::piece(data.begin() + offset[i - 1], data.begin() + offset[i]);
      }

      void open(const path_type& path) { read(path); }
      
      void read(const path_type& path)
      {
	typedef utils::repository repository_type;
	
	repository_type rep(path, repository_type::read);

	data.open(rep.path("data"));
	offset.open(rep.path("offset"));
      }

      void write(const path_type& file) const
      {
	typedef utils::repository repository_type;
	
	if (path() == file) return;
	
	repository_type rep(file, repository_type::write);
	
	data.write(rep.path("data"));
	offset.write(rep.path("offset"));
      }
      
    private:
      data_type   data;
      offset_type offset;
    };

    typedef KeyVocab<feature_type>   feature_vocab_type;
    typedef KeyVocab<attribute_type> attribute_vocab_type;    
    typedef PackedData               feature_data_type;
    typedef PackedData               attribute_data_type;

    typedef std::vector<feature_type, std::allocator<feature_type> >     feature_name_set_type;
    typedef std::vector<attribute_type, std::allocator<attribute_type> > attribute_name_set_type;

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
    typedef utils::arc_list<size_type, rule_pair_set_type, 8,
			    std::equal_to<size_type>,
			    std::allocator<std::pair<size_type, rule_pair_set_type> > > cache_rule_set_type;

    struct cache_phrase_type
    {
      rule_ptr_type rule;
      size_type     pos;
      cache_phrase_type() : rule(), pos(size_type(-1)) {}
    };
    
    typedef utils::array_power2<cache_rule_set_type, 1024 * 4, std::allocator<cache_rule_set_type> > cache_rule_map_type;
    typedef utils::array_power2<cache_phrase_type,   1024 * 2, std::allocator<cache_phrase_type> >   cache_phrase_set_type;

    typedef std::pair<word_type, size_type> word_node_type;

    typedef google::dense_hash_map<word_node_type, size_type, utils::hashmurmur<size_t>, std::equal_to<word_node_type> > cache_node_type;

  public:
    GrammarStaticImpl(const std::string& parameter)
      : max_span(15),
	caching(false)
    {
      cache_node.set_empty_key(word_node_type());
      read(parameter);
    }

    GrammarStaticImpl(const GrammarStaticImpl& x)
      : rule_db(x.rule_db),
	source_db(x.source_db),
	target_db(x.target_db),
	score_db(x.score_db),
	attr_db(x.attr_db),
	feature_data(x.feature_data),
	attribute_data(x.attribute_data),
	feature_vocab(x.feature_vocab),
	attribute_vocab(x.attribute_vocab),
	vocab(x.vocab),
	feature_names(x.feature_names),
	attribute_names(x.attribute_names),
	max_span(x.max_span),
	caching(x.caching)
    {
      cache_node.set_empty_key(word_node_type());
    }

    GrammarStaticImpl& operator=(const GrammarStaticImpl& x)
    {
      clear();
      
      rule_db         = x.rule_db;
      source_db       = x.source_db;
      target_db       = x.target_db;
      score_db        = x.score_db;
      attr_db         = x.attr_db;
      feature_data   = x.feature_data;
      attribute_data = x.attribute_data;
      feature_vocab   = x.feature_vocab;
      attribute_vocab = x.attribute_vocab;
      vocab           = x.vocab;
      feature_names   = x.feature_names;
      attribute_names = x.attribute_names;
      max_span        = x.max_span;
      caching         = x.caching;
      
      return *this;
    }
    
  public:

    void clear()
    {
      rule_db.clear();
      source_db.clear();
      target_db.clear();
      score_db.clear();
      attr_db.clear();
      feature_data.clear();
      attribute_data.clear();
      feature_vocab.clear();
      attribute_vocab.clear();
      vocab.clear();
      feature_names.clear();
      attribute_names.clear();

      cache_rule_sets.clear();
      cache_sources.clear();
      cache_targets.clear();

      cache_node.clear();

      max_span = 15;
      caching = false;
    }
    
    size_type find(const word_type& word, size_type node) const
    {
      if (caching && word.is_non_terminal()) {
	cache_node_type& __cache_node = const_cast<cache_node_type&>(cache_node);
	
	std::pair<cache_node_type::iterator, bool> result = __cache_node.insert(std::make_pair(std::make_pair(word, node), 0));
	if (result.second) {
	  const word_type::id_type id = vocab[word];
	  result.first->second = rule_db.find(&id, 1, node);
	}
	return result.first->second;
      } else {
	const word_type::id_type id = vocab[word];
	return rule_db.find(&id, 1, node);
      }
    }
    
    // valid implies that you can continue searching from node...
    bool is_valid(size_type node) const { return rule_db.is_valid(node); }
    bool has_children(size_type node) const { return rule_db.has_children(node); }
    
    // exists implies data associated with the node exists...
    bool exists(size_type node) const { return rule_db.exists(node); }
    
    const rule_pair_set_type& read_rule_set(size_type node) const
    {
      FeatureVectorCODEC   feature_codec;
      AttributeVectorCODEC attribute_codec;
      
      const size_type cache_pos = hasher_type::operator()(node) & (cache_rule_sets.size() - 1);
      
      cache_rule_set_type& cache = const_cast<cache_rule_set_type&>(cache_rule_sets[cache_pos]);
      
      std::pair<cache_rule_set_type::iterator, bool> result = cache.find(node);
      if (! result.second) {
	typedef utils::piece code_set_type;
	
	rule_pair_set_type& options = result.first->second;
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
	    
	    const rule_ptr_type rule_source = read_phrase(lhs, pos_source, cache_sources, source_db);
	    const rule_ptr_type rule_target = read_phrase(lhs, pos_target, cache_targets, target_db);
	    
	    if (rule_target->rhs.empty())
	      options.push_back(rule_pair_type(rule_source, rule_target));
	    else {
	      rule_type rule_sorted_source(*rule_source);
	      rule_type rule_sorted_target(*rule_target);
	      
	      cicada::sort(rule_sorted_source, rule_sorted_target);
	      
	      options.push_back(rule_pair_type(rule_type::create(rule_sorted_source), rule_type::create(rule_sorted_target)));
	    }
	    
	    if (! feature_data.empty())
	      feature_codec.decode(feature_data[pos_feature].begin(),
				   feature_data[pos_feature].end(),
				   feature_vocab,
				   options.back().features);
	    
	    if (! attribute_data.empty())
	      attribute_codec.decode(attribute_data[pos_feature].begin(),
				     attribute_data[pos_feature].end(),
				     attribute_vocab,
				     options.back().attributes);

	    
	    for (size_t feature = 0; feature < score_db.size(); ++ feature) {
	      const score_type score = score_db[feature][pos_feature];
	      
	      // ignore zero score...
	      if (score == 0.0) continue;
	      
	      // when zero, we will use inifinity...
	      options.back().features[feature_names[feature]] = (score <= boost::numeric::bounds<score_type>::lowest()
								 ? - std::numeric_limits<feature_set_type::mapped_type>::infinity()
								 : (score >= boost::numeric::bounds<score_type>::highest()
								    ? std::numeric_limits<feature_set_type::mapped_type>::infinity()
								    : feature_set_type::mapped_type(score)));
	    }

	    for (size_t attr = 0; attr < attr_db.size(); ++ attr) {
	      const score_type score = attr_db[attr][pos_feature];
	      
	      options.back().attributes[attribute_names[attr]] = (score <= boost::numeric::bounds<score_type>::lowest()
								  ? - std::numeric_limits<feature_set_type::mapped_type>::infinity()
								  : (score >= boost::numeric::bounds<score_type>::highest()
								     ? std::numeric_limits<feature_set_type::mapped_type>::infinity()
								     : feature_set_type::mapped_type(score)));
	    }
	    
	    ++ pos_feature;
	  }
	}
      }
      
      return result.first->second;
    }

  private:    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
    typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > id_set_type;

    sequence_type phrase_impl;
    id_set_type   phrase_id_impl;

    const rule_ptr_type& read_phrase(const symbol_type& lhs,
				     size_type pos,
				     const cache_phrase_set_type& cache_phrases,
				     const phrase_db_type& phrase_db) const
    {
      const size_type cache_pos = hasher_type::operator()(pos, lhs.id()) & (cache_phrases.size() - 1);
      
      cache_phrase_type& cache = const_cast<cache_phrase_type&>(cache_phrases[cache_pos]);
      if (cache.pos != pos || ! cache.rule || cache.rule->lhs != lhs) {
	typedef utils::piece code_set_type;
	
	code_set_type codes(phrase_db[pos].begin(), phrase_db[pos].end());
	
	code_set_type::const_iterator hiter = codes.begin();
	code_set_type::const_iterator citer = codes.begin();
	code_set_type::const_iterator citer_end = codes.end();
	
	id_set_type& phrase_id = const_cast<id_set_type&>(phrase_id_impl);
	phrase_id.clear();
	
	for (size_type code_pos = 0; citer != citer_end; ++ code_pos) {
	  word_type::id_type id;
	  const size_type offset = utils::group_aligned_decode(id, &(*hiter), code_pos);
	  
	  citer = hiter + offset;
	  hiter += offset & (- size_type((code_pos & 0x03) == 0x03));
	  
	  phrase_id.push_back(id);
	}
	
	sequence_type& phrase = const_cast<sequence_type&>(phrase_impl);
	phrase.resize(phrase_id.size());
	
	sequence_type::iterator piter = phrase.begin();
	id_set_type::const_iterator iiter_end = phrase_id.end();
	for (id_set_type::const_iterator iiter = phrase_id.begin(); iiter != iiter_end; ++ iiter, ++ piter)
	  *piter = vocab[*iiter];
	
	cache.pos = pos;
	cache.rule = rule_type::create(rule_type(lhs, phrase.begin(), phrase.end()));
      }
      return cache.rule;
    }

  public:
    size_type feature_size() const { return score_db.size(); }
    bool empty() const { return score_db.empty(); }
    path_type path() const { return rule_db.path().parent_path(); }
    bool is_open() const { return ! score_db.empty(); }


    void quantize();
    void read(const std::string& parameter);
    void write(const path_type& path) const;

    void read_keyed_text(const std::string& path);
    void read_text(const std::string& path);
    void read_binary(const std::string& path);

  private:
    
    rule_db_type    rule_db;
    
    phrase_db_type  source_db;
    phrase_db_type  target_db;
    
    score_db_type   score_db;
    score_db_type   attr_db;

    feature_data_type   feature_data;
    attribute_data_type attribute_data;
    
    feature_vocab_type   feature_vocab;
    attribute_vocab_type attribute_vocab;
    
    vocab_type      vocab;
    
    
    feature_name_set_type   feature_names;
    attribute_name_set_type attribute_names;
    
    // caching..
    cache_rule_map_type   cache_rule_sets;
    
    cache_phrase_set_type cache_sources;
    cache_phrase_set_type cache_targets;

    cache_node_type cache_node;

  public:
    int max_span;
    bool caching;
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
	throw std::runtime_error(std::string("no map file? ") + score_map_file.string());
      
      std::ifstream is(score_map_file.string().c_str());
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
      
      std::ofstream os(rep.path("score-map").string().c_str());
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
	::sync();
	
	while (! score_set_type::quantized_set_type::exists(path))
	  boost::thread::yield();

	utils::tempfile::permission(path);
	
	score_db[feature].quantized.open(path);
	score_db[feature].score.clear();
      }

    for (size_t attr = 0; attr < attr_db.size(); ++ attr)
      if (attr_db[attr].score.is_open()) {
	
	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.attr.quantized.XXXXXX");
	utils::tempfile::insert(path);
	
	boost::iostreams::filtering_ostream os;
	os.push(utils::packed_sink<quantized_type, std::allocator<quantized_type> >(path));
	
	counts.clear();
	codemap.clear();
	std::fill(codebook.begin(), codebook.end(), 0.0);
	
	score_set_type::score_set_type::const_iterator liter_end = attr_db[attr].score.end();
	for (score_set_type::score_set_type::const_iterator liter = attr_db[attr].score.begin(); liter != liter_end; ++ liter)
	  ++ counts[*liter];
	
	Quantizer::quantize(counts, codebook, codemap);
	
	for (score_set_type::score_set_type::const_iterator liter = attr_db[attr].score.begin(); liter != liter_end; ++ liter) {
	  codemap_type::const_iterator citer = codemap.find(*liter);
	  if (citer == codemap.end())
	    throw std::runtime_error("no codemap?");
	  
	  os.write((char*) &(citer->second), sizeof(quantized_type));
	}
	
	for (int i = 0; i < 256; ++ i)
	  attr_db[attr].maps[i] = codebook[i];
	
	os.pop();
	::sync();
	
	while (! score_set_type::quantized_set_type::exists(path))
	  boost::thread::yield();
	
	utils::tempfile::permission(path);
	
	attr_db[attr].quantized.open(path);
	attr_db[attr].score.clear();
      }
  }
  
  void GrammarStaticImpl::read(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
	  
    const parameter_type param(parameter);
    
    const path_type path = param.name();
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no grammar file") + param.name());
    
    bool key_value = false;
    parameter_type::const_iterator kiter = param.find("key-value");
    if (kiter != param.end())
      key_value = utils::lexical_cast<bool>(kiter->second);

    if (boost::filesystem::is_directory(path))
      read_binary(parameter);
    else if (key_value)
      read_keyed_text(parameter);
    else
      read_text(parameter);
    
    parameter_type::const_iterator siter = param.find("max-span");
    if (siter != param.end())
      max_span = utils::lexical_cast<int>(siter->second);

    parameter_type::const_iterator citer = param.find("cache");
    if (citer != param.end())
      caching = utils::lexical_cast<bool>(citer->second);
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
    
    if (! feature_data.empty())
      feature_data.write(rep.path("feature-data"));
    if (! feature_vocab.empty())
      feature_vocab.write(rep.path("feature-vocab"));
    
    if (! attribute_data.empty())
      attribute_data.write(rep.path("attribute-data"));
    if (! attribute_vocab.empty())
      attribute_vocab.write(rep.path("attribute-vocab"));
    
    const size_type feature_size = score_db.size();
    const size_type attribute_size = attr_db.size();
    
    for (size_t feature = 0; feature < feature_size; ++ feature) {
      std::ostringstream stream_score;
      stream_score << "score-" << std::setfill('0') << std::setw(6) << feature;
      
      score_db[feature].write(rep.path(stream_score.str()));

      const std::string name(std::string("feature") + utils::lexical_cast<std::string>(feature));
      
      rep[name] = feature_names[feature];
    }
    
    for (size_t attribute = 0; attribute < attribute_size; ++ attribute) {
      std::ostringstream stream_score;
      stream_score << "attribute-" << std::setfill('0') << std::setw(6) << attribute;
      
      attr_db[attribute].write(rep.path(stream_score.str()));

      const std::string name(std::string("attribute") + utils::lexical_cast<std::string>(attribute));
      
      rep[name] = attribute_names[attribute];
    }
    
    rep["feature-size"]   = utils::lexical_cast<std::string>(feature_size);
    rep["attribute-size"] = utils::lexical_cast<std::string>(attribute_size);
  }
  
  void GrammarStaticImpl::read_binary(const std::string& parameter)
  {
    typedef utils::repository repository_type;
    typedef cicada::Parameter parameter_type;
	  
    const parameter_type param(parameter);
    
    const path_type path = param.name();
    repository_type rep(path, repository_type::read);
    
    rule_db.open(rep.path("rule"));
    source_db.open(rep.path("source"));
    target_db.open(rep.path("target"));
    
    vocab.open(rep.path("vocab"));
    
    if (boost::filesystem::exists(rep.path("feature-data")))
      feature_data.open(rep.path("feature-data"));
    if (boost::filesystem::exists(rep.path("feature-vocab")))
      feature_vocab.open(rep.path("feature-vocab"));

    if (boost::filesystem::exists(rep.path("attribute-data")))
      attribute_data.open(rep.path("attribute-data"));
    if (boost::filesystem::exists(rep.path("attribute-vocab")))
      attribute_vocab.open(rep.path("attribute-vocab"));
    
    repository_type::const_iterator iter = rep.find("feature-size");
    if (iter == rep.end())
      throw std::runtime_error("no feature size?");

    const size_type feature_size = utils::lexical_cast<size_type>(iter->second);
    
    feature_names.reserve(feature_size);
    feature_names.resize(feature_size);
    score_db.reserve(feature_size);
    score_db.resize(feature_size);
    
    for (size_t feature = 0; feature < feature_size; ++ feature) {
      std::ostringstream stream_score;
      stream_score << "score-" << std::setfill('0') << std::setw(6) << feature;
      
      score_db[feature].read(rep.path(stream_score.str()));
      
      const std::string name(std::string("feature") + utils::lexical_cast<std::string>(feature));

      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	feature_names[feature] = piter->second;
      else {
	repository_type::const_iterator iter = rep.find(name);
	if (iter != rep.end())
	  feature_names[feature] = iter->second;
	else
	  throw std::runtime_error(std::string("no feature name?: ") + name);
      }
    }
    
    repository_type::const_iterator aiter = rep.find("attribute-size");
    if (aiter == rep.end()) return;
    
    const size_type attribute_size = utils::lexical_cast<size_type>(aiter->second);
    
    attribute_names.reserve(attribute_size);
    attribute_names.resize(attribute_size);
    attr_db.reserve(attribute_size);
    attr_db.resize(attribute_size);
    
    for (size_t attribute = 0; attribute < attribute_size; ++ attribute) {
      std::ostringstream stream_score;
      stream_score << "attribute-" << std::setfill('0') << std::setw(6) << attribute;
      
      attr_db[attribute].read(rep.path(stream_score.str()));
      
      const std::string name(std::string("attribute") + utils::lexical_cast<std::string>(attribute));
      
      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	attribute_names[attribute] = piter->second;
      else {
	repository_type::const_iterator iter = rep.find(name);
	if (iter != rep.end())
	  attribute_names[attribute] = iter->second;
	else
	  throw std::runtime_error(std::string("no attribute name?: ") + name);
      }
    }
  }
  
  typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
  typedef std::vector<float, std::allocator<float> > scores_parsed_type;
  typedef boost::fusion::tuple<std::string, phrase_parsed_type, phrase_parsed_type, scores_parsed_type, scores_parsed_type> rule_parsed_type;

  template <typename Iterator>
  struct rule_grammar_parser_static : boost::spirit::qi::grammar<Iterator, rule_parsed_type(), boost::spirit::standard::space_type>
  {
    
    rule_grammar_parser_static() : rule_grammar_parser_static::base_type(rule_grammar)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      lhs %= (qi::lexeme[standard::char_('[') >> +(standard::char_ - standard::space - ']') >> standard::char_(']')]);
      phrase %= *qi::lexeme[+(standard::char_ - standard::space) - "|||"];
      rule_grammar %= (qi::hold[lhs >> "|||"] | qi::attr("")) >> phrase >> "|||" >> phrase >> "|||" >> +qi::float_ >> -("|||" >> *qi::float_);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> lhs;
    boost::spirit::qi::rule<Iterator, phrase_parsed_type(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, rule_parsed_type(), boost::spirit::standard::space_type> rule_grammar;
  };

  struct GrammarParser
  {
    typedef boost::filesystem::path path_type;

    typedef std::vector<std::string, std::allocator<std::string> > phrase_parsed_type;
    typedef std::pair<std::string, double> score_parsed_type;
    typedef std::vector<score_parsed_type, std::allocator<score_parsed_type> > scores_parsed_type;
    
    typedef std::pair<std::string, AttributeVector::data_type> attr_parsed_type;
    typedef std::vector<attr_parsed_type, std::allocator<attr_parsed_type> > attrs_parsed_type;
    
    typedef boost::fusion::tuple<std::string, phrase_parsed_type, phrase_parsed_type, scores_parsed_type, attrs_parsed_type> rule_parsed_type;
    
    template <typename Iterator>
    struct parser : boost::spirit::qi::grammar<Iterator, rule_parsed_type(), boost::spirit::standard::space_type>
    {
    
      parser() : parser::base_type(rule_grammar)
      {
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	lhs %= (qi::lexeme[standard::char_('[') >> +(standard::char_ - standard::space - ']') >> standard::char_(']')]);
	phrase %= *(qi::lexeme[+(standard::char_ - standard::space) - "|||"]);
	
	score %= qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi)) >> (standard::char_ - standard::space))] >> '=' >> qi::double_;
	scores %= -(score % (+standard::space));
	
	data %= data_string | double_dot | int64_;
	
	attribute %= qi::lexeme[+(standard::char_ - standard::space - '=')] >> '=' >> data;
	attributes %= *attribute;
      
	rule_grammar %= (qi::hold[lhs >> "|||"] | qi::attr("")) >> phrase >> "|||" >> phrase >> -("|||" >> scores) >> -("|||" >> attributes);
      }
  
      typedef boost::spirit::standard::space_type space_type;

      boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1> int64_;
      boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > double_dot;
    
      boost::spirit::qi::rule<Iterator, std::string(), space_type> lhs;
      boost::spirit::qi::rule<Iterator, phrase_parsed_type(), space_type> phrase;
      
      boost::spirit::qi::rule<Iterator, score_parsed_type()>  score;
      boost::spirit::qi::rule<Iterator, scores_parsed_type()> scores;
      
      utils::json_string_parser<Iterator> data_string;
      
      boost::spirit::qi::rule<Iterator, AttributeVector::data_type(), space_type> data;
      boost::spirit::qi::rule<Iterator, attr_parsed_type(), space_type>           attribute;
      boost::spirit::qi::rule<Iterator, attrs_parsed_type(), space_type>          attributes;
      
      boost::spirit::qi::rule<Iterator, rule_parsed_type(), space_type> rule_grammar;
    };

    template <typename Data, typename Codec>
    struct DataStream
    {
      typedef uint64_t off_type;
      typedef std::vector<char, std::allocator<char> > buffer_type;
      
      DataStream(const path_type& path)
	: offset(0)
      {
	typedef utils::repository repository_type;
	
	repository_type rep(path, repository_type::write);
	
	os_data.push(boost::iostreams::file_sink(rep.path("data").string(), std::ios_base::out | std::ios_base::trunc), 1024 * 1024);
	os_offset.push(utils::vertical_coded_sink<off_type, std::allocator<off_type> >(rep.path("offset")));
      }
      
      template <typename Iterator>
      void insert(Iterator first, Iterator last)
      {
	buffer.clear();
	codec.encode(Data(first, last), std::back_inserter(buffer));
	offset += buffer.size();
	
	os_data.write((char*)&(*buffer.begin()), buffer.size());
	os_offset.write((char*)&offset, sizeof(off_type));
      }
      
      void clear()
      {
	offset = 0;
	os_data.pop();
	os_offset.pop();
      }

      bool empty() const { return ! offset; }
      
      Codec codec;
      buffer_type buffer;
      
      off_type offset;
      boost::iostreams::filtering_ostream os_data;
      boost::iostreams::filtering_ostream os_offset;
    };
    
    typedef DataStream<Transducer::feature_set_type, FeatureVectorCODEC>     feature_stream_type;
    typedef DataStream<Transducer::attribute_set_type, AttributeVectorCODEC> attribute_stream_type;
    
    template <typename Data>
    struct VocabStream : public utils::hashmurmur<uint64_t>
    {
      typedef utils::hashmurmur<uint64_t> hasher_type;

      VocabStream(const path_type& path)
      {
	typedef succinctdb::succinct_hash_stream<char, std::allocator<char> > succinct_hash_stream_type;
	
	const size_t size = Data::allocated();
	const size_t bin_size = utils::bithack::max(size, size_t(32));
	
	succinct_hash_stream_type data(path, bin_size);
	for (size_t i = 0; i != size; ++ i) {
	  const Data key(i);
	  data.insert(&(*key.begin()), key.size(), hasher_type::operator()(key.begin(), key.end(), 0));
	}

      }
    };

    typedef VocabStream<Transducer::feature_set_type::feature_type>     feature_vocab_type;
    typedef VocabStream<Transducer::attribute_set_type::attribute_type> attribute_vocab_type;
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

  template <typename Options, typename Codes, typename Id>
  inline
  void encode_options(const Options& options, Codes& codes, const Id& id_source)
  {
    typedef Id id_type;
    
    codes.clear();
    codes.resize(options.size() * 16 + 16, 0);
    
    typename Codes::iterator hiter = codes.begin();
    typename Codes::iterator citer = codes.begin();
    size_t pos = 0;

    const id_type id_feature = boost::get<0>(options.front());
	     
    const size_t offset_feature = utils::group_aligned_encode(id_feature, &(*hiter), pos);
    citer = hiter + offset_feature;
    hiter += offset_feature & (- size_t((pos & 0x03) == 0x03));
    ++ pos;
	     
    const size_t offset_source = utils::group_aligned_encode(id_source, &(*hiter), pos);
    citer = hiter + offset_source;
    hiter += offset_source & (- size_t((pos & 0x03) == 0x03));
    ++ pos;

    typename Options::const_iterator piter_end = options.end();
    for (typename Options::const_iterator piter = options.begin(); piter != piter_end; ++ piter) {
      const cicada::Symbol::id_type id_lhs = boost::get<1>(*piter);
      const id_type                 id_target = boost::get<2>(*piter);

      const size_t offset_lhs = utils::group_aligned_encode(id_lhs, &(*hiter), pos);
      citer = hiter + offset_lhs;
      hiter += offset_lhs & (- size_t((pos & 0x03) == 0x03));
      ++ pos;
	       
      const size_t offset_target = utils::group_aligned_encode(id_target, &(*hiter), pos);
      citer = hiter + offset_target;
      hiter += offset_target & (- size_t((pos & 0x03) == 0x03));
      ++ pos;
    }
    
    codes.resize(citer - codes.begin());
  }
  
  void GrammarStaticImpl::read_keyed_text(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const path_type path = param.name();
    
    typedef succinctdb::succinct_hash<byte_type, std::allocator<byte_type> > phrase_map_type;
    
    typedef boost::tuple<id_type, symbol_type::id_type, id_type> rule_option_type;
    typedef std::vector<rule_option_type, std::allocator<rule_option_type> > rule_option_set_type;
    
    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;

    typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > id_set_type;
    typedef std::vector<byte_type, std::allocator<byte_type> >  code_set_type;
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no file? ") + path.string());
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    const path_type path_rule   = utils::tempfile::directory_name(tmp_dir / "cicada.rule.XXXXXX");
    const path_type path_source = utils::tempfile::directory_name(tmp_dir / "cicada.source.XXXXXX");
    const path_type path_target = utils::tempfile::directory_name(tmp_dir / "cicada.target.XXXXXX");
    const path_type path_vocab  = utils::tempfile::directory_name(tmp_dir / "cicada.vocab.XXXXXX");
    
    const path_type path_feature_data  = utils::tempfile::directory_name(tmp_dir / "cicada.feature-data.XXXXXX");
    const path_type path_feature_vocab = utils::tempfile::directory_name(tmp_dir / "cicada.feature-vocab.XXXXXX");
    const path_type path_attribute_data  = utils::tempfile::directory_name(tmp_dir / "cicada.attribute-data.XXXXXX");
    const path_type path_attribute_vocab = utils::tempfile::directory_name(tmp_dir / "cicada.attribute-vocab.XXXXXX");
    
    utils::tempfile::insert(path_rule);
    utils::tempfile::insert(path_source);
    utils::tempfile::insert(path_target);
    utils::tempfile::insert(path_vocab);
    utils::tempfile::insert(path_feature_data);
    utils::tempfile::insert(path_feature_vocab);
    utils::tempfile::insert(path_attribute_data);
    utils::tempfile::insert(path_attribute_vocab);
    
    rule_db.open(path_rule, rule_db_type::WRITE);
    std::auto_ptr<phrase_map_type> source_map(new phrase_map_type(1024 * 1024 * 128));
    std::auto_ptr<phrase_map_type> target_map(new phrase_map_type(1024 * 1024 * 128));

    GrammarParser::feature_stream_type   feature_stream(path_feature_data);
    GrammarParser::attribute_stream_type attribute_stream(path_attribute_data);
    
    id_type id_rule = 0;
    sequence_type source_prev;
    sequence_type source;
    sequence_type target;

    GrammarParser::rule_parsed_type rule;
    
    id_set_type source_index;
    
    code_set_type codes_source;
    code_set_type codes_target;
    code_set_type codes_option;

    rule_option_set_type rule_options;

    utils::compress_istream is(path, 1024 * 1024);
    
    std::string line;
    
    GrammarParser::parser<std::string::const_iterator> rule_parser;
    
    size_type arity_source = 0;
    
    for (size_t num_line = 0; std::getline(is, line); ++ num_line) {
      if (line.empty()) continue;

      boost::fusion::get<0>(rule).clear();
      boost::fusion::get<1>(rule).clear();
      boost::fusion::get<2>(rule).clear();
      boost::fusion::get<3>(rule).clear();
      boost::fusion::get<4>(rule).clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);
      
      if (! result || iter != iter_end) {
	std::cerr << "invalid line: " << num_line << ": " << line << std::endl;
	continue;
      }

      // skip empty source...
      if (boost::fusion::get<1>(rule).empty()) continue;
      
      // lhs...
      const std::string& lhs = boost::fusion::get<0>(rule);
      const word_type::id_type id_lhs = word_type(lhs.empty() ? vocab_type::X : word_type(lhs)).id();
      
      // source
      source.clear();
      source.insert(source.end(), boost::fusion::get<1>(rule).begin(), boost::fusion::get<1>(rule).end());
      
      // target
      target.clear();
      target.insert(target.end(), boost::fusion::get<2>(rule).begin(), boost::fusion::get<2>(rule).end());

      if (source != source_prev) {
	if (! rule_options.empty()) {
	  encode_phrase(source_prev, codes_source);
	     
	  const id_type id_source = source_map->insert(&(*codes_source.begin()), codes_source.size(),
						       hasher_type::operator()(codes_source.begin(), codes_source.end(), 0));
	  
	  // encode options
	  encode_options(rule_options, codes_option, id_source);
	  
	  // encode source.. we will use index-stripped indexing!
	  encode_index(source_prev.begin(), source_prev.end(), source_index);
	  
	  // insert...
	  rule_db.insert(&(*source_index.begin()), source_index.size(), &(*codes_option.begin()), codes_option.size());
	}
	
	rule_options.clear();
	
	source_prev = source;
	
	arity_source = 0;
	for (sequence_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter)
	  arity_source += siter->is_non_terminal();
      }
      
      if (! target.empty()) {
	size_type arity_target = 0;
	for (sequence_type::const_iterator titer = target.begin(); titer != target.end(); ++ titer)
	  arity_target += titer->is_non_terminal();
	
	if (arity_source != arity_target)
	  throw std::runtime_error("# of non-terminals do not match...");
      }

      // encode features/attributes
      feature_stream.insert(boost::fusion::get<3>(rule).begin(), boost::fusion::get<3>(rule).end());
      attribute_stream.insert(boost::fusion::get<4>(rule).begin(), boost::fusion::get<4>(rule).end());
      
      // encode target...
      encode_phrase(target, codes_target);
      
      const id_type id_target = target_map->insert(&(*codes_target.begin()), codes_target.size(),
						   hasher_type::operator()(codes_target.begin(), codes_target.end(), 0));
      
      // put into rule_options...
      rule_options.push_back(boost::make_tuple(id_rule ++, id_lhs, id_target));
    }
    
    if (! rule_options.empty()) {
      encode_phrase(source_prev, codes_source);
      
      const id_type id_source = source_map->insert(&(*codes_source.begin()), codes_source.size(),
						   hasher_type::operator()(codes_source.begin(), codes_source.end(), 0));
      
      // encode options
      encode_options(rule_options, codes_option, id_source);
      
      // encode source.. we will use index-stripped indexing!
      encode_index(source_prev.begin(), source_prev.end(), source_index);
      
      // insert...
      rule_db.insert(&(*source_index.begin()), source_index.size(), &(*codes_option.begin()), codes_option.size()); 
    }
    
    source_map->write(path_source);
    source_map.reset();
    
    target_map->write(path_target);
    target_map.reset();

    const bool has_features   = ! feature_stream.empty();
    const bool has_attributes = ! attribute_stream.empty();
    
    feature_stream.clear();
    attribute_stream.clear();
    
    if (has_features) {
      GrammarParser::feature_vocab_type   __feature_vocab(path_feature_vocab);
    }
    
    if (has_attributes) {
      GrammarParser::attribute_vocab_type __attribute_vocab(path_attribute_vocab);
    }
    
    rule_db.close();
    
    word_type::write(path_vocab);
    
    ::sync();

    while (! phrase_db_type::exists(path_source))
      boost::thread::yield();
    while (! phrase_db_type::exists(path_target))
      boost::thread::yield();
    while (! vocab_type::exists(path_vocab))
      boost::thread::yield();
    while (! rule_db_type::exists(path_rule))
      boost::thread::yield();
    
    source_db.open(path_source);
    target_db.open(path_target);
    rule_db.open(path_rule);
    vocab.open(path_vocab);
    
    if (has_features) {
      feature_data.open(path_feature_data);
      feature_vocab.open(path_feature_vocab);
    }
    
    if (has_attributes) {
      attribute_data.open(path_attribute_data);
      attribute_vocab.open(path_attribute_vocab);
    }
  }

  void GrammarStaticImpl::read_text(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const path_type path = param.name();
    
    typedef succinctdb::succinct_hash<byte_type, std::allocator<byte_type> > phrase_map_type;
    
    typedef ScoreSetStream score_stream_type;
    typedef std::vector<score_stream_type, std::allocator<score_stream_type> > score_stream_set_type;
    
    // feature-id, lhs-id, target-id
    typedef boost::tuple<id_type, symbol_type::id_type, id_type> rule_option_type;
    typedef std::vector<rule_option_type, std::allocator<rule_option_type> > rule_option_set_type;

    typedef std::vector<symbol_type, std::allocator<symbol_type> > sequence_type;
    
    typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > id_set_type;
    
    typedef std::vector<byte_type, std::allocator<byte_type> >  code_set_type;
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no file? ") + path.string());
    
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
    std::auto_ptr<phrase_map_type> source_map(new phrase_map_type(1024 * 1024 * 128));
    std::auto_ptr<phrase_map_type> target_map(new phrase_map_type(1024 * 1024 * 128));
    
    score_stream_set_type score_streams;
    score_stream_set_type attr_streams;
    
    id_type id_rule = 0;
    
    int feature_size = -1;
    int attribute_size = -1;

    sequence_type source_prev;
    sequence_type source;
    sequence_type target;
    
    rule_parsed_type rule;
    
    id_set_type source_index;
    
    code_set_type codes_source;
    code_set_type codes_target;
    code_set_type codes_option;
    
    rule_option_set_type rule_options;
    
    utils::compress_istream is(path, 1024 * 1024);
    
    std::string line;
    
    // we will construct this parser everytimt...
    rule_grammar_parser_static<std::string::const_iterator> rule_parser;

    size_type arity_source = 0;

    for (size_t num_line = 0; std::getline(is, line); ++ num_line) {
      if (line.empty()) continue;

      boost::fusion::get<0>(rule).clear();
      boost::fusion::get<1>(rule).clear();
      boost::fusion::get<2>(rule).clear();
      boost::fusion::get<3>(rule).clear();
      boost::fusion::get<4>(rule).clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);
      
      if (! result || iter != iter_end) {
	std::cerr << "invalid line: " << num_line << ": " << line << std::endl;
	continue;
      }

      // skip empty source...
      if (boost::fusion::get<1>(rule).empty()) continue;
      
      // lhs...
      const std::string& lhs = boost::fusion::get<0>(rule);
      const word_type::id_type id_lhs = word_type(lhs.empty() ? vocab_type::X : word_type(lhs)).id();
      
      // source
      source.clear();
      source.insert(source.end(), boost::fusion::get<1>(rule).begin(), boost::fusion::get<1>(rule).end());
      
      // target
      target.clear();
      target.insert(target.end(), boost::fusion::get<2>(rule).begin(), boost::fusion::get<2>(rule).end());
      
      if (source != source_prev) {
	
	if (! rule_options.empty()) {
	  // encode options...
	   
	  encode_phrase(source_prev, codes_source);
	     
	  const id_type id_source = source_map->insert(&(*codes_source.begin()), codes_source.size(),
						       hasher_type::operator()(codes_source.begin(), codes_source.end(), 0));

	  // encode options
	  encode_options(rule_options, codes_option, id_source);
	   
	  // encode source.....
	  encode_index(source_prev.begin(), source_prev.end(), source_index);

	  // insert...
	  rule_db.insert(&(*source_index.begin()), source_index.size(), &(*codes_option.begin()), codes_option.size());
	}
	
	rule_options.clear();
	
	source_prev = source;

	arity_source = 0;
	for (sequence_type::const_iterator siter = source.begin(); siter != source.end(); ++ siter)
	  arity_source += siter->is_non_terminal();
      }
       
      if (! target.empty()) {
	size_type arity_target = 0;
	for (sequence_type::const_iterator titer = target.begin(); titer != target.end(); ++ titer)
	  arity_target += titer->is_non_terminal();
	
	if (arity_source != arity_target)
	  throw std::runtime_error("# of non-terminals do not match...");
      }
       
#if 0
      std::cerr << "lhs: " << lhs;
      std::cerr << " source: ";
      std::copy(source_prev.begin(), source_prev.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << "target: ";
      std::copy(target.begin(), target.end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << "features: ";
      std::copy(boost::fusion::get<3>(rule).begin(), boost::fusion::get<3>(rule).end(), std::ostream_iterator<score_type>(std::cerr, " "));
      std::cerr << "attribute: ";
      std::copy(boost::fusion::get<4>(rule).begin(), boost::fusion::get<4>(rule).end(), std::ostream_iterator<score_type>(std::cerr, " "));
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

      if (attribute_size < 0) {
	attribute_size = boost::fusion::get<4>(rule).size();
	
	attr_streams.reserve(attribute_size);
	attr_streams.resize(attribute_size);
	
	for (int attribute = 0; attribute < attribute_size; ++ attribute) {
	  attr_streams[attribute].path = utils::tempfile::file_name(tmp_dir / "cicada.attribute.XXXXXX");
	  utils::tempfile::insert(attr_streams[attribute].path);
	  
	  attr_streams[attribute].ostream.reset(new utils::compress_ostream(attr_streams[attribute].path, 1024 * 1024));
	}
      } else if (attribute_size != static_cast<int>(boost::fusion::get<4>(rule).size()))
	throw std::runtime_error("invalid # of attributes...");
      
      for (int feature = 0; feature < feature_size; ++ feature)
	score_streams[feature].ostream->write((char*) &boost::fusion::get<3>(rule)[feature], sizeof(score_type));
      
      for (int attribute = 0; attribute < attribute_size; ++ attribute)
	attr_streams[attribute].ostream->write((char*) &boost::fusion::get<4>(rule)[attribute], sizeof(score_type));
      
      // encode target...
      encode_phrase(target, codes_target);
      
      const id_type id_target = target_map->insert(&(*codes_target.begin()), codes_target.size(),
						   hasher_type::operator()(codes_target.begin(), codes_target.end(), 0));
      
      // put into rule_options...
      rule_options.push_back(boost::make_tuple(id_rule ++, id_lhs, id_target));
    }
     
    if (! rule_options.empty()) {
      // encode options...
	
      encode_phrase(source_prev, codes_source);
	     
      const id_type id_source = source_map->insert(&(*codes_source.begin()), codes_source.size(),
						   hasher_type::operator()(codes_source.begin(), codes_source.end(), 0));
      
      // encode options
      encode_options(rule_options, codes_option, id_source);
	   
      // encode source.. we will use index-stripped indexing!
      encode_index(source_prev.begin(), source_prev.end(), source_index);

      // insert...
      rule_db.insert(&(*source_index.begin()), source_index.size(), &(*codes_option.begin()), codes_option.size());
    }

    source_map->write(path_source);
    source_map.reset();
    
    target_map->write(path_target);
    target_map.reset();
    
    rule_db.close();
    
    word_type::write(path_vocab);
    
    ::sync();
    
    while (! phrase_db_type::exists(path_source))
      boost::thread::yield();
    while (! phrase_db_type::exists(path_target))
      boost::thread::yield();
    while (! vocab_type::exists(path_vocab))
      boost::thread::yield();
    while (! rule_db_type::exists(path_rule))
      boost::thread::yield();
    
    source_db.open(path_source);
    target_db.open(path_target);
    rule_db.open(path_rule);
    vocab.open(path_vocab);

    if (feature_size < 0)
      feature_size = 0;
    
    // scores...
    score_db.reserve(feature_size);
    score_db.resize(feature_size);
    
    feature_names.clear();
    feature_names.reserve(feature_size);
    feature_names.resize(feature_size, feature_type());
    
    for (int feature = 0; feature < feature_size; ++ feature) {
      score_streams[feature].ostream->reset();
      ::sync();
      
      while (! score_set_type::score_set_type::exists(score_streams[feature].path))
	boost::thread::yield();
      
      utils::tempfile::permission(score_streams[feature].path);
      score_db[feature].score.open(score_streams[feature].path);

      const std::string name("feature" + utils::lexical_cast<std::string>(feature));

      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	feature_names[feature] = feature_type(piter->second);
      
      // default name...!
      if (feature_names[feature] == feature_type())
	feature_names[feature] = "rule-table-" + utils::lexical_cast<std::string>(feature);
    }

    if (attribute_size < 0)
      attribute_size = 0;
    
    // attributes
    attr_db.reserve(attribute_size);
    attr_db.resize(attribute_size);
    
    attribute_names.clear();
    attribute_names.reserve(attribute_size);
    attribute_names.resize(attribute_size, attribute_type());
    
    for (int attribute = 0; attribute < attribute_size; ++ attribute) {
      attr_streams[attribute].ostream->reset();
      ::sync();
      
      while (! score_set_type::score_set_type::exists(attr_streams[attribute].path))
	boost::thread::yield();
      
      utils::tempfile::permission(attr_streams[attribute].path);
      attr_db[attribute].score.open(attr_streams[attribute].path);

      const std::string name(std::string("attribute") + utils::lexical_cast<std::string>(attribute));

      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	attribute_names[attribute] = attribute_type(piter->second);
      
      // default name...!
      if (attribute_names[attribute] == attribute_type())
	attribute_names[attribute] = std::string("rule-table-") + utils::lexical_cast<std::string>(attribute);
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
    // max-span checking + distance checking
    // we need this last - first == 1 when intersecting with lattice...
    return pimpl->max_span <= 0 || distance <= pimpl->max_span || last - first == 1;
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
  
  const GrammarStatic::rule_pair_set_type& GrammarStatic::rules(const id_type& node) const
  {
    static const rule_pair_set_type __empty;
    
    return (pimpl->is_valid(node) && pimpl->exists(node) ? pimpl->read_rule_set(node) : __empty);
  }
  
  void GrammarStatic::quantize()
  {
    pimpl->quantize();
  }
  
  void GrammarStatic::write(const path_type& path) const
  {
    pimpl->write(path);
  }
};
