//
//  Copyright(C) 2010-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_DISABLE_ASSERTS
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
#include "utils/compact_map.hpp"
#include "utils/succinct_vector.hpp"
#include "utils/resource.hpp"
#include "utils/unordered_map.hpp"
#include "utils/getline.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/hashmurmur3.hpp"

#include <boost/lexical_cast.hpp>
#include <boost/bind.hpp>

#include <boost/filesystem.hpp>

#include <boost/spirit/include/qi.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <boost/thread.hpp>

namespace cicada
{
  
  // TODO: we need to index by index-stripped non-terminals!

  static const size_t DEBUG_DOT = 1000000;
  static const size_t DEBUG_WRAP = 100;
  static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

  class GrammarStaticImpl : public utils::hashmurmur<uint64_t>
  {
  public:
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
      void populate() { data.populate(); }
      
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

      void populate()
      {
	data.populate();
	offset.populate();
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

      static
      bool exists(const path_type& path)
      {
	return (utils::repository::exists(path)
		&& data_type::exists(path / "data")
		&& offset_type::exists(path / "offset"));
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
      typedef utils::succinct_vector_mapped<>                                              binarized_set_type;
      typedef boost::array<score_type, 256>                                                score_map_type;
      
      ScoreSet() {}
      ScoreSet(const path_type& path) { read(path); }
      
      void read(const path_type& path);
      void write(const path_type& file) const;

      void populate()
      {
	score.populate();
	quantized.populate();
	binarized.populate();
      }
      
      void clear()
      {
	score.clear();
	quantized.clear();
	binarized.clear();
      }
      void close() { clear(); }
      
      score_type operator[](size_type pos) const
      {
	return (binarized.is_open()
		? maps[binarized[pos]]
		: (quantized.is_open()
		   ? maps[quantized[pos]]
		   : score[pos]));
      }
      
      path_type path() const
      {
	return (binarized.is_open()
		? binarized.path().parent_path()
		: (quantized.is_open()
		   ? quantized.path().parent_path()
		   : score.path().parent_path()));
      }
      bool empty() const { return quantized.empty() && binarized.empty() && score.empty(); }
      size_type size() const
      {
	return (binarized.is_open()
		? binarized.size()
		: (quantized.is_open()
		   ? quantized.size()
		   : score.size()));
      }
      
      score_set_type     score;
      quantized_set_type quantized;
      binarized_set_type binarized;
      score_map_type     maps;
    };
    
    typedef ScoreSet score_set_type;
    typedef std::vector<score_set_type, std::allocator<score_set_type> > score_db_type;
    
    // caching...
    typedef utils::arc_list<size_type, rule_pair_set_type, 4,
			    std::equal_to<size_type>,
			    std::allocator<std::pair<size_type, rule_pair_set_type> > > cache_rule_set_type;

    struct cache_phrase_type
    {
      rule_ptr_type rule;
      size_type     pos;
      cache_phrase_type() : rule(), pos(size_type(-1)) {}
    };
    
    struct cache_node_type
    {
      size_type node;
      size_type next;
      uint32_t  id;

      cache_node_type() : node(size_type(-1)), next(size_type(-1)), id(uint32_t(-1)) {}
    };

    typedef utils::array_power2<cache_rule_set_type, 1024 * 2, std::allocator<cache_rule_set_type> > cache_rule_map_type;
    typedef utils::array_power2<cache_phrase_type,   1024 * 2, std::allocator<cache_phrase_type> >   cache_phrase_set_type;
    typedef utils::array_power2<cache_node_type,     1024 * 8, std::allocator<cache_node_type> >     cache_node_set_type;
    
    typedef std::vector<size_type, std::allocator<size_type> > cache_root_type;
    
    typedef std::pair<word_type, size_type> word_node_type;

    struct unassigned_cache
    {
      word_node_type operator()() const
      {
	utils::unassigned<word_type> __unassigned;
	return word_node_type(__unassigned(), size_type(-1));
      }
    };

    struct deleted_cache
    {
      word_node_type operator()() const
      {
	utils::deleted<word_type> __deleted;
	return word_node_type(__deleted(), size_type(-1));
      }
    };

  public:
    GrammarStaticImpl(const std::string& parameter)
      : max_span(0),
	debug(0)
    {
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
	debug(x.debug)
    { }

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
      debug           = x.debug;
      
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
      cache_nodes.clear();
      cache_root.clear();

      max_span = 0;
    }
    
    size_type find(const word_type& word, size_type node) const
    {
      if (node == 0 && word.id() < 1024 * 16) {
	cache_root_type& cache = const_cast<cache_root_type&>(cache_root);
	
	if (word.id() >= cache.size()) {
	  cache.resize(word.id() + 1, 0);
	  
	  if (cache.capacity() > 1024 * 16) {
	    cache.resize(1024 * 16, 0);
	    cache_root_type(cache).swap(cache);
	  }
	}
	
	if (cache[word.id()] == 0) {
	  const word_type::id_type id = vocab[word];
	  
	  cache[word.id()] = (id != word_type::id_type(-1) ? rule_db.find(&id, 1, 0) : rule_db_type::out_of_range());
	}
	
	return cache[word.id()];
      } else {
	typedef utils::hashmurmur3<size_t> hasher_type;
	
	const size_type cache_pos = hasher_type()(word.id(), node) & (cache_nodes.size() - 1);
	cache_node_type& cache = const_cast<cache_node_type&>(cache_nodes[cache_pos]);
	
	if (cache.node != node || cache.id != word.id()) {
	  const word_type::id_type id = vocab[word];
	  
	  cache.next = (id != word_type::id_type(-1) ? rule_db.find(&id, 1, node) : rule_db_type::out_of_range());
	  cache.node = node;
	  cache.id   = word.id();
	}
	
	return cache.next;
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
      
      const size_type cache_pos = node & (cache_rule_sets.size() - 1);
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
	    
	    // rehash...
	    options.back().features.rehash(score_db.size());
	    
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
	
	rule_pair_set_type(options).swap(options);
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
      typedef utils::hashmurmur3<size_t> hasher_type;
      
      const size_type cache_pos = hasher_type()(pos, lhs.id()) & (cache_phrases.size() - 1);
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

    void binarize();
    void quantize();
    void read(const std::string& parameter);
    void write(const path_type& path) const;

    void populate()
    {
      //boost::thread_group workers;

      rule_db.populate();
      source_db.populate();
      target_db.populate();
      //workers.add_thread(new boost::thread(boost::bind(&rule_db_type::populate, boost::ref(rule_db))));
      //workers.add_thread(new boost::thread(boost::bind(&phrase_db_type::populate, boost::ref(source_db))));
      //workers.add_thread(new boost::thread(boost::bind(&phrase_db_type::populate, boost::ref(target_db))));
      
      for (size_t feature = 0; feature < score_db.size(); ++ feature) {
	score_db[feature].populate();
	//workers.add_thread(new boost::thread(boost::bind(&score_set_type::populate, boost::ref(score_db[feature]))));
      }
      
      for (size_t attr = 0; attr < attr_db.size(); ++ attr) {
	attr_db[attr].populate();
	//workers.add_thread(new boost::thread(boost::bind(&score_set_type::populate, boost::ref(attr_db[attr]))));
      }
      
      feature_data.populate();
      attribute_data.populate();
      //workers.add_thread(new boost::thread(boost::bind(&feature_data_type::populate, boost::ref(feature_data))));
      //workers.add_thread(new boost::thread(boost::bind(&attribute_data_type::populate, boost::ref(attribute_data))));
      
      
      feature_vocab.populate();
      attribute_vocab.populate();
      //workers.add_thread(new boost::thread(boost::bind(&feature_vocab_type::populate, boost::ref(feature_vocab))));
      //workers.add_thread(new boost::thread(boost::bind(&attribute_vocab_type::populate, boost::ref(attribute_vocab))));
      
      vocab.populate();
      
      //workers.join_all();
    }
    
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
    
    cache_node_set_type   cache_nodes;
    cache_root_type       cache_root;

  public:
    int max_span;
    int debug;
  };
  
  void GrammarStaticImpl::ScoreSet::read(const path_type& path)
  {
    typedef utils::repository repository_type;
    
    clear();

    repository_type rep(path, repository_type::read);

    if (boost::filesystem::exists(rep.path("binarized"))) {
      binarized.open(rep.path("binarized"));
      
      const path_type score_map_file = rep.path("score-map");
      
      if (! boost::filesystem::exists(score_map_file))
	throw std::runtime_error(std::string("no map file? ") + score_map_file.string());
      
      std::ifstream is(score_map_file.string().c_str());
      is.read((char*) &(*maps.begin()), sizeof(score_type) * maps.size());
    } else if (boost::filesystem::exists(rep.path("quantized"))) {
      quantized.open(rep.path("quantized"));
      
      const path_type score_map_file = rep.path("score-map");
      
      if (! boost::filesystem::exists(score_map_file))
	throw std::runtime_error(std::string("no map file? ") + score_map_file.string());
      
      std::ifstream is(score_map_file.string().c_str());
      is.read((char*) &(*maps.begin()), sizeof(score_type) * maps.size());
    } else if (boost::filesystem::exists(rep.path("score")))
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
    }
    
    if (binarized.is_open()) {
      binarized.write(rep.path("binarized"));
      
      std::ofstream os(rep.path("score-map").string().c_str());
      os.write((char*) &(*maps.begin()), sizeof(score_type) * maps.size());
    }
    
    if (score.is_open())
      score.write(rep.path("score"));
  }
  
  void GrammarStaticImpl::binarize()
  {
    typedef std::map<score_type, size_type, std::less<score_type>,
		     std::allocator<std::pair<const score_type, size_type> > > counts_type;
    
    typedef utils::succinct_vector<> bit_vector_type;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    utils::resource start;

    for (size_t feature = 0; feature < score_db.size(); ++ feature)
      if (score_db[feature].score.is_open()) {
	
	counts_type counts;
	
	score_set_type::score_set_type::const_iterator liter_end = score_db[feature].score.end();
	for (score_set_type::score_set_type::const_iterator liter = score_db[feature].score.begin(); liter != liter_end && counts.size() <= 2; ++ liter)
	  ++ counts[*liter];
	
	if (counts.size() > 2 || counts.empty()) continue;
	
	// we perform binarization!

	if (debug)
	  std::cerr << "binarize feature: " << feature << std::endl;

	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.score.binary.XXXXXX");
	utils::tempfile::insert(path);

	bit_vector_type bits;
	
	bits.set(score_db[feature].size() - 1, false);
	
	if (counts.size() == 1) // we have a single value only!
	  score_db[feature].maps[0] = counts.begin()->first;
	else {
	  counts_type::const_iterator citer1 = counts.begin();
	  counts_type::const_iterator citer2 = counts.begin();
	  ++ citer2;
	  
	  if (citer1->second > citer2->second) {
	    score_db[feature].maps[0] = citer1->first;
	    score_db[feature].maps[1] = citer2->first;
	  } else {
	    score_db[feature].maps[0] = citer2->first;
	    score_db[feature].maps[1] = citer1->first;
	  }
	  
	  size_t pos = 0;
	  score_set_type::score_set_type::const_iterator liter_end = score_db[feature].score.end();
	  for (score_set_type::score_set_type::const_iterator liter = score_db[feature].score.begin(); liter != liter_end; ++ liter, ++ pos)
	    if (*liter == score_db[feature].maps[1])
	      bits.set(pos, true);
	}
	
	bits.write(path);
	
	while (! score_set_type::binarized_set_type::exists(path)) {
	  ::sync();
	  boost::thread::yield();
	}
	
	utils::tempfile::permission(path);
	
	score_db[feature].binarized.open(path);
	score_db[feature].score.clear();
      }
    
    for (size_t attr = 0; attr < attr_db.size(); ++ attr)
      if (attr_db[attr].score.is_open()) {

	counts_type counts;
	
	score_set_type::score_set_type::const_iterator liter_end = attr_db[attr].score.end();
	for (score_set_type::score_set_type::const_iterator liter = attr_db[attr].score.begin(); liter != liter_end && counts.size() <= 2; ++ liter)
	  ++ counts[*liter];
	
	if (counts.size() > 2 || counts.empty()) continue;
	
	// we perform binarization!
	if (debug)
	  std::cerr << "binarize attribute: " << attr << std::endl;
	
	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.attr.binary.XXXXXX");
	utils::tempfile::insert(path);
	
	bit_vector_type bits;
	
	bits.set(attr_db[attr].size() - 1, false);
	
	if (counts.size() == 1) // we have a single value only!
	  attr_db[attr].maps[0] = counts.begin()->first;
	else {
	  counts_type::const_iterator citer1 = counts.begin();
	  counts_type::const_iterator citer2 = counts.begin();
	  ++ citer2;
	  
	  if (citer1->second > citer2->second) {
	    attr_db[attr].maps[0] = citer1->first;
	    attr_db[attr].maps[1] = citer2->first;
	  } else {
	    attr_db[attr].maps[0] = citer2->first;
	    attr_db[attr].maps[1] = citer1->first;
	  }
	  
	  size_t pos = 0;
	  score_set_type::score_set_type::const_iterator liter_end = attr_db[attr].score.end();
	  for (score_set_type::score_set_type::const_iterator liter = attr_db[attr].score.begin(); liter != liter_end; ++ liter, ++ pos)
	    if (*liter == attr_db[attr].maps[1])
	      bits.set(pos, true);
	}
	
	bits.write(path);
	
	while (! score_set_type::binarized_set_type::exists(path)) {
	  ::sync();
	  boost::thread::yield();
	}
	
	utils::tempfile::permission(path);
	
	attr_db[attr].binarized.open(path);
	attr_db[attr].score.clear();
      }

    utils::resource end;
    
    if (debug)
      std::cerr << "binarization:"
		<< " cpu time: " << end.cpu_time() - start.cpu_time()
		<< " user time: " << end.user_time() - start.user_time()
		<< std::endl;
  }

  void GrammarStaticImpl::quantize()
  {
    typedef score_type base_type;

    typedef utils::unordered_map<base_type, size_type, boost::hash<base_type>, std::equal_to<base_type>,
				 std::allocator<std::pair<const base_type, size_type> > >::type hashed_type;

    typedef std::map<base_type, size_type, std::less<base_type>,
		     std::allocator<std::pair<const base_type, size_type> > > counts_type;
    typedef utils::unordered_map<base_type, quantized_type, boost::hash<base_type>, std::equal_to<base_type>,
				 std::allocator<std::pair<const base_type, quantized_type> > >::type codemap_type;
    typedef boost::array<base_type, 256> codebook_type;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    utils::resource start;

    hashed_type      hashed;
    counts_type      counts;
    codebook_type    codebook;
    codemap_type     codemap;
    
    for (size_t feature = 0; feature < score_db.size(); ++ feature)
      if (score_db[feature].score.is_open()) {

	if (debug)
	  std::cerr << "quantize feature: " << feature << std::endl;
	
	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.score.quantized.XXXXXX");
	utils::tempfile::insert(path);
	
	boost::iostreams::filtering_ostream os;
	os.push(utils::packed_sink<quantized_type, std::allocator<quantized_type> >(path));
	
	hashed.clear();
	counts.clear();
	codemap.clear();
	std::fill(codebook.begin(), codebook.end(), 0.0);
	
	score_set_type::score_set_type::const_iterator liter_end = score_db[feature].score.end();
	for (score_set_type::score_set_type::const_iterator liter = score_db[feature].score.begin(); liter != liter_end; ++ liter)
	  ++ hashed[*liter];
	
	counts.insert(hashed.begin(), hashed.end());
	hashed.clear();
	
	Quantizer::quantize(counts, codebook, codemap);
	
	for (score_set_type::score_set_type::const_iterator liter = score_db[feature].score.begin(); liter != liter_end; ++ liter) {
	  codemap_type::const_iterator citer = codemap.find(*liter);
	  if (citer == codemap.end())
	    throw std::runtime_error("no codemap?");
	  
	  os.write((char*) &(citer->second), sizeof(quantized_type));
	}
	
	for (int i = 0; i < 256; ++ i)
	  score_db[feature].maps[i] = codebook[i];
	
	os.reset();
	
	while (! score_set_type::quantized_set_type::exists(path)) {
	  ::sync();
	  boost::thread::yield();
	}

	utils::tempfile::permission(path);
	
	score_db[feature].quantized.open(path);
	score_db[feature].score.clear();
      }

    for (size_t attr = 0; attr < attr_db.size(); ++ attr)
      if (attr_db[attr].score.is_open()) {

	if (debug)
	  std::cerr << "quantize attribute: " << attr << std::endl;
	
	const path_type path = utils::tempfile::directory_name(tmp_dir / "cicada.attr.quantized.XXXXXX");
	utils::tempfile::insert(path);
	
	boost::iostreams::filtering_ostream os;
	os.push(utils::packed_sink<quantized_type, std::allocator<quantized_type> >(path));
	
	hashed.clear();
	counts.clear();
	codemap.clear();
	std::fill(codebook.begin(), codebook.end(), 0.0);
	
	score_set_type::score_set_type::const_iterator liter_end = attr_db[attr].score.end();
	for (score_set_type::score_set_type::const_iterator liter = attr_db[attr].score.begin(); liter != liter_end; ++ liter)
	  ++ hashed[*liter];
	
	counts.insert(hashed.begin(), hashed.end());
	hashed.clear();
	
	Quantizer::quantize(counts, codebook, codemap);
	
	for (score_set_type::score_set_type::const_iterator liter = attr_db[attr].score.begin(); liter != liter_end; ++ liter) {
	  codemap_type::const_iterator citer = codemap.find(*liter);
	  if (citer == codemap.end())
	    throw std::runtime_error("no codemap?");
	  
	  os.write((char*) &(citer->second), sizeof(quantized_type));
	}
	
	for (int i = 0; i < 256; ++ i)
	  attr_db[attr].maps[i] = codebook[i];
	
	os.reset();
	
	while (! score_set_type::quantized_set_type::exists(path)) {
	  ::sync();
	  boost::thread::yield();
	}
	
	utils::tempfile::permission(path);
	
	attr_db[attr].quantized.open(path);
	attr_db[attr].score.clear();
      }

    utils::resource end;
    
    if (debug)
      std::cerr << "quantization:"
		<< " cpu time: " << end.cpu_time() - start.cpu_time()
		<< " user time: " << end.user_time() - start.user_time()
		<< std::endl;
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
    
    parameter_type::const_iterator diter = param.find("debug");
    if (diter != param.end())
      debug = utils::lexical_cast<int>(diter->second);
    
    if (boost::filesystem::is_directory(path))
      read_binary(parameter);
    else if (key_value)
      read_keyed_text(parameter);
    else
      read_text(parameter);
    
    parameter_type::const_iterator siter = param.find("max-span");
    if (siter != param.end())
      max_span = utils::lexical_cast<int>(siter->second);

    parameter_type::const_iterator piter = param.find("populate");
    if (piter != param.end() && utils::lexical_cast<bool>(piter->second))
      populate();
  }
  
  void GrammarStaticImpl::write(const path_type& file) const
  {
    typedef utils::repository repository_type;
    
    if (file == path()) return;
    
    if (debug)
      std::cerr << "grammar writing at: " << file.string() << std::endl;
    
    utils::resource start;

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

    utils::resource end;
    
    if (debug)
      std::cerr << "writing:"
		<< " cpu time: " << end.cpu_time() - start.cpu_time()
		<< " user time: " << end.user_time() - start.user_time()
		<< std::endl;
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
    if (aiter != rep.end()) {
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
      word %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];
      phrase %= *word;
      rule_grammar %= (qi::hold[lhs >> "|||"] | qi::attr("")) >> phrase >> "|||" >> phrase >> "|||" >> +qi::float_ >> -("|||" >> *qi::float_);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> lhs;
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> word;
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
	
	word %= qi::lexeme[+(standard::char_ - standard::space) - ("|||" >> (standard::space | qi::eoi))];
	phrase %= *word;
	
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
      boost::spirit::qi::rule<Iterator, std::string(), space_type> word;
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
	os_data.reset();
	os_offset.reset();
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

    typedef GrammarStaticImpl::hasher_type hasher_type;
    
    template <typename Data>
    struct VocabStream : public hasher_type
    {
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
    std::unique_ptr<phrase_map_type> source_map(new phrase_map_type(1024 * 1024 * 128));
    std::unique_ptr<phrase_map_type> target_map(new phrase_map_type(1024 * 1024 * 128));

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

    utils::resource read_start;
    
    size_t num_line = 0;
    for (/**/; utils::getline(is, line); ++ num_line) {
      if (line.empty()) continue;

      boost::fusion::get<0>(rule).clear();
      boost::fusion::get<1>(rule).clear();
      boost::fusion::get<2>(rule).clear();
      boost::fusion::get<3>(rule).clear();
      boost::fusion::get<4>(rule).clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);
      
      if (debug) {
	if ((num_line + 1) % DEBUG_DOT == 0)
	  std::cerr << '.';
	if ((num_line + 1) % DEBUG_LINE == 0)
	  std::cerr << std::endl;
      } 

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
    
    utils::resource read_end;
    
    if (debug) {
      if ((num_line / DEBUG_DOT) % DEBUG_WRAP)
	std::cerr << std::endl;
      
      std::cerr << "# of rules: " << num_line
		<< " cpu time: " << read_end.cpu_time() - read_start.cpu_time()
		<< " user time: " << read_end.user_time() - read_start.user_time()
		<< std::endl;
    }

    utils::resource index_start;
    
    source_map->prune(static_cast<const hasher_type&>(*this));
    source_map->write(path_source);
    source_map.reset();
    
    target_map->prune(static_cast<const hasher_type&>(*this));
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

    while (! phrase_db_type::exists(path_source)) {
      ::sync();
      boost::thread::yield();
    }
    while (! phrase_db_type::exists(path_target)) {
      ::sync();
      boost::thread::yield();
    }
    while (! vocab_type::exists(path_vocab)) {
      ::sync();
      boost::thread::yield();
    }
    while (! rule_db_type::exists(path_rule)) {
      ::sync();
      boost::thread::yield();
    }
    
    source_db.open(path_source);
    target_db.open(path_target);
    rule_db.open(path_rule);
    vocab.open(path_vocab);
    
    if (has_features) {
      while (! feature_data_type::exists(path_feature_data)) {
	::sync();
	boost::thread::yield();
      }
      
      feature_data.open(path_feature_data);
      feature_vocab.open(path_feature_vocab);
    }
    
    if (has_attributes) {
      while (! attribute_data_type::exists(path_attribute_data)) {
	::sync();
	boost::thread::yield();
      }
      
      attribute_data.open(path_attribute_data);
      attribute_vocab.open(path_attribute_vocab);
    }
    
    utils::resource index_end;
    
    if (debug)
      std::cerr << "indexing"
		<< " cpu time: " << index_end.cpu_time() - index_start.cpu_time()
		<< " user time: " << index_end.user_time() - index_start.user_time()
		<< std::endl;
  }

  void GrammarStaticImpl::read_text(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const path_type path = param.name();
    
    std::string feature_prefix;
    std::string attribute_prefix;
    {
      parameter_type::const_iterator piter = param.find("feature-prefix");
      if (piter != param.end())
	feature_prefix = piter->second;
      parameter_type::const_iterator aiter = param.find("attribute-prefix");
      if (aiter != param.end())
	attribute_prefix = aiter->second;
    }

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
    std::unique_ptr<phrase_map_type> source_map(new phrase_map_type(1024 * 1024 * 128));
    std::unique_ptr<phrase_map_type> target_map(new phrase_map_type(1024 * 1024 * 128));
    
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

    utils::resource read_start;

    size_t num_line = 0;
    for (/**/; utils::getline(is, line); ++ num_line) {
      if (line.empty()) continue;

      boost::fusion::get<0>(rule).clear();
      boost::fusion::get<1>(rule).clear();
      boost::fusion::get<2>(rule).clear();
      boost::fusion::get<3>(rule).clear();
      boost::fusion::get<4>(rule).clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();
      
      const bool result = boost::spirit::qi::phrase_parse(iter, iter_end, rule_parser, boost::spirit::standard::space, rule);

      if (debug) {
	if ((num_line + 1) % DEBUG_DOT == 0)
	  std::cerr << '.';
	if ((num_line + 1) % DEBUG_LINE == 0)
	  std::cerr << std::endl;
      } 
      
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
      
      for (int feature = 0; feature < feature_size; ++ feature) {
	if (boost::fusion::get<3>(rule)[feature] <= boost::numeric::bounds<score_type>::lowest()) 
	  throw std::runtime_error("we are unable to handle such a lower feature value");
	
	score_streams[feature].ostream->write((char*) &boost::fusion::get<3>(rule)[feature], sizeof(score_type));
      }
      
      for (int attribute = 0; attribute < attribute_size; ++ attribute) {
	if (boost::fusion::get<4>(rule)[attribute] <= boost::numeric::bounds<score_type>::lowest())
	  throw std::runtime_error("we are unable to handle such a lower attribute value");

	attr_streams[attribute].ostream->write((char*) &boost::fusion::get<4>(rule)[attribute], sizeof(score_type));
      }
      
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
    
    utils::resource read_end;

    if (debug) {
      if ((num_line / DEBUG_DOT) % DEBUG_WRAP)
	std::cerr << std::endl;
      
      std::cerr << "# of rules: " << num_line
		<< " cpu time: " << read_end.cpu_time() - read_start.cpu_time()
		<< " user time: " << read_end.user_time() - read_start.user_time()
		<< std::endl;
    }

    utils::resource index_start;
    
    source_map->prune(static_cast<const hasher_type&>(*this));
    source_map->write(path_source);
    source_map.reset();
    
    target_map->prune(static_cast<const hasher_type&>(*this));
    target_map->write(path_target);
    target_map.reset();
    
    rule_db.close();
    
    word_type::write(path_vocab);
    
    while (! phrase_db_type::exists(path_source)) {
      ::sync();
      boost::thread::yield();
    }
    while (! phrase_db_type::exists(path_target)) {
      ::sync();
      boost::thread::yield();
    }
    while (! vocab_type::exists(path_vocab)) {
      ::sync();
      boost::thread::yield();
    }
    while (! rule_db_type::exists(path_rule)) {
      ::sync();
      boost::thread::yield();
    }
    
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
      
      while (! score_set_type::score_set_type::exists(score_streams[feature].path)) {
	::sync();
	boost::thread::yield();
      }
      
      utils::tempfile::permission(score_streams[feature].path);
      score_db[feature].score.open(score_streams[feature].path);

      const std::string name("feature" + utils::lexical_cast<std::string>(feature));

      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	feature_names[feature] = feature_type(piter->second);
      
      // default name...!
      if (feature_names[feature] == feature_type())
	feature_names[feature] = feature_prefix + "rule-table-" + utils::lexical_cast<std::string>(feature);
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
      
      while (! score_set_type::score_set_type::exists(attr_streams[attribute].path)) {
	::sync();
	boost::thread::yield();
      }
      
      utils::tempfile::permission(attr_streams[attribute].path);
      attr_db[attribute].score.open(attr_streams[attribute].path);
      
      const std::string name(std::string("attribute") + utils::lexical_cast<std::string>(attribute));
      
      parameter_type::const_iterator piter = param.find(name);
      if (piter != param.end())
	attribute_names[attribute] = attribute_type(piter->second);
      
      // default name...!
      if (attribute_names[attribute] == attribute_type())
	attribute_names[attribute] = attribute_prefix + "rule-table-" + utils::lexical_cast<std::string>(attribute);
    }

    utils::resource index_end;
    
    if (debug)
      std::cerr << "indexing"
		<< " cpu time: " << index_end.cpu_time() - index_start.cpu_time()
		<< " user time: " << index_end.user_time() - index_start.user_time()
		<< std::endl;
    
    // perform binarization, if possible!
    binarize();
  }
  
  
  GrammarStatic::GrammarStatic(const std::string& parameter)
    : pimpl(new impl_type(parameter)) {}

  GrammarStatic::~GrammarStatic() { std::unique_ptr<impl_type> tmp(pimpl); }

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
