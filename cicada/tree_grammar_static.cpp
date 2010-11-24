
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "tree_grammar_static.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "succinct_db/succinct_trie_database.hpp"
#include "succinct_db/succinct_trie_db.hpp"
#include "succinct_db/succinct_hash.hpp"

#include "utils/bithack.hpp"
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

#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/tuple/tuple.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>

#include <google/dense_hash_set>

namespace cicada
{
  struct TreeGrammarStaticImpl : public utils::hashmurmur<uint64_t>
  {
    friend class TreeGrammarStatic;

    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol  symbol_type;
    typedef cicada::Symbol  word_type;
    typedef cicada::Feature feature_type;
    typedef cicada::Vocab   vocab_type;
    
    
    typedef TreeTransducer::rule_type          rule_type;
    typedef TreeTransducer::rule_ptr_type      rule_ptr_type;
    typedef TreeTransducer::rule_pair_type     rule_pair_type;
    typedef TreeTransducer::rule_pair_set_type rule_pair_set_type;
    
    typedef TreeTransducer::feature_set_type feature_set_type;
    
    
    typedef float          score_type;
    typedef uint8_t        quantized_type;
    
    // this id_type is the same as the succinct-hash's pos_type...
    typedef uint32_t       id_type;
    
    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef char byte_type;
    typedef char mapped_type;
    
    typedef std::allocator<std::pair<id_type, mapped_type> > rule_alloc_type;
    
    typedef succinctdb::succinct_hash_mapped<byte_type, std::allocator<byte_type> >               rule_db_type;
    
    typedef succinctdb::succinct_trie_db<word_type::id_type, id_type, std::allocator<std::pair<word_type::id_type, id_type> > > edge_db_type;
    typedef succinctdb::succinct_trie_database<id_type, mapped_type, rule_alloc_type > rule_pair_db_type;
    
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
    
    typedef utils::arc_list<size_type, rule_pair_set_type, 16,
			    std::equal_to<size_type>,
			    std::allocator<std::pair<size_type, rule_pair_set_type> > > cache_rule_pair_set_type;
    
    struct cache_rule_type
    {
      rule_ptr_type rule;
      size_type pos;
      
      cache_rule_type() : rule(), pos(size_type(-1)) {}
    };
    
    typedef utils::array_power2<cache_rule_pair_set_type, 1024 * 16, std::allocator<cache_rule_pair_set_type> > cache_rule_pair_map_type;
    typedef utils::array_power2<cache_rule_type,          1024 * 16, std::allocator<cache_rule_type> >          cache_rule_set_type;

    TreeGrammarStaticImpl(const std::string& parameter) { read(parameter); }
    TreeGrammarStaticImpl(const TreeGrammarStaticImpl& x)
      : edge_db(x.edge_db), 
	rule_db(x.rule_db),
	source_db(x.source_db),
	target_db(x.target_db),
	score_db(x.score_db),
	vocab(x.vocab),
	feature_names(x.feature_names) {}

    TreeGrammarStaticImpl& operator=(const TreeGrammarStaticImpl& x)
    {
      clear();
      
      edge_db       = x.edge_db;
      rule_db       = x.rule_db;
      source_db     = x.source_db;
      target_db     = x.target_db;
      score_db      = x.score_db;
      vocab         = x.vocab;
      feature_names = x.feature_names;
      
      return *this;
    }
    
    void clear()
    {
      edge_db.clear();
      rule_db.clear();
      source_db.clear();
      target_db.clear();
      score_db.clear();
      vocab.clear();
      feature_names.clear();

      cache_rule.clear();
      cache_source.clear();
      cache_target.clear();
    }

    size_type find_edge(const word_type& word) const
    {
      size_type node = 0;
      return find_edge(word, node);
    }
    
    size_type find_edge(const word_type& word, size_type node) const
    {
      const word_type::id_type id = vocab[word];
      return edge_db.find(&id, 1, node);
    }
    
    size_type find(const id_type& id) const
    {
      size_type node = 0;
      return find(id, node);
    }
    
    size_type find(const id_type& id, size_type node) const
    {
      return rule_db.find(&id, 1, node);
    }
    
    bool is_valid_edge(size_type node) const { return edge_db.is_valid(node); }
    bool has_children_edge(size_type node) const { return edge_db.has_children(node); }
    bool exists_edge(size_type node) const { return edge_db.exists(node); }
    
    bool is_valid(size_type node) const { return rule_db.is_valid(node); }
    bool has_children(size_type node) const { return rule_db.has_children(node); }
    bool exists(size_type node) const { return rule_db.exists(node); }
    
        struct rule_unique_hash 
    {
      size_t operator()(const rule_ptr_type& x) const
      {
	return (x.get() ? hash_value(*x) : size_t(0));
      }
    };

    struct rule_unique_equal
    {
      bool operator()(const rule_ptr_type& x, const rule_ptr_type& y) const
      {
	return x == y || (x.get() && y.get() && *x == *y);
      }
    };
    typedef google::dense_hash_set<rule_ptr_type, rule_unique_hash, rule_unique_equal > rule_unique_map_type;

    const rule_pair_set_type& read_rule_set(size_type node) const
    {
      const size_type cache_pos = hasher_type::operator()(node) & (cache_rule.size() - 1);
      
      cache_rule_pair_set_type& cache = const_cast<cache_rule_pair_set_type&>(cache_rule[cache_pos]);
      
      std::pair<cache_rule_pair_set_type::iterator, bool> result = cache.find(node);
      if (! result.second) {
	typedef std::vector<byte_type, std::allocator<byte_type> >  code_set_type;

	rule_unique_map_type rules_unique;
	rules_unique.set_empty_key(rule_ptr_type());
	
	rule_pair_set_type& options = result.first->second;
	options.clear();
	
	rule_pair_db_type::cursor cursor_end = rule_db.cend(node);
	for (rule_pair_db_type::cursor cursor = rule_db.cbegin(node); cursor != cursor_end; ++ cursor) {
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

	  const rule_ptr_type rule_source = read_rule(pos_source, cache_source, source_db);
	  
	  while (citer != citer_end) {
	     const size_type offset_target = utils::group_aligned_decode(pos_target, &(*hiter), code_pos);
	     citer = hiter + offset_target;
	     hiter += offset_target & (- size_type((code_pos & 0x03) == 0x03));
	     ++ code_pos;
	     
	     const rule_ptr_type rule_target = read_rule(pos_target, cache_target, target_db);
	     
	     rule_ptr_type rule_sorted_source(new rule_type(*rule_source));
	     rule_ptr_type rule_sorted_target(new rule_type(*rule_target));
	     
	     cicada::sort(*rule_sorted_source, *rule_sorted_target);
	     
	     if (*rule_source == *rule_sorted_source)
	       rule_sorted_source = rule_source;
	     else
	       rule_sorted_source = *(rules_unique.insert(rule_sorted_source).first);
	     
	     if (*rule_target == *rule_sorted_target)
	      rule_sorted_target = rule_target;
	    else
	      rule_sorted_target = *(rules_unique.insert(rule_sorted_target).first);
	     
	     options.push_back(rule_pair_type(rule_sorted_source, rule_sorted_target));
	     
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
	     
	     ++ pos_feature;
	  }
	}
      }
      
      return result.first->second;
    }

  private:
    
    const rule_ptr_type& read_rule(size_type pos,
				   const cache_rule_set_type& caches,
				   const rule_db_type& db) const
    {
      const size_type cache_pos = hasher_type::operator()(pos) & (caches.size() - 1);
      
      cache_rule_type& cache = const_cast<cache_rule_type&>(caches[cache_pos]);
      if (cache.pos != pos) {
	cache.pos = pos;
	cache.rule.reset(new rule_type(std::string(db[pos].begin(), db[pos].end())));
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

    void read_text(const std::string& path);
    void read_binary(const path_type& path);
    
  private:
    edge_db_type          edge_db;
    rule_pair_db_type     rule_db;
    rule_db_type          source_db;
    rule_db_type          target_db;
    score_db_type         score_db;
    vocab_type            vocab;
    feature_name_set_type feature_names;
    
    // caching..
    cache_rule_pair_map_type cache_rule;
    cache_rule_set_type      cache_source;
    cache_rule_set_type      cache_target;
  };


  void TreeGrammarStaticImpl::ScoreSet::read(const path_type& path)
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
  
  void TreeGrammarStaticImpl::ScoreSet::write(const path_type& file) const
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

  void TreeGrammarStaticImpl::quantize()
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
  
  
  void TreeGrammarStaticImpl::read(const std::string& parameter)
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
  }

  
  void TreeGrammarStaticImpl::write(const path_type& file) const
  {
    typedef utils::repository repository_type;
    
    if (file == path()) return;
    
    repository_type rep(file, repository_type::write);
    
    edge_db.write(rep.path("edge"));
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
  
  
  void TreeGrammarStaticImpl::read_binary(const path_type& path)
  {
    typedef utils::repository repository_type;
	
    repository_type rep(path, repository_type::read);
    
    edge_db.open(rep.path("edge"));
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
  }

  struct TreeScoreSetStream
  {
    typedef boost::filesystem::path path_type;
    typedef utils::compress_ostream ostream_type;
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    
    TreeScoreSetStream() : ostream() {}
    TreeScoreSetStream(const path_type& _path)
      : ostream(new ostream_type(_path, 1024 * 1024)), path(_path) {}
    
    ostream_ptr_type ostream;
    path_type        path;
  };

  inline
  void encode_rule(const TreeRule& rule, std::string& buffer)
  {
    buffer.clear();
    
    boost::iostreams::filtering_ostream os;
    os.push(boost::iostreams::back_inserter(buffer));
    
    os << rule;
  }
  
  template <typename Options, typename Codes, typename Id>
  inline
  void encode_tree_options(const Options& options, Codes& codes, const Id& id_source)
  {
    typedef Id id_type;
    
    codes.clear();
    codes.resize(options.size() * 8 + 16, 0);
    
    typename Codes::iterator hiter = codes.begin();
    typename Codes::iterator citer = codes.begin();
    size_t pos = 0;
    
    const id_type id_feature = options.front().first;
    
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
      const id_type& id_target = piter->second;
      
      const size_t offset_target = utils::group_aligned_encode(id_target, &(*hiter), pos);
      citer = hiter + offset_target;
      hiter += offset_target & (- size_t((pos & 0x03) == 0x03));
      ++ pos;
    }
    
    codes.resize(citer - codes.begin());
  }

  template <typename Path, typename Buffer, typename EdgeDb>
  inline
  void encode_path(const Path& path, Buffer& buffer, EdgeDb& edge_db)
  {
    typedef std::vector<char, std::allocator<char> > codes_type;
    typedef typename codes_type::size_type size_type;
    
    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;

    buffer.clear();
    codes_type codes;

    char buf[8];
    const size_type buf_size = utils::byte_aligned_encode(Vocab::NONE.id(), buf);
    const typename Buffer::value_type id_none = edge_db.insert(buf, buf_size, hasher_type()(buf, buf + buf_size, 0));
    
    typename Path::const_iterator piter_end = path.end();
    for (typename Path::const_iterator piter = path.begin(); piter != piter_end; ++ piter) {
      typedef typename Path::value_type node_type;
      
      codes.clear();
      codes.reserve(piter->size() * 8);
      
      typename node_type::const_iterator niter_end = piter->end();
      for (typename node_type::const_iterator niter = piter->begin(); niter != niter_end; ++ niter) {
	if (*niter == Vocab::NONE) {
	  buffer.push_back(edge_db.insert(&(*codes.begin()), codes.size(), hasher_type()(codes.begin(), codes.end(), 0)));
	  codes.clear();
	} else {
	  const size_type buf_size = utils::byte_aligned_encode(niter->id(), buf);
	  codes.insert(codes.end(), buf, buf + buf_size);
	}
      }
      
      buffer.push_back(id_none);
    }
  }
  
  
  // how do we store single-tree...? use raw string? perform compression? (quicklz etc. ... currently NO...)
  
  void TreeGrammarStaticImpl::read_text(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    const path_type path = param.name();
    
    typedef succinctdb::succinct_hash<byte_type, std::allocator<byte_type> > rule_map_type;
    typedef succinctdb::succinct_hash<byte_type, std::allocator<byte_type> > edge_map_type;
    
    typedef TreeScoreSetStream score_stream_type;
    typedef std::vector<score_stream_type, std::allocator<score_stream_type> > score_stream_set_type;
    
    // feature-id, lhs-id, target-id
    typedef std::pair<id_type, id_type > rule_pair_option_type;
    typedef std::vector<rule_pair_option_type, std::allocator<rule_pair_option_type> > rule_pair_option_set_type;
    
    typedef std::vector<byte_type, std::allocator<byte_type> >  codes_type;
    typedef std::vector<score_type, std::allocator<score_type> > scores_type;
    typedef std::vector<id_type, std::allocator<id_type> > index_type;
    typedef std::vector<word_type, std::allocator<word_type> > node_type;
    typedef std::vector<node_type, std::allocator<node_type> > hyperpath_type;
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no file? ") + path.file_string());
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    
    const path_type path_edge   = utils::tempfile::directory_name(tmp_dir / "cicada.edge.XXXXXX");
    const path_type path_rule   = utils::tempfile::directory_name(tmp_dir / "cicada.rule.XXXXXX");
    const path_type path_source = utils::tempfile::directory_name(tmp_dir / "cicada.source.XXXXXX");
    const path_type path_target = utils::tempfile::directory_name(tmp_dir / "cicada.target.XXXXXX");
    const path_type path_vocab  = utils::tempfile::directory_name(tmp_dir / "cicada.vocab.XXXXXX");
    
    utils::tempfile::insert(path_edge);
    utils::tempfile::insert(path_rule);
    utils::tempfile::insert(path_source);
    utils::tempfile::insert(path_target);
    utils::tempfile::insert(path_vocab);
    
    rule_db.open(path_rule, rule_pair_db_type::WRITE);
    edge_map_type edge_map(1024 * 1024 * 4);
    
    rule_map_type source_map(1024 * 1024 * 4);
    rule_map_type target_map(1024 * 1024 * 4);
    
    score_stream_set_type score_streams;
    
    id_type id_rule = 0;
    
    int feature_size = -1;

    rule_type   source_prev;
    rule_type   source;
    rule_type   target;
    scores_type scores;
    
    std::string buffer_source;
    std::string buffer_target;
    codes_type  buffer_options;
    index_type  buffer_index;
    hyperpath_type hyperpath;
    
    rule_pair_option_set_type rule_options;
    
    utils::compress_istream is(path, 1024 * 1024);
    std::string line;

    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;
    
    qi::rule<std::string::const_iterator, scores_type(), standard::space_type> scores_parser;
    
    scores_parser %= +qi::float_;
    
    while (std::getline(is, line)) {
      source.clear();
      target.clear();
      scores.clear();

      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();
      
      if (! source.assign(iter, iter_end)) continue;
      if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) continue;
      if (! target.assign(iter, iter_end)) continue;
      if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) continue;
      if (! qi::phrase_parse(iter, iter_end, scores_parser, standard::space, scores)) continue;
      if (iter != iter_end) continue;
      
      if (source != source_prev) {
	
	if (! rule_options.empty()) {
	  encode_rule(source_prev, buffer_source);
	  
	  const id_type id_source = source_map.insert(buffer_source.c_str(), buffer_source.size(),
						      hasher_type::operator()(buffer_source.begin(), buffer_source.end(), 0));
	  
	  encode_tree_options(rule_options, buffer_options, id_source);
	  
	  source_prev.hyperpath(hyperpath);
	  
	  encode_path(hyperpath, buffer_index, edge_map);
	  
	  rule_db.insert(&(*buffer_index.begin()), buffer_index.size(), &(*buffer_options.begin()), buffer_options.size());
	}
	
	rule_options.clear();
	source_prev = source;
      }
      
      // encode features...
      if (feature_size < 0) {
	feature_size = scores.size();
	
	score_streams.reserve(feature_size);
	score_streams.resize(feature_size);
	
	for (int feature = 0; feature < feature_size; ++ feature) {
	  score_streams[feature].path = utils::tempfile::file_name(tmp_dir / "cicada.feature.XXXXXX");
	  utils::tempfile::insert(score_streams[feature].path);
	  
	  score_streams[feature].ostream.reset(new utils::compress_ostream(score_streams[feature].path, 1024 * 1024));
	}
      } else if (feature_size != static_cast<int>(scores.size()))
	throw std::runtime_error("invalid # of features...");
      
      for (int feature = 0; feature < feature_size; ++ feature)
	score_streams[feature].ostream->write((char*) &scores[feature], sizeof(score_type));
      
      encode_rule(target, buffer_target);
      
      const id_type id_target = target_map.insert(buffer_target.c_str(), buffer_target.size(),
						  hasher_type::operator()(buffer_target.begin(), buffer_target.end(), 0));
      
      rule_options.push_back(std::make_pair(id_rule ++, id_target));
    }
    
    if (! rule_options.empty()) {
      encode_rule(source_prev, buffer_source);
	  
      const id_type id_source = source_map.insert(buffer_source.c_str(), buffer_source.size(),
						  hasher_type::operator()(buffer_source.begin(), buffer_source.end(), 0));
      
      encode_tree_options(rule_options, buffer_options, id_source);
            
      source_prev.hyperpath(hyperpath);
      
      encode_path(hyperpath, buffer_index, edge_map);
      
      rule_db.insert(&(*buffer_index.begin()), buffer_index.size(), &(*buffer_options.begin()), buffer_options.size());
    }
    
    // source trees...
    source_map.write(path_source);
    source_map.clear();
    source_db.open(path_source);

    // target trees...
    target_map.write(path_target);
    target_map.clear();
    target_db.open(path_target);
    
    // rules....
    rule_db.close();
    rule_db.open(path_rule);
    
    // uncover edge_db from edge_map...
    {
      typedef std::vector<char, std::allocator<char> > codes_type;
      typedef std::vector<word_type::id_type, std::allocator<word_type::id_type> > buffer_type;
      
      edge_db.open(path_edge, edge_db_type::WRITE);

      codes_type codes;
      buffer_type buffer;
      for (id_type id = 0; id != edge_map.size(); ++ id) {
	codes.clear();
	codes.insert(codes.end(), edge_map[id].begin(), edge_map[id].end());
	
	word_type::id_type word_id;

	buffer.clear();
	codes_type::const_iterator citer = codes.begin();
	codes_type::const_iterator citer_end = codes.end();
	
	while (citer != citer_end) {
	  citer += utils::byte_aligned_decode(word_id, &(*citer));
	  buffer.push_back(word_id);
	}
	
	edge_db.insert(&(*buffer.begin()), buffer.size(), id);
      }
      
      edge_db.close();
      edge_db.open(path_edge);
    }
    
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
  

  TreeGrammarStatic::TreeGrammarStatic(const std::string& parameter)
    : pimpl(new impl_type(parameter)) {}
  
  TreeGrammarStatic::~TreeGrammarStatic() { std::auto_ptr<impl_type> tmp(pimpl); }
  
  TreeGrammarStatic::TreeGrammarStatic(const TreeGrammarStatic& x)
    : pimpl(new impl_type(*x.pimpl)) {}

  TreeGrammarStatic& TreeGrammarStatic::operator=(const TreeGrammarStatic& x)
  {
    *pimpl = *x.pimpl;
    return *this;
  }

  TreeGrammarStatic::transducer_ptr_type TreeGrammarStatic::clone() const
  {
    return transducer_ptr_type(new TreeGrammarStatic(*this));
  }
  
  TreeGrammarStatic::edge_id_type TreeGrammarStatic::edge(const symbol_type& symbol) const
  {
    const id_type node = pimpl->find_edge(symbol.non_terminal(), 0);
    
    return (pimpl->is_valid_edge(node) && pimpl->exists_edge(node) ? pimpl->edge_db[node] : edge_id_type(-1));
  }

  
  
  TreeGrammarStatic::edge_id_type TreeGrammarStatic::edge(const symbol_set_type& symbols) const
  {
    return edge(&(*symbols.begin()), &(*symbols.end()));
  }
  
  TreeGrammarStatic::edge_id_type TreeGrammarStatic::edge(const symbol_type* first, const symbol_type* last) const
  {
    id_type node = 0;
    for (/**/; first != last; ++ first) {
      node = pimpl->find_edge(first->non_terminal(), node);
      
      if (! pimpl->is_valid_edge(node)) break;
    }
    
    return (pimpl->is_valid_edge(node) && pimpl->exists_edge(node) ? pimpl->edge_db[node] : edge_id_type(-1));
  }


  TreeGrammarStatic::id_type TreeGrammarStatic::root() const
  {
    return 0;
  }
  
  TreeGrammarStatic::id_type TreeGrammarStatic::next(const id_type& node, const edge_id_type& edge) const
  {
    const impl_type::size_type pos = pimpl->find(edge, node);
    
    return (pimpl->is_valid(pos) ? id_type(pos) : id_type(0));
  }
  
  bool TreeGrammarStatic::has_next(const id_type& node) const
  {
    return (pimpl->is_valid(node) && pimpl->has_children(node));
  }
  
  const TreeGrammarStatic::rule_pair_set_type& TreeGrammarStatic::rules(const id_type& node) const
  {
    static const rule_pair_set_type __empty;
    
    return (pimpl->is_valid(node) && pimpl->exists(node) ? pimpl->read_rule_set(node) : __empty);
  }

  
  
  void TreeGrammarStatic::quantize()
  {
    pimpl->quantize();
  }
  
  void TreeGrammarStatic::write(const path_type& path) const
  {
    pimpl->write(path);
  }

};
