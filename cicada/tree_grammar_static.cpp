
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "tree_grammar_static.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "succinct_db/succinct_trie_database.hpp"
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
  struct TreeGrammarStaticImpl : public utils::hashmurmur<uint64_t>
  {
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
    
    typedef rule_type::feature_set_type feature_set_type;
    
    
    typedef float          score_type;
    typedef uint8_t        quantized_type;
    
    typedef uint32_t       id_type;

    typedef uint64_t                           hash_value_type;
    typedef utils::hashmurmur<hash_value_type> hasher_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef char byte_type;
    typedef char mapped_type;
    
    typedef std::allocator<std::pair<word_type::id_type, mapped_type> > rule_alloc_type;
    
    typedef succinctdb::succinct_trie_database<word_type::id_type, mapped_type, rule_alloc_type > rule_pair_db_type;
    typedef succinctdb::succinct_hash_mapped<byte_type, std::allocator<byte_type> >               rule_db_type;
    
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
    
    typedef utils::array_power2<cache_rule_pair_type, 1024 * 16, std::allocator<cache_rule_pair_type> > cache_rule_pair_type;
    typedef utils::array_power2<cache_rule_type,      1024 * 16, std::allocator<cache_rule_type> >      cache_rule_type;

    TreeGrammarStaticImpl(const std::string& parameter) { read(parameter); }
    TreeGrammarStaticImpl(const TreeGrammarStaticImpl& x)
      : rule_db(x.rule_db),
	source_db(x.source_db),
	target_db(x.target_db),
	score_db(x.score_db),
	vocab(x.vocab),
	feature_names(x.feature_names) {}

    TreeGrammarStaticImpl& operator=(const TreeGrammarStaticImpl& x)
    {
      clear();
      
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
    rule_pair_db_type     rule_db;
    rule_db_type          source_db;
    rule_db_type          target_db;
    score_db_type         score_db;
    vocab_type            vocab;
    feature_name_set_type feature_names;
    
    // caching..
    cache_rule_pair_type cache_rule;
    cache_rule_type      cache_source;
    cache_rule_type      cache_target;
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
  
  inline
  int tree_depth(const TreeRule& tree, const int depth)
  {
    // pre-order traversal
    int max_depth = depth;
    for (TreeRule::const_iterator aiter = tree.begin(); aiter != tree.end(); ++ aiter)
      max_depth = utils::bithack::max(max_depth, tree_depth(*aiter, depth + 1));
    
    return max_depth;
  }

  inline
  void tree_add_epsilon(TreeRule& tree, const int max_depth, const int depth)
  {
    // pre-order traversal
    if (tree.empty()) {
      TreeRule* curr = &tree;
      
      for (int i = depth; i != max_depth; ++ i) {
	curr->antecedents = TreeRule::antecedent_set_type(1, TreeRule(Vocab::EPSILON));
	curr = &curr->front();
      }
    } else {
      for (TreeRule::const_iterator aiter = tree.begin(); aiter != tree.end(); ++ aiter)
	tree_add_epsilon(*aiter, max_depth, depth + 1);
    }
  }

  template <typename Path>
  inline
  void tree_to_hyperpath(const TreeRule& tree, Path& path, const int depth)
  {
    path[depth].push_back(tree.label);

    if (! tree.empty()) {
      for (TreeRule::const_iterator aiter = tree.begin(); aiter != tree.end(); ++ aiter)
	tree_to_hyperpath(*aiter, path, depth + 1);
      path[depth + 1].push_back(Vocab::NONE);
    }
  }
  
  template <typename Path>
  inline
  void tree_to_hyperpath(const TreeRule& rule, Path& path)
  {
    // compute max-depth by pre-order
    const int max_depth = tree_depth(tree, 0);
    
    // add epsilon annotation by pre-order
    TreeRule tree_epsilon(tree);
    tree_add_epsilon(tree_epsilon, max_depth, 0);
    
    // convert into hyperpath by pre-order
    path.clear();
    path.resize(max_depth + 1);
    
    tree_to_hyperpath(tree_epsilon, path, 0);
    path[0].push_back(Vocab::NONE);
  }
  
  // how do we store single-tree...? use raw string?
  
  void TreeGrammarStaticImpl::read_text(const std::string& parameter)
  {
    
    
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

  TreeGrammarStatic::id_type TreeGrammarStatic::root() const
  {
    return 0;
  }
  
  TreeGrammarStatic::id_type TreeGrammarStatic::next(const id_type& node, const symbol_type& symbol) const
  {
    const impl_type::size_type pos = pimpl->find(symbol.non_terminal(), node);
    
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
  
  void TreeGrammarStatic::write(const path_type& path)
  {
    pimpl->write(path);
  }

};
