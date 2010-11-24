
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "tree_grammar_mutable.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/compact_trie_dense_set.hpp"
#include "utils/compress_stream.hpp"

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
  struct TreeGrammarMutableImpl : public utils::hashmurmur<uint64_t>
  {
    friend class TreeGrammarMutable;

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
    
    typedef utils::compact_trie_dense_set<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type>,
					  std::allocator<symbol_type > > edge_trie_type;
    typedef edge_trie_type::id_type edge_id_type;
    
    typedef utils::compact_trie_dense<edge_id_type, rule_pair_set_type, boost::hash<edge_id_type>, std::equal_to<edge_id_type>,
				      std::allocator<std::pair<const edge_id_type, rule_pair_set_type> > > trie_type;
    typedef trie_type::id_type id_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
    
    TreeGrammarMutableImpl(const std::string& parameter)
      : trie(edge_id_type(-1)), edges(symbol_type()) { read(parameter); }
    
    edge_id_type edge(const symbol_type* first, const symbol_type* last) const
    {
      return edges.find(first, last);
    }
    
    id_type root() const { return trie.root(); }
    id_type next(id_type node, const edge_id_type& edge) const { return trie.find(node, edge); }
    bool has_next(id_type node) const { return ! trie.empty(node); }
    const rule_pair_set_type& rules(id_type node) { return trie[node]; }
    
    void insert(const std::string& pattern);
    void insert(const rule_pair_type& rule);
    
  public:
    void clear() { trie.clear();  edges.clear(); }
    void read(const std::string& parameter);
    
  public:
    trie_type      trie;
    edge_trie_type edges;
  };
  
  //
  // we will parse only scores part and rely on tree-rule to parse partial string...
  //
  
  typedef std::pair<std::string, double> score_parsed_type;
  typedef std::vector<score_parsed_type, std::allocator<score_parsed_type> > scores_parsed_type;
  
  template <typename Iterator>
  struct tree_rule_scores_parser : boost::spirit::qi::grammar<Iterator, scores_parsed_type(), boost::spirit::standard::space_type>
  {
    tree_rule_scores_parser() : tree_rule_scores_parser::base_type(scores)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::lexeme;
      using qi::attr;
      using qi::hold;
      using standard::char_;
      using qi::double_;
      using standard::space;
      
      score  %= (hold[lexeme[+(char_ - space - '=')] >> '='] | attr("")) >> double_;
      scores %= +score;
    }
    
    boost::spirit::qi::rule<Iterator, score_parsed_type(), boost::spirit::standard::space_type>  score;
    boost::spirit::qi::rule<Iterator, scores_parsed_type(), boost::spirit::standard::space_type> scores;
  };

  
  void TreeGrammarMutableImpl::read(const std::string& parameter)
  {
    typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
    

    typedef cicada::Parameter parameter_type;
	  
    const parameter_type param(parameter);
    const path_type path = param.name();
    
    if (path != "-" && ! boost::filesystem::exists(path))
      throw std::runtime_error(std::string("no grammar file") + param.name());
    
    feature_name_set_type feature_names;
    parameter_type::iterator piter_end = param.end();
    for (parameter_type::iterator piter = param.begin(); piter != piter_end; ++ piter) {
      
      std::string::const_iterator iter = piter->first.begin();
      std::string::const_iterator iter_end = piter->first.end();
      
      int feature_id = -1;
      
      const bool result = boost::spirit::qi::parse(iter, iter_end,
						   "feature" >> boost::spirit::qi::int_[boost::phoenix::ref(feature_id) = boost::spirit::qi::_1]);
      if (result && iter == iter_end && feature_id >= 0) {
	if (feature_id >= int(feature_names.size()))
	  feature_names.resize(feature_id + 1);
	feature_names[feature_id] = piter->second;
      } else
	throw std::runtime_error("unsupported key: " + piter->first);
    }

    typedef tree_rule_scores_parser<std::string::const_iterator> scores_parser_type;
    
#ifdef HAVE_TLS
    static __thread scores_parser_type* __scores_parser_tls = 0;
    static boost::thread_specific_ptr<scores_parser_type > __scores_parser;
    
    if (! __scores_parser_tls) {
      __scores_parser.reset(new scores_parser_type());
      __scores_parser_tls = __scores_parser.get();
    }
    
    scores_parser_type& scores_parser = *__scores_parser_tls;
#else
    static boost::thread_specific_ptr<scores_parser_type > __scores_parser;
    if (! __scores_parser.get())
      __scores_parser.reset(new scores_parser_type());
    
    scores_parser_type& scores_parser = *__scores_parser;
#endif
    
    
    utils::compress_istream is(path, 1024 * 1024);
    std::string line;

    rule_type          source;
    rule_type          target;
    scores_parsed_type scores;
    feature_set_type   features;
    
    while (std::getline(is, line)) {
      if (line.empty()) continue;
      
      source.clear();
      target.clear();
      scores.clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();

      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      if (! source.assign(iter, iter_end)) continue;
      if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) continue;
      if (! target.assign(iter, iter_end)) continue;
      
      if (iter != iter_end && qi::phrase_parse(iter, iter_end, "|||", standard::space))
	if (! qi::phrase_parse(iter, iter_end, scores_parser, standard::space, scores)) continue;
      
      if (iter != iter_end) continue;
      
      features.clear();
      int feature = 0;
      scores_parsed_type::const_iterator fiter_end = scores.end();
      for (scores_parsed_type::const_iterator fiter = scores.begin(); fiter != fiter_end; ++ fiter) {
	if (fiter->first.empty()) {
	  
	  if (feature < static_cast<int>(feature_names.size()) && ! feature_names[feature].empty())
	    features[feature_names[feature]] = fiter->second;
	  else {
	    // default name!
	    const std::string name = std::string("rule-table-") + boost::lexical_cast<std::string>(feature);
	    
	    features[name] = fiter->second;
	  }
	  
	  ++ feature;
	} else
	  features[fiter->first] = fiter->second;
      }
      
      insert(rule_pair_type(rule_ptr_type(new rule_type(source)), rule_ptr_type(new rule_type(target)), features));
    }
  }
  

  void TreeGrammarMutableImpl::insert(const std::string& line)
  {
    typedef tree_rule_scores_parser<std::string::const_iterator> scores_parser_type;
    
#ifdef HAVE_TLS
    static __thread scores_parser_type* __scores_parser_tls = 0;
    static boost::thread_specific_ptr<scores_parser_type > __scores_parser;
    
    if (! __scores_parser_tls) {
      __scores_parser.reset(new scores_parser_type());
      __scores_parser_tls = __scores_parser.get();
    }
    
    scores_parser_type& scores_parser = *__scores_parser_tls;
#else
    static boost::thread_specific_ptr<scores_parser_type > __scores_parser;
    if (! __scores_parser.get())
      __scores_parser.reset(new scores_parser_type());
    
    scores_parser_type& scores_parser = *__scores_parser;
#endif

    if (line.empty()) return;
    
    rule_ptr_type      source(new rule_type());
    rule_ptr_type      target(new rule_type());
    scores_parsed_type scores;
    feature_set_type   features;
    
    std::string::const_iterator iter_end = line.end();
    std::string::const_iterator iter = line.begin();
    
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    namespace phoenix = boost::phoenix;

    if (! source->assign(iter, iter_end)) return;
    if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) return;
    if (! target->assign(iter, iter_end)) return;
    
    if (iter != iter_end && qi::phrase_parse(iter, iter_end, "|||", standard::space))
      if (! qi::phrase_parse(iter, iter_end, scores_parser, standard::space, scores)) return;
    
    if (iter != iter_end) return;
    
    int feature = 0;
    scores_parsed_type::const_iterator fiter_end = scores.end();
    for (scores_parsed_type::const_iterator fiter = scores.begin(); fiter != fiter_end; ++ fiter)
      if (fiter->first.empty()) {
	// default name!
	const std::string name = std::string("rule-table-") + boost::lexical_cast<std::string>(feature);
	
	features[name] = fiter->second;
	
	++ feature;
      } else
	features[fiter->first] = fiter->second;
    
    insert(rule_pair_type(source, target, features));
  }
  

  template <typename Path, typename Edges, typename Trie>
  inline
  TreeGrammarMutableImpl::id_type encode_path(const Path& path, Edges& edges, Trie& trie)
  {
    typedef cicada::Symbol symbol_type;
    typedef std::vector<symbol_type, std::allocator<symbol_type> > symbol_set_type;
    
    symbol_set_type buffer;
    typename Trie::id_type id = trie.root();

    const typename Edges::id_type edge_none = edges.insert(edges.root(), Vocab::NONE);
    
    typename Path::const_iterator piter_end = path.end();
    for (typename Path::const_iterator piter = path.begin(); piter != piter_end; ++ piter) {
      typedef typename Path::value_type node_type;

#if 0
      std::cerr << "hyperpath: ";
      std::copy(piter->begin(), piter->end(), std::ostream_iterator<symbol_type>(std::cerr, " "));
      std::cerr << std::endl;
#endif

      buffer.clear();
      typename node_type::const_iterator niter_end = piter->end();
      for (typename node_type::const_iterator niter = piter->begin(); niter != niter_end; ++ niter) {
	if (*niter == Vocab::NONE) {
	  id = trie.insert(id, edges.insert(buffer.begin(), buffer.end()));
	  buffer.clear();
	} else
	  buffer.push_back(*niter);
      }
      
      id = trie.insert(id, edge_none);
    }
    
    return id;
  }

  void TreeGrammarMutableImpl::insert(const rule_pair_type& rule_pair)
  {
    typedef std::vector<symbol_type, std::allocator<symbol_type> > node_type;
    typedef std::vector<node_type, std::allocator<node_type> > hyperpath_type;
    
    hyperpath_type hyperpath;


    rule_pair.source->hyperpath(hyperpath);

    //std::cerr << "source: " << *(rule_pair.source) << std::endl;
    
    const id_type id = encode_path(hyperpath, edges, trie);
    
    trie[id].push_back(rule_pair);
    
  }

  TreeGrammarMutable::TreeGrammarMutable(const std::string& parameter)
    : pimpl(new impl_type(parameter)) {}
  
  TreeGrammarMutable::~TreeGrammarMutable() { std::auto_ptr<impl_type> tmp(pimpl); }
  
  TreeGrammarMutable::TreeGrammarMutable(const TreeGrammarMutable& x)
    : pimpl(new impl_type(*x.pimpl)) {}

  TreeGrammarMutable& TreeGrammarMutable::operator=(const TreeGrammarMutable& x)
  {
    *pimpl = *x.pimpl;
    return *this;
  }
  
  TreeGrammarMutable::transducer_ptr_type TreeGrammarMutable::clone() const
  {
    return transducer_ptr_type(new TreeGrammarMutable(*this));
  }
  
  TreeGrammarMutable::edge_id_type TreeGrammarMutable::edge(const symbol_type& symbol) const
  {
    return edge(&symbol, (&symbol) + 1);
  }
  
  TreeGrammarMutable::edge_id_type TreeGrammarMutable::edge(const symbol_set_type& symbols) const
  {
    return edge(&(*symbols.begin()), &(*symbols.end()));
  }
  
  TreeGrammarMutable::edge_id_type TreeGrammarMutable::edge(const symbol_type* first, const symbol_type* last) const
  {
    impl_type::edge_id_type id = pimpl->edge(first, last);
    
    return (id != pimpl->edges.npos() ? id : edge_id_type(-1));
  }


  TreeGrammarMutable::id_type TreeGrammarMutable::root() const
  {
    return pimpl->root();
  }
  
  TreeGrammarMutable::id_type TreeGrammarMutable::next(const id_type& node, const edge_id_type& edge) const
  {
    return pimpl->next(node, edge);
  }
  
  bool TreeGrammarMutable::has_next(const id_type& node) const
  {
    return pimpl->has_next(node);
  }
  
  const TreeGrammarMutable::rule_pair_set_type& TreeGrammarMutable::rules(const id_type& node) const
  {
    static const rule_pair_set_type __empty;
    return (node == pimpl->root() ? __empty : pimpl->rules(node));
  }
  
   
  void TreeGrammarMutable::read(const std::string& parameter)
  {
    pimpl->read(parameter);
  }

  void TreeGrammarMutable::clear()
  {
    pimpl->clear();
  }
  
  void TreeGrammarMutable::insert(const std::string& pattern)
  {
    pimpl->insert(pattern);
  }

  void TreeGrammarMutable::insert(const rule_pair_type& rule_pair)
  {
    pimpl->insert(rule_pair);
  }
  
};
