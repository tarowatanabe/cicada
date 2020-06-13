//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <memory>

#include "tree_grammar_mutable.hpp"
#include "parameter.hpp"
#include "quantizer.hpp"

#include "utils/trie_compact.hpp"
#include "utils/trie_set_compact.hpp"
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
#include "utils/thread_specific_ptr.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/json_string_parser.hpp"
#include "utils/resource.hpp"
#include "utils/getline.hpp"

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

namespace cicada
{
  static const size_t DEBUG_DOT = 1000000;
  static const size_t DEBUG_WRAP = 100;
  static const size_t DEBUG_LINE = DEBUG_DOT * DEBUG_WRAP;

  class TreeGrammarMutableImpl
  {
  public:
    friend class TreeGrammarMutable;

    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::Symbol    symbol_type;
    typedef cicada::Symbol    word_type;
    typedef cicada::Feature   feature_type;
    typedef cicada::Attribute attribute_type;
    typedef cicada::Vocab     vocab_type;
    
    
    typedef TreeTransducer::rule_type          rule_type;
    typedef TreeTransducer::rule_ptr_type      rule_ptr_type;
    typedef TreeTransducer::rule_pair_type     rule_pair_type;
    typedef TreeTransducer::rule_pair_set_type rule_pair_set_type;
    
    typedef TreeTransducer::feature_set_type   feature_set_type;
    typedef TreeTransducer::attribute_set_type attribute_set_type;
    
    typedef utils::trie_set_compact<symbol_type, utils::unassigned<symbol_type>, utils::unassigned<symbol_type>,
				    boost::hash<symbol_type>, std::equal_to<symbol_type>,
				    std::allocator<symbol_type > > edge_trie_type;
    typedef edge_trie_type::id_type edge_id_type;
    
    struct empty_key
    {
      edge_id_type operator()() const
      {
	return edge_id_type(-1);
      }
    };

    typedef utils::trie_compact<edge_id_type, rule_pair_set_type,
				empty_key, 
				boost::hash<edge_id_type>, std::equal_to<edge_id_type>,
				std::allocator<std::pair<const edge_id_type, rule_pair_set_type> > > trie_type;
    typedef trie_type::id_type id_type;
    
    typedef boost::filesystem::path path_type;
    
    typedef std::vector<feature_type, std::allocator<feature_type> >     feature_name_set_type;
    typedef std::vector<attribute_type, std::allocator<attribute_type> > attribute_name_set_type;
    
    TreeGrammarMutableImpl(const std::string& parameter)
      : trie(), edges(), cky(false), max_span(0), debug(0) { read(parameter); }

    TreeGrammarMutableImpl(const bool __cky=false, const int __max_span=0)
      : trie(), edges(), cky(__cky), max_span(__max_span), debug(0) {  }
    
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

    bool is_cky() const { return cky; }
    
  public:
    trie_type      trie;
    edge_trie_type edges;

    std::string feature_prefix;
    std::string attribute_prefix;

    feature_name_set_type   feature_names_default;
    attribute_name_set_type attribute_names_default;

    bool cky;
    int max_span;
    int debug;
  };
  
  //
  // we will parse only scores part and rely on tree-rule to parse partial string...
  //
  
  typedef std::pair<std::string, double> score_parsed_type;
  typedef std::vector<score_parsed_type, std::allocator<score_parsed_type> > scores_parsed_type;

  typedef std::pair<std::string, AttributeVector::data_type> attr_parsed_type;
  typedef std::vector<attr_parsed_type, std::allocator<attr_parsed_type> > attrs_parsed_type;

  typedef boost::fusion::tuple<scores_parsed_type, attrs_parsed_type> scores_attrs_parsed_type;
  
  template <typename Iterator>
  struct tree_rule_scores_parser_mutable : boost::spirit::qi::grammar<Iterator, scores_attrs_parsed_type(), boost::spirit::standard::space_type>
  {
    tree_rule_scores_parser_mutable() : tree_rule_scores_parser_mutable::base_type(scores_attrs)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      score  %= (qi::hold[qi::lexeme[+(!(qi::lit('=') >> qi::double_ >> (standard::space | qi::eoi)) >> (standard::char_ - standard::space))] >> '='] | qi::attr("")) >> qi::double_;
      scores %= -(score % (+standard::space));
      
      data %= data_string | double_dot | int64_;
      
      attribute %= (qi::hold[qi::lexeme[+(standard::char_ - standard::space - '=')] >> '='] | qi::attr("")) >> data;
      attributes %= *attribute;
      
      scores_attrs %= -("|||" >> scores) >> -("|||" >> attribute);
    }

    typedef boost::spirit::standard::space_type space_type;

    boost::spirit::qi::int_parser<AttributeVector::int_type, 10, 1, -1> int64_;
    boost::spirit::qi::real_parser<double, boost::spirit::qi::strict_real_policies<double> > double_dot;
    
    boost::spirit::qi::rule<Iterator, score_parsed_type()>  score;
    boost::spirit::qi::rule<Iterator, scores_parsed_type()> scores;

    utils::json_string_parser<Iterator> data_string;
    boost::spirit::qi::rule<Iterator, AttributeVector::data_type(), space_type> data;
    boost::spirit::qi::rule<Iterator, attr_parsed_type(), space_type>           attribute;
    boost::spirit::qi::rule<Iterator, attrs_parsed_type(), space_type>          attributes;

    boost::spirit::qi::rule<Iterator, scores_attrs_parsed_type(), space_type>   scores_attrs;
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
    attribute_name_set_type attribute_names;
    parameter_type::iterator piter_end = param.end();
    for (parameter_type::iterator piter = param.begin(); piter != piter_end; ++ piter) {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      if (utils::ipiece(piter->first) == "cky" || utils::ipiece(piter->first) == "cyk") {
	cky = true;
	continue;
      } else if (utils::ipiece(piter->first) == "max-span") {
	max_span = utils::lexical_cast<int>(piter->second);
	continue;
      } else if (utils::ipiece(piter->first) == "debug") {
	debug = utils::lexical_cast<int>(piter->second);
	continue;
      } else if (utils::ipiece(piter->first) == "feature-prefix") {
	feature_prefix = piter->second;
	continue;
      } else if (utils::ipiece(piter->first) == "attribute-prefix") {
	attribute_prefix = piter->second;
	continue;
      } else if (utils::ipiece(piter->first) == "populate")
	continue;

      {
	std::string::const_iterator iter = piter->first.begin();
	std::string::const_iterator iter_end = piter->first.end();
	
	int feature_id = -1;
	
	const bool result = qi::parse(iter, iter_end, "feature" >> qi::int_[phoenix::ref(feature_id) = qi::_1]);
	if (result && iter == iter_end && feature_id >= 0) {
	  if (feature_id >= int(feature_names.size()))
	    feature_names.resize(feature_id + 1);
	  
	  feature_names[feature_id] = piter->second;
	  continue;
	}
      }
      
      {
	std::string::const_iterator iter = piter->first.begin();
	std::string::const_iterator iter_end = piter->first.end();
	
	int attribute_id = -1;
	const bool result = qi::parse(iter, iter_end, "attribute" >> qi::int_[phoenix::ref(attribute_id) = qi::_1]);
	if (result && iter == iter_end && attribute_id >= 0) {
	  if (attribute_id >= int(attribute_names.size()))
	    attribute_names.resize(attribute_id + 1);
	  
	  attribute_names[attribute_id] = piter->second;
	  continue;
	}
      }
      
      throw std::runtime_error("unsupported key: " + piter->first);
    }

    typedef tree_rule_scores_parser_mutable<std::string::const_iterator> scores_parser_type;
    
#ifdef HAVE_TLS
    static __thread scores_parser_type* __scores_parser_tls = 0;
    static utils::thread_specific_ptr<scores_parser_type > __scores_parser;
    
    if (! __scores_parser_tls) {
      __scores_parser.reset(new scores_parser_type());
      __scores_parser_tls = __scores_parser.get();
    }
    
    scores_parser_type& scores_parser = *__scores_parser_tls;
#else
    static utils::thread_specific_ptr<scores_parser_type > __scores_parser;
    if (! __scores_parser.get())
      __scores_parser.reset(new scores_parser_type());
    
    scores_parser_type& scores_parser = *__scores_parser;
#endif
    
    
    utils::compress_istream is(path, 1024 * 1024);
    std::string line;

    rule_type                source;
    rule_type                target;
    scores_attrs_parsed_type scores_attrs;
    feature_set_type         features;
    attribute_set_type       attributes;

    utils::resource start;
    
    size_t num_line = 0;
    for (/**/; utils::getline(is, line); ++ num_line) {
      if (line.empty()) continue;
      
      source.clear();
      target.clear();
      boost::fusion::get<0>(scores_attrs).clear();
      boost::fusion::get<1>(scores_attrs).clear();
      
      std::string::const_iterator iter_end = line.end();
      std::string::const_iterator iter = line.begin();

      if (debug) {
	if ((num_line + 1) % DEBUG_DOT == 0)
	  std::cerr << '.';
	if ((num_line + 1) % DEBUG_LINE == 0)
	  std::cerr << std::endl;
      } 
      
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      if ((! source.assign(iter, iter_end))
	  || (! qi::phrase_parse(iter, iter_end, "|||", standard::space))
	  || (! target.assign(iter, iter_end))
	  || (! qi::phrase_parse(iter, iter_end, scores_parser, standard::space, scores_attrs))
	  || (iter != iter_end)) {
	std::cerr << "invalid line: " << num_line << ": " << line << std::endl;
	continue;
      }
      
      features.clear();
      int feature = 0;
      scores_parsed_type::const_iterator fiter_end = boost::fusion::get<0>(scores_attrs).end();
      for (scores_parsed_type::const_iterator fiter = boost::fusion::get<0>(scores_attrs).begin(); fiter != fiter_end; ++ fiter) {
	if (fiter->first.empty()) {
	  
	  if (feature < static_cast<int>(feature_names.size()) && ! feature_names[feature].empty())
	    features[feature_names[feature]] = fiter->second;
	  else {
	    // default name!
	    if (feature >= static_cast<int>(feature_names_default.size()))
	      feature_names_default.resize(feature + 1);
	    if (feature_names_default[feature].empty())
	      feature_names_default[feature] = feature_prefix + "tree-rule-table-" + utils::lexical_cast<std::string>(feature);
	    
	    features[feature_names_default[feature]] = fiter->second;
	  }
	  
	  ++ feature;
	} else
	  features[fiter->first] = fiter->second;
      }
      
      attributes.clear();
      int attribute = 0;
      attrs_parsed_type::const_iterator aiter_end = boost::fusion::get<1>(scores_attrs).end();
      for (attrs_parsed_type::const_iterator aiter = boost::fusion::get<1>(scores_attrs).begin(); aiter != aiter_end; ++ aiter) {
	if (aiter->first.empty()) {
	  
	  if (attribute < static_cast<int>(attribute_names.size()) && ! attribute_names[attribute].empty())
	    attributes[attribute_names[attribute]] = aiter->second;
	  else {
	    // default name!
	    if (attribute >= static_cast<int>(attribute_names_default.size()))
	      attribute_names_default.resize(attribute + 1);
	    if (attribute_names_default[attribute].empty())
	      attribute_names_default[attribute] = attribute_prefix + "tree-rule-table-" + utils::lexical_cast<std::string>(attribute);
	    
	    attributes[attribute_names_default[attribute]] = aiter->second;
	  }
	  
	  ++ attribute;
	} else
	  attributes[aiter->first] = aiter->second;
      }
      
      insert(rule_pair_type(rule_type::create(rule_type(source)), rule_type::create(rule_type(target)), features, attributes));
    }

    utils::resource end;

    if (debug) {
      if ((num_line / DEBUG_DOT) % DEBUG_WRAP)
	std::cerr << std::endl;

      std::cerr << "# of rules: " << num_line << std::endl;
      
      std::cerr << "indexing:"
		<< " cpu time: " << end.cpu_time() - start.cpu_time()
		<< " user time: " << end.user_time() - start.user_time()
		<< std::endl;
    }
  }
  

  void TreeGrammarMutableImpl::insert(const std::string& line)
  {
    typedef tree_rule_scores_parser_mutable<std::string::const_iterator> scores_parser_type;
    
#ifdef HAVE_TLS
    static __thread scores_parser_type* __scores_parser_tls = 0;
    static utils::thread_specific_ptr<scores_parser_type > __scores_parser;
    
    if (! __scores_parser_tls) {
      __scores_parser.reset(new scores_parser_type());
      __scores_parser_tls = __scores_parser.get();
    }
    
    scores_parser_type& scores_parser = *__scores_parser_tls;
#else
    static utils::thread_specific_ptr<scores_parser_type > __scores_parser;
    if (! __scores_parser.get())
      __scores_parser.reset(new scores_parser_type());
    
    scores_parser_type& scores_parser = *__scores_parser;
#endif

    if (line.empty()) return;
    
    rule_type            source;
    rule_type            target;
    scores_attrs_parsed_type scores_attrs;
    feature_set_type         features;
    attribute_set_type       attributes;
    
    std::string::const_iterator iter_end = line.end();
    std::string::const_iterator iter = line.begin();
    
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;

    if (! source.assign(iter, iter_end)) return;
    if (! qi::phrase_parse(iter, iter_end, "|||", standard::space)) return;
    if (! target.assign(iter, iter_end)) return;
    
    if (! qi::phrase_parse(iter, iter_end, scores_parser, standard::space, scores_attrs)) return;
    
    if (iter != iter_end) return;
    
    int feature = 0;
    scores_parsed_type::const_iterator fiter_end = boost::fusion::get<0>(scores_attrs).end();
    for (scores_parsed_type::const_iterator fiter = boost::fusion::get<0>(scores_attrs).begin(); fiter != fiter_end; ++ fiter)
      if (fiter->first.empty()) {
	// default name!
	if (feature >= static_cast<int>(feature_names_default.size()))
	  feature_names_default.resize(feature + 1);
	if (feature_names_default[feature].empty())
	  feature_names_default[feature] = feature_prefix + "tree-rule-table-" + utils::lexical_cast<std::string>(feature);
	
	features[feature_names_default[feature]] = fiter->second;
	
	++ feature;
      } else
	features[fiter->first] = fiter->second;
    
    int attribute = 0;
    attrs_parsed_type::const_iterator aiter_end = boost::fusion::get<1>(scores_attrs).end();
    for (attrs_parsed_type::const_iterator aiter = boost::fusion::get<1>(scores_attrs).begin(); aiter != aiter_end; ++ aiter)
      if (aiter->first.empty()) {
	// default name!
	if (attribute >= static_cast<int>(attribute_names_default.size()))
	  attribute_names_default.resize(attribute + 1);
	if (attribute_names_default[attribute].empty())
	  attribute_names_default[attribute] = attribute_prefix + "tree-rule-table-" + utils::lexical_cast<std::string>(attribute);
	
	attributes[attribute_names_default[attribute]] = aiter->second;
	
	++ attribute;
      } else
	attributes[aiter->first] = aiter->second;
    
    insert(rule_pair_type(rule_type::create(source), rule_type::create(target), features, attributes));
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
	  buffer.push_back(niter->non_terminal());
      }
      
      id = trie.insert(id, edge_none);
    }
    
    return id;
  }

  template <typename Trie>
  struct MutableFrontierIterator
  {
    typedef Trie trie_type;
    typedef typename trie_type::id_type id_type;
    
    MutableFrontierIterator(trie_type& __trie, id_type& __id) : trie(__trie), id(__id) {}
    
    template <typename Value>
    MutableFrontierIterator& operator=(const Value& value)
    {
      id = trie.insert(id, value.non_terminal().id());
      return *this;
    }
    
    MutableFrontierIterator& operator*() { return *this; }
    MutableFrontierIterator& operator++() { return *this; }
    
    trie_type& trie;
    id_type& id;
  };

  void TreeGrammarMutableImpl::insert(const rule_pair_type& rule_pair)
  {
    if (cky) {
      id_type id = trie.root();
      MutableFrontierIterator<trie_type> iter(trie, id);
      rule_pair.source->frontier(iter);
      
      trie[id].push_back(rule_pair);
    } else {
      typedef std::vector<symbol_type, std::allocator<symbol_type> > node_type;
      typedef std::vector<node_type, std::allocator<node_type> > hyperpath_type;
      
      hyperpath_type hyperpath;
      
      rule_pair.source->hyperpath(hyperpath);
      
      //std::cerr << "source: " << *(rule_pair.source) << std::endl;
      
      const id_type id = encode_path(hyperpath, edges, trie);
      
      trie[id].push_back(rule_pair);
    }
  }

  TreeGrammarMutable::TreeGrammarMutable(const bool __cky, const int __max_span)
    : pimpl(new impl_type(__cky, __max_span)) {}
  
  TreeGrammarMutable::TreeGrammarMutable(const std::string& parameter)
    : pimpl(new impl_type(parameter)) {}
  
  TreeGrammarMutable::~TreeGrammarMutable() { std::unique_ptr<impl_type> tmp(pimpl); }
  
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

  bool TreeGrammarMutable::valid_span(int first, int last, int distance) const
  {
    return pimpl->max_span <= 0 || distance <= pimpl->max_span || last - first == 1;
  }

  bool TreeGrammarMutable::is_cky() const
  {
    return pimpl->is_cky();
  }
  
  TreeGrammarMutable::edge_type TreeGrammarMutable::edge(const symbol_type& symbol) const
  {
    impl_type::edge_id_type node = pimpl->edges.find(pimpl->edges.root(), symbol.non_terminal());

    return (node != pimpl->edges.npos() ? node : edge_type::id_type(-1));
  }
  
  TreeGrammarMutable::edge_type TreeGrammarMutable::edge(const symbol_set_type& symbols) const
  {
    return edge(&(*symbols.begin()), &(*symbols.end()));
  }
  
  TreeGrammarMutable::edge_type TreeGrammarMutable::edge(const symbol_type* first, const symbol_type* last) const
  {
    impl_type::edge_id_type node = pimpl->edges.root();
    for (/**/; first != last; ++ first) {
      node = pimpl->edges.find(node, first->non_terminal());
      
      if (node == pimpl->edges.npos())
	return edge_type();
    }
    
    return (node != pimpl->edges.npos() ? node : edge_type::id_type(-1));
  }


  TreeGrammarMutable::id_type TreeGrammarMutable::root() const
  {
    return pimpl->root();
  }

  TreeGrammarMutable::id_type TreeGrammarMutable::next(const id_type& node, const symbol_type& symbol) const
  {
    if (! pimpl->is_cky())
      throw std::runtime_error("the tree grammar is not indexed for CKY");

    return pimpl->next(node, symbol.non_terminal().id());
  }
  
  
  TreeGrammarMutable::id_type TreeGrammarMutable::next(const id_type& node, const edge_type& edge) const
  {
    if (pimpl->is_cky())
      throw std::runtime_error("the tree grammar is indexed for CKY");

    return pimpl->next(node, edge.id);
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
