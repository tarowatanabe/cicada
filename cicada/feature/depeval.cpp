//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>

#include "utils/lexical_cast.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/indexed_set.hpp"
#include "utils/bithack.hpp"
#include "utils/piece.hpp"

#include "cicada/feature/depeval.hpp"
#include "cicada/parameter.hpp"
#include "cicada/eval/depeval.hpp"
#include "cicada/semiring.hpp"
#include "cicada/sentence_vector.hpp"
#include "cicada/dependency.hpp"

#include "cicada/tokenizer.hpp"

#include <boost/numeric/conversion/bounds.hpp>

namespace cicada
{
  namespace feature
  {
    class DepevalImpl
    {
    public:
      typedef cicada::Symbol         symbol_type;
      typedef cicada::Vocab          vocab_type;
      typedef cicada::Sentence       sentence_type;
      typedef cicada::SentenceVector sentence_set_type;
      typedef cicada::Dependency     dependency_type;

      typedef std::vector<dependency_type, std::allocator<dependency_type> > dependency_set_type;

      typedef std::vector<dependency_set_type, std::allocator<dependency_set_type> > dependency_document_type;

      typedef cicada::eval::Score score_type;
      typedef score_type::score_ptr_type score_ptr_type;
      
      typedef cicada::FeatureFunction feature_function_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      typedef std::pair<int, int> count_type;
      typedef utils::simple_vector<count_type, std::allocator<count_type> > count_set_type;
      
      struct count_set_hash : public utils::hashmurmur<size_t>
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	size_t operator()(const count_set_type& x) const
	{
	  return hasher_type::operator()(x.begin(), x.end(), 0);
	}
      };
      
      typedef utils::indexed_set<count_set_type, count_set_hash, std::equal_to<count_set_type>, std::allocator<count_set_type> > counts_index_type;
      

      DepevalImpl()
	: attr_dependency_pos("dependency-pos"),
	  attr_dependency_head("dependency-head"),
	  attr_dependency_dependent("dependency-dependent") {}
      
      struct __attribute_integer : public boost::static_visitor<attribute_set_type::int_type>
      {
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -1; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -1; }
      };
      
      double depeval_score(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge)
      {
	if (deprefs.empty()) return 0.0;
	
	int pos_head = -1;
	int pos_dep  = -1;
	
	attribute_set_type::const_iterator hiter = edge.attributes.find(attr_dependency_head);
	attribute_set_type::const_iterator diter = edge.attributes.find(attr_dependency_dependent);
	
	if (hiter != edge.attributes.end() && diter != edge.attributes.end()) {
	  pos_head = boost::apply_visitor(__attribute_integer(), hiter->second);
	  pos_dep  = boost::apply_visitor(__attribute_integer(), diter->second);
	}
	
	count_set_type counts(deprefs.size(), count_type(0, 0));
	
	double score_antecedent = 0.0;
	
	if (pos_head >= 0 && pos_dep > 0)
	  for (size_t i = 0; i != deprefs.size(); ++ i)
	    if (pos_dep - 1 < static_cast<int>(deprefs[i].size())) {
	      counts[i].first += (pos_head == deprefs[i][pos_dep - 1]);
	      ++ counts[i].second;
	    }
	
	for (size_t j = 0; j != states.size(); ++ j) {
	  const count_set_type& counts_antecedent = counts_index[*reinterpret_cast<const counts_index_type::index_type*>(states[j])];
	  
	  // compute current score...
	  double score_max = 0.0;
	  for (size_t i = 0; i != deprefs.size(); ++ i)
	    if (counts_antecedent[i].first)
	      score_max = std::max(score_max, double(counts_antecedent[i].first) / double(counts_antecedent[i].second));
	  
	  score_antecedent += score_max;
	  
	  for (size_t i = 0; i != deprefs.size(); ++ i) {
	    counts[i].first  += counts_antecedent[i].first;
	    counts[i].second += counts_antecedent[i].second;
	  }
	}
	
	counts_index_type::iterator citer = counts_index.insert(counts).first;
	*reinterpret_cast<counts_index_type::index_type*>(state) = citer - counts_index.begin();
	
	// compute current score...
	double score_max = 0.0;
	for (size_t i = 0; i != deprefs.size(); ++ i)
	  if (counts[i].first)
	    score_max = std::max(score_max, double(counts[i].first) / double(counts[i].second));
	
	return score_max - score_antecedent;
      }
      
      void initialize() { counts_index.clear(); }
      void clear() { counts_index.clear(); deprefs.clear(); score.reset(); }
      
      void insert(const dependency_type& dep)
      {
	deprefs.push_back(dep);
      }
      
      void insert(const sentence_type& sent)
      {
	dependency_type dep;
	sentence_type::const_iterator siter_end = sent.end();
	for (sentence_type::const_iterator siter = sent.begin(); siter != siter_end; ++ siter)
	  dep.push_back(utils::lexical_cast<dependency_type::index_type>(static_cast<const std::string&>(*siter)));
	
	deprefs.push_back(dep);
      }
      
      void insert(const score_ptr_type& __score)
      {
	score = __score;
	
	if (score && ! dynamic_cast<const cicada::eval::Depeval*>(score.get()))
	  throw std::runtime_error("this is not a depeval-score!");
      }

      
      counts_index_type counts_index;
      
      dependency_set_type      deprefs;
      dependency_document_type refset;
      score_ptr_type score;
      
      attribute_type attr_dependency_pos;
      attribute_type attr_dependency_head;
      attribute_type attr_dependency_dependent;
    };
    
    
    Depeval::Depeval(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "depeval")
	throw std::runtime_error("this is not Depeval feature: " + parameter);

      // dummy tag/tokenizer
      bool skip_sgml_tag = false;
      const cicada::Tokenizer* tokenizer = 0;
      
      std::string name;
      path_type   refset_file;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "tokenizer")
	  tokenizer = &cicada::Tokenizer::create(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "refset")
	  refset_file = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for depeval: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! refset_file.empty() && ! boost::filesystem::exists(refset_file))
	throw std::runtime_error("no refset file?: " + refset_file.string());
      
      std::auto_ptr<impl_type> depeval_impl(new impl_type());
      
      // matched count + total count
      base_type::__state_size = sizeof(impl_type::counts_index_type::index_type);
      base_type::__feature_name = (name.empty() ? std::string("depeval") : name);
      
      pimpl = depeval_impl.release();
      
      pimpl->refset.clear();
      
      if (! refset_file.empty()) {
	typedef Dependency dependency_type;
	
	namespace qi = boost::spirit::qi;
	namespace standard = boost::spirit::standard;
	
	utils::compress_istream is(refset_file, 1024 * 1024);
	std::string line;
	
	size_t id;
	dependency_type dependency;
	
	qi::uint_parser<size_t> id_parser;
	
	while (std::getline(is, line)) {
	  std::string::const_iterator iter     = line.begin();
	  std::string::const_iterator iter_end = line.end();
	  
	  if (! qi::phrase_parse(iter, iter_end, id_parser >> "|||", standard::space, id)) continue;
	  
	  if (! dependency.assign(iter, iter_end)) continue;
	  
	  if (id >= pimpl->refset.size())
	    pimpl->refset.resize(id + 1);
	  
	  pimpl->refset[id].push_back(dependency);
	}
      }
    }
    
    Depeval::~Depeval() { std::auto_ptr<impl_type> tmp(pimpl); }

    Depeval::Depeval(const Depeval& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    Depeval& Depeval::operator=(const Depeval& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Depeval::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
    {
      double score = pimpl->depeval_score(state, states, edge);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
    }

    void Depeval::apply_coarse(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features,
			    feature_set_type& estimates,
			    const bool final) const
    {
      
    }
    void Depeval::apply_predict(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Depeval::apply_scan(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  const int dot,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {}
    void Depeval::apply_complete(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void Depeval::initialize()
    {
      pimpl->initialize();
    }
    
    void Depeval::clear()
    {
      pimpl->clear();
    }

    void Depeval::assign(const size_type& id,
			 const hypergraph_type& hypergraph,
			 const lattice_type& lattice,
			 const span_set_type& spans,
			 const sentence_set_type& targets,
			 const ngram_count_set_type& ngram_counts)
    {
      pimpl->clear();
      
      if (! targets.empty()) {
	sentence_set_type::const_iterator titer_end = targets.end();
	for (sentence_set_type::const_iterator titer = targets.begin(); titer != titer_end; ++ titer)
	  pimpl->insert(*titer);
      } else if (! pimpl->refset.empty()) {
	if (id < pimpl->refset.size()) {
	  impl_type::dependency_set_type::const_iterator titer_end = pimpl->refset[id].end();
	  for (impl_type::dependency_set_type::const_iterator titer = pimpl->refset[id].begin(); titer != titer_end; ++ titer)
	    pimpl->insert(*titer);
	}	
      } else
	throw std::runtime_error("no reference set?");
    }
    
    
    void Depeval::assign(const score_ptr_type& score)
    {
      pimpl->insert(score);
    }

  };
};
