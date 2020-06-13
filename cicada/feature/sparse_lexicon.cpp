//
//  Copyright(C) 2011-2013 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <set>

#include "sparse_lexicon.hpp"
#include "feature_builder.hpp"

#include "cicada/lexicon.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"
#include "cicada/feature_vector_linear.hpp"

#include "utils/piece.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/alloc_vector.hpp"
#include "utils/array_power2.hpp"
#include "utils/small_vector.hpp"
#include "utils/compact_set.hpp"
#include "utils/hashmurmur3.hpp"

namespace cicada
{
  namespace feature
  {


    class SparseLexiconImpl : public utils::hashmurmur3<size_t>
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Lattice       lattice_type;
      typedef cicada::HyperGraph    hypergraph_type;

      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;
      
      typedef cicada::ClusterStemmer normalizer_type;
      typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;

      typedef utils::hashmurmur3<size_t> hasher_type;

      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

      typedef FeatureVectorLinear<feature_set_type::mapped_type>    feature_linear_set_type;

      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::rule_type rule_type;
      
      typedef rule_type::symbol_set_type phrase_type;
      
      typedef symbol_type word_type;
      typedef std::pair<word_type, word_type> word_pair_type;
      
      typedef utils::compact_set<word_type,
				 utils::unassigned<word_type>, utils::unassigned<word_type>,
				 boost::hash<word_type>, std::equal_to<word_type>,
				 std::allocator<word_type> > word_unique_type;
      
      typedef std::vector<word_pair_type, std::allocator<word_pair_type> > word_pair_set_type;
      
      struct unassigned_key : utils::unassigned<word_type>
      {
	word_pair_type operator()() const 
	{
	  return word_pair_type(utils::unassigned<word_type>::operator()(),
				utils::unassigned<word_type>::operator()());
	}
      };

      typedef utils::compact_set<word_pair_type,
				 unassigned_key, unassigned_key,
				 utils::hashmurmur3<size_t>, std::equal_to<word_pair_type>,
				 std::allocator<word_pair_type> > word_pair_unique_type;

      typedef std::vector<word_type, std::allocator<word_type> > word_set_type;
      typedef std::vector<word_set_type, std::allocator<word_set_type> > word_map_type;
      typedef std::set<size_type, std::less<size_type>, std::allocator<size_type> > pos_set_type;
      
      typedef utils::alloc_vector<feature_linear_set_type, std::allocator<feature_linear_set_type> > cache_set_type;

      struct CacheNormalize
      {
	typedef utils::small_vector<word_type, std::allocator<word_type> > word_set_type;
	
	word_type::id_type word;
	word_set_type      normalized;
	
	CacheNormalize() : word(word_type::id_type(-1)), normalized() {}
      };
      typedef CacheNormalize cache_normalize_type;
      typedef utils::array_power2<cache_normalize_type, 1024 * 4, std::allocator<cache_normalize_type> > cache_normalize_set_type;

      typedef cicada::Lexicon lexicon_type;
      
      typedef FeatureBuilder feature_builder_type;
      
      SparseLexiconImpl()
	: lexicon(0), lexicon_prefix(0), lexicon_suffix(0), caches(),
	  skip_sgml_tag(false), unique_source(false), prefix("sparse-lexicon"), forced_feature(false),
	  pair_mode(false),
	  prefix_mode(false),
	  suffix_mode(false),
	  fertility_mode(false)
      { }

      struct skipper_epsilon
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS;
	}
      };

      struct skipper_sgml
      {
	bool operator()(const symbol_type& word) const
	{
	  return word == vocab_type::EPSILON || word == vocab_type::BOS || word == vocab_type::EOS || word.is_sgml_tag();
	}
      };
      
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features)
      {
	if (skip_sgml_tag)
	  lexicon_score(edge, features, skipper_sgml());
	else
	  lexicon_score(edge, features, skipper_epsilon());
      }
      
      template <typename Skipper>
      void lexicon_score(const edge_type& edge,
			 feature_set_type& features,
			 Skipper skipper)
      {
	const phrase_type& phrase = edge.rule->rhs;
	
	phrase_type::const_iterator titer_end = phrase.end();
	for (phrase_type::const_iterator titer = phrase.begin(); titer != titer_end; ++ titer) 
	  if (titer->is_terminal() && ! skipper(*titer)) {
	    const word_type& target = *titer;
	    
	    if (! caches.exists(titer->id())) {
	      feature_set_type features;
	      
	      size_t fertility_pair = 0;
	      size_t fertility_prefix = 0;
	      size_t fertility_suffix = 0;
	      word_set_type::const_iterator witer_end = words.end();
	      for (word_set_type::const_iterator witer = words.begin(); witer != witer_end; ++ witer) 
		if (exists(*witer, target)) {
		  apply(*witer, target, features);
		  ++ fertility_pair;
		}
	      
	      word_pair_set_type::const_iterator piter_end = words_prefix.end();
	      for (word_pair_set_type::const_iterator piter = words_prefix.begin(); piter != piter_end; ++ piter) 
		if (! lexicon_prefix || exists(*lexicon_prefix, *piter, target)) {
		  apply("-", *piter, target, features);
		  ++ fertility_prefix;
		}
	      
	      word_pair_set_type::const_iterator siter_end = words_suffix.end();
	      for (word_pair_set_type::const_iterator siter = words_suffix.begin(); siter != siter_end; ++ siter) 
		if (! lexicon_suffix || exists(*lexicon_suffix, *siter, target)) {
		  apply("+", *siter, target, features);
		  ++ fertility_suffix;
		}
	      
	      if (fertility_mode) {
		if (fertility_pair)
		  apply(":", target, fertility_pair, features);
		if (fertility_prefix)
		  apply("-", target, fertility_prefix, features);
		if (fertility_suffix)
		  apply("+", target, fertility_prefix, features);
	      }
	      
	      caches[titer->id()] = features;
	    }
	    
	    features += caches[titer->id()];
	  }
      }


      template <typename Features>
      void apply(const word_type& source, const word_type& target, Features& features)
      {
	feature_builder.clear();
	feature_builder << prefix << ":" << source << "_" << target;
	
	if (forced_feature || feature_builder.exists())
	  features[feature_builder] += 1.0;
	
	if (! normalizers_source.empty()) {
	  const cache_normalize_type::word_set_type& normalized_source = normalize(source,  normalizers_source, cache_source);
	  
	  cache_normalize_type::word_set_type::const_iterator siter_end = normalized_source.end();
	  for (cache_normalize_type::word_set_type::const_iterator siter = normalized_source.begin(); siter != siter_end; ++ siter) {
	    
	    feature_builder.clear();
	    feature_builder << prefix << ":" << *siter << "_" << target;
	    
	    if (forced_feature || feature_builder.exists())
	      features[feature_builder] += 1.0;
	    
	    if (! normalizers_target.empty()) {
	      const cache_normalize_type::word_set_type& normalized_target = normalize(target, normalizers_target, cache_target);
	      
	      cache_normalize_type::word_set_type::const_iterator titer_end = normalized_target.end();
	      for (cache_normalize_type::word_set_type::const_iterator titer = normalized_target.begin(); titer != titer_end; ++ titer) {
		
		feature_builder.clear();
		feature_builder << prefix << ":" << *siter << "_" << *titer;
		
		if (forced_feature || feature_builder.exists())
		  features[feature_builder] += 1.0;
	      }
	    }
	  }
	}
	
	if (! normalizers_target.empty()) {
	  const cache_normalize_type::word_set_type& normalized_target = normalize(target, normalizers_target, cache_target);
	  
	  cache_normalize_type::word_set_type::const_iterator titer_end = normalized_target.end();
	  for (cache_normalize_type::word_set_type::const_iterator titer = normalized_target.begin(); titer != titer_end; ++ titer) {
	    feature_builder.clear();
	    feature_builder << prefix << ":" << source << "_" << *titer;
	    
	    if (forced_feature || feature_builder.exists())
	      features[feature_builder] += 1.0;
	  }
	}
      }
      
      template <typename Features>
      void apply(const char* tag, const word_pair_type& source, const word_type& target, Features& features)
      {
	feature_builder.clear();
	feature_builder << prefix << ":" << tag << "|" << source.first << ":" << source.second << "_" << target;
	
	if (forced_feature || feature_builder.exists())
	  features[feature_builder] += 1.0;
	
	if (! normalizers_source.empty()) {
	  const cache_normalize_type::word_set_type normalized_source_prev = normalize(source.first,  normalizers_source, cache_source);
	  const cache_normalize_type::word_set_type normalized_source_next = normalize(source.second, normalizers_source, cache_source);

	  cache_normalize_type::word_set_type::const_iterator piter_end = normalized_source_prev.end();
	  for (cache_normalize_type::word_set_type::const_iterator piter = normalized_source_prev.begin(); piter != piter_end; ++ piter) {
	    cache_normalize_type::word_set_type::const_iterator niter_end = normalized_source_next.end();
	    for (cache_normalize_type::word_set_type::const_iterator niter = normalized_source_next.begin(); niter != niter_end; ++ niter) {
	      feature_builder.clear();
	      feature_builder << prefix << ":" << tag << "|" << *piter << ":" << *niter << "_" << target;
	      
	      if (forced_feature || feature_builder.exists())
		features[feature_builder] += 1.0;
	      
	      if (! normalizers_target.empty()) {
		const cache_normalize_type::word_set_type& normalized_target = normalize(target, normalizers_target, cache_target);
		
		cache_normalize_type::word_set_type::const_iterator titer_end = normalized_target.end();
		for (cache_normalize_type::word_set_type::const_iterator titer = normalized_target.begin(); titer != titer_end; ++ titer) {
		  feature_builder.clear();
		  feature_builder << prefix << ":" << tag << "|" << *piter << ":" << *niter << "_" << *titer;
		  
		  if (forced_feature || feature_builder.exists())
		    features[feature_builder] += 1.0;
		}
	      }
	    }
	  }
	}
	
	if (! normalizers_target.empty()) {
	  const cache_normalize_type::word_set_type& normalized_target = normalize(target, normalizers_target, cache_target);
	  
	  cache_normalize_type::word_set_type::const_iterator titer_end = normalized_target.end();
	  for (cache_normalize_type::word_set_type::const_iterator titer = normalized_target.begin(); titer != titer_end; ++ titer) {
	    feature_builder.clear();
	    feature_builder << prefix << ":" << tag << "|" << source.first << ":" << source.second << "_" << *titer;
	    
	    if (forced_feature || feature_builder.exists())
	      features[feature_builder] += 1.0;
	  }
	}
      }
      
      template <typename Features>
      void apply(const char* tag, const word_type& target, const int fertility, Features& features)
      {
	const int fertility_power2 = utils::bithack::branch(utils::bithack::is_power2(fertility),
							    fertility,
							    static_cast<int>(utils::bithack::next_largest_power2(fertility)));
	
	feature_builder.clear();
	feature_builder << prefix << ":fertility" << tag << "|" << target << "|" << fertility_power2;
	
	if (forced_feature || feature_builder.exists())
	  features[feature_builder] += 1.0;
      }
      
      void assign(const lattice_type& lattice)
      {
	if (skip_sgml_tag)
	  assign(lattice, skipper_sgml());
	else
	  assign(lattice, skipper_epsilon());
      }
      
      template <typename Skipper>
      void assign(const lattice_type& lattice, Skipper skipper)
      {
	clear();
	
	pos_set_type  positions;
	
	lattice_prev.clear();
	lattice_prev.resize(lattice.size() + 1);
	lattice_prev.front().push_back(vocab_type::BOS);
	
	for (size_t pos = 0; pos != lattice.size(); ++ pos) {
	  positions.clear();
          
	  lattice_type::arc_set_type::const_iterator aiter_end = lattice[pos].end();
	  for (lattice_type::arc_set_type::const_iterator aiter = lattice[pos].begin(); aiter != aiter_end; ++ aiter) {
	    if (! skipper(aiter->label)) {
	      
	      if (pair_mode && exists(aiter->label))
		words.push_back(aiter->label);
	      
	      lattice_prev[pos + aiter->distance].push_back(aiter->label);
	      
	      // we will compute pair of lattice_prev[pos] and words
	      if (prefix_mode || suffix_mode) {
		word_set_type::const_iterator piter_end = lattice_prev[pos].end();
		for (word_set_type::const_iterator piter = lattice_prev[pos].begin(); piter != piter_end; ++ piter) {
		  
		  if (prefix_mode)
		    if (! lexicon_prefix || exists(*lexicon_prefix, *piter, aiter->label))
		      words_prefix.push_back(std::make_pair(*piter, aiter->label));
		  
		  if (suffix_mode)
		    if (*piter != vocab_type::BOS && (! lexicon_suffix || exists(*lexicon_suffix, *piter, aiter->label)))
		      words_suffix.push_back(std::make_pair(*piter, aiter->label));
		}
	      }
	    } else
	      positions.insert(pos + aiter->distance);
	  }
	  
	  // copy lattice_prev[pos] into  positons.
	  pos_set_type::const_iterator piter_end = positions.end();
	  for (pos_set_type::const_iterator piter = positions.begin(); piter != piter_end; ++ piter)
	    lattice_prev[*piter].insert(lattice_prev[*piter].end(), lattice_prev[pos].begin(), lattice_prev[pos].end());
	}
	
	if (suffix_mode) {
	  // we will compute pair of lattice_prev[lattice.size()] and EOS
	  word_set_type::const_iterator piter_end = lattice_prev[lattice.size()].end();
	  for (word_set_type::const_iterator piter = lattice_prev[lattice.size()].begin(); piter != piter_end; ++ piter)
	    if (! lexicon_suffix || exists(*lexicon_suffix, *piter, vocab_type::EOS))
	      words_suffix.push_back(std::make_pair(*piter, vocab_type::EOS));
	}
	
	if (pair_mode) {
	  if (unique_source) {
	    uniques.clear();
	    uniques.insert(words.begin(), words.end());
	    
	    words.clear();
	    words.insert(words.end(), uniques.begin(), uniques.end());
	  }
	  
	  std::sort(words.begin(), words.end());
	} else
	  words.clear();
	
	if (prefix_mode) {
	  if (unique_source) {
	    uniques_prefix.clear();
	    uniques_prefix.insert(words_prefix.begin(), words_prefix.end());
	    
	    words_prefix.clear();
	    words_prefix.insert(words_prefix.end(), uniques_prefix.begin(), uniques_prefix.end());
	  }
	  
	  std::sort(words_prefix.begin(), words_prefix.end());
	} else
	  words_prefix.clear();
	
	if (suffix_mode) {
	  if (unique_source) {
	    uniques_suffix.clear();
	    uniques_suffix.insert(words_suffix.begin(), words_suffix.end());
	    
	    words_suffix.clear();
	    words_suffix.insert(words_suffix.end(), uniques_suffix.begin(), uniques_suffix.end());
	  }
	  
	  std::sort(words_suffix.begin(), words_suffix.end());
	} else
	  words_suffix.clear();
      }
      
      const cache_normalize_type::word_set_type& normalize(const word_type& word,
							   normalizer_set_type& normalizers,
							   cache_normalize_set_type& caches)
      {
	cache_normalize_type& cache = caches[hasher_type::operator()(word.id()) & (caches.size() - 1)];
	if (cache.word != word.id()) {
	  cache.word = word.id();
	  cache.normalized.clear();
	  
	  for (size_t i = 0; i != normalizers.size(); ++ i) {
	    const word_type normalized = normalizers[i](word);
	    if (word != normalized)
	      cache.normalized.push_back(normalized);
	  }
	}
	
	return cache.normalized;
      }
      
      void clear()
      {
	words.clear();
	caches.clear();

	words_prefix.clear();
	words_suffix.clear();
      }
      
      void clear_cache()
      {
	cache_source.clear();
	cache_target.clear();
      }

      bool exists(const word_type& source, const word_type& target)
      {
	return (! lexicon) || (lexicon->exists(&source, (&source) + 1, target));
      }
      
      bool exists(const word_type& source)
      {
	return (! lexicon) || (lexicon->exists(&source, (&source) + 1));
      }
      
      bool exists(const lexicon_type& lexicon, const word_pair_type& prev, const word_type& next)
      {
        word_type codes[2];
        codes[0] = prev.first;
        codes[1] = prev.second;
        return lexicon.exists(codes, codes + 2, next);
      }
      
      bool exists(const lexicon_type& lexicon, const word_type& prev, const word_type& next)
      {
        word_type codes[2];
        codes[0] = prev;
        codes[1] = next;
        return lexicon.exists(codes, codes + 2);
      }

      lexicon_type* lexicon;
      lexicon_type* lexicon_prefix;
      lexicon_type* lexicon_suffix;
      
      normalizer_set_type normalizers_source;
      normalizer_set_type normalizers_target;
      
      cache_normalize_set_type cache_source;
      cache_normalize_set_type cache_target;

      cache_set_type   caches;
      
      word_unique_type uniques;
      word_set_type    words;
      
      word_pair_set_type words_prefix;
      word_pair_set_type words_suffix;
      
      word_pair_unique_type  uniques_prefix;
      word_pair_unique_type  uniques_suffix;
      word_map_type          lattice_prev;
      
      feature_builder_type feature_builder;
      
      bool skip_sgml_tag;
      bool unique_source;
      
      std::string prefix;
      bool forced_feature;

      bool pair_mode;
      bool prefix_mode;
      bool suffix_mode;
      bool fertility_mode;
    };
    
    SparseLexicon::SparseLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "sparse-lexicon")
	throw std::runtime_error("this is not sparse lexicon feature: " + parameter);

      impl_type::normalizer_set_type normalizers_source;
      impl_type::normalizer_set_type normalizers_target;
      
      bool skip_sgml_tag = false;
      bool unique_source = false;
      
      std::string name;

      std::string lexicon;
      std::string lexicon_prefix;
      std::string lexicon_suffix;
      
      bool pair_mode = false;
      bool prefix_mode = false;
      bool suffix_mode = false;
      bool fertility_mode = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "cluster-source") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers_source.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "cluster-target") {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers_target.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (utils::ipiece(piter->first) == "stemmer-source")
	  normalizers_source.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "stemmer-target")
	  normalizers_target.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "unique-source")
	  unique_source = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "lexicon")
          lexicon = piter->second;
	else if (utils::ipiece(piter->first) == "lexicon-prefix")
          lexicon_prefix = piter->second;
	else if (utils::ipiece(piter->first) == "lexicon-suffix")
          lexicon_suffix = piter->second;
	else if (utils::ipiece(piter->first) == "pair")
	  pair_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "prefix")
	  prefix_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "suffix")
	  suffix_mode = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "fertility")
	  fertility_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for sparse lexicon: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! lexicon.empty())
	pair_mode = true;
      
      if (! lexicon_prefix.empty())
	prefix_mode = true;
      
      if (! lexicon_suffix.empty())
	suffix_mode = true;
      
      if (int(pair_mode) + prefix_mode + suffix_mode == 0)
	pair_mode = true;
      
      std::unique_ptr<impl_type> lexicon_impl(new impl_type());

      lexicon_impl->normalizers_source.swap(normalizers_source);
      lexicon_impl->normalizers_target.swap(normalizers_target);
      
      lexicon_impl->skip_sgml_tag = skip_sgml_tag;
      lexicon_impl->unique_source = unique_source;
      lexicon_impl->prefix = (name.empty() ? std::string("sparse-lexicon") : name);

      if (! lexicon.empty())
	lexicon_impl->lexicon = &cicada::Lexicon::create(lexicon);
      if (! lexicon_prefix.empty())
	lexicon_impl->lexicon_prefix = &cicada::Lexicon::create(lexicon_prefix);
      if (! lexicon_suffix.empty())
	lexicon_impl->lexicon_suffix = &cicada::Lexicon::create(lexicon_suffix);

      lexicon_impl->pair_mode   = pair_mode;
      lexicon_impl->prefix_mode = prefix_mode;
      lexicon_impl->suffix_mode = suffix_mode;
      lexicon_impl->fertility_mode = fertility_mode;
      
      base_type::__state_size = 0;
      base_type::__feature_name = (name.empty() ? std::string("sparse-lexicon") : name);
      base_type::__sparse_feature = true;
      
      pimpl = lexicon_impl.release();
    }
    
    SparseLexicon::~SparseLexicon() { std::unique_ptr<impl_type> tmp(pimpl); }
    
    SparseLexicon::SparseLexicon(const SparseLexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {
      pimpl->clear_cache();

      pimpl->lexicon = 0;
      pimpl->lexicon_prefix = 0;
      pimpl->lexicon_suffix = 0;
      if (x.pimpl->lexicon)
        pimpl->lexicon = &cicada::Lexicon::create(x.pimpl->lexicon->path().string());
      if (x.pimpl->lexicon_prefix)
        pimpl->lexicon_prefix = &cicada::Lexicon::create(x.pimpl->lexicon_prefix->path().string());
      if (x.pimpl->lexicon_suffix)
        pimpl->lexicon_suffix = &cicada::Lexicon::create(x.pimpl->lexicon_suffix->path().string());
    }

    SparseLexicon& SparseLexicon::operator=(const SparseLexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      pimpl->clear_cache();
      
      pimpl->lexicon = 0;
      pimpl->lexicon_prefix = 0;
      pimpl->lexicon_suffix = 0;
      if (x.pimpl->lexicon)
        pimpl->lexicon = &cicada::Lexicon::create(x.pimpl->lexicon->path().string());
      if (x.pimpl->lexicon_prefix)
        pimpl->lexicon_prefix = &cicada::Lexicon::create(x.pimpl->lexicon_prefix->path().string());
      if (x.pimpl->lexicon_suffix)
        pimpl->lexicon_suffix = &cicada::Lexicon::create(x.pimpl->lexicon_suffix->path().string());
      
      return *this;
    }
    
    void SparseLexicon::apply(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      const bool final) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      feature_set_type feats;
      
      pimpl->lexicon_score(edge, feats);

      features.update(feats, static_cast<const std::string&>(base_type::feature_name()));
    }

    void SparseLexicon::apply_coarse(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void SparseLexicon::apply_predict(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      const bool final) const
    {
      apply(state, states, edge, features, final);
    }

    
    void SparseLexicon::apply_scan(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   const int dot,
				   feature_set_type& features,
				   const bool final) const
    {}
    void SparseLexicon::apply_complete(state_ptr_type& state,
				       const state_ptr_set_type& states,
				       const edge_type& edge,
				       feature_set_type& features,
				       const bool final) const
    {}
    
    
    void SparseLexicon::assign(const size_type& id,
			       const hypergraph_type& hypergraph,
			       const lattice_type& lattice,
			       const span_set_type& spans,
			       const sentence_set_type& targets,
			       const ngram_count_set_type& ngram_counts)
    {
      //
      // how do we assign lexion feature from hypergraph...???
      // we assume that the lattice is always filled with the source-word...!
      //
      
      pimpl->clear();
      
      if (! lattice.empty())
	pimpl->assign(lattice);
    }
  };
};
