
#include "global_lexicon.hpp"

namespace cicada
{

  namespace feature
  {
  
    class GlobalLexiconImpl
    {
    
    };
  
  
    GlobalLexicon::GlobalLexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "ngram-tree")
	throw std::runtime_error("is this really ngram tree feature function? " + parameter);

      bool source = false;
      bool target = false;
      
      int stemmer_prefix_size = 0;
      int stemmer_suffix_size = 0;
      bool stemmer_digits = false;
      
      boost::filesystem::path cluster_path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	} else if (strcasecmp(piter->first.c_str(), "cluster") == 0)
	  cluster_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "prefix") == 0)
	  stemmer_prefix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "suffix") == 0)
	  stemmer_suffix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "digits") == 0)
	  stemmer_digits = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for ngram-tree: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");

      if (stemmer_prefix_size < 0)
	throw std::runtime_error("negative prefix size?");
      if (stemmer_suffix_size < 0)
	throw std::runtime_error("negative suffix size?");
      
      std::auto_ptr<impl_type> ngram_tree_impl(source
					       ? dynamic_cast<impl_type*>(new __GlobalLexiconImpl<__ngram_tree_extract_source>())
					       : dynamic_cast<impl_type*>(new __GlobalLexiconImpl<__ngram_tree_extract_target>()));

      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.file_string());
	
	ngram_tree_impl->cluster = &cicada::Cluster::create(cluster_path);
      }
      
      if (stemmer_prefix_size > 0)
	ngram_tree_impl->stemmer_prefix = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_size));
      
      if (stemmer_suffix_size > 0)
	ngram_tree_impl->stemmer_suffix = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_size));

      if (stemmer_digits)
	ngram_tree_impl->stemmer_digits = &cicada::Stemmer::create("digits");

      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(impl_type::id_type) * 2;
      base_type::__feature_name = std::string("ngram-tree-") + (source ? "source" : "target");
      base_type::__sparse_feature = true;
      
      pimpl = ngram_tree_impl.release();
    }
    
    GlobalLexicon::~GlobalLexicon() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    GlobalLexicon::GlobalLexicon(const GlobalLexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(0)
    {
      typedef __GlobalLexiconImpl<__ngram_tree_extract_source> ngram_tree_source_type;
      typedef __GlobalLexiconImpl<__ngram_tree_extract_target> ngram_tree_target_type;
      
      if (dynamic_cast<const ngram_tree_source_type*>(x.pimpl))
	pimpl = new ngram_tree_source_type();
      else
	pimpl = new ngram_tree_target_type();
    }
    
    GlobalLexicon& GlobalLexicon::operator=(const GlobalLexicon& x)
    {
      typedef __GlobalLexiconImpl<__ngram_tree_extract_source> ngram_tree_source_type;
      typedef __GlobalLexiconImpl<__ngram_tree_extract_target> ngram_tree_target_type;
      
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      std::auto_ptr<impl_type> tmp(pimpl);
      
      if (dynamic_cast<const ngram_tree_source_type*>(x.pimpl))
	pimpl = new ngram_tree_source_type();
      else
	pimpl = new ngram_tree_target_type();
      
      return *this;
    }
    
    void GlobalLexicon::apply(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {
      
      
    }

    void GlobalLexicon::apply_coarse(state_ptr_type& state,
				     const state_ptr_set_type& states,
				     const edge_type& edge,
				     feature_set_type& features,
				     feature_set_type& estimates,
				     const bool final) const
    {
      
    }

    void GlobalLexicon::initialize()
    {
      pimpl->clear();
    }
  };
  
};
