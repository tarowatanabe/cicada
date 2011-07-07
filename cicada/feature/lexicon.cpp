//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "lexicon.hpp"

namespace cicada
{
  namespace feature
  {

    class LexiconImpl
    {
    public:
      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      
    };
    
    Lexicon::Lexicon(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "bleu")
	throw std::runtime_error("this is not Lexicon feature: " + parameter);

      int order = 4;
      bool exact = false;
      bool skip_sgml_tag = false;

      const cicada::Tokenizer* tokenizer = 0;
      
      std::string name;
      path_type   refset_file;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "order")
	  order = utils::lexical_cast<int>(piter->second);
	else if (utils::ipiece(piter->first) == "exact")
	  exact = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	  skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "tokenizer")
	  tokenizer = &cicada::Tokenizer::create(piter->second);
	else if (utils::ipiece(piter->first) == "name")
	  name = piter->second;
	else if (utils::ipiece(piter->first) == "refset")
	  refset_file = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (! refset_file.empty() && ! boost::filesystem::exists(refset_file))
	throw std::runtime_error("no refset file?: " + refset_file.string());
      
      std::auto_ptr<impl_type> bleu_impl(new impl_type(order, exact, skip_sgml_tag, tokenizer));
      
      // two-side context + length + counts-id 
      base_type::__state_size = sizeof(symbol_type) * order * 2 + sizeof(int) + sizeof(impl_type::id_type);
      base_type::__feature_name = (name.empty() ? std::string("bleu") : name);
      
      pimpl = bleu_impl.release();
      
      pimpl->refset.clear();
      
      if (! refset_file.empty()) {
	typedef boost::tokenizer<utils::space_separator, utils::piece::const_iterator, utils::piece> tokenizer_type;
	
	utils::compress_istream is(refset_file, 1024 * 1024);
	std::string line;
	
	while (std::getline(is, line)) {
	  utils::piece line_piece(line);
	  tokenizer_type tokenizer(line_piece);
	  
	  tokenizer_type::iterator iter = tokenizer.begin();
	  if (iter == tokenizer.end()) continue;
	  
	  const int id = utils::lexical_cast<int>(*iter);
	  ++ iter;
	  
	  if (iter == tokenizer.end()) continue;
	  if (*iter != "|||") continue;
	  ++ iter;
	  
	  if (id >= static_cast<int>(pimpl->refset.size()))
	    pimpl->refset.resize(id + 1);
	  
	  pimpl->refset[id].push_back(sentence_type(iter, tokenizer.end()));
	}
      }
    }
    
    Lexicon::~Lexicon() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    Lexicon::Lexicon(const Lexicon& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}

    Lexicon& Lexicon::operator=(const Lexicon& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Lexicon::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
    {
      
    }

    void Lexicon::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {
      
    }
    void Lexicon::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    void Lexicon::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Lexicon::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    
    void Lexicon::assign(const size_type& id,
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
      
    }
    
  };
};
