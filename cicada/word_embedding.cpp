// -*- mode: c++ -*-
//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#define BOOST_DISABLE_ASSERTS
#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>

#include <boost/fusion/tuple.hpp>

#include "word_embedding.hpp"

#include "utils/map_file.hpp"
#include "utils/unordered_map.hpp"
#include "utils/indexed_set.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/repository.hpp"
#include "utils/spinlock.hpp"

namespace cicada
{
  template <typename Value>
  inline
  Value repository_value(const utils::repository& rep, const std::string& key)
  {
    utils::repository::const_iterator iter = rep.find(key);
    if (iter == rep.end())
      throw std::runtime_error("no " + key + "?");
    return utils::lexical_cast<Value>(iter->second);
  }

  
  class WordEmbeddingMapped : public WordEmbedding
  {
  private:
    typedef utils::map_file<parameter_type, std::allocator<parameter_type> > mapped_type;
    
  public:
    WordEmbeddingMapped() : rows_(0) {}
    WordEmbeddingMapped(const path_type& path)
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::read);
      
      rows_ = repository_value<size_type>(rep, "embedding");
      
      mapped_.open(rep.path("input.bin"));
      vocab_.open(rep.path("vocab"));
    }
    
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      
      rep["embedding"] = utils::lexical_cast<std::string>(rows_);
      
      mapped_.write(rep.path("input.bin"));
      vocab_.write(rep.path("vocab"));
    }

    matrix_type operator()(const word_type& word) const
    {
      return matrix_type(const_cast<parameter_type*>(mapped_.begin() + rows_ * vocab_[word]), rows_, 1);
    }

    size_type dimension() const { return rows_; }
    
  private:
    vocab_type  vocab_;
    mapped_type mapped_;
    size_type   rows_;
  };
  
  class WordEmbeddingMemory : public WordEmbedding
  {
  private:
    typedef utils::indexed_set<word_type, boost::hash<word_type>, std::equal_to<word_type>,
			       std::allocator<word_type> > word_set_type;
  public:
    WordEmbeddingMemory() : matrix_(), words_(), oov_(0) {}
    WordEmbeddingMemory(const path_type& path)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      typedef std::vector<parameter_type, std::allocator<parameter_type> > parameter_set_type;
typedef boost::fusion::tuple<std::string, parameter_set_type > embedding_parsed_type;
      typedef boost::spirit::istream_iterator iterator_type;
      
      qi::rule<iterator_type, std::string(), standard::blank_type>           word;
      qi::rule<iterator_type, embedding_parsed_type(), standard::blank_type> parser; 
      
      word   %= qi::lexeme[+(standard::char_ - standard::space)];
      parser %= word >> *qi::double_ >> (qi::eol | qi::eoi);
      
      utils::compress_istream is(path, 1024 * 1024);
      is.unsetf(std::ios::skipws);
      
      iterator_type iter(is);
      iterator_type iter_end;
      
      embedding_parsed_type parsed;

      size_type rows = size_type(-1);
      
      oov_ = word_type::id_type(-1);

      while (iter != iter_end) {
	boost::fusion::get<0>(parsed).clear();
	boost::fusion::get<1>(parsed).clear();
	
	if (! boost::spirit::qi::phrase_parse(iter, iter_end, parser, standard::blank, parsed))
	  if (iter != iter_end)
	    throw std::runtime_error("embedding parsing failed");
	
	if (rows == size_type(-1))
	  rows = boost::fusion::get<1>(parsed).size();
	
	if (boost::fusion::get<1>(parsed).size() != rows)
	  throw std::runtime_error("invalid embedding size");
	
	const word_type word = boost::fusion::get<0>(parsed);
	
	word_set_type::const_iterator witer = words_.insert(word).first;
	const word_type::id_type id = witer - words_.begin();

	if (word == vocab_type::UNK)
	  oov_ = id;
	
	if (id >= matrix_.cols())
	  matrix_.conservativeResize(rows, id + 1);
	
	matrix_.col(id) = matrix_type(&(*boost::fusion::get<1>(parsed).begin()), rows, 1);
      }
      
      // estimate OOV by performing averaging...
      if (oov_ == word_type::id_type(-1)) {
	word_set_type::const_iterator witer = words_.insert(vocab_type::UNK).first;
	const word_type::id_type oov_ = witer - words_.begin();
	
	if (oov_ >= matrix_.cols())
	  matrix_.conservativeResize(rows, oov_ + 1);
	
	matrix_.col(oov_).setZero();
	
	const word_type::id_type vocabulary_size = matrix_.cols();
	const double factor = 1.0 / (vocabulary_size - 1);
	for (word_type::id_type id = 0; id != vocabulary_size; ++ id)
	  if (id != oov_)
	    matrix_.col(oov_) += matrix_.col(id) * factor;
      }
    }
    
    void write(const path_type& path) const
    {
      typedef utils::repository repository_type;
      
      repository_type rep(path, repository_type::write);
      
      rep["embedding"] = utils::lexical_cast<std::string>(matrix_.rows());
      
      utils::compress_ostream os(rep.path("input.bin"), 1024 * 1024);
      
      os.write((char*) matrix_.data(), sizeof(parameter_type) * matrix_.rows() * matrix_.cols());
      
      vocab_type vocab;
      
      vocab.open(rep.path("vocab"), words_.size());
      
      word_set_type::const_iterator witer_end = words_.end();
      for (word_set_type::const_iterator witer = words_.begin(); witer != witer_end; ++ witer)
	vocab.insert(*witer);
    }
    
    matrix_type operator()(const word_type& word) const
    {
      word_set_type::const_iterator witer = words_.find(word);
      
      word_type::id_type id = oov_;
      if (witer != words_.end())
	id = witer - words_.begin();
      
      const size_type rows = matrix_.rows();
      
      return matrix_type(const_cast<parameter_type*>(matrix_.data() + rows * id), rows, 1);
    }

    size_type dimension() const { return matrix_.rows(); }
    
  private:
    tensor_type        matrix_;
    word_set_type      words_;
    word_type::id_type oov_;
  };


  const WordEmbedding& WordEmbedding::create(const path_type& path)
  {
    typedef boost::shared_ptr<WordEmbedding> shared_type;
    typedef std::string string_type;
    typedef utils::unordered_map<string_type, shared_type,
				 boost::hash<string_type>, std::equal_to<string_type>,
				 std::allocator<std::pair<const string_type, shared_type> > >::type mapped_type;
    
    typedef utils::spinlock         mutex_type;
    typedef mutex_type::scoped_lock lock_type;
    
    static mutex_type  mutex;
    static mapped_type mapped;
    
    lock_type lock(mutex);
    
    const string_type rep = path.string();
    
    mapped_type::iterator miter = mapped.find(rep);
    if (miter == mapped.end())
      miter = mapped.insert(std::make_pair(rep,
					   boost::filesystem::is_directory(path)
					   ? shared_type(new WordEmbeddingMapped(path))
					   : shared_type(new WordEmbeddingMemory(path)))).first;

    return *(miter->second);
  }
};
