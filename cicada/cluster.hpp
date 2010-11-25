// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__CLUSTER__HPP__
#define __CICADA__CLUSTER__HPP__ 1

#include <string>
#include <algorithm>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/packed_vector.hpp>

#include <boost/filesystem.hpp>

namespace cicada
{
  class Cluster
  {
  public:
    typedef Symbol    symbol_type;
    typedef Vocab     vocab_type;

    typedef symbol_type          word_type;
    typedef symbol_type::id_type id_type;
    
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    typedef uint64_t  hash_value_type;
    
    typedef boost::filesystem::path path_type;
    
  private:
    typedef utils::packed_vector_mapped<id_type, std::allocator<id_type> > cluster_set_type;
    
  public:
    Cluster() {}
    Cluster(const path_type& path) { open(path); }

    bool empty() const { return clusters.empty(); }
    
    void clear()
    {
      vocab.clear();
      clusters.clear();
      file.clear();
    }
    void close() { clear(); }
    
    void open(const path_type& path);
    void write(const path_type& path) const;
    const path_type& path() const { return file; }
    
    symbol_type operator()(const symbol_type& word) const { return operator[](word); }
    symbol_type operator[](const symbol_type& word) const
    {
      // empty word, non-terminals are not clusters...
      if (word == vocab_type::EMPTY || word.is_non_terminal() || empty())
	return word;
      
      const size_type word_size = word.size();
      
      // SGML-like symbols are not clusters...
      if (word_size >= 3 && word[0] == '<' && word[word_size - 1] == '>')
	return word;

      const id_type id = vocab[word];
      
      if (id < clusters.size()){
	const id_type cluster = clusters[id];
	return (cluster == 0 ? vocab_type::UNK : vocab[cluster - 1]);
      } else
	return vocab_type::UNK;
    }

  public:
    static Cluster& create(const path_type& path);
    
  private:
    vocab_type       vocab;
    cluster_set_type clusters;
    path_type        file;
  };
};

#endif
