// -*- mode: c++ -*-
//
//  Copyright(C) 2014 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__WORD_EMBEDDING__HPP__
#define __CICADA__WORD_EMBEDDING__HPP__ 1

#include <stdexcept>

#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <Eigen/Core>


namespace cicada
{
  class WordEmbedding
  {
  public:
    typedef Symbol                  word_type;
    typedef Vocab                   vocab_type;
    
    typedef size_t             size_type;
    typedef ptrdiff_t          difference_type;
    
    typedef float              parameter_type;

    typedef boost::filesystem::path path_type;

    typedef Eigen::Matrix<parameter_type, Eigen::Dynamic, Eigen::Dynamic>    tensor_type;
    typedef Eigen::Map<tensor_type>                                          matrix_type;
    
  public:
    virtual
    ~WordEmbedding() {}
    
  public:
    virtual
    void write(const path_type& path) const = 0;
    
    virtual
    matrix_type operator()(const word_type& word) const = 0;
    
    virtual
    size_type dimension() const = 0;

  public:
    static
    const WordEmbedding& create(const path_type& path);
  };
  
};

#endif
