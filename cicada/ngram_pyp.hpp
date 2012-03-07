//
//  Copyright(C) 2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__NGRAM_PYP__HPP__
#define __CICADA__NGRAM_PYP__HPP__ 1

// PYPLM!

#include <stdint.h>

#include <stdexcept>
#include <vector>

#include <boost/filesystem.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>

#include <utils/packed_vector.hpp>
#include <utils/succinct_vector.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/array_power2.hpp>
#include <utils/spinlock.hpp>
#include <utils/bithack.hpp>

namespace cicada
{
  class NGramPYP : public utils::hashmurmur<size_t>
  {
  public:
    typedef Symbol                  word_type;
    typedef Vocab                   vocab_type;

    typedef size_t                  size_type;
    typedef ptrdiff_t               difference_type;
    typedef word_type::id_type      id_type;
    typedef uint64_t                count_type;
    
    typedef boost::filesystem::path   path_type;
    typedef utils::hashmurmur<size_t> hasher_type;
    
    typedef utils::packed_vector_mapped<id_type, std::allocator<id_type> >       id_set_type;
    typedef utils::packed_vector_mapped<count_type, std::allocator<count_type> > count_set_type;
    typedef utils::succinct_vector_mapped<std::allocator<int32_t> >              position_set_type;
    typedef std::vector<size_type, std::allocator<size_type> >                   off_set_type;
    
    typedef utils::spinlock spinlock_type;
    typedef spinlock_type::scoped_lock     lock_type;
    typedef spinlock_type::scoped_try_lock trylock_type;
    
    
  private:
    id_set_type       index_;
    count_set_type    count_;
    position_set_type position_;
    
    vocab_type        vocab_;
    
    int               order_;
    path_type         path_;
  };
  
};

#endif
