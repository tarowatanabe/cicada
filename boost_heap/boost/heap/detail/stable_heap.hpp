// boost heap: helper classes for stable priority queues
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_DETAIL_STABLE_HEAP_HPP
#define BOOST_HEAP_DETAIL_STABLE_HEAP_HPP

#include <cassert>
#include <limits>
#include <stdexcept>
#include <utility>

#include <boost/throw_exception.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include <boost/parameter.hpp>

namespace boost {
namespace heap {
namespace detail {


template<bool ConstantSize, class SizeType>
struct size_holder
{
   static const bool constant_time_size = ConstantSize;
   typedef SizeType  size_type;

    size_holder(void):
        size_(0)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    size_holder(size_holder && rhs):
        size_(rhs.size_)
    {
        rhs.size_ = 0;
    }

    size_holder & operator=(size_holder && rhs)
    {
        size_ = rhs.size_;
        rhs.size_ = 0;
    }
#endif

    SizeType get_size() const
    {  return size_;  }

    void set_size(SizeType size)
    {  size_ = size; }

    void decrement()
    {  --size_; }

    void increment()
    {  ++size_; }

    void add(SizeType value)
    {  size_ += value; }

    void sub(SizeType value)
    {  size_ -= value; }

    void swap(size_holder & rhs)
    {  std::swap(size_, rhs.size_); }

    SizeType size_;
};

template<class SizeType>
struct size_holder<false, SizeType>
{
    static const bool constant_time_size = false;
    typedef SizeType  size_type;

    size_holder(void)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    size_holder(size_holder && rhs)
    {}

    size_holder & operator=(size_holder && rhs)
    {}
#endif

    size_type get_size() const
    {  return 0;  }

    void set_size(size_type)
    {}

    void decrement()
    {}

    void increment()
    {}

    void add(SizeType value)
    {}

    void sub(SizeType value)
    {}

    void swap(size_holder & rhs)
    {}
};


template <typename T,
          typename Cmp,
          bool constant_time_size,
          bool stable = false
         >
struct heap_base:
    Cmp,
    size_holder<constant_time_size, size_t>
{
    typedef T value_type;
    typedef T internal_type;
    typedef size_holder<constant_time_size, size_t> size_holder_type;
    typedef Cmp compare_type;

    heap_base (Cmp const & cmp = Cmp()):
        Cmp(cmp)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    heap_base(heap_base && rhs):
        Cmp(std::move(static_cast<Cmp&>(rhs))),
        size_holder_type(std::move(static_cast<size_holder_type&>(rhs)))
    {}

    heap_base & operator=(heap_base && rhs)
    {
        Cmp::operator=(std::move(static_cast<Cmp&>(rhs)));
        size_holder_type::operator=(std::move(static_cast<size_holder_type&>(rhs)));
    }
#endif

    bool operator()(internal_type const & lhs, internal_type const & rhs) const
    {
        return Cmp::operator()(lhs, rhs);
    }

    internal_type make_node(T const & val)
    {
        return val;
    }

    static T & get_value(internal_type & val)
    {
        return val;
    }

    static T const & get_value(internal_type const & val)
    {
        return val;
    }

    Cmp const & get_cmp(void)
    {
        return *this;
    }

    void swap(heap_base & rhs)
    {
        std::swap(static_cast<Cmp&>(*this), static_cast<Cmp&>(rhs));
        size_holder<constant_time_size, size_t>::swap(rhs);
    }
};

template <typename T,
          typename Cmp,
          bool constant_time_size>
struct heap_base<T, Cmp, constant_time_size, true>:
    Cmp,
    size_holder<constant_time_size, size_t>
{
    typedef T value_type;
    typedef std::pair<T, size_t> internal_type;
    typedef size_holder<constant_time_size, size_t> size_holder_type;
    typedef Cmp compare_type;

    heap_base (Cmp const & cmp = Cmp()):
        Cmp(cmp), counter_(0)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    heap_base(heap_base && rhs):
        Cmp(std::move(static_cast<Cmp&>(rhs))),
        size_holder_type(std::move(static_cast<size_holder_type&>(rhs))), counter_(rhs.counter_)
    {
        rhs.counter_ = 0;
    }

    heap_base & operator=(heap_base && rhs)
    {
        Cmp::operator=(std::move(static_cast<Cmp&>(rhs)));
        size_holder_type::operator=(std::move(static_cast<size_holder_type&>(rhs)));

        counter_ = rhs.counter_;
        rhs.counter_ = 0;
    }
#endif

    bool operator()(internal_type const & lhs, internal_type const & rhs) const
    {
        if (!Cmp::operator()(lhs.first, rhs.first) && !Cmp::operator()(rhs.first, lhs.first))
            /* lhs and rhs are equal */
            return lhs.second > rhs.second;

        return Cmp::operator()(lhs.first, rhs.first);
    }

    bool operator()(T const & lhs, T const & rhs) const
    {
        return Cmp::operator()(lhs, rhs);
    }

    internal_type make_node(T const & val)
    {
        std::size_t count = ++counter_;
        if (counter_ == std::numeric_limits<std::size_t>::max())
            BOOST_THROW_EXCEPTION(std::runtime_error("boost::heap counter overflow"));
        return std::make_pair(val, count);
    }

    static T & get_value(internal_type & val)
    {
        return val.first;
    }

    static T const & get_value(internal_type const & val)
    {
        return val.first;
    }

    Cmp const & get_cmp(void)
    {
        return *this;
    }

    void swap(heap_base & rhs)
    {
        std::swap(static_cast<Cmp&>(*this), static_cast<Cmp&>(rhs));
        std::swap(counter_, rhs.counter_);
        size_holder<constant_time_size, size_t>::swap(rhs);
    }

private:
    std::size_t counter_;
};

template <typename node_pointer,
          typename extractor,
          typename reference
         >
struct node_handle
{
    explicit node_handle(node_pointer n = 0):
        node_(n)
    {}

    reference operator*() const
    {
        return extractor::get_value(node_->value);
    }

    node_pointer node_;
};

template <typename value_type,
          typename internal_type,
          typename extractor
         >
struct value_extractor
{
    value_type const & operator()(internal_type const & data) const
    {
        return extractor::get_value(data);
    }
};

template <typename T,
          typename ContainerIterator,
          typename Extractor>
class stable_heap_iterator:
    public boost::iterator_adaptor<stable_heap_iterator<T, ContainerIterator, Extractor>,
                                   ContainerIterator,
                                   T const,
                                   boost::random_access_traversal_tag>
{
    typedef boost::iterator_adaptor<stable_heap_iterator,
                                    ContainerIterator,
                                    T const,
                                    boost::random_access_traversal_tag> super_t;

public:
    stable_heap_iterator(void):
        super_t(0)
    {}

    explicit stable_heap_iterator(ContainerIterator const & it):
        super_t(it)
    {}

private:
    friend class boost::iterator_core_access;

    T const & dereference() const
    {
        return Extractor::get_value(*super_t::base());
    }
};

template <typename T, typename Parspec, bool constant_time_size>
struct make_heap_base
{
    typedef typename parameter::binding<Parspec, tag::compare, std::less<T> >::type compare_argument;
    typedef typename parameter::binding<Parspec, tag::allocator, std::allocator<T> >::type allocator_argument;

    static const bool stable = parameter::binding<Parspec, tag::stable, boost::mpl::false_>::type::value;

    typedef heap_base<T, compare_argument, constant_time_size, stable> type;
};

} /* namespace detail */
} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_DETAIL_STABLE_HEAP_HPP */
