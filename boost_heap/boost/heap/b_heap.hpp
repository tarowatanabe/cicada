// boost heap: d-ary heap as containter adaptor
//
// Copyright (C) 2010 Tim Blechmann, adapted from Poul-Henning Kamp
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_B_HEAP_HPP
#define BOOST_HEAP_B_HEAP_HPP

#include <algorithm>
#include <cmath>
#include <vector>

#include <boost/assert.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/integer/static_log2.hpp>

#include "detail/heap_comparison.hpp"
#include "detail/stable_heap.hpp"
#include "detail/mutable_heap.hpp"

#include "d_ary_heap.hpp"

namespace boost
{
namespace heap
{
namespace detail
{

template <unsigned int n>
struct is_power_of_two
{
    static const bool value = (n%2==0) && (is_power_of_two<(n>>1)>::value);
};

template <>
struct is_power_of_two<2>
{
    static const bool value = true;
};

typedef parameter::parameters<optional<tag::allocator>,
                              optional<tag::compare>,
                              optional<tag::stable>,
                              optional<tag::objects_per_page>
                             > b_heap_signature;


/* base class for b heap */
template <typename T,
          class BoundArgs,
          class IndexUpdater>
class b_heap:
    private make_heap_base<T, BoundArgs, false>::type
{
    typedef make_heap_base<T, BoundArgs, false> heap_base_maker;

    typedef typename heap_base_maker::type super_t;
    typedef typename super_t::internal_type internal_type;

    typedef std::vector<internal_type, typename heap_base_maker::allocator_argument> container_type;
    typedef typename container_type::const_iterator container_iterator;

    typedef typename IndexUpdater::template rebind<internal_type>::other index_updater;

    container_type q_;

public:
    typedef T value_type;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename heap_base_maker::compare_argument compare_type;
    typedef typename heap_base_maker::allocator_argument allocator_type;

    /// \copydoc boost::heap::priority_queue::iterator
    typedef detail::stable_heap_iterator<T, container_iterator, super_t> iterator;

    typedef iterator const_iterator;

    typedef typename allocator_type::reference reference;
    typedef typename allocator_type::const_reference const_reference;
    typedef typename allocator_type::pointer pointer;
    typedef typename allocator_type::const_pointer const_pointer;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit b_heap(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {
        q_.push_back(internal_type()); //padding first element
    }

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    b_heap(b_heap const & rhs):
        super_t(rhs), q_(rhs.q_)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    b_heap(b_heap && rhs):
        super_t(std::move(rhs)), q_(std::move(rhs.q_))
    {}

    b_heap & operator=(b_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        q_ = std::move(rhs.q_);
        return *this;
    }
#endif

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const &)
    b_heap & operator=(b_heap const & rhs)
    {
        static_cast<super_t&>(*this) = static_cast<super_t const &>(rhs);
        q_ = rhs.q_;
        return *this;
    }

    /// \copydoc boost::heap::priority_queue::empty
    bool empty(void) const
    {
        return q_.size() == 1;
    }

    /// \copydoc boost::heap::priority_queue::size
    size_type size(void) const
    {
        return q_.size() - 1;
    }

    /// \copydoc boost::heap::priority_queue::max_size
    size_type max_size(void) const
    {
        return q_.max_size();
    }

    /// \copydoc boost::heap::priority_queue::clear
    void clear(void)
    {
        q_.clear();
        q_.push_back(internal_type());
    }

    /// \copydoc boost::heap::priority_queue::get_allocator
    allocator_type get_allocator(void) const
    {
        return q_.get_allocator();
    }

    /// \copydoc boost::heap::priority_queue::top
    value_type const & top(void) const
    {
        BOOST_ASSERT(!empty());
        return super_t::get_value(q_[1]);
    }

    /// \copydoc boost::heap::priority_queue::push
    void push(value_type const & v)
    {
        q_.push_back(super_t::make_node(v));
        reset_index(q_.size() - 1, q_.size() - 1);
        siftup(q_.size() - 1);
    }

    /// \copydoc boost::heap::priority_queue::pop
    void pop(void)
    {
        BOOST_ASSERT(!empty());
        std::swap(q_[1], q_.back());
        q_.pop_back();

        if (empty())
            return;

        reset_index(1, 1);
        siftdown(1);
    }

    /// \copydoc boost::heap::priority_queue::swap
    void swap(b_heap & rhs)
    {
        super_t::swap(rhs);
        q_.swap(rhs.q_);
    }

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: N log(N)
     *
     * */
    void merge(b_heap const & rhs)
    {
        q_.reserve(q_.size() + rhs.q_.size());

        for (iterator it = rhs.begin(); it != rhs.end(); ++it)
            push(*it);
    }

    /**
     * \b Effects: Merges all elements from rhs to this. Rhs is cleared.
     *
     * \b Complexity: N log(N)
     *
     * */
    void merge_and_clear(b_heap & rhs)
    {
        merge(rhs);
        rhs.q_.clear();
    }

    /// \copydoc boost::heap::priority_queue::begin
    iterator begin(void) const
    {
        iterator ret = iterator(q_.begin());
        ++ret;
        return ret;
    }

    /// \copydoc boost::heap::priority_queue::end
    iterator end(void) const
    {
        return iterator(q_.end());
    }

    /// \copydoc boost::heap::priority_queue::reserve
    void reserve (size_type element_count)
    {
        q_.reserve(element_count);
    }

private:
    void reset_index(size_type index, size_type new_index)
    {
        assert(index < q_.size());
        index_updater()(q_[index], new_index);
    }

    void siftdown(size_type index)
    {
        while (not_leaf(index))
        {
            size_type max_child_index = max_child(index);
            if (!super_t::operator()(q_[max_child_index], q_[index])) {
                reset_index(index, max_child_index);
                reset_index(max_child_index, index);
                std::swap(q_[max_child_index], q_[index]);
                index = max_child_index;
            }
            else
                return;
        }
    }

    /* returns new index */
    void siftup(size_type index)
    {
        while (index != 1)
        {
            size_type parent = parent_index(index);

            if (super_t::operator()(q_[parent], q_[index])) {
                reset_index(index, parent);
                reset_index(parent, index);
                std::swap(q_[parent], q_[index]);
                index = parent;
            }
            else
                return;
        }
    }

    bool not_leaf(size_type index) const
    {
        const size_t first_child = first_child_index(index);
        return first_child < q_.size();
    }

    size_type max_child(size_type index) const
    {
        const size_t first_index = first_child_index(index);
        const size_t end_index = first_index + 2;
        typename container_type::const_iterator min_element = std::max_element(q_.begin() + first_index,
                                                                               std::min(q_.begin() + end_index,
                                                                                        q_.end()),
                                                                               static_cast<super_t const &>(*this));
        return min_element - q_.begin();
    }

    static const size_type objects_per_page = parameter::binding<BoundArgs,
                                                                 tag::objects_per_page,
                                                                 boost::mpl::int_<8>
                                                                >::type::value;

    BOOST_STATIC_ASSERT(is_power_of_two<objects_per_page>::value);
    BOOST_STATIC_ASSERT(objects_per_page >= 8);

    /* directly adapted from Poul-Henning Kamp */
    static const size_type bh_psize = objects_per_page;
    static const size_type bh_shift = boost::static_log2< bh_psize >::value;
    static const size_type bh_mask = bh_psize - 1;

    static const size_type bh_half = bh_psize / 2;
    static const size_type bh_hshift = bh_shift - 2;
    static const size_type bh_hmask = bh_mask >> 1;

    static size_type bh_pg(size_type index)
    {
        return index >> bh_shift;
    }

    static size_type bh_po(size_type index)
    {
        return index & bh_mask;
    }

    static size_type parent_index(size_type idx)
    {
        size_type po = bh_po(idx);
        if (idx < bh_psize || po > 3)
            return (idx & ~bh_mask) | (po >> 1);

        if (po < 2)
        {
            size_type ip = (idx - bh_psize) >> bh_shift;
            ip += (ip & ~bh_hmask);
            ip |= bh_psize / 2;
            return ip;
        }
        else
            return idx - 2;
    }

    static size_type first_child_index(size_type idx)
    {
        if (idx > bh_mask && !(idx & (bh_mask - 1)))
            return idx + 2;
        if (idx & (bh_psize >> 1))
        {
            size_type i1 = (idx & ~bh_mask) >> 1;
            i1 |= idx & (bh_mask >> 1);
            i1 += 1;
            i1 <<= bh_shift;
            return i1;
        }
        else
            return idx + (idx & bh_mask);
    }

    template<typename U,
             typename V,
             typename W,
             typename X>
    struct rebind {
        typedef b_heap<U, typename b_heap_signature::bind<boost::heap::stable<heap_base_maker::stable>,
                                                                  boost::heap::compare<V>,
                                                                  boost::heap::allocator<W>
                                                                 >::type,
                           X
                          > other;
    };

    template <class U> friend class priority_queue_mutable_wrapper;
    void update(size_type index)
    {
        if (index == 1) {
            siftdown(index);
            return;
        }
        size_type parent = parent_index(index);

        if (super_t::operator()(q_[parent], q_[index]))
            siftup(index);
        else
            siftdown(index);
    }

    void erase(size_type index)
    {
        while (index != 1)
        {
            size_type parent = parent_index(index);

            reset_index(index, parent);
            reset_index(parent, index);
            std::swap(q_[parent], q_[index]);
            index = parent;
        }
        pop();
    }

    void increase(size_type index)
    {
        siftup(index);
    }

    void decrease(size_type index)
    {
        siftdown(index);
    }
};

} /* namespace detail */

/**
 * \class b_heap
 * \brief b-heap class
 *
 * This class implements an immutable priority queue. Internally, the b-ary heap is represented
 * as dynamically sized array (std::vector), that directly stores the values.
 *
 * The b-heap is a cache-aware algorithm, that was proposed by Poul-Henning Kamp in the article
 * "You're Doing It Wrong", published in ACM Queue, June 2010.
 *
 * The template parameter T is the type to be managed by the container. For best efficiency, its size
 * should be a divisor of the cache line size. The user can specify additional options and if no
 * options are provided default options are used. The memory layout of the data structure can be adapted
 * with the \c objects_per_page<> parameter, that specifies the number of objects, that are grouped to a
 * continuous memory region.
 *
 * The container supports the following options:
 * - \c stable<>, defaults to \c stable<false>
 * - \c compare<>, defaults to \c compare<std::less<T> >
 * - \c allocator<>, defaults to \c allocator<std::allocator<T> >
 * - \c objects_per_page<>, defaults to \c 8, has to be a power of two greater or equal 8
 *
 *
 */
#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
#else
template <typename T,
          class A0 = boost::parameter::void_,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_
         >
#endif
class b_heap:
    public detail::b_heap<T, typename detail::b_heap_signature::bind<A0, A1, A2, A3>::type, detail::nop_index_updater<T> >
{
    typedef detail::b_heap<T, typename detail::b_heap_signature::bind<A0, A1, A2, A3>::type, detail::nop_index_updater<T> > super_t;

public:
    static const bool constant_time_size = true;

    b_heap(void)
    {}

    b_heap(b_heap const & rhs):
        super_t(rhs)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    b_heap(b_heap && rhs):
        super_t(std::move(rhs))
    {}

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    b_heap & operator=(b_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        return *this;
    }
#endif
};

/**
 * \class b_heap_mutable
 * \brief mutable b_heap class
 *
 * This class implements a mutable priority queue. The mutability is realized by storing the objects
 * inside a list and maintaining the heap by using pointers to the list. While this introduces some
 * indirection, it may be preferred to use this instead of the b_heap class, if copying the values
 * is not a lightweight operation.
 *
 * The template parameter T is the type to be managed by the container.
 * The user can specify additional options and if no options are provided default options are used.
 *
 * The container supports the following options:
 * - \c stable<>, defaults to \c stable<false>
 * - \c compare<>, defaults to \c compare<std::less<T> >
 * - \c allocator<>, defaults to \c allocator<std::allocator<T> >
 * - \c objects_per_page<> (actually denotes the number of pointers per page), defaults to \c 8
 */
#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
#else
template <typename T,
          class A0 = boost::parameter::void_,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_
         >
#endif
class b_heap_mutable:
    public detail::priority_queue_mutable_wrapper<detail::b_heap<T,
                                                                     typename detail::b_heap_signature::bind<A0, A1, A2, A3>::type,
                                                                     detail::nop_index_updater<T>
                                                                    >
                                                 >
{
    typedef detail::priority_queue_mutable_wrapper<detail::b_heap<T,
                                                              typename detail::b_heap_signature::bind<A0, A1, A2, A3>::type,
                                                              detail::nop_index_updater<T>
                                                             >
                                          > super_t;
public:
    static const bool constant_time_size = true;

    typedef typename super_t::compare_type compare_type;

    typedef typename super_t::value_type value_type;
    typedef typename super_t::size_type size_type;
    typedef typename super_t::allocator_type allocator_type;

    /// \copydoc boost::heap::priority_queue::iterator
    typedef typename super_t::iterator iterator;
    typedef typename super_t::const_iterator const_iterator;
    typedef typename super_t::const_reference const_reference;
    typedef typename super_t::const_pointer const_pointer;

    typedef typename super_t::handle_type handle_type;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit b_heap_mutable(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    b_heap_mutable(b_heap_mutable const & rhs):
        super_t(rhs)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    b_heap_mutable(b_heap_mutable && rhs):
        super_t(std::move(rhs))
    {}

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    b_heap_mutable & operator=(b_heap_mutable && rhs)
    {
        super_t::operator=(std::move(rhs));
        return *this;
    }
#endif

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const &)
    b_heap_mutable & operator=(b_heap_mutable const & rhs)
    {
        super_t::operator=(rhs);
        return *this;
    }

    /// \copydoc boost::heap::priority_queue::reserve
    void reserve(typename super_t::size_type element_count)
    {
        super_t::q_.reserve(element_count);
    }

    /**
     * \b Effects: Casts an iterator to a node handle.
     *
     * \b Complexity: Constant.
     * */
    static handle_type s_handle_from_iterator(iterator const & it)
    {
        return super_t::s_handle_from_iterator(&*it);
    }
};

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(b_heap<T, Options...> const & lhs, b_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(b_heap<T, A0, A1, A2, A3> const & lhs, b_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(b_heap<T, Options...> const & lhs, b_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(b_heap<T, A0, A1, A2, A3> const & lhs, b_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(b_heap<T, Options...> const & lhs, b_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(b_heap<T, A0, A1, A2, A3> const & lhs, b_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(b_heap<T, Options...> const & lhs, b_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(b_heap<T, A0, A1, A2, A3> const & lhs, b_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(b_heap<T, Options...> const & lhs, b_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(b_heap<T, A0, A1, A2, A3> const & lhs, b_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(b_heap<T, Options...> const & lhs, b_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(b_heap<T, A0, A1, A2, A3> const & lhs, b_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}


#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(b_heap_mutable<T, Options...> const & lhs, b_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(b_heap_mutable<T, A0, A1, A2, A3> const & lhs, b_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(b_heap_mutable<T, Options...> const & lhs, b_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(b_heap_mutable<T, A0, A1, A2, A3> const & lhs, b_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(b_heap_mutable<T, Options...> const & lhs, b_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(b_heap_mutable<T, A0, A1, A2, A3> const & lhs, b_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(b_heap_mutable<T, Options...> const & lhs, b_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(b_heap_mutable<T, A0, A1, A2, A3> const & lhs, b_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(b_heap_mutable<T, Options...> const & lhs, b_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(b_heap_mutable<T, A0, A1, A2, A3> const & lhs, b_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(b_heap_mutable<T, Options...> const & lhs, b_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(b_heap_mutable<T, A0, A1, A2, A3> const & lhs, b_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}

} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_B_HEAP_HPP */
