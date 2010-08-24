// boost heap: d-ary heap as containter adaptor
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_D_ARY_HEAP_HPP
#define BOOST_HEAP_D_ARY_HEAP_HPP

#include <algorithm>
#include <vector>

#include <boost/assert.hpp>

#include "detail/heap_comparison.hpp"
#include "detail/stable_heap.hpp"
#include "detail/mutable_heap.hpp"


namespace boost
{
namespace heap
{
namespace detail
{

template <typename T>
struct nop_index_updater
{
    void operator()(T &, std::size_t) const
    {}

    template <typename U>
    struct rebind {
        typedef nop_index_updater<U> other;
    };
};


typedef parameter::parameters<parameter::required<tag::arity>,
                              optional<tag::allocator>,
                              optional<tag::compare>,
                              optional<tag::stable>
                             > d_ary_heap_signature;


/* base class for d-ary heap */
template <typename T,
          class BoundArgs,
          class IndexUpdater>
class d_ary_heap:
    private make_heap_base<T, BoundArgs, false>::type
{
    typedef make_heap_base<T, BoundArgs, false> heap_base_maker;

    typedef typename heap_base_maker::type super_t;
    typedef typename super_t::internal_type internal_type;

    typedef std::vector<internal_type, typename heap_base_maker::allocator_argument> container_type;
    typedef typename container_type::const_iterator container_iterator;

    typedef typename IndexUpdater::template rebind<internal_type>::other index_updater;

    container_type q_;

    static const unsigned int D = parameter::binding<BoundArgs, tag::arity>::type::value;

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
    explicit d_ary_heap(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    d_ary_heap(d_ary_heap const & rhs):
        super_t(rhs), q_(rhs.q_)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    d_ary_heap(d_ary_heap && rhs):
        super_t(std::move(rhs)), q_(std::move(rhs.q_))
    {}

    d_ary_heap & operator=(d_ary_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        q_ = std::move(rhs.q_);
        return *this;
    }
#endif

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const &)
    d_ary_heap & operator=(d_ary_heap const & rhs)
    {
        static_cast<super_t&>(*this) = static_cast<super_t const &>(rhs);
        q_ = rhs.q_;
        return *this;
    }

    /// \copydoc boost::heap::priority_queue::empty
    bool empty(void) const
    {
        return q_.empty();
    }

    /// \copydoc boost::heap::priority_queue::size
    size_type size(void) const
    {
        return q_.size();
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
        return super_t::get_value(q_.front());
    }

    /// \copydoc boost::heap::priority_queue::push
    void push(value_type const & v)
    {
        q_.push_back(super_t::make_node(v));
        reset_index(size() - 1, size() - 1);
        siftup(q_.size() - 1);
    }

    /// \copydoc boost::heap::priority_queue::pop
    void pop(void)
    {
        BOOST_ASSERT(!empty());
        std::swap(q_.front(), q_.back());
        q_.pop_back();

        if (q_.empty())
            return;

        reset_index(0, 0);
        siftdown(0);
    }

    /// \copydoc boost::heap::priority_queue::swap
    void swap(d_ary_heap & rhs)
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
    void merge(d_ary_heap const & rhs)
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
    void merge_and_clear(d_ary_heap & rhs)
    {
        merge(rhs);
        rhs.q_.clear();
    }

    /// \copydoc boost::heap::priority_queue::begin
    iterator begin(void) const
    {
        return iterator(q_.begin());
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
        while (index != 0)
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
        const size_t end_index = first_index + D;
        typename container_type::const_iterator min_element = std::max_element(q_.begin() + first_index,
                                                                               std::min(q_.begin() + end_index,
                                                                                        q_.end()),
                                                                               static_cast<super_t const &>(*this));
        return min_element - q_.begin();
    }

    static size_type parent_index(size_type index)
    {
        return (index - 1) / D;
    }

    static size_type first_child_index(size_type index)
    {
        return index * D + 1;
    }

    template<typename U,
             typename V,
             typename W,
             typename X>
    struct rebind {
        typedef d_ary_heap<U, typename d_ary_heap_signature::bind<boost::heap::stable<heap_base_maker::stable>,
                                                                  boost::heap::arity<D>,
                                                                  boost::heap::compare<V>,
                                                                  boost::heap::allocator<W>
                                                                 >::type,
                           X
                          > other;
    };

    template <class U> friend class priority_queue_mutable_wrapper;

    void update(size_type index)
    {
        if (index == 0) {
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
        while (index != 0)
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
 * \class d_ary_heap
 * \brief d-ary heap class
 *
 * This class implements an immutable priority queue. Internally, the d-ary heap is represented
 * as dynamically sized array (std::vector), that directly stores the values.
 *
 * The template parameter T is the type to be managed by the container.
 * The user can specify additional options and if no options are provided default options are used.
 *
 * The container supports the following options:
 * - \c arity<>, required
 * - \c stable<>, defaults to \c stable<false>
 * - \c compare<>, defaults to \c compare<std::less<T> >
 * - \c allocator<>, defaults to \c allocator<std::allocator<T> >
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
class d_ary_heap:
    public detail::d_ary_heap<T, typename detail::d_ary_heap_signature::bind<A0, A1, A2, A3>::type, detail::nop_index_updater<T> >
{
    typedef detail::d_ary_heap<T, typename detail::d_ary_heap_signature::bind<A0, A1, A2, A3>::type, detail::nop_index_updater<T> > super_t;

public:
    static const bool constant_time_size = true;

    d_ary_heap(void)
    {}

    d_ary_heap(d_ary_heap const & rhs):
        super_t(rhs)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    d_ary_heap(d_ary_heap && rhs):
        super_t(std::move(rhs))
    {}

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    d_ary_heap & operator=(d_ary_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        return *this;
    }
#endif
};

/**
 * \class d_ary_heap_mutable
 * \brief mutable d-ary heap class
 *
 * This class implements a mutable priority queue. The mutability is realized by storing the objects
 * inside a list and maintaining the heap by using pointers to the list. While this introduces some
 * indirection, it may be preferred to use this instead of the d_ary_heap class, if copying the values
 * is not a lightweight operation.
 *
 * The template parameter T is the type to be managed by the container.
 * The user can specify additional options and if no options are provided default options are used.
 *
 * The container supports the following options:
 * - \c arity<>, required
 * - \c stable<>, defaults to \c stable<false>
 * - \c compare<>, defaults to \c compare<std::less<T> >
 * - \c allocator<>, defaults to \c allocator<std::allocator<T> >
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
class d_ary_heap_mutable:
    public detail::priority_queue_mutable_wrapper<detail::d_ary_heap<T,
                                                                     typename detail::d_ary_heap_signature::bind<A0, A1, A2, A3>::type,
                                                                     detail::nop_index_updater<T>
                                                                    >
                                                 >
{
    typedef detail::priority_queue_mutable_wrapper<detail::d_ary_heap<T,
                                                              typename detail::d_ary_heap_signature::bind<A0, A1, A2, A3>::type,
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
    explicit d_ary_heap_mutable(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    d_ary_heap_mutable(d_ary_heap_mutable const & rhs):
        super_t(rhs)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    d_ary_heap_mutable(d_ary_heap_mutable && rhs):
        super_t(std::move(rhs))
    {}

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    d_ary_heap_mutable & operator=(d_ary_heap_mutable && rhs)
    {
        super_t::operator=(std::move(rhs));
        return *this;
    }
#endif

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const &)
    d_ary_heap_mutable & operator=(d_ary_heap_mutable const & rhs)
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
bool operator==(d_ary_heap<T, Options...> const & lhs, d_ary_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(d_ary_heap<T, A0, A1, A2, A3> const & lhs, d_ary_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(d_ary_heap<T, Options...> const & lhs, d_ary_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(d_ary_heap<T, A0, A1, A2, A3> const & lhs, d_ary_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(d_ary_heap<T, Options...> const & lhs, d_ary_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(d_ary_heap<T, A0, A1, A2, A3> const & lhs, d_ary_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(d_ary_heap<T, Options...> const & lhs, d_ary_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(d_ary_heap<T, A0, A1, A2, A3> const & lhs, d_ary_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(d_ary_heap<T, Options...> const & lhs, d_ary_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(d_ary_heap<T, A0, A1, A2, A3> const & lhs, d_ary_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(d_ary_heap<T, Options...> const & lhs, d_ary_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(d_ary_heap<T, A0, A1, A2, A3> const & lhs, d_ary_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}


#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(d_ary_heap_mutable<T, Options...> const & lhs, d_ary_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(d_ary_heap_mutable<T, A0, A1, A2, A3> const & lhs, d_ary_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(d_ary_heap_mutable<T, Options...> const & lhs, d_ary_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(d_ary_heap_mutable<T, A0, A1, A2, A3> const & lhs, d_ary_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(d_ary_heap_mutable<T, Options...> const & lhs, d_ary_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(d_ary_heap_mutable<T, A0, A1, A2, A3> const & lhs, d_ary_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(d_ary_heap_mutable<T, Options...> const & lhs, d_ary_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(d_ary_heap_mutable<T, A0, A1, A2, A3> const & lhs, d_ary_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(d_ary_heap_mutable<T, Options...> const & lhs, d_ary_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(d_ary_heap_mutable<T, A0, A1, A2, A3> const & lhs, d_ary_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(d_ary_heap_mutable<T, Options...> const & lhs, d_ary_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(d_ary_heap_mutable<T, A0, A1, A2, A3> const & lhs, d_ary_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}


} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_D_ARY_HEAP_HPP */
