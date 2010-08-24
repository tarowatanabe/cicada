// boost heap: wrapper for stl heap
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_PRIORITY_QUEUE_HPP
#define BOOST_HEAP_PRIORITY_QUEUE_HPP

#include <algorithm>
#include <queue>
#include <vector>

#include <boost/assert.hpp>

#include "detail/heap_comparison.hpp"
#include "detail/stable_heap.hpp"

namespace boost
{

namespace heap
{

namespace detail
{

typedef parameter::parameters<optional<tag::allocator>,
                              optional<tag::compare>,
                              optional<tag::stable>
                             > priority_queue_signature;
}

/**
 * \class priority_queue
 * \brief priority queue, based on stl heap functions
 *
 * The priority_queue class is a wrapper for the stl heap functions.<br>
 * The template parameter T is the type to be managed by the container.
 * The user can specify additional options and if no options are provided default options are used.
 *
 * The container supports the following options:
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
class priority_queue:
    private detail::make_heap_base<T, typename detail::priority_queue_signature::bind<A0, A1, A2>::type, false>::type
{
    typedef detail::make_heap_base<T, typename detail::priority_queue_signature::bind<A0, A1, A2>::type, false> heap_base_maker;

    typedef typename heap_base_maker::type super_t;
    typedef typename super_t::internal_type internal_type;
    typedef std::vector<internal_type, typename heap_base_maker::allocator_argument> container_type;

    container_type q_;

public:
    typedef T value_type;
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::const_iterator container_iterator;
    typedef typename container_type::allocator_type allocator_type;
    typedef typename heap_base_maker::compare_argument compare_type;

    static const bool constant_time_size = true;

    /**
     * \b Note: The iterator does not traverse the priority queue in order of the priorities.
     * */
    typedef detail::stable_heap_iterator<T, container_iterator, super_t> iterator;
    typedef iterator const_iterator;

    typedef typename allocator_type::reference reference;
    typedef typename allocator_type::const_reference const_reference;
    typedef typename allocator_type::pointer pointer;
    typedef typename allocator_type::const_pointer const_pointer;

    /**
     * \b Effects: constructs an empty priority queue.
     *
     * \b Complexity: Constant.
     *
     * */
    explicit priority_queue(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {}

    /**
     * \b Effects: copy-constructs priority queue from rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    priority_queue (priority_queue const & rhs):
        super_t(rhs), q_(rhs.q_)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /**
     * \b Effects: c++0x-style move constructor.
     *
     * \b Complexity: Constant.
     *
     * \b Note: Only available, if BOOST_HAS_RVALUE_REFS is defined
     * */
    priority_queue(priority_queue && rhs):
        super_t(std::move(rhs)), q_(std::move(rhs.q_))
    {}

    /**
     * \b Effects: c++0x-style move assignment.
     *
     * \b Complexity: Constant.
     *
     * \b Note: Only available, if BOOST_HAS_RVALUE_REFS is defined
     * */
    priority_queue & operator=(priority_queue && rhs)
    {
        super_t::operator=(std::move(rhs));
        q_ = std::move(rhs.q_);
        return *this;
    }
#endif

    /**
     * \b Effects: Assigns priority queue from rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    priority_queue & operator=(priority_queue const & rhs)
    {
        static_cast<super_t&>(*this) = static_cast<super_t const &>(rhs);
        q_ = rhs.q_;
        return *this;
    }

    /**
     * \b Effects: Returns true, if the priority queue contains no elements.
     *
     * \b Complexity: Constant.
     *
     * */
    bool empty(void) const
    {
        return q_.empty();
    }

    /**
     * \b Effects: Returns the number of elements contained in the priority queue.
     *
     * \b Complexity: Constant.
     *
     * */
    size_type size(void) const
    {
        return q_.size();
    }

    /**
     * \b Effects: Returns the maximum number of elements the priority queue can contain.
     *
     * \b Complexity: Constant.
     *
     * */
    size_type max_size(void) const
    {
        return q_.max_size();
    }

    /**
     * \b Effects: Removes all elements from the priority queue.
     *
     * \b Complexity: Linear.
     *
     * */
    void clear(void)
    {
        q_.clear();
    }

    /**
     * \b Effects: Returns allocator.
     *
     * \b Complexity: Constant.
     *
     * */
    allocator_type get_allocator(void) const
    {
        return q_.get_allocator();
    }

    /**
     * \b Effects: Returns a const_reference to the maximum element.
     *
     * \b Complexity: Constant.
     *
     * */
    const_reference top(void) const
    {
        BOOST_ASSERT(!empty());
        return super_t::get_value(q_.front());
    }

    /**
     * \b Effects: Adds a new element to the priority queue.
     *
     * \b Complexity: Logarithmic (amortized). Linear (worst case).
     *
     * */
    void push(const_reference v)
    {
        q_.push_back(super_t::make_node(v));
        std::push_heap(q_.begin(), q_.end(), static_cast<super_t const &>(*this));
    }

    /**
     * \b Effects: Removes the top element from the priority queue.
     *
     * \b Complexity: Logarithmic (amortized). Linear (worst case).
     *
     * */
    void pop(void)
    {
        BOOST_ASSERT(!empty());
        std::pop_heap(q_.begin(), q_.end(), static_cast<super_t const &>(*this));
        q_.pop_back();
    }

    /**
     * \b Effects: Swaps two priority queues.
     *
     * \b Complexity: Constant.
     *
     * */
    void swap(priority_queue & rhs)
    {
        super_t::swap(rhs);
        q_.swap(rhs.q_);
    }

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    void merge(priority_queue const & rhs)
    {
        q_.reserve(q_.size() + rhs.q_.size());

        for (iterator it = rhs.begin(); it != rhs.end(); ++it)
            q_.push_back(super_t::make_node(*it));

        std::make_heap(q_.begin(), q_.end(), static_cast<super_t const &>(*this));
    }

    /**
     * \b Effects: Merges all elements from rhs to this. Rhs is cleared.
     *
     * \b Complexity: Linear.
     *
     * */
    void merge_and_clear(priority_queue & rhs)
    {
        merge(rhs);
        rhs.clear();
    }

    /**
     * \b Effects: Returns an iterator to the first element contained in the priority queue.
     *
     * \b Complexity: Constant.
     *
     * */
    iterator begin(void) const
    {
        return iterator(q_.begin());
    }

    /**
     * \b Effects: Returns an iterator to the end of the priority queue.
     *
     * \b Complexity: Constant.
     *
     * */
    iterator end(void) const
    {
        return iterator(q_.end());
    }

    /**
     * \b Effects: Reserves memory for element_count elements
     *
     * \b Complexity: Linear.
     *
     * \b Node: Invalidates iterators
     *
     * */
    void reserve(size_type element_count)
    {
        q_.reserve(element_count);
    }
};

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(priority_queue<T, Options...> const & lhs, priority_queue<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(priority_queue<T, A0, A1, A2, A3> const & lhs, priority_queue<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(priority_queue<T, Options...> const & lhs, priority_queue<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(priority_queue<T, A0, A1, A2, A3> const & lhs, priority_queue<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(priority_queue<T, Options...> const & lhs, priority_queue<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(priority_queue<T, A0, A1, A2, A3> const & lhs, priority_queue<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(priority_queue<T, Options...> const & lhs, priority_queue<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(priority_queue<T, A0, A1, A2, A3> const & lhs, priority_queue<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(priority_queue<T, Options...> const & lhs, priority_queue<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(priority_queue<T, A0, A1, A2, A3> const & lhs, priority_queue<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(priority_queue<T, Options...> const & lhs, priority_queue<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(priority_queue<T, A0, A1, A2, A3> const & lhs, priority_queue<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}

} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_PRIORITY_QUEUE_HPP */
