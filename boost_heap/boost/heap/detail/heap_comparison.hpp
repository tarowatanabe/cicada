// boost heap: heap node helper classes
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_DETAIL_HEAP_COMPARISON_HPP
#define BOOST_HEAP_DETAIL_HEAP_COMPARISON_HPP

#include <boost/concept/assert.hpp>

#include "../heap_concepts.hpp"

namespace boost
{
namespace heap
{
namespace detail
{

template <typename Heap1,
          typename Heap2
         >
bool heap_equality(Heap1 const & lhs, Heap2 const & rhs)
{
    BOOST_CONCEPT_ASSERT((boost::heap::PriorityQueue<Heap1>));
    BOOST_CONCEPT_ASSERT((boost::heap::PriorityQueue<Heap2>));

    if (Heap1::constant_time_size && Heap1::constant_time_size)
        if (lhs.size() != rhs.size())
            return false;

    if (lhs.empty() && rhs.empty())
        return true;

    Heap1 lhs_copy(lhs);
    Heap2 rhs_copy(rhs);

    for (;;)
    {
        if (lhs_copy.top() != rhs_copy.top())
            return false;

        lhs_copy.pop();
        rhs_copy.pop();

        if (lhs_copy.empty() && rhs_copy.empty())
            return true;

        if (lhs_copy.empty())
            return false;

        if (rhs_copy.empty())
            return false;
    }
}

template <typename Heap1,
          typename Heap2
         >
bool heap_compare(Heap1 const & lhs, Heap2 const & rhs)
{
    typename Heap1::size_type left_size = lhs.size();
    typename Heap2::size_type right_size = rhs.size();
    if (left_size < right_size)
        return true;

    if (left_size > right_size)
        return false;

    Heap1 lhs_copy(lhs);
    Heap2 rhs_copy(rhs);

    for (;;)
    {
        if (lhs_copy.top() < rhs_copy.top())
            return true;

        if (lhs_copy.top() > rhs_copy.top())
            return false;

        lhs_copy.pop();
        rhs_copy.pop();

        if (lhs_copy.empty() && rhs_copy.empty())
            return false;
    }
}


}
}
}

#endif // BOOST_HEAP_DETAIL_HEAP_COMPARISON_HPP
