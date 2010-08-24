// boost heap: concepts
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_CONCEPTS_HPP
#define BOOST_HEAP_CONCEPTS_HPP

#include <boost/concept_check.hpp>

namespace boost {
namespace heap {

template <class C>
struct PriorityQueue:
    boost::ForwardContainer<C>
{
    typedef typename C::iterator iterator;
    typedef typename C::allocator_type allocator_type;
    typedef typename C::compare_type compare_type;
    typedef typename C::value_type value_type;


    BOOST_CONCEPT_USAGE(PriorityQueue)
    {
        BOOST_CONCEPT_ASSERT((boost::Assignable<value_type>));
        BOOST_CONCEPT_ASSERT((boost::InputIterator<iterator>));

        BOOST_CONCEPT_ASSERT((boost::Const_BinaryPredicate<compare_type, value_type, value_type>));

        i = c.begin();
        i = c.end();
        c.swap(c2);
        c.clear();
        a = c.get_allocator();

        v = c.top();
        c.pop();
        c.merge(c2);
        c.merge_and_clear(c2);
    }

private:
    iterator i;
    C c, c2;
    allocator_type a;
    typename C::value_type v;
};

template <class C>
struct ImmutablePriorityQueue:
    PriorityQueue<C>
{
    BOOST_CONCEPT_USAGE(ImmutablePriorityQueue)
    {
        typename ImmutablePriorityQueue::value_type v;
        c.push(v);
    }

    C c;
};

template <class C>
struct MutablePriorityQueue:
    PriorityQueue<C>
{
    typedef typename C::handle_type handle_type;

    BOOST_CONCEPT_USAGE(MutablePriorityQueue)
    {
        BOOST_CONCEPT_ASSERT((boost::Assignable<typename MutablePriorityQueue::handle_type>));

        typename MutablePriorityQueue::value_type v;
        typename MutablePriorityQueue::handle_type h = c.push(v);
        c.update(h, v);
        c.increase(h, v);
        c.decrease(h, v);

        c.update(h);
        c.increase(h);
        c.decrease(h);
    }

    C c;
};

}}

#endif /* BOOST_HEAP_CONCEPTS_HPP */
