// boost heap: pairing heap
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_PAIRING_HEAP_HPP
#define BOOST_HEAP_PAIRING_HEAP_HPP

#include <algorithm>
#include <vector>

#include <boost/assert.hpp>

#include "detail/heap_comparison.hpp"
#include "detail/heap_node.hpp"
#include "detail/parameter.hpp"
#include "detail/stable_heap.hpp"
#include "detail/tree_iterator.hpp"


namespace boost
{

namespace heap
{

namespace detail
{

typedef parameter::parameters<optional<tag::allocator>,
                              optional<tag::compare>,
                              optional<tag::stable>,
                              optional<tag::constant_time_size>
                             > pairing_heap_signature;

template <typename T, typename Parspec>
struct make_pairing_heap_base
{
    static const bool constant_time_size = parameter::binding<Parspec,
                                                              tag::constant_time_size,
                                                              boost::mpl::true_
                                                             >::type::value;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::type base_type;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::allocator_argument allocator_argument;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::compare_argument compare_argument;

    typedef heap_node<typename base_type::internal_type, false> node_type;

    typedef typename allocator_argument::template rebind<node_type>::other allocator_type;

    struct type:
        base_type,
        allocator_type
    {
        type(compare_argument const & arg):
            base_type(arg)
        {}

#ifdef BOOST_HAS_RVALUE_REFS
        type(type && rhs):
            base_type(std::move(static_cast<base_type&>(rhs))),
            allocator_type(std::move(static_cast<allocator_type&>(rhs)))
        {}

        type & operator=(type && rhs)
        {
            base_type::operator=(std::move(static_cast<base_type&>(rhs)));
            allocator_type::operator=(std::move(static_cast<allocator_type&>(rhs)));
        }
#endif
    };
};

}

/**
 * \class pairing_heap
 * \brief pairing heap
 *
 * Pairing heaps are self-adjusting binary heaps. Although design and implementation are rather simple,
 * the complexity analysis is yet unsolved. For details, consult:
 *
 * Pettie, Seth (2005), "Towards a final analysis of pairing heaps",
 * Proc. 46th Annual IEEE Symposium on Foundations of Computer Science, pp. 174–183
 *
 * The template parameter T is the type to be managed by the container.
 * The user can specify additional options and if no options are provided default options are used.
 *
 * The container supports the following options:
 * - \c stable<>, defaults to \c stable<false>
 * - \c compare<>, defaults to \c compare<std::less<T> >
 * - \c allocator<>, defaults to \c allocator<std::allocator<T> >
 * - \c constant_time_size<>, defaults to \c constant_time_size<true>
 *
 * \b warning: due to the recursive nature of this data structure, it may not be suited for hundred thousands
 *             elements. beware of stack overflows
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
class pairing_heap:
    private detail::make_pairing_heap_base<T,
                                           typename detail::pairing_heap_signature::bind<A0, A1, A2, A3>::type
                                          >::type
{
    typedef detail::pairing_heap_signature::bind<A0, A1, A2, A3> bound_args;
    typedef detail::make_pairing_heap_base<T, typename bound_args::type> base_maker;
    typedef typename base_maker::type super_t;

    typedef typename super_t::internal_type internal_type;
    typedef typename super_t::size_holder_type size_holder;
    typedef typename base_maker::allocator_argument allocator_argument;

public:
    typedef typename base_maker::compare_argument compare_type;
    typedef typename base_maker::allocator_type allocator_type;

    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::size_t difference_type;

    typedef typename allocator_argument::reference reference;
    typedef typename allocator_argument::const_reference const_reference;
    typedef typename allocator_argument::pointer pointer;
    typedef typename allocator_argument::const_pointer const_pointer;

private:
    typedef typename allocator_type::pointer node_pointer;
    typedef typename allocator_type::const_pointer const_node_pointer;

    typedef typename base_maker::node_type node;
    typedef detail::value_extractor<value_type, internal_type, super_t> value_extractor;

public:
    typedef detail::node_handle<node_pointer, super_t, reference> handle_type;
    static const bool constant_time_size = super_t::constant_time_size;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit pairing_heap(compare_type const & cmp = compare_type()):
        super_t(cmp), root(NULL)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    pairing_heap(pairing_heap const & rhs):
        super_t(rhs), root(NULL)
    {
        if (rhs.empty())
            return;

        clone_tree(rhs);
    }

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    pairing_heap(pairing_heap && rhs):
        super_t(std::move(rhs)), root(rhs.root)
    {
        rhs.root = NULL;
    }

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    pairing_heap & operator=(pairing_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        root = rhs.root;
        rhs.root = NULL;
        return *this;
    }
#endif

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const & rhs)
    pairing_heap & operator=(pairing_heap const & rhs)
    {
        clear();
        size_holder::set_size(rhs.get_size());
        static_cast<super_t&>(*this) = rhs;

        clone_tree(rhs);
        return *this;
    }

    ~pairing_heap(void)
    {
        while (!empty())
            pop();
    }

    /// \copydoc boost::heap::priority_queue::empty
    bool empty(void) const
    {
        return root == NULL;
    }

    /// \copydoc boost::heap::binomial_heap::size
    size_type size(void) const
    {
        if (constant_time_size)
            return size_holder::get_size();

        if (root == NULL)
            return 0;
        else
            return detail::count_nodes(root);
    }

    /// \copydoc boost::heap::priority_queue::max_size
    size_type max_size(void) const
    {
        return allocator_type::max_size();
    }

    /// \copydoc boost::heap::priority_queue::clear
    void clear(void)
    {
        if (empty())
            return;

        root->template clear_subtree<allocator_type>(*this);
        root->~node();
        allocator_type::deallocate(root, 1);
        root = NULL;
        size_holder::set_size(0);
    }

    /// \copydoc boost::heap::priority_queue::get_allocator
    allocator_type get_allocator(void) const
    {
        return *this;
    }

    /// \copydoc boost::heap::priority_queue::swap
    void swap(pairing_heap & rhs)
    {
        super_t::swap(rhs);
        std::swap(root, rhs.root);
    }


    /// \copydoc boost::heap::priority_queue::top
    const_reference top(void) const
    {
        BOOST_ASSERT(!empty());

        return super_t::get_value(root->value);
    }

    /**
     * \b Effects: Adds a new element to the priority queue. Returns handle to element
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * */
    handle_type push(const_reference v)
    {
        size_holder::increment();

        node_pointer n = allocator_type::allocate(1);

        new(n) node(super_t::make_node(v));

        merge_node(n);
        return handle_type(n);
    }

    /**
     * \b Effects: Removes the top element from the priority queue.
     *
     * \b Complexity: Logarithmic (amortized).
     *
     * */
    void pop(void)
    {
        BOOST_ASSERT(!empty());

        erase(handle_type(root));
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * */
    void update (handle_type handle, const_reference v)
    {
        handle.node_->value = super_t::make_node(v);
        update(handle);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void update (handle_type handle)
    {
        node_pointer n = handle.node_;

        n->unlink();
        if (!n->children.empty())
            n = merge_nodes(n, merge_children(n));

        if (n != root)
            merge_node(n);
    }

     /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * \b Note: The new value is expected to be greater than the current one
     * */
    void increase (handle_type handle, const_reference v)
    {
        update(handle, v);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void increase (handle_type handle)
    {
        update(handle);
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * \b Note: The new value is expected to be less than the current one
     * */
    void decrease (handle_type handle, const_reference v)
    {
        update(handle, v);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * \b Note: The new value is expected to be less than the current one. If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void decrease (handle_type handle)
    {
        update(handle);
    }

    /**
     * \b Effects: Removes the element handled by \c handle from the priority_queue.
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     * */
    void erase(handle_type handle)
    {
        node_pointer n = handle.node_;
        if (n != root)
        {
            n->unlink();
            if (!n->children.empty())
                merge_node(merge_children(n));
        }
        else
        {
            if (!n->children.empty())
                root = merge_children(n);
            else
                root = NULL;
        }

        size_holder::decrement();
        n->~node();
        allocator_type::deallocate(n, 1);
    }

    /// \copydoc boost::heap::priority_queue::iterator
    typedef detail::tree_iterator<node,
                                  const value_type,
                                  allocator_type,
                                  value_extractor,
                                  detail::pointer_to_reference<node>,
                                  false
                                 > iterator;

    typedef iterator const_iterator;

    /// \copydoc boost::heap::priority_queue::begin
    iterator begin(void) const
    {
        return iterator(root);
    }

    /// \copydoc boost::heap::priority_queue::end
    iterator end(void) const
    {
        return iterator();
    }

    /// \copydoc boost::heap::d_ary_heap_mutable::s_handle_from_iterator
    static handle_type s_handle_from_iterator(iterator const & it)
    {
        return super_t::s_handle_from_iterator(&*it);
    }

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    void merge(pairing_heap const & rhs)
    {
        pairing_heap copy(rhs);
        merge_and_clear(copy);
    }

    /**
     * \b Effects: Merge all elements from rhs into this
     *
     * \cond
     * \b Complexity: \f$2^2log(log(N))\f$ (amortized).
     * \endcond
     *
     * \b Complexity: 2**2*log(log(N)) (amortized).
     *
     * */
    void merge_and_clear(pairing_heap & rhs)
    {
        if (rhs.empty())
            return;

        merge_node(rhs.root);

        size_holder::add(rhs.get_size());
        rhs.set_size(0);
        rhs.root = NULL;
    }

private:
#if !defined(BOOST_DOXYGEN_INVOKED)
    void clone_tree(pairing_heap const & rhs)
    {
        assert(root == NULL);
        if (rhs.empty())
            return;

        root = allocator_type::allocate(1);

        new(root) node(*rhs.root, static_cast<allocator_type&>(*this));
    }

    void merge_node(node_pointer other)
    {
        assert(other);
        if (root != NULL)
            root = merge_nodes(root, other);
        else
            root = other;
    }

    node_pointer merge_children(node_pointer n)
    {
        assert(!n->children.empty());

        node_pointer first_child = static_cast<node_pointer>(&n->children.front());
        n->children.pop_front();
        if (n->children.empty())
            return first_child;

        node_pointer second_child = static_cast<node_pointer>(&n->children.front());
        n->children.pop_front();

        node_pointer merged = merge_nodes(first_child, second_child);

        if (n->children.empty())
            return merged;
        else
            return merge_nodes(merged, merge_children(n));
    }

    node_pointer merge_nodes(node_pointer node1, node_pointer node2)
    {
        if (super_t::operator()(node1->value, node2->value))
            std::swap(node1, node2);

        node2->unlink();
        node1->children.push_front(*node2);
        return node1;
    }

    node_pointer root;
#endif
};

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(pairing_heap<T, Options...> const & lhs, pairing_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(pairing_heap<T, A0, A1, A2, A3> const & lhs, pairing_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(pairing_heap<T, Options...> const & lhs, pairing_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(pairing_heap<T, A0, A1, A2, A3> const & lhs, pairing_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(pairing_heap<T, Options...> const & lhs, pairing_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(pairing_heap<T, A0, A1, A2, A3> const & lhs, pairing_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(pairing_heap<T, Options...> const & lhs, pairing_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(pairing_heap<T, A0, A1, A2, A3> const & lhs, pairing_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(pairing_heap<T, Options...> const & lhs, pairing_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(pairing_heap<T, A0, A1, A2, A3> const & lhs, pairing_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(pairing_heap<T, Options...> const & lhs, pairing_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(pairing_heap<T, A0, A1, A2, A3> const & lhs, pairing_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}

} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_PAIRING_HEAP_HPP */
