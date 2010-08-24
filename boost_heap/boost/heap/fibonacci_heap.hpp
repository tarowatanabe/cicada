// boost heap: fibonacci heap
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_FIBONACCI_HEAP_HPP
#define BOOST_HEAP_FIBONACCI_HEAP_HPP

#include <algorithm>
#include <vector>

#include <boost/array.hpp>
#include <boost/assert.hpp>

#include "detail/heap_comparison.hpp"
#include "detail/heap_node.hpp"
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
                             > fibonacci_heap_signature;

template <typename T, typename Parspec>
struct make_fibonacci_heap_base
{
    static const bool constant_time_size = parameter::binding<Parspec,
                                                              tag::constant_time_size,
                                                              boost::mpl::true_
                                                             >::type::value;

    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::type base_type;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::allocator_argument allocator_argument;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::compare_argument compare_argument;
    typedef marked_heap_node<typename base_type::internal_type> node_type;

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
 * \class fibonacci_heap
 * \brief fibonacci heap
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
class fibonacci_heap:
    private detail::make_fibonacci_heap_base<T,
                                             typename detail::fibonacci_heap_signature::bind<A0, A1, A2, A3>::type
                                            >::type
{
    typedef detail::fibonacci_heap_signature::bind<A0, A1, A2, A3> bound_args;
    typedef detail::make_fibonacci_heap_base<T, typename bound_args::type> base_maker;
    typedef typename base_maker::type super_t;

    typedef typename super_t::size_holder_type size_holder;
    typedef typename super_t::internal_type internal_type;
    typedef typename base_maker::allocator_argument allocator_argument;

public:
    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::size_t difference_type;

    typedef typename base_maker::compare_argument compare_type;
    typedef typename base_maker::allocator_type allocator_type;

    typedef typename allocator_argument::reference reference;
    typedef typename allocator_argument::const_reference const_reference;
    typedef typename allocator_argument::pointer pointer;
    typedef typename allocator_argument::const_pointer const_pointer;

    static const bool constant_time_size = base_maker::constant_time_size;

private:
    typedef typename allocator_type::pointer node_pointer;
    typedef typename allocator_type::const_pointer const_node_pointer;

    typedef detail::heap_node_list node_list_type;
    typedef typename node_list_type::iterator node_list_iterator;
    typedef typename node_list_type::const_iterator node_list_const_iterator;

    typedef typename base_maker::node_type node;

    typedef detail::value_extractor<value_type, internal_type, super_t> value_extractor;

public:
    typedef detail::node_handle<node_pointer, super_t, reference> handle_type;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit fibonacci_heap(compare_type const & cmp = compare_type()):
        super_t(cmp), top_element(0)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    fibonacci_heap(fibonacci_heap const & rhs):
        super_t(rhs), top_element(0)
    {
        if (rhs.empty())
            return;

        clone_forest(rhs);
    }

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    fibonacci_heap(fibonacci_heap && rhs):
        super_t(std::move(rhs)), top_element(rhs.top_element)
    {
        roots.splice(roots.begin(), rhs.roots);
        rhs.top_element = NULL;
    }

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    fibonacci_heap & operator=(fibonacci_heap && rhs)
    {
        clear();

        super_t::operator=(std::move(rhs));
        roots.splice(roots.begin(), rhs.roots);
        top_element = rhs.top_element;
        rhs.top_element = NULL;
        return *this;
    }
#endif

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const &)
    fibonacci_heap & operator=(fibonacci_heap const & rhs)
    {
        clear();
        size_holder::set_size(rhs.size());
        static_cast<super_t&>(*this) = rhs;

        if (rhs.empty())
            top_element = NULL;
        else
            clone_forest(rhs);
        return *this;
    }

    ~fibonacci_heap(void)
    {
        clear();
    }

    /// \copydoc boost::heap::priority_queue::empty
    bool empty(void) const
    {
        if (constant_time_size)
            return size() == 0;
        else
            return roots.empty();
    }

    /// \copydoc boost::heap::priority_queue::size
    size_type size(void) const
    {
        if (constant_time_size)
            return size_holder::get_size();

        if (empty())
            return 0;
        else
            return detail::count_list_nodes<node, node_list_type>(roots);
    }

    /// \copydoc boost::heap::priority_queue::max_size
    size_type max_size(void) const
    {
        return allocator_type::max_size();
    }

    /// \copydoc boost::heap::priority_queue::clear
    void clear(void)
    {
        typedef detail::node_disposer<node, typename node_list_type::value_type, allocator_type> disposer;
        roots.clear_and_dispose(disposer(*this));

        size_holder::set_size(0);
        top_element = NULL;
    }

    /// \copydoc boost::heap::priority_queue::get_allocator
    allocator_type get_allocator(void) const
    {
        return *this;
    }

    /// \copydoc boost::heap::priority_queue::swap
    void swap(fibonacci_heap & rhs)
    {
        super_t::swap(rhs);
        std::swap(top_element, rhs.top_element);
        roots.swap(rhs.roots);
    }


    /// \copydoc boost::heap::priority_queue::top
    value_type const & top(void) const
    {
        BOOST_ASSERT(!empty());

        return super_t::get_value(top_element->value);
    }

    /**
     * \b Effects: Adds a new element to the priority queue. Returns handle to element
     *
     * \b Complexity: Constant.
     *
     * \b Note: Does not invalidate iterators.
     *
     * */
    handle_type push(const_reference v)
    {
        size_holder::increment();

        node_pointer n = allocator_type::allocate(1);

        new(n) node(super_t::make_node(v));
        roots.push_front(*n);

        if (!top_element || super_t::operator()(top_element->value, n->value))
            top_element = n;
        return handle_type(n);
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

        node_pointer element = top_element;
        roots.erase(node_list_type::s_iterator_to(*element));

        add_children_to_root(element);

        element->~node();
        allocator_type::deallocate(element, 1);

        size_holder::decrement();
        if (!empty())
            consolidate();
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic if current value < v, Constant otherwise.
     *
     * */
    void update (handle_type handle, const_reference v)
    {
        if (super_t::operator()(super_t::get_value(handle.node_->value), v))
            increase(handle, v);
        else
            decrease(handle, v);
    }

    /** \copydoc boost::heap::fibonacci_heap::update(handle_type, const_reference)
     *
     * \b Rationale: The lazy update function is a modification of the traditional update, that just invalidates
     *               the iterator the the object referred to by the handle.
     * */
    void update_lazy(handle_type handle, const_reference v)
    {
        handle.node_->value = super_t::make_node(v);
        update_lazy(handle);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void update (handle_type handle)
    {
        node_pointer n = handle.node_;
        node_pointer parent = n->get_parent();

        if (parent) {
            n->parent = NULL;
            roots.splice(roots.begin(), parent->children, node_list_type::s_iterator_to(*n));
        }
        add_children_to_root(n);
        consolidate();
    }

    /** \copydoc boost::heap::fibonacci_heap::update (handle_type handle)
     *
     * \b Rationale: The lazy update function is a modification of the traditional update, that just invalidates
     *               the iterator the the object referred to by the handle.
     * */
    void update_lazy (handle_type handle)
    {
        node_pointer n = handle.node_;
        node_pointer parent = n->get_parent();

        if (parent) {
            n->parent = NULL;
            roots.splice(roots.begin(), parent->children, node_list_type::s_iterator_to(*n));
        }
        add_children_to_root(n);
    }


     /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Constant.
     *
     * \b Note: The new value is expected to be greater than the current one
     * */
    void increase (handle_type handle, const_reference v)
    {
        handle.node_->value = super_t::make_node(v);
        increase(handle);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \b Complexity: Constant.
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void increase (handle_type handle)
    {
        node_pointer n = handle.node_;

        if (n->parent) {
            if (super_t::operator()(n->get_parent()->value, n->value)) {
                node_pointer parent = n->get_parent();
                cut(n);
                cascading_cut(parent);
            }
        }

        if (super_t::operator()(top_element->value, n->value)) {
            top_element = n;
            return;
        }
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: The new value is expected to be less than the current one
     * */
    void decrease (handle_type handle, const_reference v)
    {
        handle.node_->value = super_t::make_node(v);
        decrease(handle);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \b Complexity: Logarithmic.
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
     * \b Complexity: Logarithmic.
     * */
    void erase(handle_type const & handle)
    {
        node_pointer n = handle.node_;
        node_pointer parent = n->get_parent();

        if (parent)
            parent->children.erase(node_list_type::s_iterator_to(*n));
        else
            roots.erase(node_list_type::s_iterator_to(*n));

        add_children_to_root(n);
        consolidate();

        n->~node();
        allocator_type::deallocate(n, 1);

        size_holder::decrement();
    }

    /// \copydoc boost::heap::priority_queue::iterator
    typedef detail::recursive_tree_iterator<node,
                                            node_list_const_iterator,
                                            value_type,
                                            value_extractor,
                                            detail::list_iterator_converter<node, node_list_type>
                                           > iterator;

    typedef iterator const_iterator;

    /// \copydoc boost::heap::priority_queue::begin
    iterator begin(void) const
    {
        return iterator(roots.begin());
    }

    /// \copydoc boost::heap::priority_queue::end
    iterator end(void) const
    {
        return iterator(roots.end());
    }

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    void merge(fibonacci_heap const & rhs)
    {
        for (iterator it = rhs.begin(); it != rhs.end(); ++it)
            push(*it);
    }

    /**
     * \b Effects: Merge all elements from rhs into this
     *
     * \b Complexity: Constant.
     *
     * */
    void merge_and_clear(fibonacci_heap & rhs)
    {
        size_holder::add(rhs.get_size());

        if (!top_element ||
            (rhs.top_element && super_t::operator()(top_element->value, rhs.top_element->value)))
            top_element = rhs.top_element;

        roots.splice(roots.end(), rhs.roots);

        rhs.set_size(0);
    }

    /// \copydoc boost::heap::d_ary_heap_mutable::s_handle_from_iterator
    static handle_type s_handle_from_iterator(iterator const & it)
    {
        return super_t::s_handle_from_iterator(&*it);
    }

private:
#if !defined(BOOST_DOXYGEN_INVOKED)
    void clone_forest(fibonacci_heap const & rhs)
    {
        assert(roots.empty());
        typedef typename node::template node_cloner<allocator_type> node_cloner;
        roots.clone_from(rhs.roots, node_cloner(*this, NULL), detail::nop_disposer());

        top_element = detail::find_max_child<node_list_type, node, super_t, value_extractor>(roots, super_t::get_cmp());
    }

    void cut(node_pointer n)
    {
        node_pointer parent = n->get_parent();
        roots.splice(roots.begin(), parent->children, node_list_type::s_iterator_to(*n));
        n->parent = 0;
        n->mark = false;
    }

    void cascading_cut(node_pointer n)
    {
        node_pointer parent = n->get_parent();

        if (parent) {
            if (!parent->mark)
                parent->mark = true;
            else {
                cut(n);
                cascading_cut(parent);
            }
        }
    }

    void add_children_to_root(node_pointer n)
    {
        for (node_list_iterator it = n->children.begin(); it != n->children.end(); ++it) {
            node_pointer child = static_cast<node_pointer>(&*it);
            child->parent = 0;
        }

        roots.splice(roots.end(), n->children);
    }

    void consolidate(void)
    {
        static const size_type max_log2 = sizeof(size_type) * 8;
        boost::array<node_pointer, max_log2> aux;
        aux.assign(NULL);

        node_list_iterator it = roots.begin();
        top_element = static_cast<node_pointer>(&*it);

        do {
            node_pointer n = static_cast<node_pointer>(&*it);
            ++it;
            size_type node_rank = n->child_count();

            if (aux[node_rank] == NULL)
                aux[node_rank] = n;
            else {
                do {
                    node_pointer other = aux[node_rank];
                    if (super_t::operator()(n->value, other->value))
                        std::swap(n, other);

                    if (other->parent)
                        n->children.splice(n->children.end(), other->parent->children, node_list_type::s_iterator_to(*other));
                    else
                        n->children.splice(n->children.end(), roots, node_list_type::s_iterator_to(*other));

                    other->parent = n;

                    aux[node_rank] = NULL;
                    node_rank = n->child_count();
                } while (aux[node_rank] != NULL);
                aux[node_rank] = n;
            }

            if (super_t::operator()(top_element->value, n->value))
                top_element = n;
        }
        while (it != roots.end());
    }

    mutable node_pointer top_element;
    node_list_type roots;
#endif
};

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(fibonacci_heap<T, Options...> const & lhs, fibonacci_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(fibonacci_heap<T, A0, A1, A2, A3> const & lhs, fibonacci_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(fibonacci_heap<T, Options...> const & lhs, fibonacci_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(fibonacci_heap<T, A0, A1, A2, A3> const & lhs, fibonacci_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(fibonacci_heap<T, Options...> const & lhs, fibonacci_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(fibonacci_heap<T, A0, A1, A2, A3> const & lhs, fibonacci_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(fibonacci_heap<T, Options...> const & lhs, fibonacci_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(fibonacci_heap<T, A0, A1, A2, A3> const & lhs, fibonacci_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(fibonacci_heap<T, Options...> const & lhs, fibonacci_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(fibonacci_heap<T, A0, A1, A2, A3> const & lhs, fibonacci_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(fibonacci_heap<T, Options...> const & lhs, fibonacci_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(fibonacci_heap<T, A0, A1, A2, A3> const & lhs, fibonacci_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}

} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_FIBONACCI_HEAP_HPP */
