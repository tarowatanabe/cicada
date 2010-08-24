// boost heap: binomial heap
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_BINOMIAL_HEAP_HPP
#define BOOST_HEAP_BINOMIAL_HEAP_HPP

#include <algorithm>
#include <vector>

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
                             > binomial_heap_signature;

template <typename T, typename Parspec>
struct make_binomial_heap_base
{
    static const bool constant_time_size = parameter::binding<Parspec,
                                                              tag::constant_time_size,
                                                              boost::mpl::true_
                                                             >::type::value;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::type base_type;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::allocator_argument allocator_argument;
    typedef typename detail::make_heap_base<T, Parspec, constant_time_size>::compare_argument compare_argument;

    typedef parent_pointing_heap_node<typename base_type::internal_type> node_type;

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
 * \class binomial_heap
 * \brief binomial heap
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
class binomial_heap:
    private detail::make_binomial_heap_base<T,
                                            typename detail::binomial_heap_signature::bind<A0, A1, A2, A3>::type
                                           >::type
{
    typedef detail::binomial_heap_signature::bind<A0, A1, A2, A3> bound_args;
    typedef detail::make_binomial_heap_base<T, typename bound_args::type> base_maker;
    typedef typename base_maker::type super_t;

    typedef typename super_t::internal_type internal_type;
    typedef typename super_t::size_holder_type size_holder;
    typedef typename base_maker::allocator_argument allocator_argument;

public:
    static const bool constant_time_size = super_t::constant_time_size;

    typedef T value_type;
    typedef std::size_t size_type;
    typedef std::size_t difference_type;

    typedef typename base_maker::compare_argument compare_type;
    typedef typename base_maker::allocator_type allocator_type;

    typedef typename allocator_argument::reference reference;
    typedef typename allocator_argument::const_reference const_reference;
    typedef typename allocator_argument::pointer pointer;
    typedef typename allocator_argument::const_pointer const_pointer;

private:
    typedef typename allocator_type::pointer node_pointer;
    typedef typename allocator_type::const_pointer const_node_pointer;

    typedef typename base_maker::node_type node;

    typedef typename node::child_list node_list_type;
    typedef typename node_list_type::iterator node_list_iterator;
    typedef typename node_list_type::const_iterator node_list_const_iterator;

    typedef detail::value_extractor<value_type, internal_type, super_t> value_extractor;

public:
    typedef detail::node_handle<node_pointer, super_t, reference> handle_type;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit binomial_heap(compare_type const & cmp = compare_type()):
        super_t(cmp), top_element(0)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    binomial_heap(binomial_heap const & rhs):
        super_t(rhs), top_element(0)
    {
        if (rhs.empty())
            return;

        clone_forest(rhs);
    }

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const &)
    binomial_heap & operator=(binomial_heap const & rhs)
    {
        clear();
        size_holder::set_size(rhs.get_size());
        static_cast<super_t&>(*this) = rhs;

        if (rhs.empty())
            top_element = NULL;
        else
            clone_forest(rhs);
        return *this;
    }

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    binomial_heap(binomial_heap && rhs):
        super_t(std::move(rhs)), top_element(rhs.top_element)
    {
        trees.splice(trees.begin(), rhs.trees);
        rhs.top_element = NULL;
    }

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    binomial_heap & operator=(binomial_heap && rhs)
    {
        clear();
        super_t::operator=(std::move(rhs));
        trees.splice(trees.begin(), rhs.trees);
        top_element = rhs.top_element;
        rhs.top_element = NULL;
        return *this;
    }
#endif

    ~binomial_heap(void)
    {
        clear();
    }

    /// \copydoc boost::heap::priority_queue::empty
    bool empty(void) const
    {
        return top_element == NULL;
    }

    /**
     * \b Effects: Returns the number of elements contained in the priority queue.
     *
     * \b Complexity: Constant, if configured with constant_time_size<true>, otherwise linear.
     *
     * */
    size_type size(void) const
    {
        if (constant_time_size)
            return size_holder::get_size();

        if (empty())
            return 0;
        else
            return detail::count_list_nodes<node, node_list_type>(trees);
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
        trees.clear_and_dispose(disposer(*this));

        size_holder::set_size(0);
        top_element = NULL;
    }

    /// \copydoc boost::heap::priority_queue::get_allocator
    allocator_type get_allocator(void) const
    {
        return *this;
    }

    /// \copydoc boost::heap::priority_queue::swap
    void swap(binomial_heap & rhs)
    {
        super_t::swap(rhs);
        std::swap(top_element, rhs.top_element);
        trees.swap(rhs.trees);
    }

    /// \copydoc boost::heap::priority_queue::top
    const_reference top(void) const
    {
        BOOST_ASSERT(!empty());

        return super_t::get_value(top_element->value);
    }

    /**
     * \b Effects: Adds a new element to the priority queue. Returns handle to element
     *
     * \b Complexity: Logarithmic.
     *
     * */
    handle_type push(const_reference v)
    {
        node_pointer n = allocator_type::allocate(1);

        new(n) node(super_t::make_node(v));

        insert_node(trees.begin(), n);

        if (!top_element || super_t::operator()(top_element->value, n->value))
            top_element = n;

        size_holder::increment();
        sanity_check();
        return handle_type(n);
    }

    /**
     * \b Effects: Removes the top element from the priority queue.
     *
     * \b Complexity: Logarithmic.
     *
     * */
    void pop(void)
    {
        BOOST_ASSERT(!empty());

        node_pointer element = top_element;

        trees.erase(node_list_type::s_iterator_to(*element));
        size_holder::decrement();

        if (element->child_count())
        {
            size_type sz = (1 << element->child_count()) - 1;
            binomial_heap children(element->children, sz);
            if (trees.empty())
                swap(children);
            else
                merge_and_clear_nodes(children);
        }

        sanity_check();

        if (trees.empty())
            top_element = NULL;
        else
            top_element = detail::find_max_child<node_list_type, node, super_t, value_extractor>(trees, super_t::get_cmp());

        element->~node();
        allocator_type::deallocate(element, 1);
        sanity_check();
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic.
     *
     * */
    void update (handle_type handle, const_reference v)
    {
        if (super_t::operator()(super_t::get_value(handle.node_->value), v))
            increase(handle, v);
        else
            decrease(handle, v);
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
        node_pointer this_node = handle.node_;

        if (this_node->parent)
        {
            if (super_t::operator()(super_t::get_value(this_node->parent->value), super_t::get_value(this_node->value)))
                increase(handle);
            else
                decrease(handle);
        }
        else
            decrease(handle);
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic.
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
     * \b Complexity: Logarithmic.
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void increase (handle_type handle)
    {
        node_pointer n = handle.node_;
        siftup(n, *this);

        top_element = detail::find_max_child<node_list_type, node, super_t, value_extractor>(trees, super_t::get_cmp());
        sanity_check();
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
        node_pointer n = handle.node_;

        siftdown(n);

        if (n == top_element)
            top_element = detail::find_max_child<node_list_type, node, super_t, value_extractor>(trees, super_t::get_cmp());
    }

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    void merge(binomial_heap const & rhs)
    {
        binomial_heap copy(rhs);
        merge_and_clear(copy);
    }

    /**
     * \b Effects: Merge all elements from rhs into this
     *
     * \b Complexity: Logarithmic.
     *
     * */
    void merge_and_clear(binomial_heap & rhs)
    {
        if (rhs.empty())
            return;

        if (empty()) {
            swap(rhs);
            return;
        }

        size_type new_size = size_holder::get_size() + rhs.get_size();
        merge_and_clear_nodes(rhs);

        size_holder::set_size(new_size);
        rhs.set_size(0);
        rhs.top_element = NULL;
    }

public:
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
        return iterator(trees.begin());
    }

    /// \copydoc boost::heap::priority_queue::end
    iterator end(void) const
    {
        return iterator(trees.end());
    }

    /**
     * \b Effects: Removes the element handled by \c handle from the priority_queue.
     *
     * \b Complexity: Logarithmic.
     * */
    void erase(handle_type handle)
    {
        node_pointer n = handle.node_;
        siftup(n, force_inf());
        top_element = n;
        pop();
    }

    /// \copydoc boost::heap::d_ary_heap_mutable::s_handle_from_iterator
    static handle_type s_handle_from_iterator(iterator const & it)
    {
        return handle_type(&*it);
    }

private:
#if !defined(BOOST_DOXYGEN_INVOKED)
    void merge_and_clear_nodes(binomial_heap & rhs)
    {
        assert (!empty());
        assert (!rhs.empty());

        node_list_iterator this_iterator = trees.begin();
        node_list_iterator rhs_iterator = rhs.trees.begin();
        node_pointer carry_node = NULL;

        while (!rhs.trees.empty())
        {
            node_pointer rhs_node = static_cast<node_pointer>(&rhs.trees.front());
            size_type rhs_degree = rhs_node->child_count();

            if (super_t::operator()(top_element->value, rhs_node->value))
                top_element = rhs_node;

        try_again:
            node_pointer this_node = static_cast<node_pointer>(&*this_iterator);
            size_type this_degree = this_node->child_count();
            sorted_by_degree();
            rhs.sorted_by_degree();

            if (this_degree == rhs_degree)
            {
                if (carry_node)
                {
                    if (carry_node->child_count() < this_degree) {
                        trees.insert(this_iterator, *carry_node);
                        carry_node = NULL;
                    }
                    else
                    {
                        rhs.trees.pop_front();
                        carry_node = merge_trees(carry_node, rhs_node);
                    }
                    ++this_iterator;
                }
                else
                {
                    this_iterator = trees.erase(this_iterator);
                    rhs.trees.pop_front();
                    carry_node = merge_trees(this_node, rhs_node);
                }

                if (this_iterator == trees.end())
                    break;
                else
                    continue;
            }

            if (this_degree < rhs_degree)
            {
                if (carry_node)
                {
                    if (carry_node->child_count() < this_degree)
                    {
                        trees.insert(this_iterator, *carry_node);
                        carry_node = NULL;
                        ++this_iterator;
                    }
                    else if (carry_node->child_count() == rhs_degree)
                    {
                        rhs.trees.pop_front();
                        carry_node = merge_trees(carry_node, rhs_node);
                        continue;
                    }
                    else
                    {
                        this_iterator = trees.erase(this_iterator);
                        carry_node = merge_trees(this_node, carry_node);
                    }
                    goto try_again;
                }
                else
                {
                    ++this_iterator;
                    if (this_iterator == trees.end())
                        break;
                    goto try_again;
                }

                if (this_iterator == trees.end())
                    break;
                else
                    continue;
            }

            if (this_degree > rhs_degree)
            {
                rhs.trees.pop_front();
                if (carry_node)
                {
                    if (carry_node->child_count() < rhs_degree) {
                        trees.insert(this_iterator, *carry_node);
                        trees.insert(this_iterator, *rhs_node);
                        carry_node = NULL;
                    }
                    else
                        carry_node = merge_trees(rhs_node, carry_node);
                }
                else
                    trees.insert(this_iterator, *rhs_node);
            }
        }

        if (!rhs.trees.empty()) {
            if (carry_node) {
                node_list_iterator rhs_it = rhs.trees.begin();
                while (static_cast<node_pointer>(&*rhs_it)->child_count() < carry_node->child_count())
                    ++rhs_it;
                rhs.insert_node(rhs_it, carry_node);
                rhs.increment();
                merge_and_clear_nodes(rhs);
            }
            else
                trees.splice(trees.end(), rhs.trees, rhs.trees.begin(), rhs.trees.end());
            return;
        }

        if (carry_node)
            insert_node(this_iterator, carry_node);
    }

    void clone_forest(binomial_heap const & rhs)
    {
        assert(trees.empty());
        typedef typename node::template node_cloner<allocator_type> node_cloner;
        trees.clone_from(rhs.trees, node_cloner(*this, NULL), detail::nop_disposer());

        top_element = detail::find_max_child<node_list_type, node, super_t, value_extractor>(trees, super_t::get_cmp());
    }

    struct force_inf
    {
        template <typename X>
        bool operator()(X const &, X const &) const
        {
            return false;
        }
    };

    template <typename Compare>
    void siftup(node_pointer n, Compare const & cmp)
    {
        while (n->parent)
        {
            node_pointer parent = n->parent;
            int parent_children = parent->child_count();
            node_pointer grand_parent = parent->parent;
            if (cmp(super_t::get_value(n->value), super_t::get_value(parent->value)))
                return;

            n->remove_from_parent();

            n->swap_children(parent);
            n->update_children();
            parent->update_children();

            if (grand_parent) {
                parent->remove_from_parent();
                grand_parent->add_child(n);
            }
            else
            {
                node_list_iterator it = trees.erase(node_list_type::s_iterator_to(*parent));
                trees.insert(it, *n);
            }
            n->add_child(parent);
            int n_children = n->child_count();
            assert(parent_children == n_children);
        }
    }

    void siftdown(node_pointer n)
    {
        while (n->child_count()) {
            node_pointer max_child = detail::find_max_child<node_list_type, node, super_t, value_extractor>(n->children, super_t::get_cmp());

            if (super_t::get_cmp()(super_t::get_value(max_child->value), super_t::get_value(n->value)))
                return;

            max_child->remove_from_parent();

            n->swap_children(max_child);
            n->update_children();
            max_child->update_children();

            node_pointer parent = n->parent;
            if (parent)
            {
                n->remove_from_parent();
                max_child->add_child(n);
                parent->add_child(max_child);
            }
            else
            {
                node_list_iterator position = trees.erase(node_list_type::s_iterator_to(*n));
                max_child->add_child(n);
                trees.insert(position, *max_child);
            }
        }
    }

    void insert_node(node_list_iterator it, node_pointer n)
    {
        if (it != trees.end())
            assert(static_cast<node_pointer>(&*it)->child_count() >= n->child_count());

        for (;;) {
            assert(!n->is_linked());
            if (it == trees.end())
                break;

            node_pointer this_node = static_cast<node_pointer>(&*it);
            size_type this_degree = this_node->child_count();
            size_type n_degree = n->child_count();
            if (this_degree == n_degree)
            {
                assert(it->is_linked());
                it = trees.erase(it);

                n = merge_trees(this_node, n);
            }
            else
                break;
        }
        trees.insert(it, *n);
    }

    explicit binomial_heap(node_list_type & child_list, size_type size):
        super_t(compare_type())
    {
        size_holder::set_size(size);
        top_element = static_cast<node_pointer>(&*child_list.begin());
        for (node_list_iterator it = child_list.begin(); it != child_list.end(); ++it)
        {
            node_pointer n = static_cast<node_pointer>(&*it);
            n->parent = NULL;
        }

        trees.splice(trees.end(), child_list, child_list.begin(), child_list.end());

        trees.sort(detail::cmp_by_degree<node>());
        sanity_check();
    }

    node_pointer merge_trees (node_pointer node1, node_pointer node2)
    {
        assert(node1->child_count() == node2->child_count());

        if (super_t::operator()(node1->value, node2->value))
            std::swap(node1, node2);

        if (node2->parent)
            node2->remove_from_parent();

        node1->add_child(node2);
        return node1;
    }

    void sorted_by_degree(void) const
    {
#ifdef BOOST_HEAP_SANITYCHECKS
        int degree = -1;

        for (node_list_const_iterator it = trees.begin(); it != trees.end(); ++it)
        {
            const_node_pointer n = static_cast<const_node_pointer>(&*it);
            assert(int(n->child_count()) > degree);
            degree = n->child_count();

            assert((detail::is_heap<node, super_t>(n, *this)));

            size_type child_nodes = detail::count_nodes<node>(n);
            assert(child_nodes == size_type(1 << static_cast<const_node_pointer>(&*it)->child_count()));
        }
#endif
    }

    void sanity_check(void) const
    {
#ifdef BOOST_HEAP_SANITYCHECKS
        sorted_by_degree();

        if (constant_time_size)
            assert((size_holder::get_size() == detail::count_list_nodes<node, node_list_type>(trees)));
#endif
    }

    node_pointer top_element;
    node_list_type trees;
#endif // BOOST_DOXYGEN_INVOKED
};

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(binomial_heap<T, Options...> const & lhs, binomial_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(binomial_heap<T, A0, A1, A2, A3> const & lhs, binomial_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(binomial_heap<T, Options...> const & lhs, binomial_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(binomial_heap<T, A0, A1, A2, A3> const & lhs, binomial_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(binomial_heap<T, Options...> const & lhs, binomial_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(binomial_heap<T, A0, A1, A2, A3> const & lhs, binomial_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(binomial_heap<T, Options...> const & lhs, binomial_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(binomial_heap<T, A0, A1, A2, A3> const & lhs, binomial_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(binomial_heap<T, Options...> const & lhs, binomial_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(binomial_heap<T, A0, A1, A2, A3> const & lhs, binomial_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(binomial_heap<T, Options...> const & lhs, binomial_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(binomial_heap<T, A0, A1, A2, A3> const & lhs, binomial_heap<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}

} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_D_ARY_HEAP_HPP */
