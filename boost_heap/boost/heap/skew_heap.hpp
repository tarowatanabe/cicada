// boost heap: skew heap
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_SKEW_HEAP_HPP
#define BOOST_HEAP_SKEW_HEAP_HPP

#include <algorithm>
#include <vector>

#include <boost/assert.hpp>
#include <boost/array.hpp>

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

template <typename node_pointer, bool store_parent_pointer>
struct parent_holder
{
    parent_holder(void):
        parent_(NULL)
    {}

    void set_parent(node_pointer parent)
    {
        assert(static_cast<node_pointer>(this) != parent);
        parent_ = parent;
    }

    node_pointer get_parent(void) const
    {
        return parent_;
    }

    node_pointer parent_;
};

template <typename node_pointer>
struct parent_holder<node_pointer, false>
{
    void set_parent(node_pointer parent)
    {}

    node_pointer get_parent(void) const
    {
        return NULL;
    }
};


template <typename value_type, bool store_parent_pointer>
struct skew_heap_node:
    parent_holder<skew_heap_node<value_type, store_parent_pointer>*, store_parent_pointer>
{
    typedef parent_holder<skew_heap_node<value_type, store_parent_pointer>*, store_parent_pointer> super_t;

    typedef boost::array<skew_heap_node*, 2> child_list_type;
    typedef typename child_list_type::iterator child_iterator;
    typedef typename child_list_type::const_iterator const_child_iterator;

    skew_heap_node(value_type const & v):
        value(v)
    {
        children.assign(0);
    }

    template <typename Alloc>
    skew_heap_node (skew_heap_node const & rhs, Alloc & allocator, skew_heap_node * parent):
        value(rhs.value)
    {
        super_t::set_parent(parent);
        node_cloner<skew_heap_node, skew_heap_node, Alloc> cloner(allocator);
        clone_child(0, rhs, cloner);
        clone_child(1, rhs, cloner);
    }

    template <typename Cloner>
    void clone_child(int index, skew_heap_node const & rhs, Cloner & cloner)
    {
        if (rhs.children[index])
            children[index] = cloner(*rhs.children[index], this);
        else
            children[index] = NULL;
    }

    template <typename Alloc>
    void clear_subtree(Alloc & alloc)
    {
        node_disposer<skew_heap_node, skew_heap_node, Alloc> disposer(alloc);
        dispose_child(children[0], disposer);
        dispose_child(children[1], disposer);
    }

    template <typename Disposer>
    void dispose_child(skew_heap_node * node, Disposer & disposer)
    {
        if (node)
            disposer(node);
    }

    std::size_t count_children(void) const
    {
        size_t ret = 1;
        if (children[0])
            ret += children[0]->count_children();
        if (children[1])
            ret += children[1]->count_children();

        return ret;
    }

    template <typename HeapBase>
    bool is_heap(typename HeapBase::compare_type const & cmp) const
    {
        for (const_child_iterator it = children.begin(); it != children.end(); ++it)
        {
            const skew_heap_node * child = *it;

            if (child == NULL)
                continue;

            if (store_parent_pointer)
                assert(child->get_parent() == this);

            if (cmp(HeapBase::get_value(value), HeapBase::get_value(child->value)) ||
                !child->is_heap<HeapBase>(cmp))
                return false;
        }
        return true;
    }


    value_type value;
    boost::array<skew_heap_node*, 2> children;
};


typedef parameter::parameters<optional<tag::allocator>,
                              optional<tag::compare>,
                              optional<tag::stable>,
                              optional<tag::store_parent_pointer>
                             > skew_heap_signature;

template <typename T, typename Parspec>
struct make_skew_heap_base
{
    static const bool constant_time_size = parameter::binding<Parspec,
                                                              tag::constant_time_size,
                                                              boost::mpl::true_
                                                             >::type::value;

    typedef typename make_heap_base<T, Parspec, constant_time_size>::type base_type;
    typedef typename make_heap_base<T, Parspec, constant_time_size>::allocator_argument allocator_argument;
    typedef typename make_heap_base<T, Parspec, constant_time_size>::compare_argument compare_argument;

    static const bool store_parent_pointer = parameter::binding<Parspec,
                                                              tag::store_parent_pointer,
                                                              boost::mpl::false_>::type::value;

    typedef skew_heap_node<typename base_type::internal_type, store_parent_pointer> node_type;

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



template <typename T,
          class A0 = boost::parameter::void_,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_,
          class A4 = boost::parameter::void_
         >
class skew_heap:
    protected make_skew_heap_base<T,
                                  typename skew_heap_signature::bind<A0, A1, A2, A3, A4>::type
                                 >::type
{
protected:
    typedef skew_heap_signature::bind<A0, A1, A2, A3, A4> bound_args;
    typedef make_skew_heap_base<T, typename bound_args::type> base_maker;
    typedef typename base_maker::type super_t;

    typedef typename super_t::internal_type internal_type;
    typedef typename super_t::size_holder_type size_holder;
    typedef typename base_maker::allocator_argument allocator_argument;

    static const bool store_parent_pointer = base_maker::store_parent_pointer;

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

    static const bool constant_time_size = super_t::constant_time_size;

protected:
    typedef typename base_maker::node_type node;
    typedef typename allocator_type::pointer node_pointer;
    typedef typename allocator_type::const_pointer const_node_pointer;

    typedef detail::value_extractor<value_type, internal_type, super_t> value_extractor;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit skew_heap(compare_type const & cmp = compare_type()):
        super_t(cmp), root(NULL)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    skew_heap(skew_heap const & rhs):
        super_t(rhs), root(0)
    {
        if (rhs.empty())
            return;

        clone_tree(rhs);
    }

public:
    /// \copydoc boost::heap::priority_queue::operator=(priority_queue const & rhs)
    skew_heap & operator=(skew_heap const & rhs)
    {
        clear();
        size_holder::set_size(rhs.get_size());
        static_cast<super_t&>(*this) = rhs;

        clone_tree(rhs);
        return *this;
    }

protected:
#ifdef BOOST_HAS_RVALUE_REFS
    skew_heap(skew_heap && rhs):
        super_t(std::move(rhs)), root(rhs.root)
    {
        rhs.root = NULL;
    }

    skew_heap & operator=(skew_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        root = rhs.root;
        rhs.root = NULL;
        return *this;
    }
#endif

    ~skew_heap(void)
    {
        clear();
    }

public:
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
            return root->count_children();
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
    void swap(skew_heap & rhs)
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
     * \b Effects: Removes the top element from the priority queue.
     *
     * \b Complexity: Logarithmic (amortized).
     *
     * */
    void pop(void)
    {
        BOOST_ASSERT(!empty());

        node_pointer top = root;

        root = merge_children(root);
        size_holder::decrement();

        if (root)
            assert(root->get_parent() == NULL);
        else
            assert(size_holder::get_size() == 0);

        top->~node();
        allocator_type::deallocate(top, 1);
        sanity_check();
    }

    /// \copydoc boost::heap::priority_queue::iterator
    typedef typename boost::mpl::if_c<false,
                                      recursive_tree_iterator<node,
                                                              typename node::child_iterator,
                                                              value_type,
                                                              value_extractor,
                                                              list_iterator_converter<node,
                                                                                      typename node::child_list_type
                                                                                     >
                                                             >,
                                      tree_iterator<node,
                                                    const value_type,
                                                    allocator_type,
                                                    value_extractor,
                                                    dereferencer<node>,
                                                    true
                                                   >
                                      >::type iterator;

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

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: Linear.
     *
     * */
    void merge(skew_heap const & rhs)
    {
        skew_heap copy(rhs);
        merge_and_clear(copy);
    }

    /**
     * \b Effects: Merge all elements from rhs into this
     *
     * \b Complexity: Logarithmic (amortized).
     *
     * */
    void merge_and_clear(skew_heap & rhs)
    {
        if (rhs.empty())
            return;

        merge_node(rhs.root);

        size_holder::add(rhs.get_size());
        rhs.set_size(0);
        rhs.root = NULL;
    }

protected:
#if !defined(BOOST_DOXYGEN_INVOKED)
    node_pointer push_internal(const_reference v)
    {
        size_holder::increment();

        node_pointer n = super_t::allocate(1);
        new(n) node(super_t::make_node(v));

        merge_node(n);
        return n;
    }

    void clone_tree(skew_heap const & rhs)
    {
        assert(root == NULL);
        if (rhs.empty())
            return;

        root = allocator_type::allocate(1);

        new(root) node(*rhs.root, static_cast<allocator_type&>(*this), NULL);
    }

    void merge_node(node_pointer other)
    {
        assert(other);
        if (root != NULL)
            root = merge_nodes(root, other, NULL);
        else
            root = other;

        sanity_check();
    }

    node_pointer merge_nodes(node_pointer node1, node_pointer node2, node_pointer new_parent)
    {
        if (node1 == NULL) {
            if (node2)
                node2->set_parent(new_parent);
            return node2;
        }
        if (node2 == NULL) {
            node1->set_parent(new_parent);
            return node1;
        }

        node_pointer merged = merge_nodes_recursive(node1, node2, new_parent);
        return merged;
    }

    node_pointer merge_children(node_pointer node)
    {
        node_pointer parent = node->get_parent();
        node_pointer merged_children = merge_nodes(node->children[0], node->children[1], parent);

        return merged_children;
    }

    node_pointer merge_nodes_recursive(node_pointer node1, node_pointer node2, node_pointer new_parent)
    {
        if (super_t::operator()(node1->value, node2->value))
            std::swap(node1, node2);

        node * parent = node1;
        node * child = node2;

        if (parent->children[1]) {
            node * merged = merge_nodes(parent->children[1], child, parent);
            parent->children[1] = merged;
            merged->set_parent(parent);
        } else {
            parent->children[1] = child;
            child->set_parent(parent);
        }


        std::swap(parent->children[0], parent->children[1]);
        parent->set_parent(new_parent);
        return parent;
    }

    void sanity_check(void)
    {
#ifdef BOOST_HEAP_SANITYCHECKS
        if (root)
            assert( root->template is_heap<super_t>(super_t::get_cmp()) );

        if (constant_time_size)
        {
            size_type stored_size = size_holder::get_size();

            size_type counted_size;
            if (root == NULL)
                counted_size = 0;
            else
                counted_size = root->count_children();

            assert(counted_size == stored_size);
        }
#endif
    }

    node_pointer root;
#endif
};

} /* namespace detail */

/**
 * \class skew_heap
 * \brief skew heap
 *
 *
 * The template parameter T is the type to be managed by the container.
 * The user can specify additional options and if no options are provided default options are used.
 *
 * The container supports the following options:
 * - \c stable<>, defaults to \c stable<false>
 * - \c compare<>, defaults to \c compare<std::less<T> >
 * - \c allocator<>, defaults to \c allocator<std::allocator<T> >
 * - \c constant_time_size<>, defaults to \c constant_time_size<true>
 * - \c store_parent_pointer<>, defaults to \c store_parent_pointer<true>. Maintaining a parent pointer adds some
 *   maintenance and size overhead, but iterating a heap is more efficient.
 *
 */
#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
#else
template <typename T,
          class A0 = boost::parameter::void_,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_,
          class A4 = boost::parameter::void_
         >
#endif
class skew_heap:
    public detail::skew_heap<T, A0, A1, A2, A3, A4>
{
private:
    typedef detail::skew_heap<T, A0, A1, A2, A3, A4> super_t;

public:
#ifdef BOOST_DOXYGEN_INVOKED
    typedef T value_type;
    typedef detail::unspecified size_type;
    typedef detail::unspecified difference_type;
    typedef detail::unspecified pointer;
    typedef detail::unspecified const_pointer;
    typedef detail::unspecified reference;
    typedef detail::unspecified const_reference;

    typedef detail::unspecified iterator;
    typedef detail::unspecified const_iterator;

    typedef detail::unspecified compare_type;
    typedef detail::unspecified allocator_type;
#else
    typedef typename super_t::compare_type compare_type;
    typedef typename super_t::reference reference;
    typedef typename super_t::const_reference const_reference;
#endif

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit skew_heap(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    skew_heap(skew_heap const & rhs):
        super_t(rhs)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    skew_heap(skew_heap && rhs):
        super_t(std::move(rhs))
    {}

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    skew_heap & operator=(skew_heap && rhs)
    {
        super_t::operator=(std::move(rhs));
        return *this;
    }
#endif

    /**
     * \b Effects: Adds a new element to the priority queue. Returns handle to element
     *
     * \b Complexity: Logarithmic (amortized).
     *
     * */
    void push(const_reference v)
    {
        super_t::push_internal(v);
    }
};







/**
 * \class skew_heap_mutable
 * \brief mutable skew heap
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
class skew_heap_mutable:
    public detail::skew_heap<T, store_parent_pointer<true>, A0, A1, A2, A3>
{
private:
    typedef detail::skew_heap<T, store_parent_pointer<true>, A0, A1, A2, A3> super_t;

public:
#ifdef BOOST_DOXYGEN_INVOKED
    typedef T value_type;
    typedef detail::unspecified size_type;
    typedef detail::unspecified difference_type;
    typedef detail::unspecified pointer;
    typedef detail::unspecified const_pointer;
    typedef detail::unspecified reference;
    typedef detail::unspecified const_reference;

    typedef detail::unspecified iterator;
    typedef detail::unspecified const_iterator;

    typedef detail::unspecified compare_type;
    typedef detail::unspecified allocator_type;

#else

    typedef typename super_t::compare_type compare_type;
    typedef typename super_t::reference reference;
    typedef typename super_t::const_reference const_reference;

    typedef typename super_t::iterator iterator;

private:
    typedef typename super_t::node node;
    typedef typename super_t::node_pointer node_pointer;

public:
#endif

    typedef detail::node_handle<typename super_t::node_pointer,
                                typename super_t::super_t,
                                typename super_t::reference
                               > handle_type;

    /// \copydoc boost::heap::priority_queue::priority_queue(compare_type const &)
    explicit skew_heap_mutable(compare_type const & cmp = compare_type()):
        super_t(cmp)
    {}

    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue const &)
    skew_heap_mutable(skew_heap_mutable const & rhs):
        super_t(rhs)
    {}

#ifdef BOOST_HAS_RVALUE_REFS
    /// \copydoc boost::heap::priority_queue::priority_queue(priority_queue &&)
    skew_heap_mutable(skew_heap_mutable && rhs):
        super_t(std::move(rhs))
    {}

    /// \copydoc boost::heap::priority_queue::operator=(priority_queue &&)
    skew_heap_mutable & operator=(skew_heap_mutable && rhs)
    {
        super_t::operator=(std::move(rhs));
        return *this;
    }
#endif

    /**
     * \b Effects: Adds a new element to the priority queue. Returns handle to element
     *
     * \b Complexity: Logarithmic (amortized).
     *
     * */
    handle_type push(const_reference v)
    {
        return handle_type(super_t::push_internal(v));
    }

    /**
     * \b Effects: Removes the element handled by \c handle from the priority_queue.
     *
     * \b Complexity: Logarithmic (amortized).
     * */
    void erase (handle_type object)
    {
        node_pointer this_node = object.node_;

        unlink_node(this_node);
        super_t::size_holder::decrement();

        super_t::sanity_check();
        this_node->~node();
        super_t::allocator_type::deallocate(this_node, 1);
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic (amortized).
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
     * \b Complexity: Logarithmic (amortized).
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void update (handle_type handle)
    {
        node_pointer this_node = handle.node_;

        if (this_node->get_parent())
        {
            if (super_t::operator()(super_t::get_value(this_node->get_parent()->value),
                                    super_t::get_value(this_node->value)))
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
     * \b Complexity: Logarithmic (amortized).
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
     * \b Complexity: Logarithmic (amortized).
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void increase (handle_type handle)
    {
        node_pointer this_node = handle.node_;

        if (this_node == super_t::root)
            return;

        node_pointer parent = this_node->get_parent();

        if (this_node == parent->children[0])
            parent->children[0] = NULL;
        else
            parent->children[1] = NULL;

        this_node->set_parent(NULL);
        super_t::merge_node(this_node);
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic (amortized).
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
     * \b Complexity: Logarithmic (amortized).
     *
     * \b Note: The new value is expected to be less than the current one. If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void decrease (handle_type handle)
    {
        node_pointer this_node = handle.node_;

        unlink_node(this_node);
        this_node->children.assign(0);
        this_node->set_parent(NULL);
        super_t::merge_node(this_node);
    }

    /// \copydoc boost::heap::d_ary_heap_mutable::s_handle_from_iterator
    static handle_type s_handle_from_iterator(iterator const & it)
    {
        return handle_type(&*it);
    }

private:
    void unlink_node(node_pointer node)
    {
        node_pointer parent = node->get_parent();
        node_pointer merged_children = super_t::merge_children(node);

        if (parent) {
            if (node == parent->children[0])
                parent->children[0] = merged_children;
            else
                parent->children[1] = merged_children;
        }
        else
            super_t::root = merged_children;
    }
};



#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(skew_heap<T, Options...> const & lhs, skew_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3,
          class A4
         >
bool operator==(skew_heap<T, A0, A1, A2, A3, A4> const & lhs, skew_heap<T, A0, A1, A2, A3, A4> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(skew_heap<T, Options...> const & lhs, skew_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3,
          class A4
         >
bool operator!=(skew_heap<T, A0, A1, A2, A3, A4> const & lhs, skew_heap<T, A0, A1, A2, A3, A4> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(skew_heap<T, Options...> const & lhs, skew_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3,
          class A4
         >
bool operator<(skew_heap<T, A0, A1, A2, A3, A4> const & lhs, skew_heap<T, A0, A1, A2, A3, A4> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(skew_heap<T, Options...> const & lhs, skew_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3,
          class A4
         >
bool operator>=(skew_heap<T, A0, A1, A2, A3, A4> const & lhs, skew_heap<T, A0, A1, A2, A3, A4> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(skew_heap<T, Options...> const & lhs, skew_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3,
          class A4
         >
bool operator>(skew_heap<T, A0, A1, A2, A3, A4> const & lhs, skew_heap<T, A0, A1, A2, A3, A4> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(skew_heap<T, Options...> const & lhs, skew_heap<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3,
          class A4
         >
bool operator<=(skew_heap<T, A0, A1, A2, A3, A4> const & lhs, skew_heap<T, A0, A1, A2, A3, A4> const & rhs)
#endif
{
    return !(lhs > rhs);
}


#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator==(skew_heap_mutable<T, Options...> const & lhs, skew_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator==(skew_heap_mutable<T, A0, A1, A2, A3> const & lhs, skew_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_equality(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator!=(skew_heap_mutable<T, Options...> const & lhs, skew_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator!=(skew_heap_mutable<T, A0, A1, A2, A3> const & lhs, skew_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs == rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<(skew_heap_mutable<T, Options...> const & lhs, skew_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<(skew_heap_mutable<T, A0, A1, A2, A3> const & lhs, skew_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return detail::heap_compare(lhs, rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>=(skew_heap_mutable<T, Options...> const & lhs, skew_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>=(skew_heap_mutable<T, A0, A1, A2, A3> const & lhs, skew_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs < rhs);
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator>(skew_heap_mutable<T, Options...> const & lhs, skew_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator>(skew_heap_mutable<T, A0, A1, A2, A3> const & lhs, skew_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return rhs < lhs;
}

#ifdef BOOST_DOXYGEN_INVOKED
template<class T, class ...Options>
bool operator<=(skew_heap_mutable<T, Options...> const & lhs, skew_heap_mutable<T, Options...> const & rhs)
#else
template <typename T,
          class A0,
          class A1,
          class A2,
          class A3
         >
bool operator<=(skew_heap_mutable<T, A0, A1, A2, A3> const & lhs, skew_heap_mutable<T, A0, A1, A2, A3> const & rhs)
#endif
{
    return !(lhs > rhs);
}


} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_SKEW_HEAP_HPP */
