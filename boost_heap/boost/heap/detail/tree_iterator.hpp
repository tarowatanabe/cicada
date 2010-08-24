// boost heap: node tree iterator helper classes
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_DETAIL_TREE_ITERATOR_HPP
#define BOOST_HEAP_DETAIL_TREE_ITERATOR_HPP

#include <cassert>

#include <functional>
#include <vector>

#include <boost/iterator/iterator_adaptor.hpp>

namespace boost {
namespace heap {
namespace detail {


template<typename type>
struct identity:
    public std::unary_function<type,type>
{
    type& operator()(type& __x) const
    { return __x; }

    const type& operator()(const type& __x) const
    { return __x; }
};

template<typename type>
struct caster:
    public std::unary_function<type,type>
{
    template <typename U>
    type& operator()(U& __x) const
    { return static_cast<type&>(__x); }

    template <typename U>
    const type& operator()(const U& __x) const
    { return static_cast<const type&>(__x); }
};

template<typename Node>
struct dereferencer
{
    template <typename Iterator>
    Node * operator()(Iterator const & it)
    {
        return static_cast<Node *>(*it);
    }
};

template<typename Node>
struct pointer_to_reference
{
    template <typename Iterator>
    Node * operator()(Iterator const & it)
    {
        return static_cast<Node *>(&*it);
    }
};


/* tree iterator helper class
 *
 * Requirements:
 * Node provides child_iterator
 * ValueExtractor can convert Node->value to ValueType
 *
 * */
template <typename Node,
          typename ValueType,
          typename Alloc = std::allocator<Node>,
          typename ValueExtractor = identity<typename Node::value_type>,
          typename PointerExtractor = dereferencer<Node>,
          bool check_null_pointer = false
         >
class tree_iterator:
    public boost::iterator_adaptor<tree_iterator<Node,
                                                 ValueType,
                                                 Alloc,
                                                 ValueExtractor,
                                                 PointerExtractor,
                                                 check_null_pointer
                                                >,
                                   Node *,
                                   ValueType,
                                   boost::forward_traversal_tag
                                  >,
    ValueExtractor,
    PointerExtractor
{
    typedef boost::iterator_adaptor<tree_iterator<Node,
                                                  ValueType,
                                                  Alloc,
                                                  ValueExtractor,
                                                  PointerExtractor,
                                                  check_null_pointer
                                                 >,
                                    Node *,
                                    ValueType,
                                    boost::forward_traversal_tag
                                   > adaptor_type;

    friend class boost::iterator_core_access;

public:
    tree_iterator(void):
        adaptor_type(0)
    {}

    explicit tree_iterator(Node * it):
        adaptor_type(it)
    {
        discover_nodes(it);
    }

    bool operator!=(tree_iterator const & rhs) const
    {
        return adaptor_type::base() != rhs.base();
    }

    bool operator==(tree_iterator const & rhs) const
    {
        return !operator!=(rhs);
    }

private:
    void increment(void)
    {
        if (unvisited.empty())
            adaptor_type::base_reference() = 0;
        else
        {
            Node * next = unvisited.back();
            unvisited.pop_back();
            discover_nodes(next);
            adaptor_type::base_reference() = next;
        }
    }

    ValueType const & dereference() const
    {
        return ValueExtractor::operator()(adaptor_type::base_reference()->value);
    }

    void discover_nodes(Node * n)
    {
        for (typename Node::child_iterator it = n->children.begin(); it != n->children.end(); ++it)
        {
            Node * n = PointerExtractor::operator()(it);
            if (check_null_pointer && n == NULL)
                continue;
            unvisited.push_back(n);
        }
    }

    std::vector<Node *, Alloc> unvisited;
};

template <typename Node, typename NodeList>
struct list_iterator_converter
{
    struct NodeList::const_iterator operator()(const Node * node)
    {
        return NodeList::s_iterator_to(*node);
    }
};

template <typename Node,
          typename NodeIterator,
          typename ValueType,
          typename ValueExtractor = identity<typename Node::value_type>,
          typename IteratorCoverter = identity<NodeIterator>
         >
class recursive_tree_iterator:
    public boost::iterator_adaptor<recursive_tree_iterator<Node,
                                                           NodeIterator,
                                                           ValueType,
                                                           ValueExtractor,
                                                           IteratorCoverter
                                                          >,
                                    NodeIterator,
                                    ValueType const,
                                    boost::bidirectional_traversal_tag>,
    ValueExtractor, IteratorCoverter
{
    typedef boost::iterator_adaptor<recursive_tree_iterator<Node,
                                                            NodeIterator,
                                                            ValueType,
                                                            ValueExtractor,
                                                            IteratorCoverter
                                                           >,
                                    NodeIterator,
                                    ValueType const,
                                    boost::bidirectional_traversal_tag> adaptor_type;

    friend class boost::iterator_core_access;


public:
    recursive_tree_iterator(void):
        adaptor_type(0)
    {}

    explicit recursive_tree_iterator(NodeIterator const & it):
        adaptor_type(it)
    {}

    void increment(void)
    {
        NodeIterator next = adaptor_type::base_reference();

        const Node * n = get_node(next);
        if (n->children.empty())
        {
            const Node * parent = get_node(next)->get_parent();

            ++next;

            for (;;) {
                if (parent == NULL || next != parent->children.end())
                    break;

                next = IteratorCoverter::operator()(parent);
                parent = get_node(next)->get_parent();
                ++next;
            }
        }
        else
            next = n->children.begin();

        adaptor_type::base_reference() = next;
        return;
    }

    ValueType const & dereference() const
    {
        return ValueExtractor::operator()(get_node(adaptor_type::base_reference())->value);
    }

    static const Node * get_node(NodeIterator const & it)
    {
        return static_cast<const Node *>(&*it);
    }
};


} /* namespace detail */
} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_DETAIL_TREE_ITERATOR_HPP */
