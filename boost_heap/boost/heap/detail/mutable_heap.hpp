// boost heap
//
// Copyright (C) 2010 Tim Blechmann
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// Disclaimer: Not a Boost library.

#ifndef BOOST_HEAP_DETAIL_MUTABLE_HEAP_HPP
#define BOOST_HEAP_DETAIL_MUTABLE_HEAP_HPP

#include <list>
#include <utility>

#include <boost/noncopyable.hpp>
#include <boost/iterator/iterator_adaptor.hpp>


namespace boost
{

namespace heap
{

namespace detail
{

/* wrapper for a mutable heap container adaptors
 *
 * this wrapper introduces an additional indirection. the heap is not constructed from objects,
 * but instead from std::list iterators. this way, the mutability is achieved
 *
 */
template <typename PriorityQueueType>
class priority_queue_mutable_wrapper
{
public:
    typedef typename PriorityQueueType::value_type value_type;
    typedef typename PriorityQueueType::size_type size_type;
    typedef typename PriorityQueueType::compare_type compare_type;
    typedef typename PriorityQueueType::allocator_type allocator_type;

    typedef typename PriorityQueueType::const_reference const_reference;
    typedef typename PriorityQueueType::const_pointer const_pointer;

private:
    typedef std::pair<value_type, size_type> node_type;

    typedef std::list<node_type, allocator_type> object_list;

    typedef typename object_list::iterator list_iterator;
    typedef typename object_list::const_iterator const_list_iterator;


    template <typename value_type>
    struct index_updater
    {
        template <typename It>
        void operator()(It & it, size_type new_index)
        {
            q_type::get_value(it)->second = new_index;
        }

        template <typename U>
        struct rebind {
            typedef index_updater<U> other;
        };
    };

public:
    struct handle_type
    {
        value_type & operator*() const
        {
            return iterator->first;
        }

        handle_type (void)
        {}

    private:
        explicit handle_type(list_iterator const & it):
            iterator(it)
        {}

        list_iterator iterator;

        friend class priority_queue_mutable_wrapper;
    };

private:
    struct indirect_cmp:
        public compare_type
    {
        indirect_cmp(compare_type const & cmp = compare_type()):
            compare_type(cmp)
        {}

        bool operator()(list_iterator const & lhs, list_iterator const & rhs) const
        {
            return compare_type::operator()(lhs->first, rhs->first);
        }
    };

    typedef typename PriorityQueueType::template rebind<list_iterator,
                                                        indirect_cmp,
                                                        allocator_type, index_updater<list_iterator> >::other q_type;

protected:
    q_type q_;
    object_list objects;

protected:
    priority_queue_mutable_wrapper(compare_type const & cmp = compare_type()):
        q_(cmp)
    {}

    priority_queue_mutable_wrapper(priority_queue_mutable_wrapper const & rhs):
        objects(rhs.objects)
    {
        for (typename object_list::iterator it = objects.begin(); it != objects.end(); ++it)
            q_.push(it);
    }

    priority_queue_mutable_wrapper & operator=(priority_queue_mutable_wrapper const & rhs)
    {
        objects = rhs.objects;
        q_.clear();
        for (typename object_list::iterator it = objects.begin(); it != objects.end(); ++it)
            q_.push(it);
        return *this;
    }

#ifdef BOOST_HAS_RVALUE_REFS
    priority_queue_mutable_wrapper (priority_queue_mutable_wrapper && rhs):
        q_(std::move(rhs.q_)), objects(std::move(rhs.objects))
    {}

    priority_queue_mutable_wrapper & operator=(priority_queue_mutable_wrapper && rhs)
    {
        q_ = std::move(rhs.q_);
        objects = std::move(rhs.objects);
        return *this;
    }
#endif


public:
    class iterator:
        public boost::iterator_adaptor<iterator,
                                       const_list_iterator,
                                       value_type const,
                                       boost::bidirectional_traversal_tag>
    {
        typedef boost::iterator_adaptor<iterator,
                                       const_list_iterator,
                                       value_type const,
                                       boost::bidirectional_traversal_tag> super_t;

        friend class boost::iterator_core_access;
        friend class priority_queue_mutable_wrapper;

        iterator(void):
            super_t(0)
        {}

        explicit iterator(const_list_iterator const & it):
            super_t(it)
        {}

        value_type const & dereference() const
        {
            return super_t::base()->first;
        }
    };

    typedef iterator const_iterator;
    typedef typename object_list::difference_type difference_type;

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
        return objects.max_size();
    }

    /// \copydoc boost::heap::priority_queue::clear
    void clear(void)
    {
        q_.clear();
        objects.clear();
    }

    /// \copydoc boost::heap::priority_queue::get_allocator
    allocator_type get_allocator(void) const
    {
        return q_.get_allocator();
    }

    /// \copydoc boost::heap::priority_queue::swap
    void swap(priority_queue_mutable_wrapper & rhs)
    {
        objects.swap(rhs.objects);
        q_.swap(rhs.q_);
    }

    /// \copydoc boost::heap::priority_queue::top
    const_reference top(void) const
    {
        BOOST_ASSERT(!empty());
        return q_.top()->first;
    }

    /// \copydoc boost::heap::priority_queue::push
    handle_type push(const_reference v)
    {
        objects.push_front(std::make_pair(v, 0));
        list_iterator ret = objects.begin();
        q_.push(ret);
        return handle_type(ret);
    }

    /// \copydoc boost::heap::priority_queue::pop
    void pop(void)
    {
        BOOST_ASSERT(!empty());
        list_iterator q_top = q_.top();
        q_.pop();
        objects.erase(q_top);
    }

    /**
     * \b Effects: Merge with priority queue rhs.
     *
     * \b Complexity: N log(N)
     *
     * */
    void merge(priority_queue_mutable_wrapper const & rhs)
    {
        q_.reserve(q_.size() + rhs.q_.size());

        for (typename object_list::const_iterator it = rhs.objects.begin(); it != rhs.objects.end(); ++it)
            push(it->first);
    }

    /**
     * \b Effects: Merges all elements from rhs to this. Rhs is cleared.
     *
     * \b Complexity: N log(N)
     *
     * */
    void merge_and_clear(priority_queue_mutable_wrapper & rhs)
    {
        merge(rhs);
        rhs.clear();
    }

    /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic.
     *
     * */
    void update(handle_type handle, const_reference v)
    {
        list_iterator it = handle.iterator;
        value_type const & current_value = it->first;
        compare_type const & cmp = q_;
        if (cmp(v, current_value))
            decrease(handle, v);
        else
            increase(handle, v);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void update(handle_type handle)
    {
        list_iterator it = handle.iterator;
        size_type index = it->second;
        q_.update(index);
    }

     /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: The new value is expected to be greater than the current one
     * */
    void increase(handle_type handle, const_reference v)
    {
        BOOST_ASSERT(!compare_type()(v, handle.iterator->first));
        handle.iterator->first = v;
        increase(handle);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: The new value is expected to be greater than the current one. If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void increase(handle_type handle)
    {
        list_iterator it = handle.iterator;
        size_type index = it->second;
        q_.increase(index);
    }

     /**
     * \b Effects: Assigns \c v to the element handled by \c handle & updates the priority queue.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: The new value is expected to be less than the current one
     * */
    void decrease(handle_type handle, const_reference v)
    {
        BOOST_ASSERT(!compare_type()(handle.iterator->first, v));
        handle.iterator->first = v;
        decrease(handle);
    }

    /**
     * \b Effects: Updates the heap after the element handled by \c handle has been changed.
     *
     * \b Complexity: Logarithmic.
     *
     * \b Note: The new value is expected to be less than the current one. If this is not called, after a handle has been updated, the behavior of the data structure is undefined!
     * */
    void decrease(handle_type handle)
    {
        list_iterator it = handle.iterator;
        size_type index = it->second;
        q_.decrease(index);
    }

    /**
     * \b Effects: Removes the element handled by \c handle from the priority_queue.
     *
     * \b Complexity: Logarithmic.
     * */
    void erase(handle_type handle)
    {
        list_iterator it = handle.iterator;
        size_type index = it->second;
        q_.erase(index);
        objects.erase(it);
    }

    /// \copydoc boost::heap::priority_queue::begin
    iterator begin(void) const
    {
        return iterator(objects.begin());
    }

    /// \copydoc boost::heap::priority_queue::end
    iterator end(void) const
    {
        return iterator(objects.end());
    }

    /// \copydoc boost::heap::d_ary_heap_mutable::s_handle_from_iterator
    static handle_type s_handle_from_iterator(iterator const & it)
    {
        return handle_type(it);
    }
};


} /* namespace detail */
} /* namespace heap */
} /* namespace boost */

#endif /* BOOST_HEAP_DETAIL_MUTABLE_HEAP_HPP */
