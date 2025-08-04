#include "linalg.h"
#include <cstddef>
#include <iterator>
#include <list>
#include <memory>


namespace utility {

namespace detail {
struct fake_node {
    fake_node* _left;
    fake_node* _right;

    ~fake_node() {
        delete _left;
        delete _right;
    }
};

template<typename R>
struct node : fake_node {
    R _value;
};

}

template<typename R>
struct bool_map {
    std::list<detail::node<R>*> _leaves;
    detail::fake_node _root;


    struct iterator {
        using list_it = typename std::list<detail::node<R>*>::iterator;

        list_it _cur;

        using iterator_category = std::bidirectional_iterator_tag;
        using difference_type = ptrdiff_t;
        using value_type = R;
        using pointer = R*;
        using reference = R&;

        iterator(list_it const& i) : _cur(i) {}

        iterator& operator++() {
            ++_cur;
            return *this;
        }
        iterator operator++(int) {
            auto tmp(*this);
            ++*this;
            return tmp;
        }

        iterator& operator--() {
            --_cur;
            return *this;
        }
        iterator operator--(int) {
            auto tmp(*this);
            --*this;
            return tmp;
        }

        R& operator*() {
            return (*_cur)->_value;
        }

        R const& operator*() const {
            return (*_cur)->_value;
        }

        R* operator->() {
            return &(*_cur)->_value;
        }

        R const* operator->() const {
            return &(*_cur)->_value;
        }

        bool operator==(iterator const& other) const noexcept {
            return _cur == other._cur;
        }

        bool operator!=(iterator const& other) const noexcept {
            return _cur != other._cur;
        }

  
    };

    R& operator[](linalg::lin_vector const& v) {
        detail::fake_node* _cur = &_root;
        for (size_t i = 0; i < v.size() - 1; ++i) {
            if (v[i]) {
                if (!_cur->_right) {
                    _cur->_right = new detail::fake_node();
                } 
                _cur = _cur->_right;
            } else {
                if (!_cur->_left) {
                    _cur->_left = new detail::fake_node();
                } 
                _cur = _cur->_left;
            }
        }
        if (v.back()) {
            if (!_cur->_right) {
                _cur->_right = new detail::node<R>();
                _leaves.push_back(static_cast<detail::node<R>*>(_cur->_right));
            } 
            _cur = _cur->_right;
        } else {
            if (!_cur->_left) {
                _cur->_left = new detail::node<R>();
                _leaves.push_back(static_cast<detail::node<R>*>(_cur->_left));
            } 
            _cur = _cur->_left;
        }
        return static_cast<detail::node<R>*>(_cur)->_value;
    }

    iterator begin() noexcept {
        return iterator(_leaves.begin());
    }


    iterator end() noexcept {
        return iterator(_leaves.end());
    }
};

}