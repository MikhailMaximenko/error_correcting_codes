#include "trellis.h"
#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <ostream>
#include <set>
#include <vector>

namespace encoding {

bool trellis::check_vert_in_group(size_t comp, size_t sect, size_t vert) const noexcept {
    for (auto const& s : _groups[comp][sect]) {
        if (s.find(vert) != s.end()) {
            return true;
        }
    }
    return false;
}

bool trellis::check_comp_in_group(size_t comp) const noexcept {
    for (auto const& s : _parallel_component_groups) {
        if (s.find(comp) != s.end()) {
            return true;
        }
    }
    return false;
}
bool trellis::compare_verteces(size_t fst, size_t snd) const noexcept {
    if (_sections[0][fst]._next.size() != _sections[0][snd]._next.size()) {
        return false;
    }
    auto it1 = _sections[0][fst]._next.begin();
    auto it2 = _sections[0][snd]._next.begin();
    while (it1 != _sections[0][fst]._next.end() && it2 != _sections[0][snd]._next.end()) {
        if (it1->first != it2->first) {
            return false;
        }
        ++it1;
        ++it2;
    }
    return true;
}

bool trellis::compare_parallel_components(size_t fst, size_t snd) const noexcept {
    for (size_t i : _parallel_components[fst][0]) {
        bool f = false;
        for (size_t j : _parallel_components[snd][0]) {
            if (compare_verteces(i, j)) {
                f = true;
                break;
            }
        }
        if (!f) {
            return false;
        }
    }
    return true;
}

void trellis::partition_parallel_components() {
    for (size_t i = 0; i < _parallel_components.size(); ++i) {
        if (check_comp_in_group(i)) {
            continue;
        } else {
            _parallel_component_groups.emplace_back();
            _parallel_component_groups.back().insert(i);
        }

        for (size_t j = i + 1; j < _parallel_components.size(); ++j) {
            if (compare_parallel_components( i, j)) {
                _parallel_component_groups.back().insert(j);
            }
        }      
    }
}


void trellis::partition_group(size_t comp) {
    for (size_t i = 0; i < _parallel_components[comp][0].size(); ++i) {
        if (check_vert_in_group(comp, 0, _parallel_components[comp][0][i])) {
            continue;
        } else {
            _groups[comp][0].emplace_back();
            _groups[comp][0].back().insert(_parallel_components[comp][0][i]);
        }
        for (size_t j = i + 1; j < _parallel_components[comp][0].size(); ++j) {
            if (compare_verteces( _parallel_components[comp][0][i], _parallel_components[comp][0][j])) {
                _groups[comp][0].back().insert(_parallel_components[comp][0][j]);
            }
        }
    }

    for (auto const& br : _sections[0][*_groups[comp][0][0].begin()]._next) {
        if (check_vert_in_group(comp, 1, br.second)) {
            continue;
        } else {
            _groups[comp][1].emplace_back();
        }
        for (size_t nd : _groups[comp][0][0]) {
            assert(!check_vert_in_group(comp, 1, _sections[0][nd]._next[br.first]));
            _groups[comp][1].back().insert(_sections[0][nd]._next[br.first]);
        }
    }
}

void trellis::partition_components() {
    for (size_t i = 0; i < _parallel_components.size(); ++i) {
        _groups.emplace_back();
        _groups.back().resize(2);
        partition_group(i);
    }
}

void trellis::init_branches_arrays() {
    for (size_t i = 0; i < _parallel_component_groups.size(); ++i) {
        _branches.emplace_back();
        size_t comp = *_parallel_component_groups[i].begin();
        for (size_t k = 0; k < _groups[comp][0].size(); ++k) {
            _branches[i].emplace_back();
            for (size_t j = 0; j < _groups[comp][1].size(); ++j) {
                _branches[i][k].emplace_back();
                for (auto const& br : _sections[0][*_groups[comp][0][k].begin()]._next) {
                    if (_groups[comp][1][j].find(br.second) != _groups[comp][1][j].end()) {
                        _branches[i][k][j].push_back({std::numeric_limits<double>::infinity(), br.first});
                    }
                }
            }
        }
    }

    for (size_t i = 0; i < _parallel_components.size(); ++i) {
        _inner_branches.emplace_back();
        for (size_t j = 0; j < _groups[i][0].size(); ++j) {
            _inner_branches[i].emplace_back();

            for (size_t vert : _groups[i][0][j]) {
                _inner_branches[i][j].push_back({std::numeric_limits<double>::infinity(), vert});
            }
        }
    }

    _group_cache.resize(_parallel_components.size());
    for (auto & comp_gr : _parallel_component_groups) {
        for (size_t comp : comp_gr) {
            _group_cache[comp].resize(_groups[comp][0].size());
            for (auto & it : _group_cache[comp]) {
                it.resize(_groups[comp][1].size());
            }
        }
    }

}

}