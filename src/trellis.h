#pragma once

#include "linalg.h"
#include <cstddef>
#include <limits>
#include <map>
#include <unordered_map>
#include <vector>

namespace encoding {
struct trellis_node {
    linalg::lin_vector _incoming_coset;

    std::map<linalg::lin_vector, size_t> _next;

    trellis_node() : _next() {}
};

struct trellis {
    std::vector<std::vector<trellis_node> > _sections;

    // it consists of nodes separated by sections

    // additional info
    std::vector<std::vector<std::vector<std::set<size_t>>>> _groups; // contains vectors of ids of nodes in the same group for each section for each parallel component
    std::vector<std::vector<std::vector<size_t>>> _parallel_components; // contains vectors of ids of nodes in the same parallel component for each section
    std::vector<std::set<size_t>> _parallel_component_groups; // contains vectors of ids of parallel components in the same groups

    std::vector<std::vector<std::vector<std::vector<std::pair<double, linalg::lin_vector>>>>> _branches; // for each compgroup for each group pair store sorted ctors from z to y
    std::vector<std::vector<std::vector<std::pair<double, size_t>>>> _inner_branches; // for each compgroup for each group store sorted ctors from x to z
    std::vector<std::vector<std::vector<std::unordered_map<linalg::lin_vector, std::pair<linalg::lin_vector, double>>>>> _group_cache; // for each comp for each group pair store a map of results

    std::vector<std::pair<double, std::pair<size_t, size_t>>> _heap_storage;
    
    bool check_vert_in_group(size_t, size_t, size_t) const noexcept;
    bool check_comp_in_group(size_t) const noexcept;
    bool compare_verteces(size_t, size_t, size_t) const noexcept;
    bool compare_parallel_components(size_t fst, size_t snd) const noexcept;
    void partition_parallel_components();

    void partition_group(size_t);
    
    // for each parallel component group
    void partition_components();

    void init_branches_arrays();
};

}