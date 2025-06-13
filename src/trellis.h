#pragma once

#include "linalg.h"
#include <cstddef>
#include <limits>
#include <map>
#include <vector>

namespace encoding {
struct trellis_node {

    std::map<linalg::lin_vector, size_t> _next;

    trellis_node() : _next() {}
};

struct trellis {
    std::vector<std::vector<trellis_node> > _sections; 
    // it consists of nodes separated by sections
};
}