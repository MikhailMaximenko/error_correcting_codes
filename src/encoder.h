
#pragma once

#include "linalg.h"
#include <cstddef>

namespace encoding {

struct encoder {
    linalg::matrix generating;

    linalg::lin_vector encode(linalg::lin_vector const&) const;
};


struct decoder {
    linalg::matrix generating;
    std::size_t w_max;

    static double count_metric(linalg::lin_vector const&, std::vector<double> const&);

    linalg::lin_vector decode(std::vector<double> const&);

    linalg::lin_vector decode(linalg::lin_vector const&, std::vector<double> const&);
};

}