#include "encoder.h"
#include "linalg.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace encoding {

linalg::lin_vector encoder::encode(linalg::lin_vector const& o) const {
    return o * generating;
}

double decoder::count_metric(linalg::lin_vector const& vec, std::vector<double> const& weights) {
    double res = 0.0;
    if (vec.size() != weights.size()) {
        throw std::logic_error("code and weight are expected to have same dimentions");
    }
    for (size_t i = 0; i < vec.size(); ++i) {
        if (vec[i]) {
            res += weights[i];
        }
    }
    return res;
}

linalg::lin_vector decoder::decode(std::vector<double> const& signals) {
    std::vector<double> reliabilities;
    linalg::lin_vector hard_decisions;
    for (auto const& i : signals) {
        if (i >= 0) {
            hard_decisions.push_back(1);
        } else {
            hard_decisions.push_back(0);
        }
        reliabilities.push_back(fabs(i));
    }
    return decode(hard_decisions, reliabilities);
}

linalg::lin_vector decoder::decode(linalg::lin_vector const& hard_decisions, std::vector<double> const& reliability) {
    // save matrix to set it back in the end
    linalg::matrix tmp_gen(generating);
    
    //init state
    linalg::lin_vector res(hard_decisions.size());
    double best_score = std::numeric_limits<double>::infinity();
    
    // sort reliabilities to get permutation
    std::vector<std::pair<double, size_t>> find_permutation;
    find_permutation.reserve(reliability.size());
    for (size_t i = 0; i < reliability.size(); ++i) {
        find_permutation.emplace_back(std::make_pair(reliability[i], i));
    }
    std::sort(find_permutation.rbegin(), find_permutation.rend(), 
        [](std::pair<double, size_t> const& a, std::pair<double, size_t> const& b) { return a.first < b.first; });


    // get permutation and sorted reliabilities
    std::vector<size_t> permutation;
    std::vector<double> sorted_rels(reliability.size());
    permutation.reserve(reliability.size());
    for (size_t i = 0; i < reliability.size(); ++i) {
        sorted_rels[i] = find_permutation[i].first;
        permutation.push_back(find_permutation[i].second);
    }


    // get inversed permutation
    std::vector<size_t> back_permutation(permutation.size());
    for (size_t i = 0; i < reliability.size(); ++i) {
        back_permutation[permutation[i]] = i;
    }

    // find interesting basis and transformatin matrix
    generating.permutate(permutation);
    auto resolved = generating.resolve_basis();
    auto base_vectors = std::move(resolved.first);
    auto transformation_matrix = std::move(resolved.second);


    linalg::lin_vector permutated_decisions = hard_decisions;
    permutated_decisions.permutate(permutation);

    // get decisions from interesting positions    
    linalg::lin_vector choosen_decisions(base_vectors.size());
    for (size_t i = 0; i < base_vectors.size(); ++i) {
        choosen_decisions[i] = permutated_decisions[base_vectors[i]];
    }

    linalg::lin_vector decisions;
    for (int w = 0; w <= w_max; ++w) {
        // inverse hard decisions 

        std::vector<bool> positions(base_vectors.size(), false);
        // set inverting positions
        for (size_t i = positions.size() - 1; i > positions.size() - 1 - w; --i) {
            positions[i] = true;
        } 
        do {
            decisions = choosen_decisions;
            for (size_t i = 0; i < positions.size(); ++i) {
                if (positions[i]) {
                    // invert
                    decisions[i] = !decisions[i];
                }
            }
    
            // c := alpha * G~ 
            auto c = decisions * generating;
            // count score in permutated basis
            auto score = count_metric(c + permutated_decisions, sorted_rels);
            if (score < best_score) {

                res = c;
                best_score = score;
            }

        } while(std::next_permutation(positions.begin(), positions.end()));
    }


    res.permutate(back_permutation);
    
    linalg::lin_vector res_decision(base_vectors.size());

    // first k bits are message as generating is systemized
    for (size_t i = 0; i < base_vectors.size(); ++i) {
        res_decision[i] = res[i];
    }

    generating = std::move(tmp_gen);
    return res_decision;
};

}