#include "encoder.h"
#include "linalg.h"
#include "trellis.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace encoding {

hamming_metric::hamming_metric(linalg::lin_vector const& vect, std::vector<double> const& rels) : _given(vect) , _rels(rels) {}

double hamming_metric::count(linalg::lin_vector const& vect) const {
    assert(vect.size() == _given.size());
    return count(vect, 0, vect.size());
}

double hamming_metric::count(linalg::lin_vector const& vect, size_t begin, size_t end) const {
    
    double res = 0.0;
    for (size_t i = begin; i < end; ++i) {
        if (vect[i - begin] != _given[i]) {
            res += _rels[i];
        }
    }
    return res;
}

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




trellis_based_rml_decoder::trellis_based_rml_decoder(linalg::matrix const& mt, bool use_grey_code, bool use_uniform_optimization) : _gen(mt) , g(use_grey_code) , u(use_uniform_optimization) 
    , _special_matrices(_gen[0].size(), std::vector<linalg::matrix>(_gen[0].size()))
    , _shifts(_gen[0].size(), std::vector<linalg::matrix>(_gen[0].size()))
    , _trellises(_gen[0].size(), std::vector<trellis>(_gen[0].size()))
    , _cbt(_gen[0].size(), std::vector<std::map<linalg::lin_vector, std::pair<linalg::lin_vector, double>>>(_gen[0].size())) 
{
    init(0, _gen[0].size());
}

void trellis_based_rml_decoder::init(size_t x, size_t y) {

    build_special_matrix(x, y);

    size_t z = count_partition(x, y);
    if (y - x > 1) {
        init(x, z);
        init(z, y);
    }
    build_special_trellis(x, z - x, y - x);   
}


size_t trellis_based_rml_decoder::count_partition(size_t x, size_t y) {
    return (y + x) / 2;
}

void trellis_based_rml_decoder::build_special_matrix(size_t x, size_t y) {
    if (y <= x) {
        throw std::domain_error("bad time boundings given");
    }

    if (!_special_matrices[x][y-1].empty()) {
        return;
    }

    linalg::matrix mt = _gen;
    size_t dim = 0;
    size_t c_tr_ctors = 0;
    for (size_t i = 0; i < mt.size(); ++i) {
        if (mt[i].leading() >= x && mt[i].trailing() <= y) {
            using std::swap;
            swap(mt[i], mt[c_tr_ctors]);
            ++c_tr_ctors;
        }
    }
    
    mt = mt.puncture(x, y);
    
    auto res = mt.resolve_basis_gaussian(c_tr_ctors);

    // we assume here, that matrix is "OK", and rows in generating matrix form a k-basis, 
    // so first c_tr_ctors are linear independant

    for (size_t i = 0; i < c_tr_ctors; ++i) {
        res.first[i].resize(y - x + res.first.size() - c_tr_ctors);
        res.second[i].resize(y - x + res.first.size() - c_tr_ctors);
    }


    for (size_t i = c_tr_ctors; i < res.first.size(); ++i) {
        res.first[i].resize(y - x + res.first.size() - c_tr_ctors);
        res.second[i].resize(y - x + res.first.size() - c_tr_ctors);
        res.first[i][y - x + i - c_tr_ctors] = true;
    }

    res.first.make_tof(res.second);

    _special_matrices[x][y-1] = res.first;
    _shifts[x][y-1] = res.second;
}

linalg::lin_vector trellis_based_rml_decoder::get_coset_vect(size_t x, size_t y, linalg::lin_vector const& vect) const {
    // build_special_matrix(x, y);
    
    linalg::matrix const& ctors = _special_matrices[x][y-1];
    linalg::matrix const& shifts = _shifts[x][y-1];

    linalg::lin_vector cur = vect;

    linalg::lin_vector result(vect.size(), false);
    
    for (size_t i = 0; i < ctors.size(); ++i) {
        if (cur.leading() == ctors[i].leading()) {
            for (size_t j = 0; j < y - x; ++j) {
                cur[j] = (cur[j] != ctors[i][j]);
            }
            if (ctors[i].trailing() >= (y - x)) {
                for (size_t j = 0; j < y - x; ++j) {
                    result[j] = (result[j] != shifts[i][j]);
                    result[j] = (result[j] != ctors[i][j]);
                }
            }
            
        }
    }
    
    return result;
}


// gen matrix should be in TOF
void trellis_based_rml_decoder::build_special_trellis(size_t x, size_t z, size_t y) {

    linalg::matrix const& gen = _special_matrices[x][y + x - 1];
    linalg::matrix const& coset_shifts = _shifts[x][y + x - 1];
    
    trellis result;
    result._sections.resize(2);
    // get states for z, we assume that z corresponds to multiple start states
    // get states for y, we assume that y corresponds to multiple final states

    size_t z_dim = 0;
    size_t y_dim = 0;
    for (auto const& vect : gen) {
        size_t leading = vect.leading();
        size_t trailing = vect.trailing();
        if (z <= trailing && z > leading) { // check if z in active span
            ++z_dim;
        }
        if (y <= trailing && y > leading) { // check if y in active span
            ++y_dim;
        }
    }

    for (size_t i = 0; i < (1 << z_dim); ++i) {
        result._sections[0].emplace_back();    
    }

    for (size_t i = 0; i < (1 << y_dim); ++i) {
        result._sections[1].emplace_back();
    }

    // determine transitions
        
    auto g_s_ind = gen.get_g_s(z);

    auto g_s = gen.retrieve(g_s_ind);
    // auto g_s_shifts = coset_shifts.retrieve(g_s_ind);

    auto g_f_s_ind = gen.get_g_f_s(z, y);

    

    auto g_f_s_punctured = gen.retrieve(g_f_s_ind).puncture(z, y);
    // auto g_f_s_shifts = coset_shifts.retrieve(g_f_s_ind).puncture(z, y);


    linalg::lin_vector a;
    a.resize(g_s.size());

    size_t branches = gen.retrieve(gen.get_g_s(z)).get_g_s(y).size();


    for (size_t i = 0; i < (1 << g_s.size()); ++i) {
        linalg::lin_vector a_star, a_a_star;
        a_star.resize(y_dim - branches);
        a_a_star.insert(a_a_star.begin(), a.begin() + g_s.size() - branches, a.end());
        for (size_t j = branches; j < y_dim; ++j) {
            a_a_star.push_back(false);
        }

        linalg::lin_vector p_u = a.empty() ? linalg::lin_vector(y - z) : (a * g_s).puncture(z, y);

        for (size_t j = 0; j < (1 << (y_dim - branches)); ++j) {                
            linalg::lin_vector coset = get_coset_vect(x + z, x + y,a_star.empty() || g_f_s_punctured.empty() ? p_u :  p_u + a_star * g_f_s_punctured);
            result._sections[0][a.to_bit_mask()]._next[coset] = a_a_star.to_bit_mask();

            ++a_star;
            ++a_a_star;
        }
        ++a;
    }

    _trellises[x][y+x-1] = result;
}

void trellis_based_rml_decoder::make_cbt_i(size_t x, size_t y, hamming_metric const& metric) {
    if (x >= y) {
        return;
    }
    build_special_matrix(x, y);
    
    linalg::matrix const& mt = _special_matrices[x][y-1];
    linalg::matrix const& shifts = _shifts[x][y-1];

    linalg::matrix coset_basis;
    linalg::matrix cosets;
    linalg::matrix cosets_shifts;
    bool flag = false;
    size_t i = 0;

    for (; i < mt.size(); ++i) {
        if (mt[i].trailing() >= y - x) {
            flag = true;
        }
        if (flag) {
            cosets.emplace_back(mt[i].puncture(0, y - x));
            cosets_shifts.emplace_back(shifts[i].puncture(0, y - x));
        } else {
            coset_basis.emplace_back(mt[i].puncture(0, y - x));
        }
    }

    for (size_t k = 0; k < (1 << cosets.size()); ++k) {
        linalg::lin_vector a = cosets.empty() ? linalg::lin_vector(y - x) : linalg::lin_vector(k, cosets.size()) * cosets;
        linalg::lin_vector c = cosets.empty() ? linalg::lin_vector(y - x) : linalg::lin_vector(k, cosets.size()) * cosets_shifts;
        linalg::lin_vector coset_vect = a + c;

        for (size_t j = 0; j < (1 << coset_basis.size()); ++j) {
            linalg::lin_vector vect = !coset_basis.empty() ? linalg::lin_vector(j, coset_basis.size()) * coset_basis + a : a;

            double score = metric.count(vect, x, y);
            if (_cbt[x][y - 1].find(coset_vect) == _cbt[x][y - 1].end() || 
                    _cbt[x][y - 1][coset_vect].second > score) {
                _cbt[x][y - 1][coset_vect] = std::make_pair(vect, score);
            }
        }
    }

}

void trellis_based_rml_decoder::make_cbt_g(size_t x, size_t y, hamming_metric const& metric)  {
    throw std::logic_error("not implemented yet");
}

void trellis_based_rml_decoder::comb_cbt_v(size_t x, size_t y, hamming_metric const& metric)  {
    if (y - x <= 4) {
        if (g) {
            make_cbt_g(x, y, metric);
        } else {
            make_cbt_i(x, y, metric);
        }
        return;
    }
    size_t z = count_partition(x, y);
    comb_cbt_v(x, z, metric);
    comb_cbt_v(z, y, metric);


    linalg::matrix const& trellis_gen = _special_matrices[x][y-1];
    linalg::matrix const& shifts = _shifts[x][y-1];

    auto const& tr = _trellises[x][y-1];


    auto g_f_s_ind = trellis_gen.get_g_f_s(0, z - x);
    linalg::matrix g_f_s = trellis_gen.retrieve(g_f_s_ind).puncture(0, z - x); // it is just g_s
    // linalg::matrix g_f_s_shifts = shifts.retrieve(g_f_s_ind).puncture(0, z - x);


    auto g_f_s_y_ind = trellis_gen.get_g_f_s(0, y - x);
    linalg::matrix g_f_s_y = trellis_gen.retrieve(g_f_s_y_ind).puncture(0, y - x); // it is just g_s
    // linalg::matrix g_f_s_y_shifts = shifts.retrieve(g_f_s_y_ind).puncture(0, y - x);

    for (size_t state_mask_z = 0; state_mask_z < tr._sections[0].size(); ++state_mask_z) {

        // calculate coset shift p(u) for given set is zero, 
        // as x is the only starting state in special trellis
        linalg::lin_vector a_star(state_mask_z, g_f_s.size());
        linalg::lin_vector b = a_star.empty() ? linalg::lin_vector(z - x) : a_star * g_f_s; // + a_star * g_f_s_shifts;

        if (_cbt[x][z - 1].find(b) == _cbt[x][z - 1].end()) {
            b = get_coset_vect(x, z, b);
        }

        for(auto const& branch : tr._sections[0][state_mask_z]._next) {
            linalg::lin_vector bf = branch.first;
            linalg::lin_vector b_and_adj = get_coset_vect(x, y, b.concat(bf)); 

            // make matrix:
            // vectors of p(C_x_z))
            // then vectors of C^tr_x_z
            // b
            // apply resolve gaussian basis
            // the last vector is a true b
            // can be used always, yet this method is less effective and it is better not to apply it if we already have a coset vector
            

            if (_cbt[z][y - 1].find(bf) == _cbt[z][y - 1].end()) {
                bf = get_coset_vect(z, y, bf);
            }

            linalg::lin_vector best = _cbt[x][z - 1][b].first.concat(_cbt[z][y - 1][bf].first);
            
            double metric = _cbt[x][z - 1][b].second + _cbt[z][y - 1][bf].second;
            if (_cbt[x][y - 1].find(b_and_adj) == _cbt[x][y - 1].end() || metric < _cbt[x][y - 1][b_and_adj].second) {
                _cbt[x][y - 1][b_and_adj] = std::make_pair(best, metric);
            }
        }
    }
}

void trellis_based_rml_decoder::comb_cbt_u(size_t x, size_t y, hamming_metric const& metric)  {
    throw std::logic_error("not implemented yet");
}

linalg::lin_vector trellis_based_rml_decoder::decode(linalg::lin_vector const& decisions, std::vector<double> const& rels)  {
    hamming_metric metric(decisions, rels);

    comb_cbt_v(0, decisions.size(), metric);
    
    auto res = (*_cbt[0][decisions.size() - 1].begin()).second.first;
    for (auto& xs : _cbt) { // reset _cbt
        for (auto& ys : xs) {
            ys.clear();
        }
    }

    linalg::lin_vector decoded;
    for (auto const& vect : _gen) {
        if (vect.leading() == res.leading()) {
            decoded.push_back(true);
            res += vect;
        } else {
            decoded.push_back(false);
        }
    }

    return decoded;
}

linalg::lin_vector trellis_based_rml_decoder::decode(std::vector<double> const& signals)  {
    std::vector<double> reliabilities;
    linalg::lin_vector hard_decisions;
    for (auto const& i : signals) {
        if (i >= 0) {
            hard_decisions.push_back(1);
        } else {
            hard_decisions.push_back(0);
        }
        reliabilities.push_back(std::fabs(i));
    }
    return decode(hard_decisions, reliabilities);
}





}