#include "encoder.h"
#include "linalg.h"
#include "trellis.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <queue>
#include <stdexcept>
#include <utility>
#include <vector>

namespace encoding {

linalg::matrix build_rm_code(size_t m, size_t r) {
    linalg::matrix res;
    linalg::matrix basis;
    for (size_t i = 0; i < m; ++i) {
        linalg::lin_vector bv(1 << m, false);
        for (size_t j = 0; j < (1 << m); ++j) {
            if (j & (1 << i)) {
                bv[j] = true;
            }
        }
        basis.push_back(bv);
    }
    for (size_t layer = 0; layer <= r; ++layer) {
        linalg::lin_vector mask(m - layer, false);
        for (size_t i = 0; i < layer; ++i) {
            mask.push_back(true);
        }
        do {
            res.push_back(basis.get_and_multiply(mask));
        } while (std::next_permutation(mask.begin(), mask.end()));
    }
    res.make_tof();
    return res;
}

hamming_metric::hamming_metric(linalg::lin_vector const& vect, std::vector<double> const& rels) : _given(vect) , _rels(rels) {}

double hamming_metric::count(linalg::lin_vector const& vect) const {
    assert(vect.size() == _given.size());
    return count(vect, 0, vect.size());
}

double hamming_metric::count(linalg::lin_vector const& vect, size_t begin, size_t end) const {
    
    double res = vect[0] != _given[begin] ? _rels[begin] : -_rels[begin];
    for (size_t i = begin + 1; i < end; ++i) {
        if (vect[i - begin] != _given[i]) {
            res += _rels[i];
        } else {
            res -= _rels[i]; // now it is symmetric about zero
        }
    }
    return res;
}

double hamming_metric::count_diff(linalg::lin_vector const& vect1, linalg::lin_vector const& vect2, size_t begin, size_t end, size_t & add_cnt) const {
    double r = std::numeric_limits<double>::infinity();
    for (size_t i = begin; i < end; ++i) {
        if (vect1[i - begin] != vect2[i - begin]) {
            if (!std::isinf(r)) { 
                ++add_cnt;
            }
            if (vect1[i - begin] != _given[i]) {
                if (std::isinf(r)) {
                    r = -2 * _rels[i];
                } else {
                    r -= 2 * _rels[i];
                }
            } else {
                if (std::isinf(r)) {
                    r = 2 * _rels[i];
                } else {
                    r += 2 * _rels[i];
                }
            }
        }
    }
    if (std::isinf(r)) {
        r = 0;
    }
    return r;
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
    , _partiotions(_gen[0].size(), std::vector<std::pair<size_t, size_t>>(_gen[0].size()))
    , _gray_codes()
    , _ctors_for_make_cbt(_gen[0].size(), std::vector<std::vector<std::pair<std::pair<linalg::lin_vector, linalg::lin_vector>, linalg::lin_vector>>>(_gen[0].size()))
    , _trellises(_gen[0].size(), std::vector<trellis>(_gen[0].size()))
    , _cbt(_gen[0].size(), std::vector<std::map<linalg::lin_vector, std::pair<linalg::lin_vector, double>>>(_gen[0].size())) 
{
    if (!_gen.is_tof()) {
        throw std::logic_error("matrix should be in tof");
    }
    count_partitions();
    init(0, _gen[0].size());
}

void trellis_based_rml_decoder::init(size_t x, size_t y) {
    build_special_matrix(x, y);
    size_t z = _partiotions[x][y-1].first;

    if (x != z) {
        init(x, z);
        init(z, y);
        build_special_trellis(x, z - x, y - x);
        if (u) {
            make_uniform_decomposition(x, z, y);
        } 
    } else {
        if (g) {
            prepare_make_cbt_g(x, y);
        } else {
            prepare_make_cbt_i(x, y);
        }
    }
}


void trellis_based_rml_decoder::count_partition(size_t x, size_t y) {
    size_t metric = 0;
    build_special_matrix(x, y);

    auto const & mt = _special_matrices[x][y-1];
    size_t tr_ctors = mt.get_c_tr_ctors_number(y - x);

    size_t y_dim = 0;
    for (auto const& vect : mt) {
        size_t leading = vect.leading();
        size_t trailing = vect.trailing();
        if (y - x <= trailing && y - x > leading) { // check if y in active span
            ++y_dim;
        }
    }
    if (mt.size() >= 32 || (g && y - x >= 32)) {
        metric = -1;
    } else {
        if (g) {
            // if coset is sparse real value can be a bit smaller due to the way we count metric (via diff)
            if (get_coset_vect(x, y, linalg::lin_vector(y - x, true)).all_zeros(0, y - x)) {
                metric = (1 << (y - x - 1)) + y - x - 2 + (1 << (mt.size() - tr_ctors)) * ((1 << (tr_ctors - 1)) - 1);
            } else {
                metric = (1 << (y - x - 1)) + y - x - 2 + (1 << (mt.size() - tr_ctors)) * ((1 << tr_ctors) - 1);
            }
        } else {
            metric = (y - x - 1) * (1 << (mt.size())) + (1 << (mt.size() - tr_ctors)) * ((1 << tr_ctors) - 1);
        }
    }
    size_t best = x;
    if (y - x <= 2) {
        _partiotions[x][y - 1] = {best, metric};
        return;
    }
    for (size_t z = x + 1; z < y; ++z)  {
        size_t branches = mt.get_g_s_p(z - x, y - x).size();
        size_t add; 
        if (u) {
            auto const& mt1 = _special_matrices[x][z-1];
            auto const& mt2 = _special_matrices[z][y - 1];

            size_t tr1 = mt1.get_c_tr_ctors_number(z - x);
            size_t tr2 = mt2.get_c_tr_ctors_number(y - z);
            size_t nu = get_nu(x, z, y); 
            add = ((1 << (mt1.size() - tr1 - nu)) + (1 << (mt2.size() - tr2 - nu))) * (nu == 0 || nu == 1 ? nu : (2 * (nu - 1) * ((1 << nu) - 1) - 1)) +
                ((1 << (_special_matrices[x][y - 1].size() - tr_ctors))) * ((1.0 + 1.0 / (1 << nu)) * (1 << branches) - 1);
        } else {
            add = (1 << y_dim) * ((1 << branches) * 2 - 1); 
        }
        size_t phi = _partiotions[x][z - 1].second + _partiotions[z][y - 1].second + add;
        if (phi < metric) {
            best = z;
            metric = phi;
        }
    }

    _partiotions[x][y - 1] = {best, metric};
}

void trellis_based_rml_decoder::count_partitions() {
    for (size_t rng = 1; rng <= _gen[0].size(); ++rng) {
        for (size_t start = 0; start <= _gen[0].size() - rng; ++start) {
            count_partition(start, start + rng);
        }
    }

    // free memory
    for (size_t xs = 0; xs < _special_matrices.size(); ++xs) {
        for (size_t ys = xs; ys < _special_matrices[xs].size(); ++ys) {
            _special_matrices[xs][ys].clear();
        }
    }
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
        if (mt[i].leading() >= x && mt[i].trailing() < y) {
            using std::swap;
            swap(mt[i], mt[c_tr_ctors]);
            ++c_tr_ctors;
        }
    }
    
    mt = mt.puncture(x, y);
    
    auto res = mt.resolve_basis_gaussian();

    // we assume here, that matrix is "OK", and rows in generating matrix form a k-basis, 
    // so first c_tr_ctors are linear independant

    for (size_t i = 0; i < c_tr_ctors; ++i) {
        res[i].resize(y - x + res.size() - c_tr_ctors);
    }


    for (size_t i = c_tr_ctors; i < res.size(); ++i) {
        res[i].resize(y - x + res.size() - c_tr_ctors);
        res[i][y - x + i - c_tr_ctors] = true;
    }

    res.make_tof();

    _special_matrices[x][y-1] = res;
}

linalg::lin_vector trellis_based_rml_decoder::get_coset_vect(size_t x, size_t y, linalg::lin_vector const& vect) const {
    linalg::matrix const& ctors = _special_matrices[x][y-1];

    linalg::lin_vector cur = vect;

    for (size_t i = 0; i < ctors.size(); ++i) {
        if (cur[ctors[i].leading()] && ctors[i].trailing() < (y - x)) {
            for (size_t j = 0; j < y - x; ++j) {
                cur[j] = (cur[j] != ctors[i][j]);
            }
            
        }
    }
    return cur;
}


// gen matrix should be in TOF
void trellis_based_rml_decoder::build_special_trellis(size_t x, size_t z, size_t y) {

    linalg::matrix const& gen = _special_matrices[x][y + x - 1];
    
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
        
    auto g_s_ind = gen.get_g_s(z);
    auto g_s = gen.retrieve(g_s_ind);

    auto g_s_p = gen.retrieve(gen.get_g_s_p(z, y));
    auto g_s_s = g_s.retrieve(g_s.get_g_s(y));

    for (auto && v : g_s_s) {
        g_s_p.emplace_back(std::move(v));
    }


    auto g_f_s_ind = gen.get_g_f_s(z, y);


    auto g_f_s_punctured = gen.retrieve(g_f_s_ind).puncture(z, y);


    linalg::lin_vector a;
    a.resize(g_s.size());

    size_t branches = gen.retrieve(gen.get_g_s(z)).get_g_s(y).size();

    for (size_t i = 0; i < (1 << z_dim); ++i) {
        linalg::lin_vector a_star, a_a_star;
        a_star.resize(y_dim - branches);
        a_a_star.insert(a_a_star.begin(), a.begin() + g_s.size() - branches, a.end());
        for (size_t j = branches; j < y_dim; ++j) {
            a_a_star.push_back(false);
        }

        linalg::lin_vector p_u = a.empty() ? linalg::lin_vector(y - z) : (a * g_s_p).puncture(z, y);

        for (size_t j = 0; j < (1 << (y_dim - branches)); ++j) {                
            linalg::lin_vector coset = get_coset_vect(x + z, x + y,a_star.empty() || g_f_s_punctured.empty() ? p_u :  p_u + a_star * g_f_s_punctured);
            result._sections[0][a.to_bit_mask()]._next[coset] = a_a_star.to_bit_mask();
            result._sections[1][a_a_star.to_bit_mask()]._next[coset] = a.to_bit_mask();

            ++a_star;
            ++a_a_star;
        }
        ++a;
    }

    // find parallel components
    size_t parallel_component_mask = (1 << branches) - 1;

    result._parallel_components.resize(1 << branches);
    for (size_t i = 0; i < (1 << branches); ++i) {
        result._parallel_components[i].resize(2);
    }


    for (size_t i = 0; i < (1 << z_dim); ++i) {
        result._parallel_components[i & parallel_component_mask][0].push_back(i);
    }


    for (size_t i = 0; i < (1 << y_dim); ++i) {
        result._parallel_components[(i >> (y_dim - branches)) & parallel_component_mask][1].push_back(i);
    }

    // init incoming cosets
    for (size_t state_mask_z = 0; state_mask_z < (1 << z_dim); ++state_mask_z) {
        linalg::lin_vector a_star(state_mask_z, z_dim);

        assert(a_star.empty() || (a_star * g_s).size() == z);

        linalg::lin_vector b = get_coset_vect(x, x + z, a_star.empty() ? linalg::lin_vector(z) : (a_star * g_s_p).puncture(0, z));

        result._sections[0][state_mask_z]._incoming_coset = b;
        for (auto const& br : result._sections[0][state_mask_z]._next) {
            if (!result._sections[1][br.second]._incoming_coset.empty()) { // escape multiple time initialization
                break;
            }
            result._sections[1][br.second]._incoming_coset = get_coset_vect(x, x + y, b.concat(br.first));
            _cbt[x][x + y - 1][result._sections[1][br.second]._incoming_coset] = {linalg::lin_vector(), std::numeric_limits<double>::infinity()};
        }

    }

    _trellises[x][y+x-1] = result;
}

void trellis_based_rml_decoder::build_gray_codes(size_t n) {
    if (_gray_codes.empty()) {
        _gray_codes.push_back(std::vector<linalg::lin_vector>(1, linalg::lin_vector()));
    }
    for (size_t i = _gray_codes.size(); i <= n; ++i) {
        _gray_codes.push_back(std::vector<linalg::lin_vector>());
        for (size_t j = 0; j < _gray_codes[i - 1].size(); ++j) {
            linalg::lin_vector cur(_gray_codes[i-1][j]);
            cur.push_back(false);
            _gray_codes[i].push_back(cur);
        }
        for (ptrdiff_t j = _gray_codes[i - 1].size() - 1; j >= 0; --j) {
            linalg::lin_vector cur(_gray_codes[i-1][j]);
            cur.push_back(true);
            _gray_codes[i].push_back(cur);
        }
    }
}

void trellis_based_rml_decoder::make_uniform_decomposition(size_t x, size_t z, size_t y) {

    _trellises[x][y - 1].partition_parallel_components();

    _trellises[x][y - 1].partition_components();

    _trellises[x][y - 1].init_branches_arrays();
}

size_t trellis_based_rml_decoder::get_nu(size_t x, size_t z, size_t y) {
    linalg::matrix a;
    linalg::matrix b;
    linalg::matrix c;
    linalg::matrix c_tr;

    size_t c_tr_cnt = 0;


    for (auto const& v :_gen) {
        size_t tr = v.trailing();
        size_t le = v.leading();
        if (tr < y && le >= x) {
            a.push_back(v.puncture(x, z));
        }
        if (v.all_zeros(z, y)) {
            b.push_back(v.puncture(x, z));
        }
        if (tr < z && le >= x) {
            ++c_tr_cnt;
        }
    }
    a = a.resolve_basis_gaussian();
    b = b.resolve_basis_gaussian();

    for (auto const& v : a) {
        c.push_back(v);
    }
    for (auto const& v : b) {
        c.push_back(v);
    }
    c = c.resolve_basis_gaussian();
    return a.size() + b.size() - c.size() - c_tr_cnt;
}

void trellis_based_rml_decoder::prepare_make_cbt_i(size_t x, size_t y) {
    if (x >= y) {
        return;
    }
    
    linalg::matrix const& mt = _special_matrices[x][y-1];

    linalg::matrix coset_basis;
    linalg::matrix cosets;
    bool flag = false;
    size_t i = 0;

    for (; i < mt.size(); ++i) {
        if (mt[i].trailing() >= y - x) {
            flag = true;
        }
        if (flag) {
            cosets.emplace_back(mt[i].puncture(0, y - x));
        } else {
            coset_basis.emplace_back(mt[i].puncture(0, y - x));
        }
    }

    for (size_t k = 0; k < (1 << cosets.size()); ++k) {
        linalg::lin_vector a = cosets.empty() ? linalg::lin_vector(y - x) : linalg::lin_vector(k, cosets.size()) * cosets;
        linalg::lin_vector coset_vect = get_coset_vect(x, y, a);
        auto opposite = linalg::lin_vector();
        _cbt[x][y-1][coset_vect] = {linalg::lin_vector(y - x, false), std::numeric_limits<double>::infinity()};

        for (size_t j = 0; j < (1 << coset_basis.size()); ++j) {
            linalg::lin_vector vect = !coset_basis.empty() ? linalg::lin_vector(j, coset_basis.size()) * coset_basis + a : a;
            _ctors_for_make_cbt[x][y-1].push_back({{coset_vect, opposite}, vect});
        }
    }

}

void trellis_based_rml_decoder::prepare_make_cbt_g(size_t x, size_t y) {
    if (x >= y) {
        return;
    }

    linalg::matrix mt = _special_matrices[x][y-1];
    mt = mt.puncture(0, y - x);
    size_t basis_size = mt.size();
    build_gray_codes(y - x);

    for (size_t i = 0; i < _gray_codes[y - x].size() / 2; ++i) {
        auto const& v = _gray_codes[y - x][i];
        
        // check if v in basis
        mt.push_back(v);
        auto res = mt.resolve_basis_gaussian();
        auto coset = get_coset_vect(x, y, v);
        if (res.size() == mt.size() - 1) {
            auto opposite = get_coset_vect(x, y, -linalg::lin_vector(v));
            if (opposite == coset) {
                opposite = linalg::lin_vector();
            } else {
            }
            _cbt[x][y-1][coset] = {linalg::lin_vector(y - x, false), std::numeric_limits<double>::infinity()};
            if (!opposite.empty()) {
                _cbt[x][y-1][opposite] = {linalg::lin_vector(y - x, false), std::numeric_limits<double>::infinity()};
            }
            _ctors_for_make_cbt[x][y-1].push_back({{coset, opposite}, v});
        }
        mt.pop_back();
    }
}



void trellis_based_rml_decoder::make_cbt_i(size_t x, size_t y, hamming_metric const& metric) {
    for (auto const& p : _ctors_for_make_cbt[x][y-1]) {
        double score = metric.count(p.second, x, y);
        _additions += (y - x - 1);
        if (!std::isinf(_cbt[x][y - 1][p.first.first].second)) {
            ++_comparisons;
        }
        if (_cbt[x][y - 1][p.first.first].second > score) {
            _cbt[x][y - 1][p.first.first] = std::make_pair(p.second, score);
        }
    }
}

void trellis_based_rml_decoder::make_cbt_g(size_t x, size_t y, hamming_metric const& metric) {
    double score = metric.count(_ctors_for_make_cbt[x][y-1][0].second, x, y);
    _additions += y - x - 1;
    linalg::lin_vector const* prev = &_ctors_for_make_cbt[x][y-1][0].second;
    double prev_metric = 0;
    bool f = false;
    // i suppose that there is a small mistake in the paper and it is meant that we need 2^(k(p_x,y(c))) gray code ctors
    for (auto const& p : _ctors_for_make_cbt[x][y-1]) { 
        if (f) {
            ++_additions;
            score += metric.count_diff(*prev, p.second, x, y, _additions);
        } else {
            f = true;
        }
        prev = &p.second;
        if (p.first.second.empty()) {
            if (!std::isinf(_cbt[x][y - 1][p.first.first].second)) {
                ++_comparisons;
            }
            if (std::isinf(_cbt[x][y - 1][p.first.first].second) || fabs(_cbt[x][y - 1][p.first.first].second) < fabs(score)) {
                if (std::signbit(score)) {
                    _cbt[x][y - 1][p.first.first] = std::make_pair(p.second, score);
                } else {
                    _cbt[x][y - 1][p.first.first] = std::make_pair(-linalg::lin_vector(p.second), -score);
                }
            }
        } else {
            if (!std::isinf(_cbt[x][y - 1][p.first.first].second)) {
                ++_comparisons;
            }
            if (!std::isinf(_cbt[x][y - 1][p.first.second].second)) {
                ++_comparisons;
            }
            if (std::isinf(_cbt[x][y - 1][p.first.first].second) || _cbt[x][y - 1][p.first.first].second > score) {
                _cbt[x][y - 1][p.first.first] = std::make_pair(p.second, score);
            }
            if (std::isinf(_cbt[x][y - 1][p.first.second].second) || _cbt[x][y - 1][p.first.second].second > -score) {
                _cbt[x][y - 1][p.first.second] = std::make_pair(-linalg::lin_vector(p.second), -score);
            }
        }
        
    }
}

void trellis_based_rml_decoder::comb_cbt_v(size_t x, size_t y, hamming_metric const& metric) {
    size_t z = _partiotions[x][y-1].first;
    if (z == x) {
        if (g) {
            make_cbt_g(x, y, metric);
        } else {
            make_cbt_i(x, y, metric);
        }
        return;
    }
    comb_cbt_v(x, z, metric);
    comb_cbt_v(z, y, metric);

    auto const& tr = _trellises[x][y-1];

    for (size_t y_state = 0; y_state < tr._sections[1].size(); ++y_state) {
        linalg::lin_vector const& coset = tr._sections[1][y_state]._incoming_coset;
        auto & r = _cbt[x][y-1][coset];

        for(auto const& branch : tr._sections[1][y_state]._next) {
            linalg::lin_vector const& bf = branch.first;
            linalg::lin_vector const& b = tr._sections[0][branch.second]._incoming_coset;

            
            ++_additions;
            double metric = _cbt[x][z - 1][b].second + _cbt[z][y - 1][bf].second;
            if (!std::isinf(r.second)) {
                ++_comparisons;
            }
            if (std::isinf(r.second) || metric < r.second) {
                r = std::make_pair(_cbt[x][z - 1][b].first.concat(_cbt[z][y - 1][bf].first), metric);
            }
        }
    }
}

void trellis_based_rml_decoder::comb_cbt_u(size_t x, size_t y, hamming_metric const& metric)  {
    size_t z = _partiotions[x][y-1].first;
    if (z == x) {
        if (g) {
            make_cbt_g(x, y, metric);
        } else {
            make_cbt_i(x, y, metric);
        }
        return;
    }
    comb_cbt_u(x, z, metric);
    comb_cbt_u(z, y, metric);

    // need for each group init incoming ctors (already do) and sort them
    linalg::matrix const& trellis_gen = _special_matrices[x][y-1];

    auto & tr = _trellises[x][y-1];

    for (auto & comp : tr._branches) {
        for (auto & group : comp) {
            for (auto & gr : group) {
                for (auto & it : gr) {
                    it.first = _cbt[z][y-1][it.second].second;
                }
                // f2
                std::sort(gr.begin(), gr.end(), [&](std::pair<double, linalg::lin_vector> const& a, std::pair<double, linalg::lin_vector> const & b) {
                    ++_comparisons;
                    return a.first < b.first;
                });
            }
        }
    }

    for (auto & comp : tr._inner_branches) {
        for (auto & group : comp) {
            for (auto & it : group) {
                it.first = _cbt[x][z-1][tr._sections[0][it.second]._incoming_coset].second;
            }
            // f1
            std::sort(group.begin(), group.end(),[&](std::pair<double, size_t> const& a, std::pair<double, size_t> const & b) {
                ++_comparisons;
                return a.first < b.first;
            });
        }
    }

    for (size_t comp_group = 0; comp_group < tr._parallel_component_groups.size(); ++comp_group) {
        for (size_t component : tr._parallel_component_groups[comp_group]) {
            for (size_t i = 0; i < tr._groups[component][0].size(); ++i) {
                for (size_t j = 0; j < tr._groups[component][1].size(); ++j) {
                    size_t nu = tr._groups[comp_group][0][i].size();
                    // considering left and right groups

                    auto cmp = [&](std::pair<double, std::pair<size_t, size_t>> & a, std::pair<double, std::pair<size_t, size_t>> & b) {
                        if (a.second.first >= b.second.first && a.second.second >= b.second.second) {
                            return true;
                        }
                        if (a.second.first <= b.second.first && a.second.second <= b.second.second) {
                            return false;
                        }
                        if (std::isinf(a.first)) {
                            ++_additions;
                            a.first = tr._inner_branches[component][i][a.second.first].first + tr._branches[comp_group][i][j][a.second.second].first;
                        }
                        if (std::isinf(b.first)) {
                            ++_additions;
                            b.first = tr._inner_branches[component][i][b.second.first].first + tr._branches[comp_group][i][j][b.second.second].first;
                        }
                        ++_comparisons;
                        return a.first > b.first;
                    };

                    // use touples instead of pairs
                     std::priority_queue<std::pair<double, std::pair<size_t, size_t>>, std::vector<std::pair<double, std::pair<size_t, size_t>>>, decltype(cmp)> pq (tr._heap_storage.begin(), tr._heap_storage.end(), cmp);

                    size_t cnt = 0;

                    for (size_t b_z = 0; b_z < nu; ++b_z) {
                        pq.push({std::numeric_limits<double>::infinity(), {b_z, 0}});
                    }

                    while (cnt < nu) { // f3
                        auto const& cur = pq.top();
                        double metric = cur.first;
                        size_t b_z = cur.second.first;
                        size_t b_y = cur.second.second;
                        size_t left_node = tr._inner_branches[component][i][b_z].second;
                        linalg::lin_vector const& transition = tr._branches[comp_group][i][j][b_y].second;
                        linalg::lin_vector const& cur_coset = tr._sections[1][tr._sections[0][left_node]._next[transition]]._incoming_coset;
                        if (tr._group_cache[component][i][j].find(cur_coset) == tr._group_cache[component][i][j].end() || std::isinf(tr._group_cache[component][i][j][cur_coset].second)) {
                            if (std::isinf(metric)) {
                                ++_additions;
                                metric = tr._inner_branches[component][i][b_z].first + tr._branches[comp_group][i][j][b_y].first;
                            }
                            tr._group_cache[component][i][j][cur_coset] = {_cbt[x][z-1][tr._sections[0][left_node]._incoming_coset].first.concat(_cbt[z][y-1][transition].first), metric};
                            ++cnt;
                        }
                        pq.pop();
                        linalg::lin_vector const* trans;
                        linalg::lin_vector const* cst;
                        do { // skip b_ys if the best metric was already found
                            ++b_y;
                            if (b_y < nu) {
                                trans = &tr._branches[comp_group][i][j][b_y].second;
                                cst = &tr._sections[1][tr._sections[0][left_node]._next[*trans]]._incoming_coset;
                            }
                        } while(b_y < nu && !(tr._group_cache[component][i][j].find(*cst) == tr._group_cache[component][i][j].end() || std::isinf(tr._group_cache[component][i][j][*cst].second)));
                        if (b_y < nu) {
                            pq.push({std::numeric_limits<double>::infinity(), {b_z, b_y}});  
                        }
                    }

                    for (auto & br : tr._group_cache[component][i][j]) { // f4
                        if (!std::isinf(_cbt[x][y-1][br.first].second)) {
                            ++_comparisons;
                        }
                        if (std::isinf(_cbt[x][y-1][br.first].second) || br.second.second < _cbt[x][y-1][br.first].second) {
                            _cbt[x][y-1][br.first] = br.second;
                        }
                        br.second.second = std::numeric_limits<double>::infinity();
                    }
                }
            }
        }
    }
}

linalg::lin_vector trellis_based_rml_decoder::decode(linalg::lin_vector const& decisions, std::vector<double> const& rels)  {
    hamming_metric metric(decisions, rels);

    if (u) {
        comb_cbt_u(0, decisions.size(), metric);
    } else {
        comb_cbt_v(0, decisions.size(), metric);
    }
    
    auto res = (*_cbt[0][decisions.size() - 1].begin()).second.first;
    for (size_t xs = 0; xs < _cbt.size(); ++xs) { // reset _cbt
        for (size_t ys = 0; ys < _cbt[xs].size(); ++ys) {
            for (auto& br : _cbt[xs][ys]) {
                _cbt[xs][ys][br.first] = {linalg::lin_vector(), std::numeric_limits<double>::infinity()};
            }
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
    _additions = 0;
    _comparisons = 0;
    std::vector<double> reliabilities;
    linalg::lin_vector hard_decisions;
    for (auto const& i : signals) {
        if (i >= 0.0) {
            hard_decisions.push_back(true);
        } else {
            hard_decisions.push_back(false);
        }
        reliabilities.push_back(std::fabs(i));
    }
    return decode(hard_decisions, reliabilities);
}
}
