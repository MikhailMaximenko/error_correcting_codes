
#pragma once

#include "linalg.h"
#include "trellis.h"
#include <cmath>
#include <cstddef>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

namespace encoding {

linalg::matrix build_rm_code(size_t, size_t);

struct encoder {
    linalg::bit_matrix generating;

    linalg::bit_vector encode(linalg::bit_vector const&) const;
};


struct hamming_metric {
    linalg::bit_vector _given;
    std::vector<double> _rels;


    hamming_metric(linalg::bit_vector const& vect, std::vector<double> const& rels);

    double count(linalg::bit_vector const& vect) const;

    double count(linalg::bit_vector const& vect, size_t begin, size_t end) const;
    double count_diff(linalg::bit_vector const& vect1, linalg::bit_vector const& vect2, size_t begin, size_t end, size_t & add_cnt) const;
};

struct trellis_based_rml_decoder {
    const linalg::bit_matrix _gen;
    bool g{false};
    bool u{false};

    size_t _comparisons{0};
    size_t _additions{0};

    std::vector<std::vector<linalg::bit_matrix>> _special_matrices; // lets store all the special matrices as we need to compute it ones

    std::vector<std::vector<std::vector<std::pair<size_t, size_t>>>> _leading_trailing_cache;

    std::vector<std::vector<std::pair<size_t, size_t>>> _partiotions; // store all partitions

    std::vector<std::vector<linalg::bit_vector>> _gray_codes; // contains all the gray codes of dims 0 to max(p_x_y(C)), where x, y are bounds of make_cbt procedure

    std::vector<std::vector<std::vector<std::pair<std::pair<linalg::bit_vector, linalg::bit_vector>, std::pair<linalg::bit_vector, linalg::bit_vector>>>>> _ctors_for_make_cbt; // stores corresponding coset of vect and opposite coset and vect itself; if (g) then it stores ctors in gray order

    std::vector<std::vector<trellis>> _trellises;
    // we do not need to store shifts as we find basis vectors by gaussian method, however it might have some small speed up
    // we can also store all the trellises
    std::vector<std::vector<std::unordered_map<linalg::bit_vector, std::pair<linalg::bit_vector, double>>>> _cbt;

    hamming_metric _metric;
    linalg::bit_vector _result;


    trellis_based_rml_decoder(linalg::bit_matrix const& mt, bool use_grey_code, bool use_uniform_optimization);

    void count_partition(size_t x, size_t y);
    void count_partitions();

    void build_special_matrix(size_t x, size_t y);

    void init(size_t, size_t);

    linalg::bit_vector get_coset_vect(size_t x, size_t y, linalg::bit_vector const& vect) const;


    // gen matrix should be in TOF
    void build_special_trellis(size_t x, size_t z, size_t y);

    void make_uniform_decomposition(size_t x, size_t z, size_t y);

    size_t get_nu(size_t x, size_t z, size_t y);

    void build_gray_codes(size_t n);

    void prepare_make_cbt_i(size_t x, size_t y);
    
    void prepare_make_cbt_g(size_t x, size_t y);

    void make_cbt_i(size_t x, size_t y);

    void make_cbt_g(size_t x, size_t y);

    void comb_cbt_v(size_t x, size_t y);

    void comb_cbt_u(size_t x, size_t y);

    linalg::bit_vector decode();

    linalg::bit_vector decode(std::vector<double> const& signals);
    

};




}