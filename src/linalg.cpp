#include "linalg.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace linalg {

lin_vector::lin_vector(size_t src, size_t dim) : std::vector<bool>(dim) {
    for (size_t i = 0; i < dim; ++i) {
        (*this)[i] = (src >> (dim - i - 1)) & 1;
    }
}

lin_vector& lin_vector::operator+=(lin_vector const& o) {
    if (size() != o.size()) {
        throw std::domain_error("cannot add vectors of different size");
    }
    for (size_t i = 0; i < size(); ++i) {
        (*this)[i] = (*this)[i] != o[i];
    }
    return *this;
}

lin_vector lin_vector::operator+(lin_vector const& o) const {
    lin_vector tmp(*this);
    tmp += o;
    return tmp;
}    

int64_t lin_vector::operator*(lin_vector const& o) const {
    if (size() != o.size()) {
        throw std::domain_error("cannot multiply vectors of different size");
    }
    int64_t res = 0;
    for (size_t i = 0; i < size(); ++i) {
        res += ((*this)[i] && o[i]) ? 1 : 0;
    }
    return res;
}

lin_vector lin_vector::operator*(matrix const& o) const { // expects vector to be "horisontal"
    matrix tmp;
    tmp.push_back(*this);
    return (tmp * o)[0];
}

lin_vector& lin_vector::operator-() {
    for (size_t i = 0; i < size(); ++i) {
        at(i) = at(i) != true;
    }
    return *this;
}

lin_vector& lin_vector::multiply(lin_vector const& o) {
    if (size() != o.size()) {
        throw std::logic_error("cannot multiply vectors of different size");
    }
    for (size_t i = 0; i < size(); ++i) {
        at(i) = at(i) & o.at(i);
    }
    return *this;
}

void lin_vector::permutate(std::vector<size_t> const& permutation) {
    lin_vector res(size());

    for (size_t i = 0; i < size(); ++i) {
        res[i] = (*this)[permutation[i]];
    }
    *this = std::move(res);
}

matrix matrix::transpose() const {
    if (size() == 0) {
        return matrix();
    }
    matrix res((*this)[0].size(), lin_vector(size()));
    for (size_t i = 0; i < (*this)[0].size(); ++i) {
        for (size_t j = 0; j < size(); ++j) {
            res[i][j] = (*this)[j][i];
        }
    }
    return res;
}

lin_vector matrix::operator*(lin_vector const& o) const { // expects vector to be "vertical"
    if (size() == 0) {
        throw std::domain_error("cannot multiply with empty matrix");
    }
    if ((*this)[0].size() != o.size()) {
        throw std::domain_error("cannot multiply matrix and vector of different size");
    }
    lin_vector res(size(), false);
    for (size_t i = 0; i < size(); ++i) {
        for (size_t j = 0; j < (*this)[i].size(); ++j) {
            res[i] = res[i] != ((*this)[i][j] && o[j]);
        }
    }
    return res;
}

std::string lin_vector::to_string() const {
    std::string res;
    for (auto const& it : *this) {
        if (it) {
            res.push_back('1');
        } else {
            res.push_back('0');
        }
        // res.push_back(',');
        // res.push_back(' ');
    }
    return res;
}

lin_vector& lin_vector::operator++() noexcept {
    for (ptrdiff_t i = size() - 1; i >= 0; --i) {
        if (!(*this)[i]) {
            (*this)[i] = true;
            return *this;
        } else {
            (*this)[i] = false;
        }
    }
    return *this;
}

uint64_t lin_vector::to_bit_mask() const {
    if (size() > 64) {
        throw std::domain_error("to big vector to convert to int");
    }
    uint64_t res = 0;
    for (size_t i = 0; i < size(); ++i) {
        res <<= 1;
        if ((*this)[i]) {
            ++res;
        }
    }
    return res;
}


size_t lin_vector::leading() const noexcept {
    size_t res = 0;
    for (; res < size(); ++res) {
        if (at(res)) {
            return res;
        }
    }
    return res;
}


size_t lin_vector::trailing() const noexcept {
    size_t res = size() - 1;
    for (; res >= 1; --res) {
        if (at(res)) {
            return res;
        }
    }
    return res;
}

bool lin_vector::all_zeros(size_t from, size_t to) const {
    bool flag = false;
    for (size_t i = from; i < to; ++i) {
        flag |= at(i);
    }
    return flag;
}

lin_vector lin_vector::puncture(size_t x, size_t y) const {
    if (x > size() || y < x) {
        throw std::logic_error("bad bounds for puncturing");
    }
    lin_vector res;
    res.insert(res.begin(), begin() + x, begin() + y);
    return res;
}

lin_vector lin_vector::concat(lin_vector const & o) const {
    lin_vector res(*this);
    res.insert(res.end(), o.begin(), o.end());
    return res;
}

matrix matrix::operator*(matrix const& o) const {
    if (size() == 0 || o.size() == 0) {
        throw std::domain_error("cannot multiply with empty matrix");
    }
    if ((*this)[0].size() != o.size()) {
        throw std::domain_error("cannot multiply not matched matrices");
    }

    matrix res;
    res.resize(size());


    for (size_t i = 0; i < size(); ++i) {
        res[i].resize(o[0].size());
        for (size_t j = 0; j < o[0].size(); ++j) {
            res[i][j] = false;
            for (size_t k = 0; k < o.size(); ++k) {
                res[i][j] = res[i][j] != ((*this)[i][k] && o[k][j]);
            }
        }
    }

    return res;
}

matrix matrix::inverse() const {
    
    matrix tmp = (*this);

    // prepare for gaussian method
    for (size_t i = 0; i < size(); ++i) {
        for (size_t j = 0; j < size(); ++j) {
            if (j == i) {
                tmp[i].push_back(true);
            } else {
                tmp[i].push_back(0);
            }
        }
    }

    // the basis on first size() vectors would be found, 
    // this matrix is k * 2k, at right was id matrix, 
    // so at the right there is inversed matrix found by gaussian method
    tmp.resolve_basis(); 
    matrix res(size());

    for (size_t i = 0; i < size(); ++i) {
        res[i].insert(res[i].end(), tmp[i].begin() + size(), tmp[i].end());
    }

    return res;
}

void matrix::permutate(std::vector<size_t> const& permutation) {
    for (auto & vect : *this) {
        vect.permutate(permutation);
    }
}

matrix matrix::resolve_basis_gaussian() {
    matrix res;
    for (size_t i = 0; i < size(); ++i) {
        if ((*this)[i].leading() != (*this)[i].size()) {
            res.push_back((*this)[i]);
            for (size_t j = i + 1; j < size(); ++j) {
                if ((*this)[j][(*this)[i].leading()]) {
                    (*this)[j] += (*this)[i];
                }
            }
        }
    }
    return res;
}


// returns vector of basis pos and transformation matrix to old basis
std::pair<std::vector<size_t>, matrix> matrix::resolve_basis() {
    matrix old(*this);
    std::set<size_t> used_rows;
    size_t cnt = 0;
    std::vector<size_t> basis_pos;
    matrix transformation(size(), lin_vector(size(), false));


    basis_pos.reserve(size());
    for (size_t i = 0; i < (*this)[0].size(); ++i) {
        bool is_null_vect = true;
        size_t not_null_pos = 0;
        for (size_t j = 0; j < size(); ++j) { // find first non null vector
            if ((*this)[j][i] && (used_rows.find(j) == used_rows.end())) {
                is_null_vect = false;
                not_null_pos = j;
                break;
            }
        }
        if (is_null_vect) {
            continue;
        }
        for (size_t j = 0; j < size(); ++j) {
            transformation[j][cnt] = old[j][i];
        }

        ++cnt;

        basis_pos.push_back(i); // add basis vector
        used_rows.insert(not_null_pos);
        for (size_t j = 0; j < size(); ++j) { // subtract basis vector if coord corresponding to it is equal to 1
            
            if (j == not_null_pos) {
                continue;
            }
            if ((*this)[j][i]) {
                (*this)[j] += (*this)[not_null_pos];
            }
        }

        
        if (cnt >= size()) { // if found basis then break
            break;
        }
        
    }

    // make identity matrix on basis vectors by permutations
    for (size_t i = 0; i < basis_pos.size(); ++i) { 
        for (size_t j = 0; j < size(); ++j) {
            if ((*this)[j][basis_pos[i]]) {
                using std::swap;
                swap((*this)[j], (*this)[i]);
            }
        }
    }


    return std::make_pair(std::move(basis_pos), std::move(transformation)); 
}

void matrix::make_tof() noexcept {
    for (size_t i = 0; i < size(); ++i) {
        size_t best_leading = 10e9;
        size_t best_number;
        for (size_t j = i; j < size(); ++j) {
            size_t leading = (*this)[j].leading();
            if (best_leading > leading) {
                best_leading = leading;
                best_number = j;
            }
        }
        if (best_leading == (*this)[i].size()) {
            return;
        }
        using std::swap;
        swap((*this)[i], (*this)[best_number]);
        for (size_t j = i + 1; j < size(); ++j) {
            if ((*this)[j][best_leading]) {
                (*this)[j] += (*this)[i];
            }
        }
    }

    for (ptrdiff_t i = size() - 1; i >= 0; --i) {
        for (ptrdiff_t j = i - 1; j >= 0; --j) {
            if ((*this)[j].trailing() == (*this)[i].trailing()) {
                (*this)[j] += (*this)[i]; 
            }
        }
    }

}

matrix matrix::retrieve(std::vector<size_t> const& ind) const {
    matrix res;
    for (size_t index : ind) {
        res.push_back(at(index));
    }
    return res;
}

std::vector<size_t> matrix::get_g_s(ptrdiff_t h) const {
    std::vector<size_t> res;
    for (size_t i = 0; i < size(); ++i) {
        if (h > at(i).leading() && h <= at(i).trailing()) {
            res.push_back(i);
        }
    }
    return res;
}

std::vector<size_t> matrix::get_g_f(ptrdiff_t h) const {
    std::vector<size_t> res;
    for (size_t i = 0; i < size(); ++i) {
        if (h <= at(i).leading()) {
            res.push_back(i);
        }
    }
    return res;
}

std::vector<size_t> matrix::get_g_p(ptrdiff_t h) const {
    std::vector<size_t> res;
    for (size_t i = 0; i < size(); ++i) {
        if (h > at(i).trailing()) {
            res.push_back(i);
        }
    }
    return res;
}

std::vector<size_t> matrix::get_g_f_s(ptrdiff_t f, ptrdiff_t s) const {
    if (s <= f) {
        return {};
    }
    auto indeces = get_g_f(f);
    auto snd = retrieve(get_g_f(f)).get_g_s(s);
    std::vector<size_t> res;
    for (size_t i = 0; i < snd.size(); ++i) {
        res.push_back(indeces[snd[i]]);
    }
    return res;
}

std::vector<size_t> matrix::get_g_s_p(ptrdiff_t f, ptrdiff_t s) const {
    if (s <= f) {
        return {};
    }
    auto indeces = get_g_s(f);
    auto snd = retrieve(get_g_s(f)).get_g_p(s);
    std::vector<size_t> res;
    for (size_t i = 0; i < snd.size(); ++i) {
        res.push_back(indeces[snd[i]]);
    }
    return res;
}

matrix matrix::puncture(size_t x, size_t y) const {
    matrix res;
    for (auto const& vect : *this) {
        res.emplace_back(vect.puncture(x, y));
    }
    return res;
}

std::string matrix::to_string() const {
    std::string res;
    for (auto const& vect : *this) {
        res += vect.to_string();
        res += '\n';
    }
    return res;
}

matrix matrix::operator+(matrix const& o) const {
    matrix mt(*this);
    for (size_t i = 0; i < size(); ++i) {
        mt[i] += o[i];
    }
    return mt;
}

lin_vector matrix::get_and_multiply(const lin_vector& get) const {
    if (get.size() != size()) {
        throw std::logic_error("get and multiply: sizes are not matched");
    }
    lin_vector res(at(0).size(), true);
    for (size_t i = 0; i < size(); ++i) {
        if (get.at(i)) {
            res.multiply(at(i));
        }
    }
    return res;
}

}