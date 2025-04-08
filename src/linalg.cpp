#include "linalg.h"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace linalg {

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
    }
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
                (*this)[j] += (*this)[not_null_pos]; // - can be replaced with + in a ring of module 2 
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

}