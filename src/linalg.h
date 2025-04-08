#pragma once


#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace linalg {

struct matrix;

struct lin_vector : std::vector<bool> {
    using std::vector<bool>::vector;

    lin_vector& operator+=(lin_vector const&);
    lin_vector operator+(lin_vector const&) const;
    int64_t operator*(lin_vector const&) const;

    void permutate(std::vector<size_t> const&);
    lin_vector operator*(matrix const&) const;

    std::string to_string() const;

};

struct matrix : std::vector<lin_vector> {
    using std::vector<lin_vector>::vector;

    matrix transpose() const;

    lin_vector operator*(lin_vector const&) const;
    matrix operator*(matrix const&) const;

    matrix inverse() const;

    void permutate(std::vector<size_t> const&);

    std::pair<std::vector<size_t>, matrix> resolve_basis();
};

} // namespace linalg