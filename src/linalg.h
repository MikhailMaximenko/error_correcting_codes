#pragma once


#include <cstddef>
#include <cstdint>
#include <functional>
#include <set>
#include <string>
#include <vector>

namespace linalg {

struct matrix;

struct lin_vector : std::vector<bool> {
    using std::vector<bool>::vector;

    lin_vector() = default;
    lin_vector(lin_vector const&) = default;
    lin_vector(lin_vector &&) = default;

    lin_vector& operator=(lin_vector const&) = default;
    lin_vector& operator=(lin_vector &&) = default;

    ~lin_vector() = default;


    lin_vector(size_t, size_t);

    lin_vector& operator+=(lin_vector const&);
    lin_vector operator+(lin_vector const&) const;
    int64_t operator*(lin_vector const&) const;

    lin_vector& operator-();

    lin_vector & multiply(lin_vector const&);

    void permutate(std::vector<size_t> const&);
    lin_vector operator*(matrix const&) const;


    lin_vector& operator++() noexcept;

    size_t leading() const noexcept;
    size_t trailing() const noexcept;

    bool all_zeros(size_t, size_t) const;

    lin_vector puncture(size_t, size_t) const;
    lin_vector concat(lin_vector const &) const;

    uint64_t to_bit_mask() const;
    std::string to_string() const;

};

struct matrix : std::vector<lin_vector> {
    using std::vector<lin_vector>::vector;

    matrix transpose() const;

    lin_vector operator*(lin_vector const&) const;
    matrix operator*(matrix const&) const;

    matrix inverse() const;

    void permutate(std::vector<size_t> const&);

    void make_tof();

    std::vector<size_t> get_g_s(ptrdiff_t) const;

    std::vector<size_t> get_g_f(ptrdiff_t) const;

    std::vector<size_t> get_g_p(ptrdiff_t) const;

    std::vector<size_t> get_g_f_s(ptrdiff_t, ptrdiff_t) const;

    std::vector<size_t> get_g_s_p(ptrdiff_t, ptrdiff_t) const;

    matrix puncture(size_t, size_t) const;

    matrix resolve_basis_gaussian();

    matrix retrieve(std::vector<size_t> const&) const;

    std::pair<std::vector<size_t>, matrix> resolve_basis();

    std::string to_string() const;

    matrix operator+(matrix const&) const;

    lin_vector get_and_multiply(lin_vector const&) const;

    size_t get_c_tr_ctors_number(size_t) const;

    bool is_tof() const noexcept;
};

} // namespace linalg

template<>
struct std::hash<linalg::lin_vector> {
    std::hash<std::vector<bool>> hasher = std::hash<std::vector<bool>>();
    size_t operator()(const linalg::lin_vector& lv) const {
        return hasher(static_cast<std::vector<bool>const &>(lv));
    }
};