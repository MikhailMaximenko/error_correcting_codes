#include "linalg.h"
#include "encoder.h"
#include <cassert>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>




std::random_device device;
std::mt19937 gen(device());
std::bernoulli_distribution bern(0.5);


linalg::lin_vector string_to_vec(std::string const& s) {
    linalg::lin_vector vec;
    vec.reserve(s.size());
    for (char ch : s) {
        vec.push_back(((ch - '0') != 0));
    }
    return vec;
}
    
linalg::lin_vector gen_vect(size_t sz) {
    linalg::lin_vector res;
    res.reserve(sz);

    for (size_t i = 0; i < sz; ++i) {
        res.push_back(bern(gen));
    }

    return res;
}

std::vector<double> send(linalg::lin_vector const& message, std::normal_distribution<double> & distr) {
    std::vector<double> res;
    res.reserve(message.size());
    for (size_t i = 0; i < message.size(); ++i) {
        double noise = distr(gen);
        if (message[i]) {
            res.push_back(1.0 + noise);
        } else {
            res.push_back(-1.0 + noise);
        }
    }

    return res;
}

// returns number of correct decoded messages
size_t emulate(size_t n, size_t k, linalg::matrix const& gen_matrix, size_t w_max, 
            size_t iterations, double var) {
                
    
    encoding::encoder enc = {gen_matrix};
    encoding::decoder dec = {gen_matrix, w_max};
    std::normal_distribution<double> norm(0.0, var);
    size_t cnt = 0;
    std::vector<double> rels(n, 1.0);
    for (size_t i = 0; i < iterations; ++i) {
        auto message = gen_vect(k);
        auto encoded = enc.encode(message);
        // std::cout << "testing with: " << message.to_string() << " " << encoded.to_string() << "\n\n";
        auto decoded = dec.decode(send(encoded, norm));
        if (decoded == message) {
            ++cnt;
        }

    }
    return cnt;
}

void make_and_put_samples(std::ofstream &out, linalg::matrix const& gen_matrix, size_t n, size_t k, size_t w_max, size_t iters, double step) {
    double signal_noise = -2;
    
    while (signal_noise <= 7) {
        // E / N = signal_noise in W, E = 1, 
        // so N = 1 / signal_noise in W =  1 / 10^(0.1 signal_noise in dB) = 10^(-0.1*signal_noise)
        double noise_var = std::pow(10, -0.1*signal_noise);
        auto res = emulate(n, k, gen_matrix, w_max, iters, noise_var);
        out << w_max << " " << iters << " " << res << " " << signal_noise << '\n';

        signal_noise += step;
    }
}

int main(int argc, const char* argv[]) {
    if (argc < 6) {
        std::cout << "expected args: input_file, output_file, w_max, iters, step\n";
        return -1;
    }

    std::ifstream in(argv[1]);
    if (in.bad()) {
        std::cout << "could not open input file\n";
        return -1;
    }
    
    std::ofstream out(argv[2]);
    if (out.bad()) {
        std::cout << "could not open input file\n";
        return -1;
    }

    size_t w_max, iters;
    double step;

    try {
        w_max = std::stoull(argv[3]);
        iters = std::stoull(argv[4]);
        step = std::stod(argv[5]);
    } catch(std::invalid_argument const& e) {
        std::cout << "could not parse numeric args: " << e.what();
        return -1;
    } catch(std::out_of_range const& e) {
        std::cout << "too big numeric arg: " << e.what();
        return -1;
    }

    size_t n, k;
    in >> n >> k;

    linalg::matrix gen;
    for (size_t i = 0; i < k; ++i) {
        std::string s;
        in >> s;
        gen.push_back(string_to_vec(s));
    }

    auto basis = gen.resolve_basis();

    for (size_t i = 0; i < basis.first.size(); ++i) {
        if (basis.first[i] != i) {
            for (size_t j = 0; j < gen.size(); ++j) {
                using std::swap;
                swap(gen[j][i], gen[j][basis.first[i]]);
            }
        }
    }
    make_and_put_samples(out, gen, n, k, w_max, iters, step);
    return 0;
}