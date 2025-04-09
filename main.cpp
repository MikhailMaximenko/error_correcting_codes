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
std::pair<size_t, size_t> emulate(size_t n, size_t k, linalg::matrix const& gen_matrix, size_t w_max, 
            size_t errors, double deviation) {
                
    
    encoding::encoder enc = {gen_matrix};
    encoding::decoder dec = {gen_matrix, w_max};
    std::normal_distribution<double> norm(0.0, deviation);
    size_t iterations = 0;
    size_t cnt = 0;
    std::vector<double> rels(n, 1.0);
    while (iterations - cnt < errors) {
        auto message = gen_vect(k);
        auto encoded = enc.encode(message);
        auto decoded = dec.decode(send(encoded, norm));
        if (decoded == message) {
            ++cnt;
        }
        if (iterations % 100000 == 0) {
            std::cout << iterations << " " << deviation << '\n';
        }
        ++iterations;
    }
    return std::make_pair(iterations, cnt);
}

void make_and_put_samples(std::ofstream &out, linalg::matrix const& gen_matrix, size_t n, size_t k, size_t w_max, size_t errors, double step) {
    double signal_noise = -2;
    
    while (signal_noise <= 5.5) {
        // Eb / N = signal_noise,
        // Es = 1
        // Es / (Rc*N) = signal_noise,
        // 1 / (Rc*N) = signal_noise,
        // N = 1 / (Rc*signal_noise) = n / (k * signal_noise)
        // so N = 1 / signal_noise =  n / (k * 10^(0.1 signal_noise)) = n * 10^(-0.1*signal_noise) / k 
        double twice_noise_var = std::pow(10, -0.1*signal_noise) * n / k;
        // N = 2 * \sigma ^ 2 -> \sigma = sqrt(N / 2)
        auto res = emulate(n, k, gen_matrix, w_max, errors, sqrt(twice_noise_var / 2));
        out << w_max << " " << res.first << " " << res.second << " " << static_cast<double>(res.second) / res.first << " " << signal_noise  << '\n';

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

    size_t w_max, errors;
    double step;

    try {
        w_max = std::stoull(argv[3]);
        errors = std::stoull(argv[4]);
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
    make_and_put_samples(out, gen, n, k, w_max, errors, step);
    return 0;
}