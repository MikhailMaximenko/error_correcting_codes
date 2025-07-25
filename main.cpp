#include "linalg.h"
#include "encoder.h"
#include "src/encoder.h"
#include "src/linalg.h"
#include <algorithm>
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
std::vector<size_t> emulate(size_t n, size_t k, linalg::matrix const& gen_matrix, 
            size_t errors, double deviation, bool g, bool u) {
                
    
    encoding::encoder enc = {gen_matrix};
    // encoding::trellis_based_rml_decoder checking_dec(gen_matrix, false, false);
    std::normal_distribution<double> norm(0.0, deviation);
    size_t iterations = 0;
    size_t cnt = 0;
    size_t adds = 0;
    size_t cmps = 0;
    encoding::trellis_based_rml_decoder dec(gen_matrix, g,u);
    std::vector<double> rels(n, 1.0);
    while (iterations - cnt < errors) {
        auto message = gen_vect(k);
        auto encoded = enc.encode(message);
        auto snd = send(encoded, norm);
        auto decoded = dec.decode(snd);
        adds = std::max(adds, dec._additions);
        cmps = std::max(cmps, dec._comparisons);
        // std::cout << encoded.to_string() << " " << decoded.to_string() << " " << message.to_string() << "\n";
        // auto check_decoded = checking_dec.decode(snd);
        if (decoded == message) {
            ++cnt;
        } else {
            // std::cout << message.to_string() << " " << encoded.to_string() << " " << decoded.to_string() << " " << (decoded * gen_matrix).to_string() << "\n";
            // for (auto d : snd) {
            //     std::cout << d << ", ";
            // }
            // std::cout << "\n";
        }
        ++iterations;
    }
    return {iterations, cnt, adds, cmps};
}

void make_and_put_samples(std::ofstream &out, linalg::matrix const& gen_matrix, size_t n, size_t k, size_t errors, double step, bool g, bool u) {
    double signal_noise = -2;

    signal_noise = 1;
    
    while (signal_noise <= 2) {
        // Eb / N = signal_noise,
        // Es = 1
        // Es / (Rc*N) = signal_noise,
        // 1 / (Rc*N) = signal_noise,
        // N = 1 / (Rc*signal_noise) = n / (k * signal_noise)
        // so N = 1 / signal_noise =  n / (k * 10^(0.1 signal_noise)) = n * 10^(-0.1*signal_noise) / k 
        double twice_noise_var = std::pow(10, -0.1*signal_noise) * n / k;
        // N = 2 * \sigma ^ 2 -> \sigma = sqrt(N / 2)
        auto res = emulate(n, k, gen_matrix, errors, sqrt(twice_noise_var / 2), g, u);
        out << res[0] << " " << res[1] << " " << static_cast<double>(res[1]) / res[0] << " " << signal_noise << " " << res[2] << " " << res[3] << " " << (res[2] + res[3]) << '\n';

        signal_noise += step;
    }
}

int main(int argc, const char* argv[]) {
    if (argc < 6) {
        std::cout << "expected args: output_file, errors, step, g, u\n";
        return -1;
    }

    
    std::ofstream out(argv[1]);
    if (out.bad()) {
        std::cout << "could not open input file\n";
        return -1;
    }

    size_t errors;
    double step;
    bool g, u;

    try {
        errors = std::stoull(argv[2]);
        step = std::stod(argv[3]);
        g = std::stod(argv[4]);
        u = std::stod(argv[5]);
    } catch(std::invalid_argument const& e) {
        std::cout << "could not parse numeric args: " << e.what();
        return -1;
    } catch(std::out_of_range const& e) {
        std::cout << "too big numeric arg: " << e.what();
        return -1;
    }

    size_t n, k;
    // in >> n >> k;

    linalg::matrix gen = encoding::build_rm_code(6, 2);
    n = 64;
    k = gen.size();
    make_and_put_samples(out, gen, n, k, errors, step, g, u);
    gen = encoding::build_rm_code(6, 3);
    n = 64;
    k = gen.size();
    make_and_put_samples(out, gen, n, k, errors, step, g, u);
    gen = encoding::build_rm_code(6, 4);
    n = 64;
    k = gen.size();
    make_and_put_samples(out, gen, n, k, errors, step, g, u);
    // gen = encoding::build_rm_code(5, 2);
    // n = 32;
    // k = gen.size();
    // make_and_put_samples(out, gen, n, k, errors, step, g, u);


    // make_and_put_samples(out, gen, n, k, errors, step);
    return 0;
}