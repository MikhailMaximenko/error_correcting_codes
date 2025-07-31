#include "linalg.h"
#include "encoder.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
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
        if (decoded == message) {
            ++cnt;
        }
        ++iterations;
    }
    if (u) {
        return {iterations, cnt, adds, cmps, dec._partiotions[0][n - 1].second};
    }
    return {iterations, cnt, adds, cmps};
}

void make_and_put_samples(std::ofstream &out, linalg::matrix const& gen_matrix, size_t n, size_t k, size_t errors, double step, double max_noise, bool g, bool u) {
    double signal_noise = 0;
    
    while (signal_noise <= max_noise) {
        // Eb / N = signal_noise,
        // Es = 1
        // Es / (Rc*N) = signal_noise,
        // 1 / (Rc*N) = signal_noise,
        // N = 1 / (Rc*signal_noise) = n / (k * signal_noise)
        // so N = 1 / signal_noise =  n / (k * 10^(0.1 signal_noise)) = n * 10^(-0.1*signal_noise) / k 
        double twice_noise_var = std::pow(10, -0.1*signal_noise) * n / k;
        // N = 2 * \sigma ^ 2 -> \sigma = sqrt(N / 2)
        auto res = emulate(n, k, gen_matrix, errors, sqrt(twice_noise_var / 2), g, u);
        if (u) {
            out << res[0] << " " << res[1] << " " << static_cast<double>(res[1]) / res[0] << " " << signal_noise << " " << res[2] << " " << res[3] << " " << (res[2] + res[3]) << " " << res[4] << '\n';
        } else {
            out << res[0] << " " << res[1] << " " << static_cast<double>(res[1]) / res[0] << " " << signal_noise << " " << res[2] << " " << res[3] << " " << (res[2] + res[3]) << '\n';
        }

        signal_noise += step;
    }
}

int main(int argc, const char* argv[]) {
    if (argc < 8) {
        std::cout << "expected args: input_file, output_file, errors, step, max noise, g, u\n";
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

    size_t errors;
    double step, max_noise;
    bool g, u;

    try {
        errors = std::stoull(argv[3]);
        step = std::stod(argv[4]);
        max_noise = std::stod(argv[5]);
        g = std::stod(argv[6]);
        u = std::stod(argv[7]);
    } catch(std::invalid_argument const& e) {
        std::cout << "could not parse numeric args: " << e.what();
        return -1;
    } catch(std::out_of_range const& e) {
        std::cout << "too big numeric arg: " << e.what();
        return -1;
    }

    size_t n, k;
    if (u) {
        out << "fields: total iterations, correct answers, error rate, noise (dB), max adds, max comparisons, max total, upper bound\n";
    } else {
        out << "fields: total iterations, correct answers, error rate, noise (dB), max adds, max comparisons, max total\n";
    }
    linalg::matrix gen;
    in >> n >> k;

    for (size_t i = 0; i < k; ++i) {
        linalg::lin_vector v;
        std::string s;
        in >> s;
        for (char ch : s) {
            v.push_back(ch - '0');
        }
        gen.push_back(v);
    }
    gen.make_tof();

    make_and_put_samples(out, gen, n, k, errors, step, max_noise, g, u);
    return 0;
}