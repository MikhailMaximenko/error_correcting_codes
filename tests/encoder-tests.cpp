#include "encoder.h"
#include "linalg.h"

#include <gtest/gtest.h>
#include <vector>

TEST(encoding_correctness, simple_encode) {
    linalg::matrix gen = {{true, false}, {false, true}};
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 1};
    EXPECT_EQ(vect, e.encode(vect));
}

TEST(encoding_correctness, encode) {
    linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}};
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 1};
    linalg::lin_vector expected = {0, 1, 0, 1};
    EXPECT_EQ(expected, e.encode(vect));
}

TEST(decoding_correctness, decode_without_noise1) {
    linalg::matrix gen = {{1, 0, 0, 1}, 
                          {1, 1, 0, 0}};
    encoding::encoder e = {gen};

    
    linalg::lin_vector vect = {1, 1};
    auto encoded = e.encode(vect);
    std::vector<double> reliabilities = {1, 1, 1, 1};
    encoding::decoder dec = {gen, 0};

    linalg::matrix check = {
        {1, 0, 0, 1},
        {0, 1, 0, 1}
    };
    gen.resolve_basis();
    EXPECT_EQ(gen, check);
    auto decoded = dec.decode(encoded, reliabilities);
    EXPECT_EQ(vect, decoded);
}


TEST(decoding_correctness, decode_without_noise2) {
    linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}};
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 0};
    auto encoded = e.encode(vect);
    std::vector<double> reliabilities = {1, 1, 1, 1};
    encoding::decoder dec = {gen, 0};
    auto decoded = dec.decode(encoded, reliabilities);
    EXPECT_EQ(vect, decoded);
}

TEST(decoding_correctness, decode_without_noise_w_max) {
    linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}};
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 0, 1};
    auto encoded = e.encode(vect);
    std::vector<double> reliabilities = {1, 1, 1, 1};
    encoding::decoder dec = {gen, 2};
    auto decoded = dec.decode(encoded, reliabilities);
    EXPECT_EQ(vect, decoded);
}

TEST(decoding_correctness, decode_without_noise_full_basis) {
    linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 1}};
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 0, 1, 1};
    auto encoded = e.encode(vect);
    std::vector<double> reliabilities = {1, 1, 1, 1};
    encoding::decoder dec = {gen, 2};
    auto decoded = dec.decode(encoded, reliabilities);
    EXPECT_EQ(vect, decoded);
}

TEST(decoding_correctness, decode_not_ordered_rels) {
    linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}};
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 0, 1};
    auto encoded = e.encode(vect);
    std::vector<double> reliabilities = {0.9997, 0.9999, 0.9998, 0.99997};
    encoding::decoder dec = {gen, 2};
    auto decoded = dec.decode(encoded, reliabilities);
    EXPECT_EQ(vect, decoded);
}

TEST(decoding_correctness, decode_without_noise_given_8_4) {
    linalg::matrix gen = {
        {1, 1, 1, 1, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 1, 1, 1, 1}, 
        {0, 1, 0, 1, 1, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0}
    };
    encoding::encoder e = {gen};
    linalg::lin_vector vect = {1, 1, 1, 1};
    auto encoded = e.encode(vect);
    std::vector<double> reliabilities = {1, 1, 1, 1, 1, 1, 1, 1};
    encoding::decoder dec = {gen, 2};
    auto decoded = dec.decode(encoded, reliabilities);



    EXPECT_EQ(vect, decoded);
}
