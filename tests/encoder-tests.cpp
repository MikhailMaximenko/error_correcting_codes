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

// TEST(decoding_correctness, decode_without_noise1) {
//     linalg::matrix gen = {{1, 0, 0, 1}, 
//                           {1, 1, 0, 0}};
//     encoding::encoder e = {gen};

    
//     linalg::lin_vector vect = {1, 1};
//     auto encoded = e.encode(vect);
//     std::vector<double> reliabilities = {1, 1, 1, 1};
//     encoding::decoder dec = {gen, 0};

//     linalg::matrix check = {
//         {1, 0, 0, 1},
//         {0, 1, 0, 1}
//     };
//     gen.resolve_basis();
//     EXPECT_EQ(gen, check);
//     auto decoded = dec.decode(encoded, reliabilities);
//     EXPECT_EQ(vect, decoded);
// }


// TEST(decoding_correctness, decode_without_noise2) {
//     linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}};
//     encoding::encoder e = {gen};
//     linalg::lin_vector vect = {1, 0};
//     auto encoded = e.encode(vect);
//     std::vector<double> reliabilities = {1, 1, 1, 1};
//     encoding::decoder dec = {gen, 0};
//     auto decoded = dec.decode(encoded, reliabilities);
//     EXPECT_EQ(vect, decoded);
// }

// TEST(decoding_correctness, decode_without_noise_w_max) {
//     linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}};
//     encoding::encoder e = {gen};
//     linalg::lin_vector vect = {1, 0, 1};
//     auto encoded = e.encode(vect);
//     std::vector<double> reliabilities = {1, 1, 1, 1};
//     encoding::decoder dec = {gen, 2};
//     auto decoded = dec.decode(encoded, reliabilities);
//     EXPECT_EQ(vect, decoded);
// }

// TEST(decoding_correctness, decode_without_noise_full_basis) {
//     linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, 1}};
//     encoding::encoder e = {gen};
//     linalg::lin_vector vect = {1, 0, 1, 1};
//     auto encoded = e.encode(vect);
//     std::vector<double> reliabilities = {1, 1, 1, 1};
//     encoding::decoder dec = {gen, 2};
//     auto decoded = dec.decode(encoded, reliabilities);
//     EXPECT_EQ(vect, decoded);
// }

// TEST(decoding_correctness, decode_not_ordered_rels) {
//     linalg::matrix gen = {{1, 0, 0, 1}, {1, 1, 0, 0}, {1, 1, 0, 1}};
//     encoding::encoder e = {gen};
//     linalg::lin_vector vect = {1, 0, 1};
//     auto encoded = e.encode(vect);
//     std::vector<double> reliabilities = {0.9997, 0.9999, 0.9998, 0.99997};
//     encoding::decoder dec = {gen, 2};
//     auto decoded = dec.decode(encoded, reliabilities);
//     EXPECT_EQ(vect, decoded);
// }

// TEST(decoding_correctness, decode_without_noise_given_8_4) {
//     linalg::matrix gen = {
//         {1, 1, 1, 1, 0, 0, 0, 0}, 
//         {0, 0, 0, 0, 1, 1, 1, 1}, 
//         {0, 1, 0, 1, 1, 0, 1, 0},
//         {0, 0, 1, 1, 1, 1, 0, 0}
//     };
//     encoding::encoder e = {gen};
//     linalg::lin_vector vect = {1, 1, 1, 1};
//     auto encoded = e.encode(vect);
//     std::vector<double> reliabilities = {1, 1, 1, 1, 1, 1, 1, 1};
//     encoding::decoder dec = {gen, 2};
//     auto decoded = dec.decode(encoded, reliabilities);



//     EXPECT_EQ(vect, decoded);
// }

TEST(rmld_correctness, decode_correct) {
    linalg::matrix g = {
        {1, 1, 1, 1, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 1, 1, 1, 1}, 
        {0, 1, 0, 1, 1, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0}
    };

    linalg::matrix shifts(4, linalg::lin_vector(8));

    g.make_tof(shifts);

    encoding::encoder e = {g};
    linalg::lin_vector vect = {0, 0, 1, 1};
    auto encoded = e.encode(vect);

    linalg::lin_vector encoded_check = {0, 1, 1, 0, 0, 1, 1, 0};
// 10011001
    // EXPECT_EQ(encoded, encoded_check);

    std::vector<double> reliabilities = {-10, 10, 10, -10, -10, 10, 10, -10};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities);

    EXPECT_EQ(encoded_check, v);

}


TEST(rmld_correctness, decode_wrong) {
    linalg::matrix g = {
        {1, 1, 1, 1, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 1, 1, 1, 1}, 
        {0, 1, 0, 1, 1, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0}
    };

    linalg::matrix shifts(4, linalg::lin_vector(8));

    g.make_tof(shifts);

    encoding::encoder e = {g};
    linalg::lin_vector vect = {0, 0, 1, 1};
    auto encoded = e.encode(vect);

    linalg::lin_vector encoded_check = {0, 1, 1, 0, 0, 1, 1, 0};
// 10011001
    // EXPECT_EQ(encoded, encoded_check);

    std::vector<double> reliabilities = {-10, 10, 10, -10, -10, 10, -1, -10};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities);

    EXPECT_EQ(encoded_check, v);

}
// 1010 10010110 1110
// 2.02487 -0.911127 0.447729 0.00269423 -0.503658 0.66734 0.776549 0.138703 


TEST(rmld_correctness, decode_wrong2) {
    linalg::matrix g = {
        {1, 1, 1, 1, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 1, 1, 1, 1}, 
        {0, 1, 0, 1, 1, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0}
    };

    linalg::matrix shifts(4, linalg::lin_vector(8));

    g.make_tof(shifts);

    encoding::encoder e = {g};
    linalg::lin_vector vect = {1, 1, 1, 0};
    auto encoded = e.encode(vect);

    linalg::lin_vector encoded_check = {1, 0, 0, 1, 0, 1, 1, 0};
// 10011001
    // EXPECT_EQ(encoded, encoded_check);

    std::vector<double> reliabilities = {2.02487, -0.911127, 0.447729, 0.00269423, -0.503658, 0.66734, 0.776549, 0.138703};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities);

    EXPECT_EQ(encoded_check, v);

}

// 01101111 01101111010111001010001110010000 0110110101000010
// -1.8851 0.976583 0.528779 -1.95404 1.15196 1.09015 0.542011 1.03569 -1.94776 0.80406 -1.46708 1.46362 1.35543 1.31522 -0.615261 -1.19715 0.83817 -0.676406 1.6807 -1.4126 -0.872296 -0.834092 2.29801 0.533701 1.80917 -0.428016 -1.29598 1.18537 -0.991321 -0.763306 -0.979013 -0.789227

TEST(rmld_correctness, decode_wrong_large) {
    linalg::matrix g = 
    {
        {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, },
        {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, }
    };

    linalg::matrix shifts(g.size(), linalg::lin_vector(g[0].size(), false));

    g.make_tof(shifts);

    encoding::encoder e = {g};
    // linalg::lin_vector vect = {0, 1, 1, 0, 1, 1, 1, 1};
    // auto encoded = e.encode(vect);

    linalg::lin_vector encoded_check = {0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1};
// 10011001
    // EXPECT_EQ(encoded, encoded_check);

    std::vector<double> reliabilities = {-0.576804, 1.03504, 1.80771, 0.912645, -0.500792, 0.789529, 1.22792, 0.64886, -1.31845, 0.110547, 0.855006, 1.41033, 1.3656, 0.306943, -1.40199, -0.689252, 0.823275, -0.038073, -1.7673, -1.42825, 0.258857, -0.971259, -0.772703, -0.362956, 0.585685, -0.746481, -0.52608, -0.569451, -1.01687, 2.35313, 1.52079, 1.56831};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities);

    EXPECT_EQ(encoded_check, v);

}


TEST(rmld_correctness, decode_wrong_large2) {
    linalg::matrix g = 
    {
        {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, },
        {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, }
    };

    linalg::matrix shifts(g.size(), linalg::lin_vector(g[0].size(), false));

    g.make_tof(shifts);

    linalg::lin_vector encoded_check = {0,1,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,1,1,0,1,1,0,0,0};

    std::vector<double> reliabilities = {-0.944571, 1.3643, -1.21466, -1.80696, 0.887298, 0.465176, 0.838953, -1.126, -1.33395, -1.00815, -0.214494, -0.552195, -0.846559, 0.698426, 0.979281, 0.470859, 2.18831, -0.750444, 1.23306, 0.758309, -0.846314, -1.78091, 0.47336, 1.29588, 1.50673, -0.193073, -1.44247, 1.72633, 0.0470784, -1.78299, -1.44743, 0.103198};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities); // 01001110001001111011001010110100

    EXPECT_EQ(encoded_check, v);

}



TEST(rmld_correctness, decode_wrong_large3) {
    linalg::matrix g = 
    {
        {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0 },
        {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
        {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}
    };
    linalg::matrix shifts(g.size(), linalg::lin_vector(g[0].size(), false));

    g.make_tof(shifts);

    linalg::lin_vector encoded_check = {1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0};

    std::vector<double> reliabilities = {1.22055, 1.61721, 1.35983, 1.33653, -1.82543, -1.16671, 0.0244755, 1.42969, -0.985874, -1.13443, 1.03667, 0.634122, -1.09611, -0.922897, -0.0839821, -1.3403, 0.46734, -0.343023, -1.09859, 1.35925, -1.3316, 1.27189, -0.0234682, 1.42554, -0.856104, 0.39479, -1.10027, 1.32572, -1.20051, 0.211034, 0.498358, -0.762412};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities); // 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0
    EXPECT_EQ(encoded_check, v);

}

// TEST(trellis_correctness, rm2_2) {
//     linalg::matrix g = 
//     {
//         {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, 
//         {0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0},
//         {0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0},
//         {0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0},
//         {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
//         {0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0},
//         {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1} 
//     };

//     linalg::matrix shifts(g.size(), linalg::lin_vector(g[0].size(), false));

//     encoding::trellis_based_rml_decoder dec(g, false, false);

//     auto tr = dec.build_special_trellis(g, shifts, 4, 8);

//     std::cout << " printing special trellis: \n";
//     for (size_t zz = 0; zz < tr._sections[0].size(); ++zz) {
//         std::cout << zz << " transmissions:\n";
//         for (auto const& br : tr._sections[0][zz]._next) {
//             std::cout << "\t" << br.first.to_string() << " " << br.second << "\n";
//         }
//     }
//     std::cout << "\n";
// }


// TEST(get_coset_vect_correctness, check_ctor) {
//     linalg::matrix g = 
//     {
//         {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0 },
//         {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0 },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1}
//     };
//     linalg::matrix shifts(g.size(), linalg::lin_vector(g[0].size(), false));

//     g.make_tof(shifts);
//     encoding::trellis_based_rml_decoder dec(g, false, false);

    
//     // linalg::lin_vector v = {1,1,1,1,0,1,1,0,0,1,1,0,0,0,0,0};
//     // linalg::lin_vector v1 = {0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0};
//     // auto mt = dec.build_special_matrix(0, 16);
//     // mt.first.make_tof(mt.second);
//     // std::cout << mt.first.puncture(0, 16).to_string() << "\n" << mt.second.puncture(0, 16).to_string() << "\n";

//     // linalg::lin_vector v2 = dec.get_coset_vect(16, 32, v);
//     // linalg::lin_vector v3 = dec.get_coset_vect(16, 32, v1);
//     // std::cout << v2.to_string() << "\n" << v3.to_string() << "\n";
// }
