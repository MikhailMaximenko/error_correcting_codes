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

// TEST(rmld_correctness, decode_correct) {
//     linalg::matrix g = {
//         {1, 1, 1, 1, 0, 0, 0, 0}, 
//         {0, 0, 0, 0, 1, 1, 1, 1}, 
//         {0, 1, 0, 1, 1, 0, 1, 0},
//         {0, 0, 1, 1, 1, 1, 0, 0}
//     };

//     linalg::matrix shifts(4, linalg::lin_vector(8));

//     g.make_tof(shifts);

//     encoding::encoder e = {g};
//     linalg::lin_vector vect = {0, 0, 1, 1};
//     auto encoded = e.encode(vect);

//     linalg::lin_vector encoded_check = {0, 1, 1, 0, 0, 1, 1, 0};
// // 10011001
//     // EXPECT_EQ(encoded, encoded_check);

//     std::vector<double> reliabilities = {-10, 10, 10, -10, -10, 10, 10, -10};
//     encoding::trellis_based_rml_decoder dec(g, false, false);

//     auto v = dec.decode(reliabilities);

//     EXPECT_EQ(encoded_check, v);

// }


// TEST(rmld_correctness, decode_wrong) {
    // linalg::matrix g = {
    //     {1, 1, 1, 1, 0, 0, 0, 0}, 
    //     {0, 0, 0, 0, 1, 1, 1, 1}, 
    //     {0, 1, 0, 1, 1, 0, 1, 0},
    //     {0, 0, 1, 1, 1, 1, 0, 0}
    // };

//     linalg::matrix shifts(4, linalg::lin_vector(8));

//     g.make_tof(shifts);

//     encoding::encoder e = {g};
//     linalg::lin_vector vect = {0, 0, 1, 1};
//     auto encoded = e.encode(vect);

//     linalg::lin_vector encoded_check = {0, 1, 1, 0, 0, 1, 1, 0};
// // 10011001
//     // EXPECT_EQ(encoded, encoded_check);

//     std::vector<double> reliabilities = {-10, 10, 10, -10, -10, 10, -1, -10};
//     encoding::trellis_based_rml_decoder dec(g, false, false);

//     auto v = dec.decode(reliabilities);

//     EXPECT_EQ(encoded_check, v);

// }
// // 1010 10010110 1110
// // 2.02487 -0.911127 0.447729 0.00269423 -0.503658 0.66734 0.776549 0.138703 


// TEST(rmld_correctness, decode_wrong2) {
//     linalg::matrix g = {
//         {1, 1, 1, 1, 0, 0, 0, 0}, 
//         {0, 0, 0, 0, 1, 1, 1, 1}, 
//         {0, 1, 0, 1, 1, 0, 1, 0},
//         {0, 0, 1, 1, 1, 1, 0, 0}
//     };

//     linalg::matrix shifts(4, linalg::lin_vector(8));

//     g.make_tof(shifts);

//     encoding::encoder e = {g};
//     linalg::lin_vector vect = {1, 1, 1, 0};
//     auto encoded = e.encode(vect);

//     linalg::lin_vector encoded_check = {1, 0, 0, 1, 0, 1, 1, 0};
// // 10011001
//     // EXPECT_EQ(encoded, encoded_check);

//     std::vector<double> reliabilities = {2.02487, -0.911127, 0.447729, 0.00269423, -0.503658, 0.66734, 0.776549, 0.138703};
//     encoding::trellis_based_rml_decoder dec(g, false, false);

//     auto v = dec.decode(reliabilities);

//     EXPECT_EQ(encoded_check, v);

// }

// // 01101111 01101111010111001010001110010000 0110110101000010
// // -1.8851 0.976583 0.528779 -1.95404 1.15196 1.09015 0.542011 1.03569 -1.94776 0.80406 -1.46708 1.46362 1.35543 1.31522 -0.615261 -1.19715 0.83817 -0.676406 1.6807 -1.4126 -0.872296 -0.834092 2.29801 0.533701 1.80917 -0.428016 -1.29598 1.18537 -0.991321 -0.763306 -0.979013 -0.789227

// TEST(rmld_correctness, decode_wrong_large) {
//     linalg::matrix g = 
//     {
//         {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, },
//         {0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 1, 0, 1, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, },
//         {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, }
//     };

//     linalg::matrix shifts(g.size(), linalg::lin_vector(g[0].size(), false));

//     g.make_tof(shifts);

//     encoding::encoder e = {g};
//     // linalg::lin_vector vect = {0, 1, 1, 0, 1, 1, 1, 1};
//     // auto encoded = e.encode(vect);

//     linalg::lin_vector encoded_check = {0,1,1,1,0,1,1,1,0,1,1,1,1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,1,1,1};
// // 10011001
//     // EXPECT_EQ(encoded, encoded_check);

//     std::vector<double> reliabilities = {-0.576804, 1.03504, 1.80771, 0.912645, -0.500792, 0.789529, 1.22792, 0.64886, -1.31845, 0.110547, 0.855006, 1.41033, 1.3656, 0.306943, -1.40199, -0.689252, 0.823275, -0.038073, -1.7673, -1.42825, 0.258857, -0.971259, -0.772703, -0.362956, 0.585685, -0.746481, -0.52608, -0.569451, -1.01687, 2.35313, 1.52079, 1.56831};
//     encoding::trellis_based_rml_decoder dec(g, false, false);

//     auto v = dec.decode(reliabilities);

//     EXPECT_EQ(encoded_check, v);

// }


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


    g.make_tof();
    std::cout << g.to_string() << "\n";

    linalg::lin_vector encoded_check = {0,1,0,0,1,1,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,0,0,1,1,1,0,1,1,0,0,0};

    linalg::lin_vector decoded;
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }

    std::vector<double> reliabilities = {-0.944571, 1.3643, -1.21466, -1.80696, 0.887298, 0.465176, 0.838953, -1.126, -1.33395, -1.00815, -0.214494, -0.552195, -0.846559, 0.698426, 0.979281, 0.470859, 2.18831, -0.750444, 1.23306, 0.758309, -0.846314, -1.78091, 0.47336, 1.29588, 1.50673, -0.193073, -1.44247, 1.72633, 0.0470784, -1.78299, -1.44743, 0.103198};
    encoding::trellis_based_rml_decoder dec(g, false, true);

    auto v = dec.decode(reliabilities); // 01001110001001111011001010110100

    EXPECT_EQ(decoded, v);

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

    g.make_tof();

    linalg::lin_vector encoded_check = {1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 0};

    linalg::lin_vector decoded;
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }

    std::vector<double> reliabilities = {1.22055, 1.61721, 1.35983, 1.33653, -1.82543, -1.16671, 0.0244755, 1.42969, -0.985874, -1.13443, 1.03667, 0.634122, -1.09611, -0.922897, -0.0839821, -1.3403, 0.46734, -0.343023, -1.09859, 1.35925, -1.3316, 1.27189, -0.0234682, 1.42554, -0.856104, 0.39479, -1.10027, 1.32572, -1.20051, 0.211034, 0.498358, -0.762412};
    encoding::trellis_based_rml_decoder dec(g, true, true);

    auto v = dec.decode(reliabilities); // 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0
    EXPECT_EQ(decoded, v);

}

TEST(rmld_correctness, decode_wrong_large4) {
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

    g.make_tof();

    linalg::lin_vector encoded_check = {1,1,1,1,1,0,1,0,0,0,1,1,1,0,0,1,1,0,0,1,0,0,1,1,1,0,1,0,1,1,1,1};

    linalg::lin_vector decoded;
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }

    std::vector<double> reliabilities = {0.98163, 1.62571, 0.338157, 0.737762, 1.53411, -2.78978, 0.812684, -0.589047, 0.153512, -0.228798, 1.36106, 1.16975, -0.394087, -0.182766, -1.47972, 1.80939, 0.772481, -0.81109, 0.0497088, 1.98767, -0.615874, -0.424611, 0.951474, 1.99707, 0.264756, -0.975016, 0.245394, -2.11901, 0.70669, 1.31723, 0.443163, 1.16956};
    encoding::trellis_based_rml_decoder dec(g, true, true);

    auto v = dec.decode(reliabilities); // 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0
    EXPECT_EQ(decoded, v);

}

TEST(rmld_correctness, decode_wrong_large5) {
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

    g.make_tof();
    linalg::lin_vector encoded_check = {1,0,0,0,1,0,1,1,1,1,1,0,1,1,0,1,0,0,0,1,1,1,0,1,0,1,1,1,1,0,1,1};

    linalg::lin_vector decoded;
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }

    std::vector<double> reliabilities = {0.723866, -0.693951, -1.68221, -1.16534, 0.344415, -0.750087, 2.01703, 0.271928, 0.504275, 1.17766, 1.5882, -1.3428, 0.256575, 2.75327, -0.000607125, 0.615362, -2.46117, -1.53283, -1.03665, 1.97032, 1.88432, 1.51971, -1.10319, 0.516062, -1.92636, 1.97972, 0.214812, 0.540141, 1.08304, -1.84369, 0.189614, 0.109821};
    encoding::trellis_based_rml_decoder dec(g, false, false);

    auto v = dec.decode(reliabilities); // 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0
    EXPECT_EQ(decoded, v);

    encoded_check = {0,1,0,0,1,1,0,1,0,0,1,0,0,1,0,0,0,1,0,0,1,1,0,1,1,1,0,1,1,0,1,1};
    decoded.clear();
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }
    reliabilities = {-0.86696, 2.40561, -0.977012, -1.20359, 2.32433, 1.42514, -1.10551, 2.665, 0.582872 ,-1.59186 ,1.89384, -0.914534, -1.03906, 1.43699, -0.90853, -1.17326, -1.17153, 0.449967, -1.38027, -1.07219, -0.0838505, 0.923477, -0.456859, 1.34256 ,0.701526 ,1.8107, -0.27364, 0.756328 ,0.330184, 0.148568 ,0.699739, 1.97644};
    v = dec.decode(reliabilities);
    EXPECT_EQ(decoded, v);

    reliabilities = {-1.24651, -0.145936, 1.2715, 3.50491, 0.490972, -1.15913, -0.144501, 0.908666, -1.67285, -0.632916, 0.135353, -0.764073, -1.47697, 1.7497, 0.0517355, 1.30342, 1.19919, -1.19633, -0.34694, -0.216693, -1.09834, -2.0228, 1.7296, 1.92786, -1.42099, 0.237496, -1.17482, 1.79197, -1.31126, -0.820556, -0.731003, -0.731904};
    encoded_check = {0,0,1,1,1,0,0,1,0,0,0,0,0,1,0,1,1,0,0,1,0,0,1,1,0,1,0,1,0,0,0,0};
    decoded.clear();
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }
    // std::cout << decoded.to_string() << "\n";
    // std::cout << encoded_check.to_string() << "\n";
    // decoded = {0,0,1,0,1,1,1,0,0,1,1,1,0,0,0,0};
    v = dec.decode(reliabilities);
    // std::cout << (v * g).to_string() << "\n";
    // std::cout << (decoded * g).to_string() << "\n";
    // std::cout << v.to_string() << "\n";
    EXPECT_EQ(decoded, v);



    reliabilities = {1.5758, 1.83096, 1.72273, -0.762441, -2.86739, 0.607294, -1.39952, 2.31931, -2.65765, -0.377633, -0.902422, -1.63359, -0.667267, 1.43959, -0.569823, -1.41433, 0.157579, 1.3167, -1.67171, 2.49687, -0.596566, -1.93274, 1.5218, -1.25642, 1.2741, 0.665177, -2.16329, -0.184256, 2.42067, -1.80205, -0.483804, -1.20259};
    encoded_check = {1,1,1,0,0,0,0,1,0,1,0,0,0,1,0,0,1,1,0,1,0,0,1,0,1,0,0,0,1,0,0,0};
    decoded.clear();
    for (auto const& vect : g) {
        if (vect.leading() == encoded_check.leading()) {
            decoded.push_back(true);
            encoded_check += vect;
        } else {
            decoded.push_back(false);
        }
    }
    // decoded = {1,0,0,1,1,0,0,1,1,0,1,1,0,1,0,0};
    v = dec.decode(reliabilities);
    EXPECT_EQ(decoded, v);
}
// 0101110011110001 01001101001001000100110111011011
// -0.86696 2.40561 -0.977012 -1.20359 2.32433 1.42514 -1.10551 2.665 0.582872 -1.59186 1.89384 -0.914534 -1.03906 1.43699 -0.90853 -1.17326 -1.17153 0.449967 -1.38027 -1.07219 -0.0838505 0.923477 -0.456859 1.34256 0.701526 1.8107 -0.27364 0.756328 0.330184 0.148568 0.699739 1.97644

// 0010111001110000 00111001000001011001001101010000 0010110000100101
// -1.24651, -0.145936, 1.2715, 3.50491, 0.490972, -1.15913, -0.144501, 0.908666, -1.67285, -0.632916, 0.135353, -0.764073, -1.47697, 1.7497, 0.0517355, 1.30342, 1.19919, -1.19633, -0.34694, -0.216693, -1.09834, -2.0228, 1.7296, 1.92786, -1.42099, 0.237496, -1.17482, 1.79197, -1.31126, -0.820556, -0.731003, -0.731904


// 1001100110110100 11100001010001001101001010001000 1010110011001000 11000101010111000101001111001010
// 1.5758, 1.83096, 1.72273, -0.762441, -2.86739, 0.607294, -1.39952, 2.31931, -2.65765, -0.377633, -0.902422, -1.63359, -0.667267, 1.43959, -0.569823, -1.41433, 0.157579, 1.3167, -1.67171, 2.49687, -0.596566, -1.93274, 1.5218, -1.25642, 1.2741, 0.665177, -2.16329, -0.184256, 2.42067, -1.80205, -0.483804, -1.20259

TEST(gray_codes_correctness, simple_check) {
    linalg::matrix g = {
        {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0},
        {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1} 
    };
    std::vector<double> weights = {1, 1, 1, 1, 1, 1, 1, 1};
    encoding::trellis_based_rml_decoder dec(g, false, false);
    dec.build_gray_codes(9);
    for (size_t i = 0; i < dec._gray_codes[8].size() - 1; ++i) {
        EXPECT_TRUE(encoding::hamming_metric(dec._gray_codes[8][i], weights).count(dec._gray_codes[8][i + 1]) <= 1.5);
    }

}

TEST(trellis_correctness, rm2_2) {
    linalg::matrix g = 
    {
        {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, 
        {0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0},
        {0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0},
        {0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0},
        {0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0},
        {0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0},
        {0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1} 
    };


    encoding::trellis_based_rml_decoder dec(g, false, true);

    std::cout << g.to_string() << "\n";

    // dec._special_matrices[0][7] = g;
    // dec._special_matrices[4][7] = gg;
    // dec.build_special_trellis(0, 4, 8);

    // dec.make_uniform_decomposition(0, 4, 8);

    // std::cout << " printing special trellis: \n";
    // auto tr = dec._trellises[0][7];
    // for (size_t zz = 0; zz < tr._sections[0].size(); ++zz) {
    //     std::cout << zz << " transmissions:\n";
    //     for (auto const& br : tr._sections[0][zz]._next) {
    //         std::cout << "\t" << br.first.to_string() << " " << br.second << "\n";
    //     }
    // }
    // std::cout << "\n";

    // for (auto const& pc : tr._parallel_components) {
    //     std::cout << "left parallel component: \n";
    //     for (auto c : pc[0]) {
    //         std::cout << c << " ";
    //     }
    //     std::cout << "\n";
    //     std::cout << "right parallel component: \n";
    //     for (auto c : pc[1]) {
    //         std::cout << c << " ";
    //     }
    //     std::cout << "\n";
    // }

    // for (auto const& pc : tr._parallel_component_groups) {
    //     std::cout << "parallel component group: \n";
    //     for (auto c : pc) {
    //         std::cout << c << " ";
    //     }
    //     std::cout << "\n";
    // }

    // for (auto const& pc : tr._groups) {
    //     std::cout << "parallel component groups \n";
    //     for (auto const &g : pc[0]) {
    //         std::cout << "left group: \n";
    //         for (auto v : g) {
    //             std::cout << v << " ";
    //         }
    //         std::cout << "\n";
    //     }
    //     for (auto const &g : pc[1]) {
    //         std::cout << "right group: \n";
    //         for (auto v : g) {
    //             std::cout << v << " ";
    //         }
    //         std::cout << "\n";
    //     }
    //     std::cout << "\n";
    // }
}

