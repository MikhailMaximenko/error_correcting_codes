#include "linalg.h"

#include <cstddef>
#include <gtest/gtest.h>
#include <vector>

TEST(linalg_correctness, lin_vector_permutation) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 1, 0, 0};

    std::vector<size_t> permute = {0, 1, 2, 3};
    vect2.permutate(permute);
    EXPECT_EQ(vect1, vect2);

    linalg::lin_vector vect3 = {0, 1, 1, 0};
    permute = {2, 1, 0, 3};
    vect1.permutate(permute);

    EXPECT_EQ(vect1, vect3);
}

TEST(linalg_correctness, scalar_product) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 1, 0, 0};
    EXPECT_EQ(vect1 * vect2, 2);
    vect2 = {1, 0, 1, 0};
    EXPECT_EQ(vect1 * vect2, 1);
}

TEST(linalg_correctness, vector_add) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 1, 0, 0};
    linalg::lin_vector vect3 = {0, 0, 0, 0};
    EXPECT_EQ(vect1 + vect2, vect3);
    vect2 = {1, 0, 1, 0};
    vect3 = {0, 1, 1, 0};
    EXPECT_EQ(vect1 + vect2, vect3);
}

TEST(linalg_correctness, vector_matrix_product) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 0, 1, 0};
    linalg::matrix mt = {vect1, vect2};

    linalg::lin_vector vect3 = {1, 1};
    linalg::lin_vector vect4 = {0, 1, 1, 0};
    EXPECT_EQ(vect3 * mt, vect4);
    
    vect3 = {1, 0};
    vect4 = {1, 1, 0, 0};
    EXPECT_EQ(vect3 * mt, vect4);
}


TEST(linalg_correctness, matrix_vector_product) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 0, 1, 0};
    linalg::matrix mt = {vect1, vect2};

    linalg::lin_vector vect3 = {1, 1, 0, 0};
    linalg::lin_vector vect4 = {0, 1};
    EXPECT_EQ(mt * vect3, vect4);
    
    vect3 = {1, 0, 0, 1};
    vect4 = {1, 1};
    EXPECT_EQ(mt * vect3, vect4);
}


TEST(linalg_correctness, matrix_permutations) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 0, 1, 0};
    linalg::matrix mt = {vect1, vect2};

    std::vector<size_t> p = {2, 0, 3, 1};
    mt.permutate(p);

    linalg::lin_vector vect3 = {0, 1, 0, 1};
    linalg::lin_vector vect4 = {1, 1, 0, 0};
    linalg::matrix expected = {vect3, vect4};

    EXPECT_EQ(mt, expected);
}


TEST(linalg_correctness, resolve_basis1) {
    linalg::lin_vector vect1 = {1, 1, 0, 0};
    linalg::lin_vector vect2 = {1, 0, 1, 0};
    linalg::matrix mt = {vect1, vect2};
    auto check = mt;

    std::vector<size_t> pos = {0, 1};
    
    auto basis = mt.resolve_basis(); 

    linalg::matrix transform = {
        {1, 1},
        {1, 0}
    };

    EXPECT_EQ(basis.first, pos);
    EXPECT_EQ(basis.second, transform);


    vect1 = {1, 0, 1, 0};
    vect2 = {0, 1, 1, 0};
    linalg::matrix expected = {vect1, vect2};

    EXPECT_EQ(mt, expected);
    EXPECT_EQ(transform.inverse() * check, mt);
    
}

TEST(linalg_correctness, resolve_basis2) {
    linalg::lin_vector vect1 = {1, 1, 1, 0};
    linalg::lin_vector vect2 = {1, 0, 0, 0};
    linalg::lin_vector vect3 = {1, 1, 1, 1};
    linalg::matrix mt = {vect1, vect2, vect3};
    auto check = mt;
    std::vector<size_t> pos = {0, 1, 3};
    
    auto basis = mt.resolve_basis();

    EXPECT_EQ(basis.first, pos);
    EXPECT_EQ(basis.second.inverse() * check, mt);

    vect1 = {1, 0, 0, 0};
    vect2 = {0, 1, 1, 0};
    vect3 = {0, 0, 0, 1};
    linalg::matrix expected = {vect1, vect2, vect3};

    EXPECT_EQ(mt, expected);    
}


TEST(linalg_correctness, resolve_basis3_with_8_4) {
    linalg::matrix gen = {
        {1, 1, 1, 1, 0, 0, 0, 0}, 
        {0, 0, 0, 0, 1, 1, 1, 1}, 
        {0, 1, 0, 1, 1, 0, 1, 0},
        {0, 0, 1, 1, 1, 1, 0, 0}
    };

    auto back_transform_check = gen;

    linalg::matrix check = {
        {1, 0, 0, 1, 0, 1, 1, 0}, 
        {0, 1, 0, 1, 0, 1, 0, 1},
        {0, 0, 1, 1, 0, 0, 1, 1},
        {0, 0, 0, 0, 1, 1, 1, 1} 
    };


    linalg::matrix transformation = {
        {1, 1, 1, 0}, 
        {0, 0, 0, 1},
        {0, 1, 0, 1},
        {0, 0, 1, 1} 
    };


    std::vector<size_t> pos = {0, 1, 2, 4};


    auto basis = gen.resolve_basis();

    EXPECT_EQ(pos, basis.first);
    EXPECT_EQ(transformation, basis.second);

    EXPECT_EQ(check, gen);
    EXPECT_EQ(transformation.inverse() * back_transform_check, check);
}


