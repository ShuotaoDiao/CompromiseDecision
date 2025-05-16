//
//  DataStructure.cpp
//  CompromiseDecision_v1.0
//
//  Created by Shuotao Diao on 4/17/25.
//

#include "DataStructure.hpp"

double operator*(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    double res = 0;
    if (vec1.size() != vec2.size()) {
        std::cout << "Warning: vec1.size() = " << vec1.size() << std::endl;
        std::cout << "Warning: vec2.size() = " << vec2.size() << std::endl;
        throw std::invalid_argument("Vector sizes do not match.\n");
    }
    for (int index = 0; index < vec1.size(); ++index) {
        res += vec1[index] * vec2[index];
    }
    return res;
}

std::vector<double> operator*(double a, const std::vector<double>& vec1) {
    std::vector<double> res;
    for (int index = 0; index < vec1.size(); ++index) {
        res.push_back(a * vec1[index]);
    }
    return res;
}

std::vector<double> operator+(const std::vector<double>& vec1, const std::vector<double>& vec2) {
    if (vec1.size() != vec2.size()) {
        throw std::invalid_argument("Vector sizes are not matched.\n");
    }
    std::vector<double> res;
    for (int index = 0; index < vec1.size(); ++index) {
        res.push_back(vec1[index] + vec2[index]);
    }
    return res;
}

double max(double a, double b) {
    if (b > a) {
        return b;
    }
    else {
        return a;
    }
}

double min(double a, double b) {
    if (a < b) {
        return a;
    }
    else {
        return b;
    }
}

// Functions for testing
// sparse matrix
void SparseMatrix_demo_test1() {
    // test on the sparse matrix insert
    // case 1
    // M1 = | 1 3 0 |
    //      | 2 0 4 |
    std::cout << "Case 1" << std::endl;
    std::cout << "M1 = |1 3 0|" << std::endl;
    std::cout << "     |2 0 4|" << std::endl;
    std::cout << "Generating sparse matrix..." << std::endl;
    SparseMatrix M1;
    M1.insert(0, 0, 1);
    M1.insert(1, 0, 2);
    M1.insert(0, 1, 3);
    M1.insert(1, 2, 4);
    M1.insert_end(2, 3);
    std::vector<std::vector<double>> M1_dense = M1.scatter();
    std::cout << "Converting sparse matrix into a dense matrix" << std::endl;
    for (int row_idx = 0; row_idx < 2; ++row_idx) {
        for (int col_idx = 0; col_idx < 3; ++col_idx) {
            std::cout << M1_dense[row_idx][col_idx] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "***************" << std::endl;
    // case 2
    // M2 = |1 0 0|
    //      |2 0 4|
    std::cout << "Case 2" << std::endl;
    std::cout << "M2 = |1 0 0|" << std::endl;
    std::cout << "     |2 0 4|" << std::endl;
    std::cout << "Generating sparse matrix..." << std::endl;
    SparseMatrix M2;
    M2.insert(0, 0, 1);
    M2.insert(1, 0, 2);
    M2.insert(1, 2, 4);
    M2.insert_end(2, 3);
    std::cout << "Converting sparse matrix into a dense matrix" << std::endl;
    std::vector<std::vector<double>> M2_dense = M2.scatter();
    for (int row_idx = 0; row_idx < 2; ++row_idx) {
        for (int col_idx = 0; col_idx < 3; ++col_idx) {
            std::cout << M2_dense[row_idx][col_idx] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "***************" << std::endl;
    // case 3
    // M3 = |0 0 0|
    //      |0 0 4|
    std::cout << "Case 3" << std::endl;
    std::cout << "M3 = |0 0 0|" << std::endl;
    std::cout << "     |0 0 4|" << std::endl;
    std::cout << "Generating sparse matrix..." << std::endl;
    SparseMatrix M3;
    M3.insert(1, 2, 4);
    M3.insert_end(2, 3);
    std::cout << "Converting sparse matrix into a dense matrix" << std::endl;
    std::vector<std::vector<double>> M3_dense = M3.scatter();
    for (int row_idx = 0; row_idx < 2; ++row_idx) {
        for (int col_idx = 0; col_idx < 3; ++col_idx) {
            std::cout << M3_dense[row_idx][col_idx] << " ";
        }
        std::cout << std::endl;
    }
}
