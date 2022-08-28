//
//  main.cpp
//  Matrix
//
//  Created by Sami Hatna on 28/08/2022.
//

#include <cmath>
#include <iostream>
#include <vector>

#include "Matrix.cpp"
#include "Matrix.hpp"

int main() {
    // Initialise matrices
    // Matrix class is a template so can be used with int, double and float
    // types
    Matrix<double> testMat1{{{8, 2, 3}, {1, 3, 5}, {6, 8, 1}}};
    Matrix<double> testMat2{{{-2, -4, 2}, {-2, 1, 2}, {4, 2, 5}}};

    // Output matrix
    testMat1.printMatrix();

    // Get dimensions
    std::cout << std::endl
              << testMat1.numRows() << " x " << testMat1.numCols() << std::endl;

    // Get value at position
    std::cout << std::endl
              << "Value at pos 0x0: " << testMat1.getVal(0, 0) << std::endl;

    // Set value at position
    testMat1.setVal(0, 0, 5);
    std::cout << std::endl
              << "Value now at pos 0x0: " << testMat1.getVal(0, 0) << std::endl;

    // Perform elementary row operations
    testMat1.swapRow(0, 1);
    testMat1.multiplyRow(1, 5);
    testMat1.addScalarMultiple(0, 1, 4);
    std::cout << std::endl << "Apply elementary row operations:" << std::endl;
    testMat1.printMatrix();

    // Matrix multiplication
    auto multiply = testMat1 * testMat2;
    std::cout << std::endl << "Multiply Matrix 1 by Matrix 2:" << std::endl;
    multiply.printMatrix();

    // Matrix addition
    testMat1 += testMat2;
    std::cout << std::endl << "Add Matrix 2 to Matrix 1:" << std::endl;
    testMat1.printMatrix();

    // Transpose matrix
    std::cout << std::endl << "Transpose matrix:" << std::endl;
    auto transposed = testMat1.transpose();
    transposed.printMatrix();

    // Transform matrix to row echelon form
    std::cout << std::endl << "Carry to REF:" << std::endl;
    testMat1.ref();
    testMat1.printMatrix();

    // Get matrix determinant
    std::cout << std::endl
              << "Determinant: " << testMat2.determinant() << std::endl;

    // Find cofactor matrix
    std::cout << std::endl << "Cofactor matrix:" << std::endl;
    auto cofactor = testMat2.cofactor();
    cofactor.printMatrix();

    // Invert matrix
    std::cout << std::endl << "Invert matrix:" << std::endl;
    auto inverse = testMat2.inverse();
    inverse.printMatrix();

    // Orthogonalise matrix using gram-schmidt
    std::cout << std::endl << "Orthogonalise matrix:" << std::endl;
    auto orthogonal = testMat2.gramSchmidt(false);
    orthogonal.printMatrix();

    // Orthonormalize matrix
    std::cout << std::endl << "Orthonormalize matrix:" << std::endl;
    auto orthonormal = testMat2.gramSchmidt(true);
    orthonormal.printMatrix();

    // Perform QR factorization
    std::vector<Matrix<double>> qr = testMat2.QRFactorization();

    // Return vector of eigenvalues
    std::cout << std::endl << "Find eigenvalues: ";
    std::vector<double> eigenvalues = testMat2.eigenvalues();
    for (auto i : eigenvalues) {
        std::cout << i << ", ";
    }
    std::cout << std::endl << std::endl;

    return 0;
}
