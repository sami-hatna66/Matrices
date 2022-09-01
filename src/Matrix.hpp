//
//  Matrix.hpp
//  Matrix
//
//  Created by Sami Hatna on 28/08/2022.
//

#ifndef Matrix_hpp
#define Matrix_hpp

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <vector>

template <typename T> class Matrix {
  private:
    std::vector<std::vector<T>> content;

  public:
    Matrix(std::vector<std::vector<T>> c);
    void printMatrix();
    Matrix<T> &operator=(Matrix<T> &newMat);
    Matrix<T> &operator*=(Matrix<T> rhs);
    Matrix<T> &operator*=(double scalar);
    Matrix<T> &operator+=(Matrix<T> rhs);
    Matrix<T> &operator-=(Matrix<T> rhs);
    unsigned long numRows();
    unsigned long numCols();
    T getVal(int row, int col);
    std::vector<T> getCol(int col);
    void setVal(int row, int col, T newVal);
    void setContent(std::vector<std::vector<T>> newContent);
    void swapRow(int row1, int row2);
    void multiplyRow(int row, T scalar);
    void addScalarMultiple(int rowTarget, int rowAdding, T scalar);
    Matrix transpose();
    void ref();
    double determinant();
    Matrix cofactor();
    Matrix inverse();
    Matrix gramSchmidt(bool isOrthonormal);
    std::vector<T> projection(std::vector<T> vec, std::vector<T> dir);
    std::vector<Matrix> QRFactorization();
    T dotProduct(std::vector<T> vec1, std::vector<T> vec2);
    std::vector<T> eigenvalues();
};

#endif /* Matrix_hpp */
