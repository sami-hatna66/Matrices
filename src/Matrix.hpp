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

// Declarations
template <typename T> class Matrix {
  private:
    std::vector<std::vector<T>> content;

  public:
    Matrix(std::vector<std::vector<T>> c);
    void printMatrix();
    Matrix<T> &operator=(Matrix<T> &newMat);
    Matrix<T> &operator*=(Matrix<T> rhs);
    Matrix<T> &operator*=(T scalar);
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
    T determinant();
    Matrix cofactor();
    Matrix inverse();
    Matrix gramSchmidt(bool isOrthonormal);
    std::vector<T> projection(std::vector<T> vec, std::vector<T> dir);
    std::vector<Matrix> QRFactorization();
    T dotProduct(std::vector<T> vec1, std::vector<T> vec2);
    std::vector<T> eigenvalues();
};

template <typename T> Matrix<T> operator*(Matrix<T> lhs, Matrix<T> rhs);
template <typename T> Matrix<T> operator*(Matrix<T> lhs, T scalar);
template <typename T> Matrix<T> operator*(T scalar, Matrix<T> rhs);
template <typename T> Matrix<T> operator+(Matrix<T> lhs, Matrix<T> rhs);
template <typename T> Matrix<T> operator-(Matrix<T> lhs, Matrix<T> rhs);

// Implementations
template <typename T> Matrix<T>::Matrix(std::vector<std::vector<T>> c) {
    // Check if valid matrix
    try {
        if (c.size() == 0) {
            throw 505;
        } else {
            bool flag = true;
            unsigned long prevLength = c[0].size();
            for (int i = 1; i < c.size(); i++) {
                if (prevLength != c[i].size()) {
                    throw 505;
                }
            }
            content = c;
        }
    } catch (...) {
        std::cout << "Input error" << std::endl;
    }
}

// Pretty print matrix
template <typename T> void Matrix<T>::printMatrix() {
    int maxLen = 0;
    int row = 0;
    std::for_each(content.begin(), content.end(), [&](const std::vector<T>& col) {
        row++;
        std::for_each(col.begin(), col.end(), [&](T val) {
            unsigned long len = std::to_string(val).length();
            if (len > maxLen) {
                maxLen = len;
            }
        });
    });

    row = 0;
    std::for_each(content.begin(), content.end(), [&](const std::vector<T>& col) {
        row++;
        std::cout << "| ";
        std::for_each(col.begin(), col.end(), [&](T val) {
            std::cout << std::setw(maxLen) << std::to_string(val) << " ";
        });
        std::cout << "|" << std::endl;
    }); 
}

template <typename T> unsigned long Matrix<T>::numRows() {
    return content.size();
}

template <typename T> unsigned long Matrix<T>::numCols() {
    return content[0].size();
}

// Matrix indexing begins at 0
template <typename T> T Matrix<T>::getVal(int row, int col) {
    try {
        if (row > this->numRows() - 1) {
            throw(std::string("row"));
        } else if (col > this->numCols() - 1) {
            throw(std::string("col"));
        } else {
            return content[row][col];
        }
    } catch (std::string dim) {
        std::cout << dim << " out of range" << std::endl;
        return 0;
    }
}

// Select individual columns from a matrix (helper function for QRFactorization)
template <typename T> std::vector<T> Matrix<T>::getCol(int col) {
    try {
        if (col > this->numCols() - 1) {
            throw 505;
        } else {
            std::vector<T> column;
            for (int i = 0; i < this->numRows(); i++) {
                column.push_back(this->getVal(i, col));
            }
            return column;
        }
    } catch (...) {
        std::cout << "Column out of range";
        return {0};
    }
}

template <typename T> void Matrix<T>::setVal(int row, int col, T newVal) {
    content[row][col] = newVal;
}

template <typename T>
void Matrix<T>::setContent(std::vector<std::vector<T>> newContent) {
    content = newContent;
}

// Elementary row operations
template <typename T> void Matrix<T>::swapRow(int row1, int row2) {
    try {
        if (row1 < this->numRows() && row2 < this->numRows()) {
            swap(content[row2], content[row1]);
        } else {
            throw 505;
        }
    } catch (...) {
        std::cout << "Row index out of range" << std::endl;
    }
}

template <typename T> void Matrix<T>::multiplyRow(int row, T scalar) {
    try {
        if (row < this->numRows()) {
            for (int i = 0; i < this->numCols(); i++) {
                content[row][i] *= scalar;
            }
        } else {
            throw 505;
        }
    } catch (...) {
        std::cout << "Row index out of range" << std::endl;
    }
}

template <typename T>
void Matrix<T>::addScalarMultiple(int rowTarget, int rowAdding, T scalar) {
    try {
        if (rowAdding > this->numRows() || rowTarget > this->numRows()) {
            throw 505;
        } else {
            for (int i = 0; i < this->numCols(); i++) {
                content[rowTarget][i] += (this->getVal(rowAdding, i) * scalar);
            }
        }
    } catch (...) {
        std::cout << "Invalid dimensions" << std::endl;
    }
}

// Transpose matrix
template <typename T> Matrix<T> Matrix<T>::transpose() {
    std::vector<std::vector<T>> newContent = {};
    for (int i = 0; i < this->numCols(); i++) {
        std::vector<T> intermediary = {};
        for (int j = 0; j < this->numRows(); j++) {
            intermediary.push_back({content[j][i]});
        }
        newContent.push_back(intermediary);
    }
    Matrix<T> newMat = Matrix(newContent);
    return newMat;
}

// Convert matrix to row echelon form (unfinished)
template <typename T> void Matrix<T>::ref() {
    int col = 0;
    for (int row = 0; row < std::min(this->numRows(), this->numCols()); row++) {
        this->multiplyRow(row, 1 / this->getVal(row, col));

        for (int row2 = row + 1; row2 < this->numRows(); row2++) {
            this->addScalarMultiple(row2, row, -this->getVal(row2, col));
        }

        col++;
    }
}

// Calculate determinant using cofactor expansion
template <typename T> T Matrix<T>::determinant() {
    T det;

    try {
        if (this->numCols() != this->numRows()) {
            throw 505;
        }

        else {
            switch (this->numRows()) {
            case 1:
                det = this->getVal(0, 0);
                break;
            case 2:
                det = (this->getVal(0, 0) * this->getVal(1, 1)) -
                      (this->getVal(0, 1) * this->getVal(1, 0));
                break;
            case 3:
                det = (this->getVal(0, 0) * this->getVal(1, 1) *
                       this->getVal(2, 2)) +
                      (this->getVal(0, 1) * this->getVal(1, 2) *
                       this->getVal(2, 0)) +
                      (this->getVal(0, 2) * this->getVal(1, 0) *
                       this->getVal(2, 1)) -
                      (this->getVal(0, 2) * this->getVal(1, 1) *
                       this->getVal(2, 0)) -
                      (this->getVal(0, 0) * this->getVal(1, 2) *
                       this->getVal(2, 1)) -
                      (this->getVal(0, 1) * this->getVal(1, 0) *
                       this->getVal(2, 2));
                break;
            default:
                det = 0;
                for (int i = 0; i < this->numCols(); i++) {
                    std::vector<std::vector<T>> newMatContent;
                    for (int j = 1; j < this->numCols(); j++) {
                        std::vector<T> newRow;
                        for (int k = 0; k < this->numCols(); k++) {
                            if (k != i) {
                                newRow.push_back(this->getVal(j, k));
                            }
                        }
                        newMatContent.push_back(newRow);
                    }
                    Matrix<T> newMat = Matrix(newMatContent);
                    det += (this->getVal(0, i)) * (pow(-1, 2 + i)) *
                           newMat.determinant();
                }
            }
        }
    } catch (...) {
        std::cout << "Matrix must be square";
        return 0;
    }

    return det;
}

// Return cofactor matrix
template <typename T> Matrix<T> Matrix<T>::cofactor() {
    std::vector<std::vector<T>> cofactorContent;

    for (int i = 0; i < this->numRows(); i++) {
        std::vector<T> cofactorRow;
        for (int j = 0; j < this->numCols(); j++) {
            std::vector<std::vector<T>> prunedContent = {};
            for (int i2 = 0; i2 < this->numRows(); i2++) {
                std::vector<T> prunedRow;
                for (int j2 = 0; j2 < this->numCols(); j2++) {
                    if (i2 != i && j2 != j) {
                        prunedRow.push_back(this->getVal(i2, j2));
                    }
                }
                if (prunedRow.size() > 0) {
                    prunedContent.push_back(prunedRow);
                }
            }
            Matrix<T> prunedMat = Matrix(prunedContent);
            cofactorRow.push_back(pow(-1, i + j) * prunedMat.determinant());
        }
        cofactorContent.push_back(cofactorRow);
    }

    return Matrix<T>(cofactorContent);
}

// Return inverse of matrix
template <typename T> Matrix<T> Matrix<T>::inverse() {
    try {
        if (this->determinant() == 0) {
            throw 505;
        } else {
            Matrix adjugate = (this->cofactor()).transpose();
            return ((T)1 / this->determinant()) * adjugate;
        }
    } catch (...) {
        std::cout << "Matrix is not invertible";
        return Matrix(this->content);
    }
}

// Assignment
template <typename T> Matrix<T> &Matrix<T>::operator=(Matrix<T> &newMat) {
    this->setContent(newMat.content);
    return *this;
}

// Matrix-Matrix multiplication
template <typename T> Matrix<T> &Matrix<T>::operator*=(Matrix<T> rhs) {
    try {
        if (this->numCols() != rhs.numRows()) {
            throw 505;
        } else {
            std::vector<std::vector<T>> newContent = {};

            for (int i = 0; i < this->numRows(); i++) {
                std::vector<T> newRow;
                for (int j = 0; j < rhs.numCols(); j++) {
                    T sum = 0;
                    for (int k = 0; k < this->numCols(); k++) {
                        sum += (this->getVal(i, k) * rhs.getVal(k, j));
                    }
                    newRow.push_back(sum);
                }
                newContent.push_back(newRow);
            }

            this->setContent(newContent);
        }
    } catch (...) {
        std::cout << "Invalid dimensions" << std::endl;
    }

    return *this;
}

template <typename T> Matrix<T> operator*(Matrix<T> lhs, Matrix<T> rhs) {
    return lhs *= rhs;
}

// Matrix-Scalar multiplication
template <typename T> Matrix<T> &Matrix<T>::operator*=(T scalar) {
    for (int row = 0; row < this->numRows(); row++) {
        for (int col = 0; col < this->numCols(); col++) {
            this->setVal(row, col, scalar * this->getVal(row, col));
        }
    }
    return *this;
}

template <typename T> Matrix<T> operator*(Matrix<T> lhs, T scalar) {
    return lhs *= scalar;
}
// For commutativity
template <typename T> Matrix<T> operator*(T scalar, Matrix<T> rhs) {
    return rhs *= scalar;
}

// Addition
template <typename T> Matrix<T> &Matrix<T>::operator+=(Matrix<T> rhs) {
    try {
        if (this->numCols() != rhs.numCols() ||
            this->numRows() != rhs.numRows()) {
            throw 505;
        } else {
            for (int row = 0; row < this->numRows(); row++) {
                for (int col = 0; col < this->numCols(); col++) {
                    this->setVal(row, col,
                                 this->getVal(row, col) + rhs.getVal(row, col));
                }
            }
        }
    } catch (...) {
        std::cout << "Invalid dimensions" << std::endl;
    }

    return *this;
}

template <typename T> Matrix<T> operator+(Matrix<T> lhs, Matrix<T> rhs) {
    return lhs += rhs;
}

// Subtraction
template <typename T> Matrix<T> &Matrix<T>::operator-=(Matrix<T> rhs) {
    try {
        if (this->numCols() != rhs.numCols() ||
            this->numRows() != rhs.numRows()) {
            throw 505;
        } else {
            for (int row = 0; row < this->numRows(); row++) {
                for (int col = 0; col < this->numCols(); col++) {
                    this->setVal(row, col,
                                 this->getVal(row, col) - rhs.getVal(row, col));
                }
            }
        }
    } catch (...) {
        std::cout << "Invalid dimensions" << std::endl;
    }

    return *this;
}

template <typename T> Matrix<T> operator-(Matrix<T> lhs, Matrix<T> rhs) {
    return lhs -= rhs;
}

// Performs Gram-Schmidt process on columns of matrix
// Result can be made orthonormal by setting parameter to true
template <typename T> Matrix<T> Matrix<T>::gramSchmidt(bool isOrthonormal) {
    // Transpose matrix to get columns as individual vectors
    Matrix transposed = this->transpose();

    std::vector<std::vector<T>> qContent;
    if (isOrthonormal) {
        std::vector<T> result = transposed.content[0];
        T squareTotal = 0;
        for (int l = 0; l < result.size(); l++) {
            squareTotal += result[l] * result[l];
        }
        T magnitude = sqrt(squareTotal);
        for (int m = 0; m < result.size(); m++) {
            result[m] /= magnitude;
        }
        qContent.push_back(result);
    } else {
        qContent.push_back(transposed.content[0]);
    }

    for (int i = 1; i < this->numCols(); i++) {
        std::vector<T> result = transposed.content[i];
        for (int j = 0; j < i; j++) {
            std::vector<T> proj =
                projection(transposed.content[i], qContent[j]);
            for (int k = 0; k < proj.size(); k++) {
                result[k] -= proj[k];
            }
        }
        if (isOrthonormal) {
            T squareTotal = 0;
            for (int l = 0; l < result.size(); l++) {
                squareTotal += result[l] * result[l];
            }
            T magnitude = sqrt(squareTotal);
            for (int m = 0; m < result.size(); m++) {
                result[m] /= magnitude;
            }
        }

        qContent.push_back(result);
    }

    Matrix<T> Q = Matrix(qContent);
    return Q.transpose();
}

// Calculates projection of vec in direction of dir (helper function for
// gramSchmidt)
template <typename T>
std::vector<T> Matrix<T>::projection(std::vector<T> vec, std::vector<T> dir) {
    T numerator = 0;
    T denominator = 0;
    std::vector<T> result;

    for (int i = 0; i < vec.size(); i++) {
        numerator += (vec[i] * dir[i]);
        denominator += (dir[i] * dir[i]);
    }

    for (int j = 0; j < vec.size(); j++) {
        result.push_back(dir[j] * (numerator / denominator));
    }

    return result;
}

// Calculates the dot product of two vectors (helper function)
template <typename T>
T Matrix<T>::dotProduct(std::vector<T> vec1, std::vector<T> vec2) {
    T result = 0;
    for (int i = 0; i < vec1.size(); i++) {
        result += (vec1[i] * vec2[i]);
    }
    return result;
}

// Performs QR factorization on matrix
// Returns a vector of two matrices containing Q and R
template <typename T> std::vector<Matrix<T>> Matrix<T>::QRFactorization() {
    Matrix<T> Q = this->gramSchmidt(true);

    Matrix<T> notOrthonormal = this->gramSchmidt(false);

    std::vector<std::vector<T>> rContent = {};
    for (int i = 0; i < this->numCols(); i++) {
        std::vector<T> rRow;
        for (int j = 0; j < this->numCols(); j++) {
            if (j < i) {
                rRow.push_back(0);
            } else if (i == j) {
                std::vector<T> magCol = notOrthonormal.getCol(i);
                T result = 0;
                for (int k = 0; k < magCol.size(); k++) {
                    result += magCol[k] * magCol[k];
                }
                result = sqrt(result);
                rRow.push_back(result);
            } else {
                rRow.push_back(dotProduct(this->getCol(j), Q.getCol(i)));
            }
        }
        rContent.push_back(rRow);
    }

    Matrix<T> R = Matrix(rContent);

    std::vector<Matrix> result = {Q, R};
    return result;
}

// Find eigenvalues using QR algorithm
template <typename T> std::vector<T> Matrix<T>::eigenvalues() {
    std::vector<T> result;

    Matrix<T> A = Matrix(this->content);

    for (int i = 0; i < 10; i++) {
        std::vector<Matrix> QR = A.QRFactorization();
        Matrix<T> product = QR[1] * QR[0];
        A = product;
    }

    for (int j = 0; j < A.numCols(); j++) {
        result.push_back(round(A.getVal(j, j)));
    }

    return result;
}

#endif /* Matrix_hpp */
