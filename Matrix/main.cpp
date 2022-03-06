//
//  main.cpp
//  temp
//
//  Created by Sami Hatna on 18/12/2021.
//

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

template<typename T>
class Matrix {
private:
    vector<vector<T>> content;
public:
    Matrix(vector<vector<T>> c) {
        // Check if valid matrix
        try {
            if (c.size() == 0) {
                throw 505;
            }
            else {
                bool flag = true;
                unsigned long prevLength = c[0].size();
                for (int i = 1; i < c.size(); i++) {
                    if (prevLength != c[i].size()) {
                        flag = false;
                        break;
                    }
                }
                if (!flag) {
                    throw 505;
                }
                else {
                    content = c;
                }
            }
        } catch (...) {
            cout << "Input error" << endl;
        }
    }
    
    // Use if operator needs to access private attributes
    // friend Matrix& operator*=(Matrix& lhs, Matrix rhs);
    
    // Function declarations
    void printMatrix();
    Matrix<T>& operator=(Matrix<T>& newMat);
    unsigned long numRows();
    unsigned long numCols();
    T getVal(int row, int col);
    vector<T> getCol(int col);
    void setVal(int row, int col, T newVal);
    void setContent(vector<vector<T>> newContent);
    void swapRow(int row1, int row2);
    void multiplyRow(int row, T scalar);
    void addScalarMultiple(int rowTarget, int rowAdding, T scalar);
    Matrix transpose();
    void ref();
    double determinant();
    Matrix cofactor();
    Matrix inverse();
    Matrix gramSchmidt(bool isOrthonormal);
    vector<T> projection(vector<T> vec, vector<T> dir);
    vector<Matrix> QRFactorization();
    T dotProduct(vector<T> vec1, vector<T> vec2);
    vector<T> eigenvalues();
    
};

// Pretty print matrix
template<typename T>
void Matrix<T>::printMatrix() {
    unsigned long maxLen = 0;
    for (int i = 0; i < content.size(); i++) {
        for (int j = 0; j < content[i].size(); j++) {
            unsigned long len = to_string(content[i][j]).length();
            if (len > maxLen) {
                maxLen = len;
            }
        }
    }
    
    for (int i = 0; i < content.size(); i++) {
        string forPrint = "| ";
        for_each(content[i].begin(),
                 content[i].end(),
                 [&forPrint, &maxLen, *this](T &n){
                    for (int j = 0; j < (maxLen - to_string(n).length()); j++) {
                        forPrint += " ";
                    }
                    forPrint += to_string(n) + " ";
        });
        forPrint += "|";
        cout << forPrint << endl;
    }
}

template<typename T>
unsigned long Matrix<T>::numRows() {
    return content.size();
}

template<typename T>
unsigned long Matrix<T>::numCols() {
    return content[0].size();
}

// Matrix indexing begins at 0
template<typename T>
T Matrix<T>::getVal(int row, int col) {
    try {
        if (row > this->numRows() - 1) {
            throw (string("row"));
        }
        else if (col > this->numCols() - 1) {
            throw (string("col"));
        }
        else {
            return content[row][col];
        }
    } catch (string dim) {
        cout << dim << " out of range" << endl;
        return 0;
    }
}

// Select individual columns from a matrix (helper function for QRFactorization)
template<typename T>
vector<T> Matrix<T>::getCol(int col) {
    try {
        if (col > this->numCols() - 1) {
            throw 505;
        }
        else {
            Matrix<T> transposed = this->transpose();
            return transposed.content[col];
        }
    } catch(...) {
        cout << "Column out of range";
        return {0};
    }
}

template<typename T>
void Matrix<T>::setVal(int row, int col, T newVal) {
    content[row][col] = newVal;
}

template<typename T>
void Matrix<T>::setContent(vector<vector<T>> newContent) {
    content = newContent;
}

// Elementary row operations
template<typename T>
void Matrix<T>::swapRow(int row1, int row2) {
    try {
        if (row1 < this->numRows() && row2 < this->numRows()) {
            swap(content[row2], content[row1]);
        } else { throw 505; }
    } catch (...) {
        cout << "Row index out of range" << endl;
    }
}

template<typename T>
void Matrix<T>::multiplyRow(int row, T scalar) {
    try {
        if (row < this->numRows()) {
            for (int i = 0; i < this->numCols(); i++) {
                content[row][i] *= scalar;
            }
        } else { throw 505; }
    } catch (...) {
        cout << "Row index out of range" << endl;
    }
}

template<typename T>
void Matrix<T>::addScalarMultiple(int rowTarget, int rowAdding, T scalar) {
    try {
        if (rowAdding > this->numRows() || rowTarget > this->numRows()) {
            throw 505;
        }
        else {
            vector<T> multiplied = {};
            for (int i = 0; i < this->numCols(); i++) {
                multiplied.push_back(this->getVal(rowAdding, i) * scalar);
            }
            for (int k = 0; k < multiplied.size(); k++) {
                content[rowTarget][k] += multiplied[k];
            }
        }
    } catch (...) {
        cout << "Invalid dimensions" << endl;
    }
}

// Transpose matrix
template<typename T>
Matrix<T> Matrix<T>::transpose() {
    vector<vector<T>> newContent = {};
    for (int i = 0; i < this->numCols(); i++) {
        vector<T> intermediary = {};
        for (int j = 0; j < this->numRows(); j++) {
            intermediary.push_back({content[j][i]});
        }
        newContent.push_back(intermediary);
    }
    Matrix<T> newMat = Matrix(newContent);
    return newMat;
}

// Convert matrix to row echelon form (unfinished)
template<typename T>
void Matrix<T>::ref() {
    int col = 0;
    for (int row = 0; row < min(this->numRows(), this->numCols()); row++) {
        this->multiplyRow(row, 1/this->getVal(row, col));
        
        //this->printMatrix();
        //cout << endl;
        for (int row2 = row+1; row2 < this->numRows(); row2++) {
            this->addScalarMultiple(row2, row, -this->getVal(row2, col));
        }
        
        //this->printMatrix();
        //cout <<endl;
        col++;
    }
}

// Calculate determinant using cofactor expansion
template<typename T>
double Matrix<T>::determinant() {
    double det;
    
    try {
        if (this->numCols() != this->numRows()) {
            throw 505;
        }
        
        else {
            switch(this->numRows()) {
                case 1:
                    det = this->getVal(0, 0);
                    break;
                case 2:
                    det = (this->getVal(0, 0) * this->getVal(1, 1)) - (this->getVal(0, 1) * this->getVal(1, 0));
                    break;
                case 3:
                    det = (this->getVal(0, 0) * this->getVal(1, 1) * this->getVal(2, 2)) + (this->getVal(0, 1) * this->getVal(1, 2) * this->getVal(2, 0)) + (this->getVal(0, 2) * this->getVal(1, 0) * this->getVal(2, 1)) - (this->getVal(0, 2) * this->getVal(1, 1) * this->getVal(2, 0)) - (this->getVal(0, 0) * this->getVal(1, 2) * this->getVal(2, 1)) - (this->getVal(0, 1) * this->getVal(1, 0) * this->getVal(2, 2));
                    break;
                default:
                    det = 0;
                    for (int i = 0; i < this->numCols(); i++) {
                        vector<vector<T>> newMatContent = {};
                        for (int j = 1; j < this->numCols(); j++) {
                            vector<T> newRow = {};
                            for (int k = 0; k < this->numCols(); k++) {
                                if (k != i) {
                                    newRow.push_back(this->getVal(j, k));
                                }
                            }
                            newMatContent.push_back(newRow);
                        }
                        Matrix<T> newMat = Matrix(newMatContent);
                        det += (this->getVal(0, i)) * (pow(-1, 2+i)) * newMat.determinant();
                    }
            }
        }
    } catch(...) {
        cout << "Matrix must be square";
        return 0;
    }
    
    return det;
}

// Return cofactor matrix
template<typename T>
Matrix<T> Matrix<T>::cofactor() {
    vector<vector<T>> cofactorContent = {};
    
    for (int i = 0; i < this->numRows(); i++) {
        vector<T> cofactorRow = {};
        for (int j = 0; j < this->numCols(); j++) {
            vector<vector<T>> prunedContent = {};
            for (int i2 = 0; i2 < this->numRows(); i2++) {
                vector<T> prunedRow = {};
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
template<typename T>
Matrix<T> Matrix<T>::inverse() {
    try {
        if (this->determinant() == 0) {
            throw 505;
        }
        else {
            Matrix adjugate = (this->cofactor()).transpose();
            return (1/this->determinant()) * adjugate;
        }
    } catch(...) {
        cout << "Matrix is not invertible";
        return Matrix(this->content);
    }
}

// Assignment
template<typename T>
Matrix<T>& Matrix<T>::operator=(Matrix<T>& newMat) {
    this->setContent(newMat.content);
    return *this;
}

// Matrix-Matrix multiplication
template<typename T>
Matrix<T>& operator*=(Matrix<T>& lhs, Matrix<T> rhs) {
    try {
        if (lhs.numCols() != rhs.numRows()) { throw 505; }
        else {
            vector<vector<T>> newContent = {};
            
            for (int i = 0; i < lhs.numRows(); i++) {
                vector<T> newRow = {};
                for (int j = 0; j < rhs.numCols(); j++) {
                    T sum = 0;
                    for (int k = 0; k < lhs.numCols(); k++) {
                        sum += (lhs.getVal(i, k) * rhs.getVal(k, j));
                    }
                    newRow.push_back(sum);
                }
                newContent.push_back(newRow);
            }
            
            lhs.setContent(newContent);
        }
    } catch (...) {
        cout << "Invalid dimensions" << endl;
    }
    
    return lhs;
}

template<typename T>
Matrix<T> operator*(Matrix<T> lhs, Matrix<T> rhs) {
    return lhs *= rhs;
}

// Matrix-Scalar multiplication
template<typename T>
Matrix<T>& operator*=(Matrix<T>& lhs, double scalar) {
    for (int row = 0; row < lhs.numRows(); row++) {
        for (int col = 0; col < lhs.numCols(); col++) {
            lhs.setVal(row, col, scalar * lhs.getVal(row, col));
        }
    }
    return lhs;
}

template<typename T>
Matrix<T> operator*(Matrix<T> lhs, double scalar) {
    return lhs *= scalar;
}
// For commutativity
template<typename T>
Matrix<T> operator*(double scalar, Matrix<T> rhs) {
    return rhs *= scalar;
}

// Addition
template<typename T>
Matrix<T>& operator+=(Matrix<T>& lhs, Matrix<T> rhs) {
    try {
        if (lhs.numCols() != rhs.numCols() || lhs.numRows() != rhs.numRows()) {
            throw 505;
        }
        else {
            for (int row = 0; row < lhs.numRows(); row++) {
                for (int col = 0; col < lhs.numCols(); col++) {
                    lhs.setVal(row, col, lhs.getVal(row, col) + rhs.getVal(row, col));
                }
            }
        }
    } catch (...) {
        cout << "Invalid dimensions" << endl;
    }
    
    return lhs;
}

template<typename T>
Matrix<T> operator+(Matrix<T> lhs, Matrix<T> rhs) {
    return lhs += rhs;
}

// Subtraction
template <typename T>
Matrix<T>& operator-=(Matrix<T>& lhs, Matrix<T> rhs) {
    try {
        if (lhs.numCols() != rhs.numCols() || lhs.numRows() != rhs.numRows()) {
            throw 505;
        }
        else {
            for (int row = 0; row < lhs.numRows(); row++) {
                for (int col = 0; col < lhs.numCols(); col++) {
                    lhs.setVal(row, col, lhs.getVal(row, col) - rhs.getVal(row, col));
                }
            }
        }
    } catch (...) {
        cout << "Invalid dimensions" << endl;
    }
    
    return lhs;
}

template <typename T>
Matrix<T> operator-(Matrix<T> lhs, Matrix<T> rhs) {
    return lhs -= rhs;
}

// Performs Gram-Schmidt process on columns of matrix
// Result can be made orthonormal by setting parameter to true
template <typename T>
Matrix<T> Matrix<T>::gramSchmidt(bool isOrthonormal) {
    // Transpose matrix to get columns as individual vectors
    Matrix transposed = this->transpose();
    
    vector<vector<T>> qContent = {};
    if (isOrthonormal) {
        vector<T> result = transposed.content[0];
        T squareTotal = 0;
        for (int l = 0; l < result.size(); l++) {
            squareTotal += result[l] * result[l];
        }
        T magnitude = sqrt(squareTotal);
        for (int m = 0; m < result.size(); m++) {
            result[m] /= magnitude;
        }
        qContent.push_back(result);
    }
    else {
        qContent.push_back(transposed.content[0]);
    }
    
    for (int i = 1; i < this->numCols(); i++) {
        vector<T> result = transposed.content[i];
        for (int j = 0; j < i; j++) {
            vector<T> proj = projection(transposed.content[i], qContent[j]);
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


// Calculates projection of vec in direction of dir (helper function for gramSchmidt)
template <typename T>
vector<T> Matrix<T>::projection(vector<T> vec, vector<T> dir) {
    double numerator = 0;
    double denominator = 0;
    vector<double> result = {};
    
    for (int i = 0; i < vec.size(); i++) {
        numerator += (vec[i] * dir[i]);
        denominator += (dir[i] * dir[i]);
    }
    
    for (int j = 0; j < vec.size(); j++) {
        result.push_back(dir[j] * (numerator/denominator));
    }
    
    return result;
}

// Calculates the dot product of two vectors (helper function)
template <typename T>
T Matrix<T>::dotProduct(vector<T> vec1, vector<T> vec2) {
    T result = 0;
    for (int i = 0; i < vec1.size(); i++) {
        result += (vec1[i] * vec2[i]);
    }
    return result;
}

// Performs QR factorization on matrix
// Returns a vector of two matrices containing Q and R
template <typename T>
vector<Matrix<T>> Matrix<T>::QRFactorization() {
    Matrix<T> Q = this->gramSchmidt(true);
    
    Matrix<T> notOrthonormal = this->gramSchmidt(false);
    
    vector<vector<T>> rContent = {};
    for (int i = 0; i < this->numCols(); i++) {
        vector<T> rRow = {};
        for (int j = 0; j < this->numCols(); j++) {
            if (j < i) {
                rRow.push_back(0);
            }
            else if (i == j) {
                vector<T> magCol = notOrthonormal.getCol(i);
                T result = 0;
                for (int k = 0; k < magCol.size(); k++) {
                    result += magCol[k] * magCol[k];
                }
                result = sqrt(result);
                rRow.push_back(result);
            }
            else {
                rRow.push_back(dotProduct(this->getCol(j), Q.getCol(i)));
            }
        }
        rContent.push_back(rRow);
    }
    
    Matrix<T> R = Matrix(rContent);
    
    vector<Matrix> result = {Q, R};
    return result;
}

// Find eigenvalues using QR algorithm
template <typename T>
vector<T> Matrix<T>::eigenvalues() {
    vector<T> result = {};
    
    Matrix<T> A = Matrix(this->content);
    
    for (int i = 0; i < 10; i++) {
        vector<Matrix> QR = A.QRFactorization();
        Matrix<T> product = QR[1] * QR[0];
        A = product;
    }
    
    for (int j = 0; j < A.numCols(); j++) {
        result.push_back(round(A.getVal(j, j)));
    }
    
    return result;
}





int main() {
    // Initialise matrices
    // Matrix class is a template so can be used with int, double and float types
    Matrix testMat1 = Matrix<double>({{8,2,3}, {1,3,5}, {6,8,1}});
    Matrix testMat2 = Matrix<double>({{1,2,3}, {7,6,5}, {8,9,5}});
    
    // Output matrix
    testMat1.printMatrix();
    
    // Get dimensions
    cout << endl << testMat1.numRows() << " x " << testMat1.numCols() << endl;
    
    // Get value at position
    cout << endl << "Value at pos 0x0: " << testMat1.getVal(0, 0) << endl;
    
    // Set value at position
    testMat1.setVal(0, 0, 5);
    cout << endl << "Value now at pos 0x0: " << testMat1.getVal(0, 0) << endl;
    
    // Perform elementary row operations
    testMat1.swapRow(0, 1);
    testMat1.multiplyRow(1, 5);
    testMat1.addScalarMultiple(0, 1, 4);
    cout << endl << "Apply elementary row operations:" << endl;
    testMat1.printMatrix();
    
    // Transpose matrix
    cout << endl << "Transpose matrix:" << endl;
    Matrix transposed = testMat1.transpose();
    transposed.printMatrix();
    
    // Transform matrix to row echelon form
    cout << endl << "Carry to REF:" << endl;
    testMat1.ref();
    testMat1.printMatrix();
    
    // Get matrix determinant
    cout << endl << "Determinant: " << testMat2.determinant() << endl;
    
    // Find cofactor matrix
    cout << endl << "Cofactor matrix:" << endl;
    Matrix cofactor = testMat2.cofactor();
    cofactor.printMatrix();
    
    // Invert matrix
    cout << endl << "Invert matrix:" << endl;
    Matrix inverse = testMat2.inverse();
    inverse.printMatrix();
    
    // Orthogonalise matrix using gram-schmidt
    cout << endl << "Orthogonalise matrix:" << endl;
    Matrix orthogonal = testMat2.gramSchmidt(false);
    orthogonal.printMatrix();
    
    // Orthonormalize matrix
    cout << endl << "Orthonormalize matrix:" << endl;
    Matrix orthonormal = testMat2.gramSchmidt(true);
    orthonormal.printMatrix();
    
    // Perform QR factorization
    vector<Matrix<double>> qr = testMat2.QRFactorization();
    
    // Return vector of eigenvalues
    cout << endl << "Find eigenvalues: ";
    vector<double> eigenvalues = testMat2.eigenvalues();
    for (auto i : eigenvalues) {
        cout << i << ", ";
    }
    cout << endl << endl;
    
    return 0;
}
