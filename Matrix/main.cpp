//
//  main.cpp
//  Matrix
//
//  Created by Sami Hatna on 04/12/2021.
//

#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

class Matrix {
private:
    vector<vector<double>> content;
public:
    Matrix(vector<vector<double>> c) {
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
    
    void printMatrix();
    unsigned long numRows();
    unsigned long numCols();
    void transpose();
    
    // matrix indexing begins at 0
    double getVal(int row, int col) {
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
    
    void setVal(int row, int col, double newVal) {
        content[row][col] = newVal;
    }
    
    void setContent(vector<vector<double>> newContent) {
        content = newContent;
    }
    
    // Row operations
    void swapRow(int row1, int row2) {
        try {
            if (row1 < this->numRows() && row2 < this->numRows()) {
                swap(content[row2], content[row1]);
            } else { throw 505; }
        } catch (...) {
            cout << "Row index out of range" << endl;
        }
    }
    
    void multiplyRow(int row, double scalar) {
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
    
    void addScalarMultiple(int rowTarget, int rowAdding, int scalar) {
        try {
            if (rowAdding > this->numRows() || rowTarget > this->numRows()) {
                throw 505;
            }
            else {
                vector<double> multiplied = {};
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
    
    // Convert matrix to row echelon form
    void ref() {
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
    double determinant() {
        double det;
        
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
                this->printMatrix();
                for (int i = 0; i < this->numCols(); i++) {
                    vector<vector<double>> newMatContent = {};
                    for (int j = 1; j < this->numCols(); j++) {
                        vector<double> newRow = {};
                        for (int k = 0; k < this->numCols(); k++) {
                            if (k != i) {
                                newRow.push_back(this->getVal(j, k));
                            }
                        }
                        newMatContent.push_back(newRow);
                    }
                    cout << endl;
                    Matrix newMat = Matrix(newMatContent);
                    newMat.printMatrix();
                    det += (this->getVal(0, i)) * (pow(-1, 2+i)) * newMat.determinant();
                }
        }
        
        return det;
    }
};

unsigned long Matrix::numRows() {
    return content.size();
}

unsigned long Matrix::numCols() {
    return content[0].size();
}
 
// Transpose matrix
void Matrix::transpose() {
    vector<vector<double>> newContent = {};
    for (int i = 0; i < this->numCols(); i++) {
        vector<double> intermediary = {};
        for (int j = 0; j < this->numRows(); j++) {
            intermediary.push_back({content[j][i]});
        }
        newContent.push_back(intermediary);
    }
    this->setContent(newContent);
}

void Matrix::printMatrix() {
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
                 [&forPrint, &maxLen, *this](double &n){
                    for (int j = 0; j < (maxLen - to_string(n).length()); j++) {
                        forPrint += " ";
                    }
                    forPrint += to_string(n) + " ";
        });
        forPrint += "|";
        cout << forPrint << endl;
    }
}

// Matrix-Matrix multiplication
Matrix& operator*=(Matrix& lhs, Matrix rhs) {
    try {
        if (lhs.numCols() != rhs.numRows()) { throw 505; }
        else {
            vector<vector<double>> newContent = {};
            
            for (int i = 0; i < lhs.numRows(); i++) {
                vector<double> newRow = {};
                for (int j = 0; j < rhs.numCols(); j++) {
                    double sum = 0;
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

Matrix operator*(Matrix lhs, Matrix rhs) {
    return lhs *= rhs;
}

// Matrix-Scalar multiplication
Matrix& operator*=(Matrix& lhs, double scalar) {
    for (int row = 0; row < lhs.numRows(); row++) {
        for (int col = 0; col < lhs.numCols(); col++) {
            lhs.setVal(row, col, scalar * lhs.getVal(row, col));
        }
    }
    
    return lhs;
}

Matrix operator*(Matrix lhs, double scalar) {
    return lhs *= scalar;
}
// For commutativity
Matrix operator*(double scalar, Matrix rhs) {
    return rhs *= scalar;
}

// Addition
Matrix& operator+=(Matrix& lhs, Matrix rhs) {
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

Matrix operator+(Matrix lhs, Matrix rhs) {
    return lhs += rhs;
}

// Subtraction
Matrix& operator-=(Matrix& lhs, Matrix rhs) {
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

Matrix operator-(Matrix lhs, Matrix rhs) {
    return lhs -= rhs;
}








int main() {
    Matrix test1 = Matrix({{1, 2}, {3, 4}, {44, 55}});
    
    Matrix test2 = Matrix({{5, -7, 2, 2}, {0, 3, 0, -4}, {-5, -8, 0, 3}, {0, 5, 0, -6}});
        
    cout << test2.determinant() << endl;
    
    return 0;
}
