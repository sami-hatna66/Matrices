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
};

unsigned long Matrix::numRows() {
    return content.size();
}

unsigned long Matrix::numCols() {
    return content[0].size();
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

Matrix& operator*=(Matrix& lhs, Matrix rhs) {
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
    
    return lhs;
}

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



int main() {
    Matrix test1 = Matrix({{1, 2, 4}, {3, 4, 5}, {44, 55, 33}});
    
    Matrix test2 = Matrix({{33, 55, 44}, {55, 44, 33}, {44, 55, 33}});
    
    test2 *= test1;
    test2.printMatrix();
    
    return 0;
}
