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
    // friend Matrix operator*(Matrix&m1, int scalar);
    
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

Matrix operator+(Matrix& m1, Matrix& m2) {
    try {
        if (m1.numCols() != m2.numCols() || m1.numRows() != m2.numRows()) {
            throw 505;
        }
        else {
            vector<vector<double>> result = {};
            
            for (int row = 0; row < m1.numRows(); row++) {
                vector<double> intermediary = {};
                for (int col = 0; col < m1.numCols(); col++) {
                    intermediary.push_back(m1.getVal(row, col) + m2.getVal(row, col));
                }
                result.push_back(intermediary);
            }
            return Matrix(result);
        }
    } catch (...) {
        cout << "Invalid dimensions" << endl;
        return Matrix({{0}});
    }
}

Matrix operator*(Matrix&m, int scalar) {
    vector<vector<double>> result = {};
    
    for (int row = 0; row < m.numRows(); row++) {
        vector <double> intermediary = {};
        for (int col = 0; col < m.numCols(); col++) {
            intermediary.push_back(m.getVal(row, col) * scalar);
        }
        result.push_back(intermediary);
    }
    
    return Matrix(result);
}

// For commutativity
Matrix operator*(int scalar, Matrix &m) {
    return m * scalar;
}

Matrix operator+=(Matrix &m1, Matrix &m2) {
    return Matrix({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});
}




int main() {
    Matrix test1 = Matrix({{1, 2, 4}, {3, 4, 5}, {44, 55, 33}});
    
    Matrix test2 = Matrix({{33, 55, 44}, {55, 44, 33}, {44, 55, 33}});
    
    Matrix test3 = 7 * test1;
    test3.printMatrix();
    
    return 0;
}
