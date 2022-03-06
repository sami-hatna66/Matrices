#  C++ Matrices

A linear algebra module implemented in C++

####Features:
- Template support = working with different numerical data types
- Pretty printing
- Addition, Subtraction and Multiplication
- Get dimensions
- Assign values
- Elementary row operations
    - Swapping rows
    - Multiplying a row by a scalar
    - Adding a scalar multiple of a row to another row
- Transposition
- Row echelon form
- Calculate determinant
- Find cofactor matrix
- Find matrix inverse
- Orthogonalise or Orthonormalise using Gram-Schmidt
- Calculate vector projection
- Perform QR Factorization
- Calculate vector dot product
- Calculate Eigenvalues

###Example usage:

```cpp
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
```

###Output:

```
| 8.000000 2.000000 3.000000 |
| 1.000000 3.000000 5.000000 |
| 6.000000 8.000000 1.000000 |

3 x 3

Value at pos 0x0: 8

Value now at pos 0x0: 5

Apply elementary row operations:
| 101.000000  43.000000  65.000000 |
|  25.000000  10.000000  15.000000 |
|   6.000000   8.000000   1.000000 |

Transpose matrix:
| 101.000000  25.000000   6.000000 |
|  43.000000  10.000000   8.000000 |
|  65.000000  15.000000   1.000000 |

Carry to REF:
|  1.000000  0.425743  0.643564 |
| -0.000000  1.000000  1.692308 |
| -0.000000 -0.000000  1.000000 |

Determinant: 40

Cofactor matrix:
| -15.000000   5.000000  15.000000 |
|  17.000000 -19.000000   7.000000 |
|  -8.000000  16.000000  -8.000000 |

Invert matrix:
| -0.375000  0.425000 -0.200000 |
|  0.125000 -0.475000  0.400000 |
|  0.375000  0.175000 -0.200000 |

Orthogonalise matrix:
|  1.000000  0.982456  1.775148 |
|  7.000000 -1.122807  0.828402 |
|  8.000000  0.859649 -0.946746 |

Orthonormalize matrix:
|  0.093659  0.570568  0.815892 |
|  0.655610 -0.652077  0.380750 |
|  0.749269  0.499247 -0.435143 |

Find eigenvalues: 15, -2, -1, 

Program ended with exit code: 0
```

