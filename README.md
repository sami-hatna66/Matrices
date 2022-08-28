#  C++ Matrices

A linear algebra module implemented in C++

### Features:
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

### Example usage:

```cpp
int main() {
    // Initialise matrices
    // Matrix class is a template so can be used with int, double and float types
    Matrix testMat1 = Matrix<double>({{8,2,3}, {1,3,5}, {6,8,1}});
    Matrix testMat2 = Matrix<double>({{-2,-4,2}, {-2,1,2}, {4,2,5}});
    
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
    
    // Matrix multiplication
    Matrix multiply = testMat1 * testMat2;
    cout << endl << "Multiply Matrix 1 by Matrix 2:" << endl;
    multiply.printMatrix();
    
    // Matrix addition
    testMat1 += testMat2;
    cout << endl << "Add Matrix 2 to Matrix 1:" << endl;
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

### Output:

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

Multiply Matrix 1 by Matrix 2:
|  -28.000000 -231.000000  613.000000 |
|  -10.000000  -60.000000  145.000000 |
|  -24.000000  -14.000000   33.000000 |

Add Matrix 2 to Matrix 1:
| 99.000000 39.000000 67.000000 |
| 23.000000 11.000000 17.000000 |
| 10.000000 10.000000  6.000000 |

Transpose matrix:
| 99.000000 23.000000 10.000000 |
| 39.000000 11.000000 10.000000 |
| 67.000000 17.000000  6.000000 |

Carry to REF:
|  1.000000  0.393939  0.676768 |
|  0.000000  1.000000  0.739583 |
| -0.000000 -0.000000  1.000000 |

Determinant: -90

Cofactor matrix:
|   1.000000  18.000000  -8.000000 |
|  24.000000 -18.000000 -12.000000 |
| -10.000000  -0.000000 -10.000000 |

Invert matrix:
| -0.011111 -0.266667  0.111111 |
| -0.200000  0.200000  0.000000 |
|  0.088889  0.133333  0.111111 |

Orthogonalise matrix:
| -2.000000 -2.833333  2.337662 |
| -2.000000  2.166667  3.506494 |
|  4.000000 -0.333333  2.922078 |

Orthonormalize matrix:
| -0.408248 -0.790912  0.455842 |
| -0.408248  0.604815  0.683763 |
|  0.816497 -0.093048  0.569803 |

Find eigenvalues: 6, -5, 3, 

Program ended with exit code: 0
```
