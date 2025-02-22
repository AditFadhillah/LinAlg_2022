"""
@Project: LinalgDat2022
@File: BasicExtensions.py

@Description: Project A basic extensions.

"""
import math
import sys

sys.path.append('../Core')
from Vector import Vector
from Matrix import Matrix

__author__ = "François Lauze, University of Copenhagen"
__date__ = "3/28/22"
__version__ = "0.0.1"


def AugmentRight(A: Matrix, v: Vector) -> Matrix:
    """
    Create an augmented matrix from a matrix A and a vector v.

    See page 12 in 'Linear Algebra for Engineers and Scientists'
    by K. Hardy.

    :param A: a matrix of size M-by-N.
    :param v: a column vector of size M.

    :return: a matrix of size M-by-(N + 1)
    """
    M = A.M_Rows
    N = A.N_Cols
    B = Matrix(M, N + 1)
    for i in range(M):
        for j in range(N):
            B[i, j] = A[i, j]
    for i in range(M):
        B[i, N] = v[i]

    return B


def MatVecProduct(A: Matrix, v: Vector) -> Vector:
    """
    This function computes the matrix-vector product of a matrix A
    and a column vector v

    See page 68 in "Linear Algebra for Engineers and Scientists"
    by K. Hardy.
    :param A: an M-by-N Matrix.
    :param v: a size N Vector.

    :return: a size M Vector y such that y = A.v
    """

    M = A.M_Rows
    N = A.N_Cols
    y = Vector(M)
    for i in range(M):
        for j in range(N):
            y[i] += v[j]*A[i,j]                     # Matrix A is multiplied with vector v

    return y


def MatrixProduct(A: Matrix, B: Matrix) -> Matrix:
    """
    Compute the Matrix product of two given matrices A and B.

    See page 58 in "Linear Algebra for Engineers and Scientists"
    by K. Hardy.

    :param A: an M-by-N Matrix.
    :param B: an N-by-P Matrix.

    :returns: the M-by-P Matrix A*B.
    """
    M = A.M_Rows
    N = A.N_Cols
    P = B.N_Cols
    C = Matrix(M, P)
    for i in range(M):
        for j in range(N):
            for k in range(P):
                C[i,k] += A[i,j]*B[j,k]             # Matrix A is multiplied with matrix B
                
    return C


def Transpose(A: Matrix) -> Matrix:
    """
    Computes the transpose of a given Matrix.

    See page 69 in "Linear Algebra for Engineers and Scientists"
    by K. Hardy.

    :param A: A M-by-N Matrix.
    :returns: A N-by-M Matrix B such that B = A^T.
    """
    M = A.M_Rows
    N = A.N_Cols
    B = Matrix(N,M)
    for i in range(M):
        for j in range(N):
            B[j,i] = A[i,j]                         # The value of i and j which represent columns and rows are switched
                                                    # Matrix A rows become matrix B columns, and matrix A columns become matrix B rows
    return B


def VectorNorm(v: Vector) -> float:
    """
    Computes the Euclidean Vector norm of a given Vector.

    See page 197 in "Linear Algebra for Engineers and Scientists"
    by K.Hardy.

    :param v: An N - dimensional Vector.
    :return: The Euclidean norm of the Vector.
    """
    y = 0                                           # y is the variable for the Euclidean norm of the vector, with the start value of 0
    for i in range(len(v)):                         # For loop to go through all elements of vector v
        y += v[i]**2                                # All of the elements in vector v are 2 powered
    y = math.sqrt(y)                                # All of the 2 powered elements are then squared rooted
    
    return y