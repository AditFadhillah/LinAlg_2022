# -*- coding: utf-8 -*-
"""
@Project: LinalgDat2022
@File: AdvancedExtensions.py

@Description: Project C Determinant and Gram-Schmidt extensions.

"""

import math
import sys

sys.path.append('../Core')
from Vector import Vector
from Matrix import Matrix

Tolerance = 1e-6


def SquareSubMatrix(A: Matrix, i: int, j: int) -> Matrix:
    """
    This function creates the square submatrix given a square matrix as
    well as row and column indices to remove from it.

    Remarks:
        See page 246-247 in "Linear Algebra for Engineers and Scientists"
        by K. Hardy.

    Parameters:
        A:  N-by-N matrix
        i: int. The index of the row to remove.
        j: int. The index of the column to remove.

    Return:
        The resulting (N - 1)-by-(N - 1) submatrix.
    """
    M = A.M_Rows 
    N = A.N_Cols 
    B = Matrix[M-1,N-1]
    for a in range(M-1):
        for b in range(N-1):
            if a >= i:
                B[a,b] = A[a+1,b]
            if b >= j:
                B[a,b] = A[a,b+1]
            if a >= i and b >= j:
                B[a,b] = A[a+1,b+1]
            if a < i and b < j:
                B[a,b] = A[a,b]
    return B
    



    
    """
    M = A.M_Rows 
    N = A.N_Cols 
    row = 0
    col = 0
    B = B[row,col]
    for a in range(M):
        for b in range(N):
            if b != j:
                B[row,col] = A[M,N]
                col += 1
        row += 1
        col = 0
        if a != i:
            B[row,col] = A[M,N]
    print(B)
    return B 
    
    """
        


def Determinant(A: Matrix) -> float:    
    """
    This function computes the determinant of a given square matrix.

    Remarks:
        * See page 247 in "Linear Algebra for Engineers and Scientists"
        by K. Hardy.
        * Hint: Use SquareSubMatrix.

    Parameter:
        A: N-by-N matrix.

    Return:
        The determinant of the matrix.
    """
    raise NotImplementedError()


def VectorNorm(v: Vector) -> float:
    """
    This function computes the Euclidean norm of a Vector. This has been implemented
    in Project A and is provided here for convenience

    Parameter:
         v: Vector

    Return:
         Euclidean norm, i.e. (\sum v[i]^2)^0.5
    """
    nv = 0.0
    for i in range(len(v)):
        nv += v[i]**2
    return math.sqrt(nv)


def SetColumn(A: Matrix, v: Vector, j: int) -> Matrix:
    """
    This function copies Vector 'v' as a column of Matrix 'A'
    at column position j.

    Parameters:
        A: M-by-N Matrix.
        v: size M vector
        j: int. Column number.

    Return:
        Matrix A  after modification.

    Raise:
        ValueError if j is out of range or if len(v) != A.M_Rows.
    """
    raise NotImplementedError()


def GramSchmidt(A: Matrix) -> tuple:
    """
    This function computes the Gram-Schmidt process on a given matrix.

    Remarks:
        See page 229 in "Linear Algebra for Engineers and Scientists"
        by K. Hardy.

    Parameter:
        A: M-by-N matrix. All columns are implicitly assumed linear
        independent.

    Return:
        tuple (Q,R) where Q is a M-by-N orthonormal matrix and R is an
        N-by-N upper triangular matrix.
    """
    raise NotImplementedError()

