# -*- coding: utf-8 -*-
"""
@Project: LinalgDat2022
@File: AdvancedExtensions.py

@Description: Project C Determinant and Gram-Schmidt extensions.

"""

import math
import sys
from unittest import skip

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
    B = Matrix(M-1,N-1)
    
    for x in range(M-1):
        for y in range(N-1): # Checks all the elements in one row
            if x < i and y < j: # If A is a 3 x 3 matrix, this would add the element on the top left corner to B
                B[x,y] = A[x,y]
            if y >= j: # If A is a 3 x 3 matrix, this would add the element on the top right corner to B
                B[x,y] = A[x,y+1]
            if x >= i: # If A is a 3 x 3 matrix, this would add the element on the bottom left corner to B
                B[x,y] = A[x+1,y]
            if x >= i and y >= j: # If A is a 3 x 3 matrix, this would add the element on the bottom right corner to B
                B[x,y] = A[x+1,y+1]

    return B



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

    def RecDet(B: Matrix, j: int):
        DetM = SquareSubMatrix(B,0,j)
        if DetM.M_Rows == 1:
            return DetM[0,0]
        else:
            tempDet = 0
            for y in range(DetM.M_Rows):
                tempDet += DetM[0,y]*(-1)**(y)*RecDet(DetM, y)
            return tempDet

    M = A.M_Rows
    det = 0     

    for i in range(M):
        det += A[0,i]*(-1)**(i)*RecDet(A, i)

    return det
    


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
    M = A.M_Rows
    N = A.N_Cols
    if j > N or len(v) != M:
        raise ValueError
    for i in range(M):
        A[i,j] = v[i]  

    return A


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

    M = A.M_Rows
    N = A.N_Cols
    Q = Matrix(M, N)
    R = Matrix(N, N)
    for j in range(N):
        for a in range(M):
            Q[a, j] = A[a, j]
        for i in range(j):
            R[i, j] = 0.0
            for a in range(M):
                R[i, j] += Q[a, i] * A[a, j]
            for a in range(M):
                Q[a, j] -= R[i, j] * Q[a, i]
        for a in range(M):
            if Q[a, j] != 0:
                R[j, j] = VectorNorm(Q.Column(j))
                for b in range(M):
                    Q[b, j] /= R[j, j]
                break
    return (Q, R)