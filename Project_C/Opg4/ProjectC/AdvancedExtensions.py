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
    # Defining a function calculate determinant value of given matrix A
    def RecursiveDet(B: Matrix, j: int):
        # Calls a function to get the square submatrix 
        # after excluding i-th and j-th column
        SSM = SquareSubMatrix(B,0,j)
        if SSM.M_Rows == 1: # checks if the matrix is a vector or not
            return SSM[0,0]
        else:
            value = 0 # defining value as 0
            # Checks for all elements in the row of the submatrix
            for b in range(SSM.M_Rows):
                # Calls the function recursively to add into value
                value += SSM[0,b]*(-1)**(b)*RecursiveDet(SSM, b)
            return value

    M = A.M_Rows
    det = 0     

    # Calls the function recursively
    for x in range(M):
        det += A[0,x]*(-1)**(x)*RecursiveDet(A, x)

    return det


def VectorNorm(v: Vector) -> float:
    nv = 0.0

    for i in range(len(v)):
        nv += v[i]**2
    return math.sqrt(nv)


def SetColumn(A: Matrix, v: Vector, j: int) -> Matrix:
    M = A.M_Rows

    # Insert the row of the vector into the 
    # row of the matrix on the corresponding j column
    for i in range(M):
        A[i,j] = v[i]  

    return A


def GramSchmidt(A: Matrix) -> tuple:
    M = A.M_Rows
    N = A.N_Cols
    Q = Matrix(M, N) # Define matrix Q
    R = Matrix(N, N) # Define matrix R

    for b in range(N):
        for x in range(M):
            Q[x, b] = A[x, b]
        for a in range(b):
            for x in range(M):
                # Insert the elements for the upper triangular matrix
                R[a, b] += Q[x, a] * A[x, b]
            for x in range(M):
                # Insert the elements for the M x N orthonormal matrix
                Q[x, b] -= R[a, b] * Q[x, a]
        for x in range(M):
            if Q[x, b] != 0:
                R[b, b] = VectorNorm(Q.Column(b))
                for y in range(M):
                    Q[y, b] /= R[b, b]
                break
    return (Q, R)

