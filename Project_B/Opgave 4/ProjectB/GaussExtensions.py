# -*- coding: utf-8 -*-
"""
@Project: LinalgDat2022
@File: GaussExtensions.py

@Description: Project B Gauss extensions

"""

import math
import sys

sys.path.append('../Core')
from Vector import Vector
from Matrix import Matrix


def AugmentRight(A: Matrix, v: Vector) -> Matrix:
    
    M = A.M_Rows
    N = A.N_Cols
    if v.size() != M:
        raise ValueError("number of rows of A and length of v differ.")

    B = Matrix(M, N + 1)
    for i in range(M):
        for j in range(N):
            B[i, j] = A[i, j]
        B[i, N] = v[i]
    return B


def ElementaryRowReplacement(A: Matrix, i: int, m: float, j: int) -> Matrix:
    
    N = A.N_Cols
    for b in range(N):
        A[i,b] = A[i,b] + m*A[j,b] 
        # The numbers in row i added 
        # with m multiplied with the numbers in row j
    return A



def ElementaryRowInterchange(A: Matrix, i: int, j : int) -> Matrix:

    N = A.N_Cols
    for b in range(N):
        B = A[i,b] 
        A[i,b] = A[j,b]
        A[j,b] = B # A[i,b] becomes A[j,b] and A[j,b] becomes A[i,b]
    return A # A modified in-place after row interchange



def ElementaryRowScaling(A: Matrix, i: int, c: float) -> Matrix:
    
    N = A.N_Cols
    for b in range(N):
        A[i,b] = c*A[i,b] # The number in A[i,b] is multiplied by c
    return A # A modified in-place after row scaling.



def ForwardReduction(A: Matrix) -> Matrix:

    M = A.M_Rows # Define the rows as M
    N = A.N_Cols # Define the columns as N

    for a in range(M): # For loop that stops at the last row
        for j in range(a, N): # For loop that starts at a and stops at the last column
            pivot = 0 # Resets the pivot indicator to start a new column clean
            for i in range(a, M): # For loop that starts at a and stops at the last row
                if(A[i,j] != 0): # Checks if the number is zero or not
                    pivot = 1 # If the number is not zero it becomes a pivot point
                    ElementaryRowInterchange(A, i, a) # Calls ElementaryRowInterchange()
                                                      # to move the pivot point to the top
                    for b in range(a + 1, M): # For loop that starts at the number below a pivot point 
                        if(A[b, j] != 0): # Checks if the number below a pivot point is zero or not
                                          # If the number is not zero ElementaryRowReplacement()
                            ElementaryRowReplacement(A, b, -A[b, j] / A[a, j], a) 
                            # The row is changed so that the number below the pivot point is changed to zero
                    break 

            if(pivot == 1): # If statement to break to the next column
                break

    return A # M x N matrix which is the row-echelon form of A



def BackwardReduction(A: Matrix) -> Matrix:
    
    M = A.M_Rows # Define the rows as M
    N = A.N_Cols # Define the columns as N
    a = M - 1 # a is the last row
    b = N - 1 # b is the last column

    for j in range(b, 0, -1): # For loop that starts from the bottom of the column,
                              # ends at zero and it is going backward due to the negative increment
        if(A[M - 1, j - 1] == 0 and A[M - 1, j] != 0):  # Checks if the number the loop is on is zero or not
            b = j 
            break
    
    while a >= 0 and b >= 0: # a and b is zero or a positive number
        if(A[a,b] != 0): # Check to see if the number is zero or not
            ElementaryRowScaling(A, a, 1 / A[a,b]) 
            # ElementaryRowScaling is called to turn the number into one
            for pivot in range(a): # For loop for all numbers in the row
                if(A[pivot,b] != 0): # Checks the pivot number
                    ElementaryRowReplacement(A, pivot, -A[pivot, b], a)
                    # ElementaryRowReplacement is called to turn the number above the pivot into zero
            a -= 1 # a decrease to check the next row
            b -= 1 # b decrease to check the next column
        else: # when A[a,b] is zero 
            if(a > 1): # Checks to see if a is above one or not
                a -= 1 # a decrease to check the next row
            else:
                b -= 1 # b decrease to check the next column
    b -= 1 # b decrease to check the next column
    
    return A # M x N matrix which is the reduced form of A



def GaussElimination(A: Matrix, v: Vector) -> Vector:

    # Combine AugmentRight, ForwardReduction and BackwardReduction!
    
    B = AugmentRight(A,v) # I placed vector v to the right side of matrix A
                          # I define the augmented matrix as matrix B
    
    ForwardReduction(B) # I call ForwardReduction because I want to 
                        # fill the lower triangular matrix with zeros 
                        # and find the pivot points

    BackwardReduction(B) # I call BackwardReduction because I want to 
                         # fill the upper triangular matrix with zero
    
    M = B.M_Rows # Define the rows as M
    N = B.N_Cols # Define the columns as N
    V = Vector(M) # Define a new vector V with M rows

    for i in range(M):
        V[i] = B[i, N - 1] # Vector V is given value from matrix B 
                           # where B's column is reduced by one

    return V # This function returns M-size solution vector of the system


