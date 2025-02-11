"""
@project: LinalgDat2022 Projekt C
@file: ProjektB.py

@description: data and routines to test AdvancedExtensions module.
Do not modify, this file is automatically generated.

@author: François Lauze, University of Copenhagen
@date: Friday April 29. 2022

random state = e3759621-be13-421b-855e-99ce1206d92e
"""
from sys import path
path.append("../Core")
from Matrix import Matrix
from Vector import Vector



# Matrix size = (3, 3)
array2D_0001 = [[ -4.66000,  -3.87000,   4.95000],
                [ -0.76000,   2.43000,   1.32000],
                [ -4.80000,  -9.56000,  -5.49000]]
Matrix_0000 = Matrix.fromArray(array2D_0001)
Pair_0002 = (1, 0)
array2D_0004 = [[ -3.87000,   4.95000],
                [ -9.56000,  -5.49000]]
Matrix_0003 = Matrix.fromArray(array2D_0004)
# Matrix size = (11, 11)
array2D_0006 = [[ -4.08000,  -2.35000,  -6.17000,  -5.62000,   5.23000,  -4.88000,  -8.81000,   3.27000,  -1.28000,  -3.89000,   4.79000],
                [ -8.48000,   0.00000,   6.07000,   7.58000,  -8.23000,   3.00000,  -3.10000,   9.18000,   6.96000,  -2.24000,  -0.78000],
                [  4.46000,   5.71000,  -4.38000,   6.71000,  -7.47000,  -8.97000,  -6.09000,   8.64000,  -1.12000,   4.77000,  -5.11000],
                [ -5.79000,  -7.27000,  -5.08000,  -5.58000,   1.85000,  -3.73000,   5.65000,   7.97000,  -2.36000,  -5.54000,   1.86000],
                [ -8.37000,  -1.02000,   0.60000,   7.86000,  -4.18000,   7.85000,  -0.22000,  -0.59000,   6.79000,   9.41000,  -5.43000],
                [ -7.03000,   6.81000,  -5.04000,   5.98000,  -6.96000,   1.56000,  -2.84000,   1.64000,   7.43000,   0.44000,  -3.76000],
                [ -5.52000,   0.83000,   3.51000,  -1.95000,  -9.88000,   7.91000,  -1.45000,   3.26000,  -1.69000,   7.33000,  -1.95000],
                [  8.13000,  -7.00000,   0.22000,   7.26000,  -6.41000,  -8.00000,  -8.10000,   1.17000,  -7.86000,   9.21000,  -7.89000],
                [  8.90000,   8.00000,   3.99000,  -3.52000,  -8.77000,   6.99000,  -0.61000,  -4.86000,  -1.85000,   5.20000,  -9.16000],
                [ -6.55000,  -1.55000,   7.67000,   5.96000,   3.06000,  -4.65000,   0.59000,  -2.70000,   2.92000,   7.11000,  -2.54000],
                [ -1.23000,  -1.04000,   6.49000,   5.35000,  -5.88000,  -1.28000,  -1.36000,  -8.10000,   8.36000,   3.88000,  -1.11000]]
Matrix_0005 = Matrix.fromArray(array2D_0006)
Pair_0007 = (2, 8)
array2D_0009 = [[ -4.08000,  -2.35000,  -6.17000,  -5.62000,   5.23000,  -4.88000,  -8.81000,   3.27000,  -3.89000,   4.79000],
                [ -8.48000,   0.00000,   6.07000,   7.58000,  -8.23000,   3.00000,  -3.10000,   9.18000,  -2.24000,  -0.78000],
                [ -5.79000,  -7.27000,  -5.08000,  -5.58000,   1.85000,  -3.73000,   5.65000,   7.97000,  -5.54000,   1.86000],
                [ -8.37000,  -1.02000,   0.60000,   7.86000,  -4.18000,   7.85000,  -0.22000,  -0.59000,   9.41000,  -5.43000],
                [ -7.03000,   6.81000,  -5.04000,   5.98000,  -6.96000,   1.56000,  -2.84000,   1.64000,   0.44000,  -3.76000],
                [ -5.52000,   0.83000,   3.51000,  -1.95000,  -9.88000,   7.91000,  -1.45000,   3.26000,   7.33000,  -1.95000],
                [  8.13000,  -7.00000,   0.22000,   7.26000,  -6.41000,  -8.00000,  -8.10000,   1.17000,   9.21000,  -7.89000],
                [  8.90000,   8.00000,   3.99000,  -3.52000,  -8.77000,   6.99000,  -0.61000,  -4.86000,   5.20000,  -9.16000],
                [ -6.55000,  -1.55000,   7.67000,   5.96000,   3.06000,  -4.65000,   0.59000,  -2.70000,   7.11000,  -2.54000],
                [ -1.23000,  -1.04000,   6.49000,   5.35000,  -5.88000,  -1.28000,  -1.36000,  -8.10000,   3.88000,  -1.11000]]
Matrix_0008 = Matrix.fromArray(array2D_0009)
# Matrix size = (9, 9)
array2D_0011 = [[ -2.87000,  -7.39000,   3.53000,  -0.05000,  -2.43000,   0.31000,  -7.13000,   1.41000,  -9.80000],
                [  1.59000,  -3.37000,  -9.28000,   7.02000,   2.71000,   1.63000,  -5.78000,   8.36000,  -1.01000],
                [ -8.10000,  -3.35000,   8.65000,  -6.65000,   2.98000,   1.64000,  -0.94000,  -9.38000,  -0.77000],
                [ -3.72000,   8.15000,  -0.52000,  -3.81000,  -3.61000,   3.93000,   6.22000,  -1.27000,   1.64000],
                [ -0.07000,   4.26000,   0.28000,  -2.09000,   4.08000,   7.73000,   5.19000,  -1.61000,   9.80000],
                [  3.44000,  -5.41000,   9.15000,  -7.16000,  -9.78000,   1.10000,   0.33000,  -7.03000,   5.82000],
                [ -7.64000,  -7.69000,  -5.65000,   3.18000,   3.07000,  -3.63000,  -8.24000,   8.31000,   6.09000],
                [  4.37000,   7.19000,  -1.47000,   3.54000,   2.71000,  -9.50000,  -4.99000,  -5.95000,   1.29000],
                [ -1.05000,   4.61000,   7.67000,  -0.21000,   0.43000,   6.02000,   0.05000,   4.63000,  -1.76000]]
Matrix_0010 = Matrix.fromArray(array2D_0011)
Pair_0012 = (6, 0)
array2D_0014 = [[ -7.39000,   3.53000,  -0.05000,  -2.43000,   0.31000,  -7.13000,   1.41000,  -9.80000],
                [ -3.37000,  -9.28000,   7.02000,   2.71000,   1.63000,  -5.78000,   8.36000,  -1.01000],
                [ -3.35000,   8.65000,  -6.65000,   2.98000,   1.64000,  -0.94000,  -9.38000,  -0.77000],
                [  8.15000,  -0.52000,  -3.81000,  -3.61000,   3.93000,   6.22000,  -1.27000,   1.64000],
                [  4.26000,   0.28000,  -2.09000,   4.08000,   7.73000,   5.19000,  -1.61000,   9.80000],
                [ -5.41000,   9.15000,  -7.16000,  -9.78000,   1.10000,   0.33000,  -7.03000,   5.82000],
                [  7.19000,  -1.47000,   3.54000,   2.71000,  -9.50000,  -4.99000,  -5.95000,   1.29000],
                [  4.61000,   7.67000,  -0.21000,   0.43000,   6.02000,   0.05000,   4.63000,  -1.76000]]
Matrix_0013 = Matrix.fromArray(array2D_0014)

SSMMatrixList = [Matrix_0000, Matrix_0005, Matrix_0010]
SSMPairList = [Pair_0002, Pair_0007, Pair_0012]
SSMExpected = [Matrix_0003, Matrix_0008, Matrix_0013]
SSMArgs = [SSMMatrixList, SSMPairList, SSMExpected]



# Matrix size = (9, 9)
array2D_0016 = [[  5.57000,   2.30000,  -2.59000,   0.89000,  -2.84000,  -3.86000,   5.72000,  -3.46000,  -4.64000],
                [  2.35000,   4.37000,   7.14000,  -2.89000,   5.68000,   2.14000,  -1.48000,  -5.76000,  -8.37000],
                [  7.13000,   5.46000,  -8.63000,  -3.28000,  -3.51000,  -4.46000,   6.76000,   6.10000,   5.78000],
                [ -5.29000,  -9.08000,   0.10000,  -2.57000,   1.11000,   1.69000,   8.23000,  -3.20000,  -2.39000],
                [ -8.08000,  -7.56000,  -7.67000,  -8.25000,   7.69000,   1.11000,  -1.80000,   1.71000,  -3.03000],
                [ -7.97000,   7.81000,  -9.23000,   3.13000,   5.74000,  -4.56000,  -6.76000,  -9.06000,   6.58000],
                [  8.76000,   4.63000,   8.61000,  -2.26000,   9.41000,   2.00000,  -6.17000,   2.78000,   0.45000],
                [ -5.35000,  -1.96000,   1.47000,   2.30000,  -3.64000,   3.14000,  -7.77000,   2.25000,  -3.69000],
                [  6.49000,  -3.67000,  -8.44000,   9.68000,  -5.34000,  -1.80000,   4.15000,  -1.52000,  -5.57000]]
Matrix_0015 = Matrix.fromArray(array2D_0016)
Float_0017 = 58822159.21040
# Matrix size = (8, 8)
array2D_0019 = [[  0.49000,  -7.96000,  -8.25000,  -8.53000,  -3.58000,   3.78000,   2.49000,  -4.13000],
                [  4.92000,   9.26000,   4.35000,  -9.08000,   3.00000,  -2.54000,   4.26000,  -7.80000],
                [ -1.50000,  -4.59000,   5.34000,   4.98000,  -5.80000,   8.31000,  -0.50000,  -3.91000],
                [  4.86000,  -8.52000,  -1.49000,   8.73000,  -4.49000,   1.69000,  -0.74000,   8.74000],
                [ -2.66000,   6.83000,  -8.50000,  -9.82000,   2.76000,   6.87000,  -3.89000,   0.71000],
                [ -1.01000,   2.58000,   5.02000,   6.07000,  -2.95000,  -1.63000,  -1.13000,  -7.39000],
                [ -8.10000,   3.80000,   7.30000,   2.71000,  -7.75000,  -5.03000,  -2.34000,   2.17000],
                [ -0.35000,   1.06000,  -9.67000,   5.08000,   5.33000,   2.11000,  -7.81000,   1.45000]]
Matrix_0018 = Matrix.fromArray(array2D_0019)
Float_0020 = 48969515.61345
# Matrix size = (3, 3)
array2D_0022 = [[-2.74000,  5.56000,  5.38000],
                [ 5.20000, -1.04000,  8.95000],
                [ 0.46000, -1.15000,  1.84000]]
Matrix_0021 = Matrix.fromArray(array2D_0022)
Float_0023 = -82.86435

DeterminantMatrixList = [Matrix_0015, Matrix_0018, Matrix_0021]
DeterminantExpected = [Float_0017, Float_0020, Float_0023]
DeterminantArgs = [DeterminantMatrixList, DeterminantExpected]



# Matrix size = (6, 6)
array2D_0025 = [[ -2.11000,  -1.46000,  -9.17000,   1.21000,   4.77000,  -1.33000],
                [ -2.95000,  -1.18000,   0.35000,  -6.44000,  -7.08000,   9.37000],
                [  4.64000,   7.72000,   4.56000,  -0.45000,   3.79000,  -6.10000],
                [-10.00000,   2.55000,   3.12000,  -1.38000,  -7.76000,   5.12000],
                [  5.27000,  -1.46000,   3.36000,   3.33000,   8.32000,  -9.56000],
                [ -3.63000,   2.35000,   7.27000,   9.28000,  -6.98000,   9.37000]]
Matrix_0024 = Matrix.fromArray(array2D_0025)
Int_0026 = 1
array1D_0028 = [-1.55000, -0.98000,  2.40000, -0.99000, -5.99000,  2.55000]
Vector_0027 = Vector.fromArray(array1D_0028)
array2D_0030 = [[ -2.11000,  -1.55000,  -9.17000,   1.21000,   4.77000,  -1.33000],
                [ -2.95000,  -0.98000,   0.35000,  -6.44000,  -7.08000,   9.37000],
                [  4.64000,   2.40000,   4.56000,  -0.45000,   3.79000,  -6.10000],
                [-10.00000,  -0.99000,   3.12000,  -1.38000,  -7.76000,   5.12000],
                [  5.27000,  -5.99000,   3.36000,   3.33000,   8.32000,  -9.56000],
                [ -3.63000,   2.55000,   7.27000,   9.28000,  -6.98000,   9.37000]]
Matrix_0029 = Matrix.fromArray(array2D_0030)
# Matrix size = (11, 5)
array2D_0032 = [[  0.81000,  -0.54000,   3.08000,   9.22000,   8.76000],
                [  9.34000,  -6.62000,  -4.93000,  -7.10000,  -2.82000],
                [ -5.34000,   9.20000,   2.13000,   7.92000,  -3.65000],
                [  7.43000,  -3.50000,   0.64000,  10.00000,  -4.02000],
                [ -3.48000,  -1.26000,   1.11000,  -3.02000,  -3.88000],
                [  8.59000,   5.11000,   9.15000,  -7.02000,  -6.43000],
                [ -5.50000,   9.00000,  -9.03000,  -0.43000,   4.64000],
                [  2.31000,  -4.64000,   7.40000,   7.66000,  -6.64000],
                [ -7.68000,  -5.46000,   2.50000,   7.58000,  -9.15000],
                [ -0.89000,  10.00000,  -2.60000,   1.95000,   0.43000],
                [ -2.82000,  -9.94000,  -4.48000,  -7.73000,  -0.92000]]
Matrix_0031 = Matrix.fromArray(array2D_0032)
Int_0033 = 2
array1D_0035 = [ -9.96000,   7.47000,  -3.78000,   8.39000,  -2.69000,  -6.40000,  -6.82000,   5.56000,   3.60000,  -4.20000,  -8.67000]
Vector_0034 = Vector.fromArray(array1D_0035)
array2D_0037 = [[  0.81000,  -0.54000,  -9.96000,   9.22000,   8.76000],
                [  9.34000,  -6.62000,   7.47000,  -7.10000,  -2.82000],
                [ -5.34000,   9.20000,  -3.78000,   7.92000,  -3.65000],
                [  7.43000,  -3.50000,   8.39000,  10.00000,  -4.02000],
                [ -3.48000,  -1.26000,  -2.69000,  -3.02000,  -3.88000],
                [  8.59000,   5.11000,  -6.40000,  -7.02000,  -6.43000],
                [ -5.50000,   9.00000,  -6.82000,  -0.43000,   4.64000],
                [  2.31000,  -4.64000,   5.56000,   7.66000,  -6.64000],
                [ -7.68000,  -5.46000,   3.60000,   7.58000,  -9.15000],
                [ -0.89000,  10.00000,  -4.20000,   1.95000,   0.43000],
                [ -2.82000,  -9.94000,  -8.67000,  -7.73000,  -0.92000]]
Matrix_0036 = Matrix.fromArray(array2D_0037)
# Matrix size = (2, 11)
array2D_0039 = [[ -8.48000,  -9.34000,  -7.96000,   2.75000,  -5.43000,   2.36000,  -0.02000,  -5.52000,  -1.51000,  -7.27000,  -1.52000],
                [  9.80000,  -1.95000,  -2.26000,   5.84000,  -0.22000,  -2.50000,  -0.92000,  -8.92000,  -2.98000,   9.25000,   5.46000]]
Matrix_0038 = Matrix.fromArray(array2D_0039)
Int_0040 = 0
array1D_0042 = [-2.00000,  5.08000]
Vector_0041 = Vector.fromArray(array1D_0042)
array2D_0044 = [[ -2.00000,  -9.34000,  -7.96000,   2.75000,  -5.43000,   2.36000,  -0.02000,  -5.52000,  -1.51000,  -7.27000,  -1.52000],
                [  5.08000,  -1.95000,  -2.26000,   5.84000,  -0.22000,  -2.50000,  -0.92000,  -8.92000,  -2.98000,   9.25000,   5.46000]]
Matrix_0043 = Matrix.fromArray(array2D_0044)

SCMatrixList = [Matrix_0024, Matrix_0031, Matrix_0038]
SCIndexList = [Int_0026, Int_0033, Int_0040]
SCVectorList = [Vector_0027, Vector_0034, Vector_0041]
SCExpected = [Matrix_0029, Matrix_0036, Matrix_0043]
SCArgs = [SCMatrixList, SCVectorList, SCIndexList, SCExpected]



# Matrix size = (13, 9)
array2D_0046 = [[  7.10000,  -5.46000,   9.50000,   6.94000,  -5.07000,   6.08000,   4.56000,   0.80000,   9.97000],
                [ -9.35000,  -7.75000,  -9.54000,  -6.53000,   6.42000,  -4.04000,  -1.15000,   1.74000,  -1.23000],
                [ -4.27000,   9.41000,   4.07000,   3.51000,  -4.46000,  -5.82000,   0.73000,   7.17000,   1.60000],
                [ -1.39000,  -6.68000,  -6.49000,  -5.26000,   2.82000,   6.16000,   5.33000,   9.61000,  -4.04000],
                [  2.64000,   4.37000,   0.56000,  -4.19000,   8.21000,  -7.93000,   8.05000,   4.09000,   3.98000],
                [  7.56000,  -1.21000,   1.33000,   0.99000,  -2.77000,  -9.54000,   3.14000,   5.43000,   7.81000],
                [ -3.80000,  -6.34000,   5.76000,  -1.42000,  -7.41000,  -6.87000,   1.79000,  -8.28000,   0.89000],
                [  1.41000,  -8.79000,  -9.38000,  -0.21000,  -5.56000,   6.95000,   0.14000,   5.27000,  -4.89000],
                [ -1.66000,  -5.24000,   9.26000,   1.60000,  -9.08000,   5.47000,  -1.37000,  -6.80000,  -2.41000],
                [  2.75000,  -1.66000,   3.40000,  -0.97000,   2.36000,   9.66000,  -4.66000,  -1.32000,   5.57000],
                [ -7.60000,  -9.67000,   3.05000,   7.35000,  -4.94000,   3.36000,   8.05000,   5.18000,   6.37000],
                [  2.46000,  -6.48000,  -8.79000,  -7.16000,   4.17000,  -0.47000,  -2.24000,   8.06000,   8.58000],
                [ -8.22000,   5.29000,  -8.27000,  -4.63000,  -2.48000,  -5.90000,  -1.46000,   0.23000,   7.35000]]
Matrix_0045 = Matrix.fromArray(array2D_0046)
array2D_0048 = [[ 0.36410, -0.26039,  0.32249,  0.17916,  0.02477, -0.00336,  0.16633,  0.03902,  0.40282],
                [-0.47949, -0.29326, -0.22151, -0.19392,  0.32724, -0.22928, -0.29791,  0.10484, -0.08964],
                [-0.21898,  0.41747,  0.18997,  0.15633, -0.19286, -0.06960,  0.12432,  0.70155, -0.20569],
                [-0.07128, -0.27898, -0.21983, -0.20132,  0.01972,  0.26594,  0.61631,  0.19807, -0.26143],
                [ 0.13539,  0.17571, -0.03932, -0.37231,  0.30487, -0.26885,  0.57577, -0.09986,  0.09151],
                [ 0.38769, -0.08123, -0.05347,  0.03510, -0.24088, -0.59632, -0.02994,  0.25224,  0.01470],
                [-0.19487, -0.25503,  0.33572, -0.40904, -0.35569, -0.33420, -0.01650, -0.20279, -0.03180],
                [ 0.07231, -0.37983, -0.37473,  0.36637, -0.41262,  0.12935,  0.11820, -0.10801, -0.17632],
                [-0.08513, -0.21660,  0.44580, -0.29155, -0.35303,  0.25755,  0.00645,  0.07917, -0.14307],
                [ 0.14103, -0.08149,  0.11069, -0.23764,  0.15196,  0.46403, -0.18503,  0.27975,  0.34954],
                [-0.38975, -0.38190,  0.29619,  0.46257,  0.26907, -0.16549,  0.23224,  0.13068,  0.23257],
                [ 0.12615, -0.28559, -0.37772, -0.25962, -0.05891, -0.07708, -0.14933,  0.46411,  0.26915],
                [-0.42154,  0.25755, -0.25217, -0.03319, -0.42823,  0.02414,  0.17640, -0.11266,  0.63967]]
Matrix_0047 = Matrix.fromArray(array2D_0048)
Array2D_0050 = [[ 19.49992,   1.79945,   7.27516,   3.25321,   1.53347,   4.07025,   0.29531,   1.14910,   3.68929],
                [  0.00000,  23.48455,   2.85065,   0.44384,   4.28935, -10.69326,  -3.51674,  -2.65921,  -1.66700],
                [  0.00000,   0.00000,  23.54579,  11.87641, -11.30280,   2.47536,   3.22837, -10.82935,   2.57161],
                [  0.00000,   0.00000,   0.00000,  11.39234,  -6.77696,   5.67846,   2.37378,   5.56971,  -0.42111],
                [  0.00000,   0.00000,   0.00000,   0.00000,  14.04252,   2.52388,   3.40297,   3.12590,   0.24772],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  19.19385,  -6.71168,  -3.07046,  -7.16858],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  11.80463,   9.62107,   1.64737],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  12.70243,   7.85639],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  16.95046]]
Matrix_0049 = Matrix.fromArray(Array2D_0050)
# Matrix size = (9, 8)
array2D_0052 = [[  2.70000,  -1.30000,  -6.25000,   5.24000,   1.88000,  -4.45000,  -3.54000,  -9.66000],
                [ -6.15000,  -0.51000,  -1.67000,   4.45000,   9.16000,  -2.54000,   2.39000,  -0.56000],
                [ -7.01000,   6.26000,  -0.91000,   8.10000,   6.16000,   9.37000,  -9.65000,   8.60000],
                [ -7.97000,   2.77000,  -5.76000,   1.26000,  -0.04000,  -7.32000,   7.91000,   9.60000],
                [ -0.20000,   0.27000,  -3.93000,   4.29000,   0.30000,   5.55000,   3.40000,  -5.53000],
                [ -1.32000,   3.28000,   6.29000,   3.08000,   9.56000,   1.32000,  -2.84000,   7.79000],
                [ -5.36000,  -4.72000,  -2.71000,  -9.86000,   7.57000,  -4.73000,  -0.69000,   8.96000],
                [ -0.32000,  -4.16000,   4.87000,  -6.82000,  -5.08000,   6.48000,  -8.30000,  -3.45000],
                [  8.58000,   0.08000,  -3.17000,   7.48000,  -3.23000,  -0.52000,   8.60000,  -2.27000]]
Matrix_0051 = Matrix.fromArray(array2D_0052)
array2D_0054 = [[ 0.16681, -0.08889, -0.49650,  0.29683,  0.16871, -0.06491, -0.69978, -0.09541],
                [-0.37995, -0.15934, -0.08107,  0.76161,  0.00646, -0.30899,  0.20371, -0.05432],
                [-0.43308,  0.53146, -0.03404,  0.12665, -0.03526,  0.54645, -0.28296,  0.30512],
                [-0.49238,  0.15112, -0.39009, -0.28054, -0.28906, -0.43277,  0.13675,  0.23862],
                [-0.01236,  0.02469, -0.30092,  0.20477, -0.09843,  0.52774,  0.46863, -0.42500],
                [-0.08155,  0.31909,  0.48412,  0.22206,  0.52383, -0.16804,  0.06676,  0.06930],
                [-0.33114, -0.58452, -0.15578, -0.20292,  0.56326,  0.26519,  0.09310,  0.29489],
                [-0.01977, -0.43913,  0.38713,  0.24042, -0.53025,  0.19322, -0.15716,  0.35968],
                [ 0.53007,  0.15648, -0.30745,  0.21652,  0.06461,  0.00272,  0.33557,  0.66042]]
Matrix_0053 = Matrix.fromArray(array2D_0054)
Array2D_0056 = [[ 16.18655,  -2.68123,   1.47861,   2.11559, -10.71649,   0.75544,   3.92663, -14.51891],
                [  0.00000,   9.59392,   0.32707,  14.33789,   1.99956,   5.06998,   0.57209,   5.24086],
                [  0.00000,   0.00000,  13.02623,  -6.93360,   0.51491,   7.32608,  -9.34101,   4.20587],
                [  0.00000,   0.00000,   0.00000,   9.15997,   7.05329,   3.81961,  -2.59965,  -7.43956],
                [  0.00000,   0.00000,   0.00000,   0.00000,  11.89788,  -4.97033,   0.21775,   6.64293],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  12.06531,  -8.69703,  -1.17920],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  12.30618,   4.06835],
                [  0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   0.00000,   8.65904]]
Matrix_0055 = Matrix.fromArray(Array2D_0056)
# Matrix size = (6, 2)
array2D_0058 = [[ 0.99000,  3.86000],
                [-7.37000, -2.50000],
                [-7.30000, -3.76000],
                [-8.88000, -4.90000],
                [ 6.25000, -4.66000],
                [ 6.47000, -4.16000]]
Matrix_0057 = Matrix.fromArray(array2D_0058)
array2D_0060 = [[ 0.06043,  0.38569],
                [-0.44989, -0.15326],
                [-0.44562, -0.28480],
                [-0.54206, -0.38023],
                [ 0.38152, -0.57245],
                [ 0.39495, -0.52380]]
Matrix_0059 = Matrix.fromArray(array2D_0060)
Array2D_0062 = [[16.38184, 2.26874],
                [0.00000, 9.65252]]
Matrix_0061 = Matrix.fromArray(Array2D_0062)

GSMatrixList = [Matrix_0045, Matrix_0051, Matrix_0057]
GSExpected = [(Matrix_0047, Matrix_0049), (Matrix_0053, Matrix_0055), (Matrix_0059, Matrix_0061)]
GSArgs = [GSMatrixList, GSExpected]



