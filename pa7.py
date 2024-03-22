from helpers import gram_schmidt
from structures import Vec, Matrix
import numpy as np
import cmath

# ----------------------- PROBLEM 1 ----------------------- #
def qr_solve(A: Matrix, b: Vec):
  """
  Solves the system of equations Ax = b by using the
  QR factorization of Matrix A
  :param A: Matrix of coefficients of the system
  :param b: Vec of constants
  :return: Vec solution to the system
  """
  # Constructing U
  # U should be the set of orthonormal vectors returned
  # by applying the Gram-Schmidt Process to the columns of A
  m, n = A.dim()
  cols_of_A = []
  for i in range(m):
      cols_of_A.append(Vec(A.get_col(i + 1)))
  U = gram_schmidt(cols_of_A)  # FIXME: Replace with the appropriate line
  n = len(U)

  # Constructing Q
  # Q should be the matrix whose columns are the elements
  # of the vector in set U
  Q = Matrix([[None for j in range(n)] for i in range(n)])
  for j in range(n):
      Q.set_col(j + 1, U[j].elements)  # FIXME: Replace with the appropriate line

  # Constructing R
  R = Q.transpose() * A  # FIXME: Replace with the appropriate line

  # Constructing the solution vector x
  b_star = Vec(Q.transpose() * b)
  x = [None for i in range(n)]
  for i in range(n - 1, -1, -1):
      summation = b_star[i]
      for j in range(i + 1, n):
          summation -= R.rowsp[i][j] * x[j]
      x[i] = summation / R.rowsp[i][i]
  # FIXME: find the components of the solution vector
  #        and replace them into elements of x
  return Vec(x)

# ----------------------- PROBLEM 2 ----------------------- #
def _submatrix(A: Matrix, i: int, j: int):
  submatrix = []
  for row in (A.rowsp[:i - 1] + A.rowsp[i:]):
      modified_row = row[:j - 1] + row[j:]
      submatrix.append(modified_row)

  result_matrix = Matrix(submatrix)
  return result_matrix

# ----------------------- PROBLEM 3 ----------------------- #
def determinant(A: Matrix):
  m, n = A.dim()
  if m != n:
      return 0
  if n == 1:
      return A.get_entry(1, 1)  # FIXME: Return the correct value
  elif n == 2:
      return A.get_entry(1, 1) * A.get_entry(2, 2) - A.get_entry(1, 2) * A.get_entry(2, 1)  # FIXME: Return the correct value
  else:
      det = 0
      for j in range(1, n + 1):
          det += ((-1) ** (j + 1)) * A.get_entry(1, j) * determinant(_submatrix(A, 1, j))

      return det

# ----------------------- PROBLEM 4 ----------------------- #
def eigen_wrapper(A: Matrix):
  eigenvalues, eigenvectors = np.linalg.eig(A.rowsp)
  eigen_dict = {}
  for val, vec in zip(eigenvalues, eigenvectors.T):
      eigen_dict[val] = Vec(vec.tolist())

  return eigen_dict

# ----------------------- PROBLEM 5 ----------------------- #
def svd(A: Matrix):
  m, n = A.dim()
  aTa = A.transpose() * A
  eigen = eigen_wrapper(aTa)
  eigenvalues = np.sort_complex(list(eigen.keys())).tolist()[::-1]

  # Constructing V
  # V should be the mxm matrix whose columns
  # are the eigenvectors of matrix A.transpose() * A
  #eigenvectors = list(eigen.values())
  V = Matrix([[None for j in range(n)] for i in range(n)])
  for j in range(1, n + 1):
    eigenvec = eigen[eigenvalues[j - 1]].elements
    V.set_col(j, eigenvec)

  # Constructing Sigma
  # Sigma should be the mxn matrix of singular values.
  singular_values = [np.sqrt(eigenvalues[x]) for x in range(len(eigenvalues))]
  Sigma = Matrix([[0 for j in range(n)] for i in range(m)])
  for i in range(1, m + 1):
    if i <= n:
      Sigma.set_entry(i, i, singular_values[i - 1])
  
  #singular_values.sort(reverse=True)
  #Sigma = Matrix([[singular_values[j] if i == j else 0 for j in range(n)] for i in range(m)])

  # Constructing U
  # U should be the matrix whose j-th column is given by
  # A * vj / sj where vj is the j-th eigenvector of A.transpose() * A
  # and sj is the corresponding j-th singular value
  U = Matrix([[0 for j in range(m)] for i in range(m)])
  for j in range(1, m + 1):
    if j <= n:
      vj = eigen[eigenvalues[j - 1]]
      sj = singular_values[j - 1]
      
      uj = 1/sj * (A * vj)
      
      U.set_col(j, list(uj))
      
  return (U, Sigma, V)
  
