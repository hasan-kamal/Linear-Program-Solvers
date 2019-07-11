import numpy as np
import itertools
from numpy.linalg import matrix_rank


class BruteSolver:
	""" This class implements brute-force exhaustive search (try all possible bases) to solve LPs """

	def solve(self, c, A, b):
		"""
		This method solves the std form LP min (c.T * x) s.t. Ax = b, x >= 0 using brute-force exhaustive search (try all possible bases).

		Parameters:
			c, A, b (np arrays): specify the LP in standard form
		
		Returns:
			-1 if LP is infeasible or optimal is (+-) infinity, else
			x (np array): solution to the LP
		"""

		# ensure dimensions are okay
		assert A.shape[0] == b.shape[0], 'first dims of A and b must match, check input!'
		assert A.shape[1] == c.shape[0], 'second dim of A must match first dim of c, check input!'

		# ensure A is full rank, drop redundant rows if not
		if matrix_rank(A) < min(A.shape[0], A.shape[1]):
			print('A is not full rank, dropping redundant rows')
			_, pivots = sympy.Matrix(A).T.rref()
			A = A[list(pivots)]
			print('Shape of A after dropping redundant rows is {}'.format(A.shape))

		# try all possible basis matrices (i.e. all m-combinations of basic indices) and take best
		m = A.shape[0]
		indices = list(range(A.shape[1]))
		opt_basis = None
		opt_val = float('inf')
		opt_xb = None

		# main loop body
		iteration_number = 0
		for basic_indices in itertools.combinations(indices, m):
			iteration_number += 1

			B = A[:, list(basic_indices)]
			if matrix_rank(B) != m:
					continue
			x_b = np.dot(np.linalg.inv(B), b)

			if (x_b < 0.0).any():
				# infeasible
				continue
			
			obj = 0.0
			for i, b_i in enumerate(basic_indices):
				obj += (c[b_i] * x_b[i])

			if obj < opt_val:
				opt_val = obj
				opt_basis = basic_indices
				opt_xb = x_b

		# show how many iterations it took
		print('brute took {} iterations'.format(iteration_number))

		# infeasible LP
		if opt_basis is None:
			return -1

		# return optimal solution
		x = np.zeros(shape=(A.shape[1], ))
		for i in range(x.shape[0]):
			if i in opt_basis:
				x[i] = opt_xb[opt_basis.index(i)]
		return x