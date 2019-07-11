import numpy as np
import sympy
import itertools
from numpy.linalg import matrix_rank


class SimplexSolver:
	""" This class implements simplex algorithm to solve LPs """

	def solve(self, c, A, b):
		"""
		This method solves the std form LP min (c.T * x) s.t. Ax = b, x >= 0 using simplex algorithm.

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

		self.A = A
		self.c = c
		self.b = b

		# define some frequently used parameters
		indices = list(range(self.A.shape[1]))

		# get the initial BFS (we will keep around the list of basic indices or basis matrix B as our current solution, instead of x explicitly)
		init_bfs = self.get_initial_bfs()
		if init_bfs is None:
			return -1
		(_, basic_indices) = init_bfs
		B = self.A[:, list(basic_indices)]
		optimal = False
		opt_infinity = False
		iteration_number = 0
		obj_val = float('inf')

		# main simplex body
		while not optimal:
			# print iteration number
			print('simplex: starting iteration #{}, obj = {}'.format(iteration_number, obj_val))
			iteration_number += 1

			# compute x_b, c_b, B_inv
			B_inv = np.linalg.inv(B)
			x_b = np.dot(B_inv, self.b)
			if (x_b == 0.0).any():
				print('simplex: alert! this bfs is degenerate')
			c_b = self.c[basic_indices]

			if iteration_number == 1:
				print('initial x_b = {}, with basic_indices = {}'.format(x_b, basic_indices))

			# compute obj_val just for display purposes
			obj_val = 0.0
			for i, b_i in enumerate(basic_indices):
				obj_val += (c[b_i] * x_b[i])

			# compute reduced cost in each non-basic j-th direction
			reduced_costs = {}
			for j in indices:
				if j not in basic_indices:
					# j is a non-basic index
					A_j = self.A[:, j]
					reduced_costs[j] = self.c[j] - np.dot(c_b.T, np.dot(B_inv, A_j))


			# check if this solution is optimal
			if (np.array(reduced_costs.values()) >= 0.0).all():
				# all reduced costs are >= 0.0 so this means we are at optimal already
				optimal = True
				break


			# this solution is not optimal, go to a better neighbouring BFS
			chosen_j = None
			for j in reduced_costs.keys():
				if reduced_costs[j] < 0.0:
					chosen_j = j
					break

			d_b = -1.0 * np.dot(B_inv, self.A[:, chosen_j])
			# check if optimal is infinity
			if (d_b >= 0).all():
				# optimal is -infinity
				opt_infinity = True
				break

			# calculate theta_star and the exiting index l
			l = None
			theta_star = None
			for i, basic_index in enumerate(basic_indices):
				if d_b[i] < 0:
					if l is None:
						l = i
					if -x_b[i]/d_b[i] < -x_b[l]/d_b[l]:
						l = i
						theta_star = -x_b[i]/d_b[i]

			# form new solution by replacing basic_indices[l] with chosen_j
			basic_indices[l] = chosen_j
			basic_indices.sort()
			B = self.A[:, list(basic_indices)]

		
		if opt_infinity:
			print('Optimal is inifinity')
			return -1

		if not optimal:
			print('optimal not found')
			return -1

		# return solution
		x = np.zeros(shape=(self.A.shape[1], ))
		for i in range(x.shape[0]):
			if i in basic_indices:
				x[i] = x_b[basic_indices.index(i)]
		return x

	
	def get_initial_bfs(self):
		"""
		This is a helper method used by solve() method to compute the initial basic feasible solution (BFS) required by simplex algorithm.
		It uses the auxiliary LP technique to compute such a BFS.

		Parameters:
			None

		Returns:
			None if the original LP is infeasible, else
			Tuple (B, basic_indices) where
				-> B 			(np array): basis matrix of the BFS
				-> basic_indices    (list): list of basic indices of the BFS
		"""
		
		M = self.A.shape[0]
		N = self.A.shape[1]

		# new constraint matrix A_ and c, b vector (b >= 0 must hold, so multiply by -1 if not already)
		A_positive = np.copy(self.A)
		b = np.copy(self.b)
		for i in range(M):
			if b[i] < 0.0:
				b[i] = -1.0 * b[i]
				A_positive[i, :] = -1 * A_positive[i, :]

		A_ = np.concatenate((A_positive, np.eye(M)), axis=1)
		c = np.zeros(shape=(N + M, ))
		for i in range(N, N + M):
			c[i] = 1.0
		indices = list(range(N + M))

		# variables: [ x_0, ..., x_N-1, y_0, ..., y_M-1 ]
		# init bfs of aux problem is x = 0, y = b
		# so basic initial indices for aux problem are (N, ..., N + M - 1)
		basic_indices = list(range(N, N + M))
		B = A_[:, basic_indices]

		optimal = False
		opt_infinity = False
		iteration_number = 0
		obj_val = float('inf')

		# main simplex body
		while not optimal:

			# print iteration number
			print('get_init_bfs_aux: starting iteration #{}, obj = {}'.format(iteration_number, obj_val))
			iteration_number += 1

			# compute x_b, c_b, B_inv
			B_inv = np.linalg.inv(B)
			x_b = np.dot(B_inv, b)
			if (x_b == 0.0).any():
				print('get_init_bfs_aux: alert! this bfs is degenerate')
			c_b = c[basic_indices]

			# compute obj_val just for display purposes
			obj_val = 0.0
			for i, b_i in enumerate(basic_indices):
				obj_val += (c[b_i] * x_b[i])

			# compute reduced cost in each non-basic j-th direction
			reduced_costs = {}
			for j in indices:
				if j not in basic_indices:
					# j is a non-basic index
					A_j = A_[:, j]
					reduced_costs[j] = c[j] - np.dot(c_b.T, np.dot(B_inv, A_j))


			# check if this solution is optimal
			if (np.array(reduced_costs.values()) >= 0.0).all():
				# all reduced costs are >= 0.0 so this means we are at optimal already
				optimal = True
				break

			# this solution is not optimal, go to a better neighbouring BFS
			chosen_j = None
			for j in reduced_costs.keys():
				if reduced_costs[j] < 0.0:
					chosen_j = j
					break

			d_b = -1.0 * np.dot(B_inv, A_[:, chosen_j])
			# check if optimal is infinity
			if (d_b >= 0).all():
				# optimal is -infinity
				opt_infinity = True
				break

			# calculate theta_star and the exiting index l
			l = None
			theta_star = None
			for i, basic_index in enumerate(basic_indices):
				if d_b[i] < 0:
					if l is None:
						l = i
					if -x_b[i]/d_b[i] < -x_b[l]/d_b[l]:
						l = i
						theta_star = -x_b[i]/d_b[i]

			# form new solution by replacing basic_indices[l] with chosen_j
			basic_indices[l] = chosen_j
			basic_indices.sort()
			B = A_[:, list(basic_indices)]

		
		if obj_val != 0.0:
			print('get_init_bfs_aux: the original problem is infeasible!')
			return None

		# if basic_indices contains no artifical variables, return that
		contains_artifical = False
		for x in basic_indices:
			if x >= N:
				contains_artifical = True
				break
		if not contains_artifical:

			assert len(basic_indices) == M, 'assertion failed, please check this'
			assert matrix_rank(B) == M, 'this should have been equal, assertion failed'
			x_b = np.dot(np.linalg.inv(B), b)
			assert (x_b >= 0.0).all(), 'this does not give a feasible solution, something is wrong, assertion failed'
			print('init_bfs_aux: assertions passed, no artificial vars in basis by chance! found a valid init bfs in {} iterations'.format(iteration_number))
			basic_indices.sort()
			
			return (B, basic_indices)

		# basis contains artificial variables
		basic_indices_no_artificial = []
		for index in basic_indices:
			if index < N:
				basic_indices_no_artificial.append(index)	

		# now have to choose columns from A that are linearly independent to the current selection of basis indices
		counter = 0
		while len(basic_indices_no_artificial) < M:
			if counter in basic_indices_no_artificial:
				continue

			# check if counter-th column of A is linearly independent with current selection of indices
			B_small = self.A[:, basic_indices_no_artificial]
			B_test = np.concatenate((B_small, self.A[:, counter]), axis=1)
			if matrix_rank(B_test) == min(B_test.shape[0], B_test.shape[1]):
				# is l.i., so take this column
				basic_indices_no_artificial.append(counter)

			counter += 1

		# test if what we got is indeed a BFS to original problem
		basic_indices = basic_indices_no_artificial
		basic_indices.sort()
		B = self.A[:, basic_indices]

		assert len(basic_indices) == M, 'assertion failed, please check this'
		assert matrix_rank(B) == M, 'this should have been equal, assertion failed'
		x_b = np.dot(np.linalg.inv(B), b)
		assert (x_b >= 0.0).all(), 'this does not give a feasible solution, something is wrong, assertion failed'
		print('init_bfs_aux: assertions passed! found a valid init bfs in {} iterations'.format(iteration_number))

		return (B, basic_indices)