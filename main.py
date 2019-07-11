from solvers.simplex_solver import SimplexSolver
from solvers.interior_point_solver import InteriorPointSolver
from solvers.brute_solver import BruteSolver

import numpy as np

from utils import network_flow_to_std_LP

"""
This file contains 2 test cases for trying out implementation of LP solvers
Both tests are instances of max-flow problem
We run simplex algorithm to solve test 1
We run interior-point algorithm to solve test 2
"""

if __name__ == '__main__':

	# Test 1, max-flow problem, small graph
	print('\nTest #1')
	G = np.array([
			[0, 16, 13,  0,  0, 0 ],
			[0,  0, 10, 12,  0, 0 ],
			[0,  4,  0,  0, 14, 0 ],
			[0,  0,  9,  0,  0, 20],
			[0,  0,  0,  7,  0, 7 ],
			[0,  0,  0,  0,  0, 0 ]
		])

	# convert max-flow problem into std LP form
	c, A, b, id_to_edge = network_flow_to_std_LP(G, s=0, t=5)

	# solve LP using simplex algorithm
	solver = SimplexSolver()
	x = solver.solve(c, A, b)

	# show optimal objective value, max-flow of original problem and display optimal edge flows
	opt_obj = np.dot(c, x)
	print('Optimal objective value 	= {}'.format(opt_obj))
	print('Max flow from s -> t 		= {}'.format(-opt_obj))
	print('Optimal edge flows: \n')
	for i in id_to_edge.keys():
		a1, a2 = id_to_edge[i]
		print('Edge ({:1}, {:2}) = {:6.2f} / {:2}'.format(a1, a2, x[i], G[a1][a2]))



	# Test 2, max-flow problem, large graph
	print('\nTest #2')
	G = np.array([
			[0, 11, 15, 10, 0,  0, 0,  0,  0,  0, 0,  0],
			[0,  0,  0,  0, 0, 18, 4,  0,  0,  0, 0,  0],
			[0,  3,  8,  5, 0,  0, 0,  0,  0,  0, 0,  0],
			[0,  0,  0,  0, 6,  0, 0,  3, 11,  0, 0,  0],
			[0,  0,  0,  4, 0,  0, 0, 17,  6,  0, 0,  0],
			[0,  0,  0,  0, 3, 16, 0,  0,  0, 13, 0,  0],
			[0, 12,  0,  0, 4,  0, 0,  0,  0,  0, 0, 21],
			[0,  0,  0,  0, 0,  0, 0,  0,  4,  9, 4,  3],
			[0,  0,  0,  0, 0,  0, 0,  4,  0,  0, 5,  4],
			[0,  0,  0,  0, 0,  0, 0,  0,  0,  0, 7,  9],
			[0,  0,  0,  0, 0,  0, 0,  0,  2,  0, 0, 15],
			[0,  0,  0,  0, 0,  0, 0,  0,  0,  0, 0,  0]
		])

	# convert max-flow problem into std LP form
	c, A, b, id_to_edge = network_flow_to_std_LP(G, s=0, t=11)

	# solve LP using interior-point method
	solver = InteriorPointSolver()
	epsilon = 1e-4
	x = solver.solve(c, A, b, epsilon)

	# show optimal objective value, max-flow of original problem and display optimal edge flows
	opt_obj = np.dot(c, x)
	print('Optimal objective value 	= {}'.format(opt_obj))
	print('Max flow from s -> t 		= {}'.format(-opt_obj))
	print('Optimal edge flows: \n')
	for i in id_to_edge.keys():
		a1, a2 = id_to_edge[i]
		print('Edge ({:2}, {:2}) = {:6.2f} / {:2}'.format(a1, a2, x[i], G[a1][a2]))