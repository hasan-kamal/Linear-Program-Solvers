### Following is what each file contains:

1. `simplex_solver.py` contains the actual code that solves standard form LP using simplex algorithm
2. `test.py` contains some random test cases which checks the validity of simplex solver by comparing against brute force search
3. `network_flow.py` contains code that given an adjacency matrix of a graph creates both the primal (flow) LP and its dual (which corresponds to min-cut) LP in standard form and runs simplex solver on them