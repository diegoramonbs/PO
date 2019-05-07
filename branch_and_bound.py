#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego

from simplex import Simplex, format_log 
import numpy as np

def branch_and_bound(A, b, c, s, type, max_node=500, verbose=True):
	node_count = 0
	
	#TODO: add support to minimization problem
	best_feasible = -np.inf
	
	q = []
	q.append(Simplex(A, b, c, s, type))

	while len(q) > 0: 
		node_count += 1
		node = q.pop(0)
		
		# Execute simplex
		node.solve()

		if node_count >= max_node:
			q = []
			print('Maximum number of nodes reached.')
			return

		if verbose:
			print("Node <id={}, objective={}, solution={}>"
				.format(node_count, 
					    node.get_primal_objetive(),
					    node.get_primal_solution()))

		if node.is_infeasible():
			if verbose:
				print("Node <id={}, solution=infeasible>".format(node_count))
			continue 

		if node.get_primal_objetive() < best_feasible:
			if verbose:
				print("Node <id={}, solution=worst>".format(node_count))
			continue 

		index, value = None, None
		for i, v in enumerate(node.get_primal_solution()):
			if v.is_integer():
				continue 
			else:
				index, value = i, v

		if index is None:
			if node.get_primal_objetive() > best_feasible:
				best_feasible = node.get_primal_objetive()

		else:
			a, b, c, s, t = node.get_model()
			an = np.append(a, np.zeros((1, a.shape[1]), dtype=a.dtype), axis=0)
			an[-1, index] = 1

			rb = np.append(b, np.floor(value))
			rs = np.append(s, "<=")

			lb =  np.append(b, np.ceil(value))
			ls =  np.append(s, "=>")

			right_node = Simplex(an, lb, c, ls, t)
			right_node.set_primal_objetive(node.get_primal_objetive())
			q.append(right_node)

			left_node = Simplex(an, rb, c, rs, t)
			left_node.set_primal_objetive(node.get_primal_objetive())
			q.append(left_node)

	print("Optimal solution: {}".format(best_feasible))



if __name__ == '__main__':


	A = np.array([[1,1],[2,0]])
	b = np.array([4,5])
	c = np.array([2,1])
	s = np.array(["<=", "<="])

	"""A = np.array([[10,7],[1,1]])
	b = np.array([40,5])
	c = np.array([17,12])
	s = np.array(["<=", "<="])

	A = np.array([[8000,4000],[15,30]])
	b = np.array([40000,200])
	c = np.array([100, 150])
	s = np.array(["<=", "<="])

	A = np.array([[9, 5],[-4,5]])
	b = np.array([45, 5])
	c = np.array([10, 6])
	s = np.array(["<=", "<="])

	A = np.array([[1, 2],[5, 2]])
	b = np.array([10, 20])
	c = np.array([13, 8])
	s = np.array(["<=", "<="])

	A = np.array([[2, 6],[2, -2]])   #optimal 8
	b = np.array([15, 3])
	c = np.array([2, 3])
	s = np.array(["<=", "<="])

	A = np.array([[1, 2],[5, 2]])
	b = np.array([10, 20])
	c = np.array([13, 8])
	s = np.array(["<=", "<="])"""


	#format_log(Simplex(A, b, c, s, "max").solve())
	branch_and_bound(A, b, c, s, "max")
