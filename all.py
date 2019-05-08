#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego

import time
import numpy as np
import copy


class Model(object):

	def __init__(self, A=None, b=None, c=None, constraints=None, direction="max", M=10**10):
		self._A = A 
		self._b = b 
		self._c = c 
		self._constraints = constraints 
		self._direction = direction 
		self._M = M
		self._primal_objective = None 
		self._dual_objective = None 
		self._primal_variables = None 
		self._dual_variables = None
		self._is_primal_infeasible = False
		self._is_dual_infeasible = False
		self._basis = None 
		self._iterations = None 
		self._time_execution = None 
		self._history = []
		self._solved = False
		self._unbounded = False

	@property
	def unbounded(self):
		return self._unbounded

	@unbounded.setter 
	def unbounded(self, x):
		self._unbounded = x

	@property
	def num_variables(self):
		if self._A is not None:
			return self._A.shape[1]

	@property
	def num_constraints(self):
		if self._A is not None:
			return self._A.shape[0]

	@property
	def A(self):
		return self._A 

	@A.setter 
	def A(self, x):
		self._A = x 

	@property
	def b(self):
		return self._b 

	@b.setter 
	def b(self, x):
		self._b = x 

	@property
	def c(self):
		return self._c

	@c.setter 
	def c(self, x):
		self._c = x

	@property
	def constraints(self):
		return self._constraints

	@constraints.setter 
	def constraints(self, x):
		self._constraints = x

	@property 
	def M(self):
		return self._M 

	@M.setter 
	def M(self, M):
		self._M = M 

	@property
	def direction(self):
		return self._direction 

	@direction.setter 
	def direction(self, x):
		self._direction = x

	@property
	def primal_objective(self):
		return self._primal_objective 

	@primal_objective.setter 
	def primal_objective(self, x):
		self._primal_objective = x

	@property
	def dual_objective(self):
		return self._dual_objective

	@dual_objective.setter 
	def dual_objective(self, x):
		self._dual_objective = x

	@property
	def primal_variables(self):
		return self._primal_variables 

	@primal_variables.setter 
	def primal_variables(self, x):
		self._primal_variables = x 

	@property
	def dual_variables(self):
		return self._dual_variables 

	@dual_variables.setter 
	def dual_variables(self, x):
		self._dual_variables = x 

	@property
	def is_primal_infeasible(self):
		return self._is_primal_infeasible

	@is_primal_infeasible.setter 
	def is_primal_infeasible(self, x):
		self._is_primal_infeasible = x 

	@property
	def is_dual_infeasible(self):
		return self._is_dual_infeasible

	@is_dual_infeasible.setter 
	def is_dual_infeasible(self, x):
		self._is_dual_infeasible = x 

	@property
	def basis(self):
		return self._basis

	@basis.setter 
	def basis(self, x):
		self._basis = x

	@property
	def iterations(self):
		return self._iterations

	@iterations.setter 
	def iterations(self, x):
		self._iterations = x

	@property
	def time_execution(self):
		return self._time_execution

	@time_execution.setter 
	def time_execution(self, x):
		self._time_execution = x

	@property
	def history(self):
		return self._history

	@history.setter 
	def history(self, x):
		self._history = x

	@property
	def solved(self):
		return self._solved 

	@solved.setter
	def solved(self, x):
		self._solved = x

	def write(self, filename):
		pass 

	def read(self, filename):
		pass

	def __str__(self):
		return "Model<primal_objective={}, primal_variables={}> ".format(
			self._primal_objective,
			self._primal_variables
			)

	def __repr__(self):
		pass




def simplex(model):

	start = time.time()

	m = model.num_constraints
	n = model.num_variables 

	basis = [n+i for i in range(m)]

	# Build the table for Ax <= b
	t = np.zeros(shape=(m+1, m+n+1), dtype=np.float64)
	t[0, 0:n] = -model.c if model.direction == "max" else -model.c * (-1)
	t[1:m+1,0:n] = model.A
	t[1:m+1, n:m+n] = np.eye(m)
	t[1:m+1, -1] = model.b

	n_extra = 0

	# Build the tableau for any type of constraint
	for index, const in enumerate(model.constraints):
		if const == "<=" or const == "=<":
			continue
		elif const == "=":
			t[0, index+n] = model.M
			t[0,:] = t[0,:] - model.M*t[index+1]
		elif const == ">=" or const == "=>":
			t[0, index+n] = model.M
			t[0,:] = t[0,:] - model.M*t[index+1]
			t[index+1, index+n] = -1
			z = np.zeros((m+1, 1), dtype=np.float64)
			z[-1] = 1
			t  = np.column_stack((t , z))
			t[:,-2],  t[:,-1] = t[:,-1],t[:,-2].copy()
			t[0, index+n] = model.M
			n += 1
			n_extra += 1
		else:
			raise Exception("Invalid symbol constraint!")

	count_iteration = 0
	history = {}

	# Put current tableau and basic variables in history, log only
	history[count_iteration] = (t, basis)

	# Optimality condition
	while np.any(t[0, 0:n] < 0):

    	# Bland's rule
		q, l = -1, -1 
		for i in range(0, m+n):
			z = t[0][i]
			if z > 0: 
				continue 
			z = abs(z)
			if z > l:
				l = z
				q = i 

		# Minimum Ratio Test
		p = -1
		for i in range(1, m+1):
			# Avoid zero division
			if t[i][q] <= 0:
				continue 
			elif p == -1:
				p = i 
			elif t[i][-1] / t[i][q] < t[p][-1] / t[p][q]:
				p = i 

		if p == -1:
			raise Exception("Unbounded problem.")
		
		# Pivot - gauss jordan elimination
		for i in range(0, m+1):
			for j in range(0, m+n+1 ):
				if i != p and j != q:
					t[i][j] -= t[p][j] * t[i][q] / t[p][q]
		for i in range(0, m+1):
			if i != p:
				t [i][q] = 0.0
		for j in range(0, m+n+1 ):
			if j != q:
				t[p][j] /= t[p][q]
		t[p][q] = 1.0

		# Update basic variables
		basis[p-1] = q 

		# Put current tableau and basic variables in history, log only
		count_iteration += 1
		history[count_iteration] = (t, basis[:])

	# Get the primal solution
	x = [0] * (n-n_extra) 
	for i in range(m):
		if basis[i] < n-n_extra:
			x[basis[i]] = t[i+1][-1]

	# Get the dual solution 
	y = t[0, n:-1]


	is_primal_feasible = False

	# Check if the primal solution is feasible, Ax <=> b and x > 0 
	for ax, b, k in zip(np.dot(model.A, x), model.b, model.constraints):
		if k == "=":
			if ax != b: 
				is_primal_feasible = True
				break
		elif k == "<=" or k == "=<":
			if ax > b:
				is_primal_feasible = True
				break
		elif k == "=>" or k == ">=":
			if ax < b:
				is_primal_feasible = True
				break

	is_primal_feasible = np.all(x > 0) and not is_primal_feasible

	# Check if the dual solution is feasible, Ay <=> b and y > 0 
	is_dual_feasible = False
	for ax, b, k in zip(np.dot(model.A.T, y), model.b, model.constraints):
		if k == "=":
			if ax != b: 
				is_dual_feasible = True
				break
		elif k == "<=" or k == "=<":
			if ax > b:
				is_dual_feasible = True
				break
		elif k == "=>" or k == ">=":
			if ax < b:
				is_dual_feasible = True
				break

		is_dual_feasible = np.all(y > 0) and not is_dual_feasible

	total_time = time.time() - start


	model.solved = True 
	model.primal_objective = t[0][m+n]
	model.dual_objective = t[0][m+n]
	model.primal_variables = x
	model.dual_variables = y
	model.basis = basis
	model.is_primal_infeasible = not is_primal_feasible
	model.is_dual_infeasible = not is_dual_feasible
	model.iterations = count_iteration
	model.history = history 
	model.time_execution = total_time

	return model
		

def branch_and_bound(model, max_node=500, verbose=True):
	count_node = 0
	best_feasible = -np.inf if model.direction == "max" else np.inf
	best_model = None

	queue = []
	queue.append(model)

	while len(queue) > 0:
		count_node += 1

		# Get a node from the queue and run simplex on it and 
		# Find an optimal solution to this relaxed problem
		node = simplex(queue.pop(0))

		if count_node >= max_node:
			print('Maximum number of nodes reached.')
			break

		# If current node has infeasible solution, the it does not branch from it
		if node.is_primal_infeasible:
			if verbose:
				print("Node <id={}, solution=infeasible>".format(count_node))
			continue 

		if verbose:
			print("Node <id={}, objective={}, solution={}>"
				.format(count_node, 
					    node.primal_objective,
					    node.primal_variables))

		if node.direction == "max":
			if node.primal_objective < best_feasible:
				if verbose:
					print("Node <id={}, solution=worst>".format(count_node))
				continue 
		else:
			if node.primal_objective > best_feasible:
				if verbose:
					print("Node <id={}, solution=worst>".format(count_node))
				continue 

		# If the current solution has only integer elements
		index, value = None, None
		for i, v in enumerate(node.primal_variables):
			if not v.is_integer():
				index, value = i, v
			else:
				continue

		# Then verify if it's better than current best solution.
		if index is None:
			if node.direction == "max":
				if node.primal_objective > best_feasible:
					best_feasible = node.primal_objective
					best_model = node
			else:
				if node.primal_objective > best_feasible:
					best_feasible = node.primal_objective
					best_model = node

		# If No, it branches (relaxes) the current node and puts it in the queue.
		else:

			An = np.append(node.A, np.zeros((1, node.A.shape[1]), dtype=node.A.dtype), axis=0)
			An[-1, index] = 1

			left = copy.deepcopy(node)
			left.A = An
			left.b = np.append(left.b, np.floor(value))
			left.constraints =  np.append(left.constraints, "<=")
			queue.append(left)

			right = copy.deepcopy(node)
			right.A  = An
			right.b = np.append(right.b, np.ceil(value))
			right.constraints =  np.append(right.constraints, "=>")
			queue.append(right)

	return best_model

def print_matrix(g):
	print('\n'.join(['\t'.join(['{:+7.3e}'.format(item) for item in row])
        for row in g]))


def format_log(model):
	print("Primal infeasible: \t\t {}".format(str(model.is_primal_infeasible)))
	print("Primal objective: \t\t {}".format(model.primal_objective))
	print("Primal solution: \t\t {}".format(model.primal_variables))
	print("Basis:          \t\t {}".format(map(lambda x: x+1, model.basis)))
	print("Time execution:  \t\t {:.10f} seconds".format(model.time_execution))


	for iteration, (tableau, basis) in model.history.items():
		print("\nIteration: {}\n".format(iteration))
		print_matrix(tableau)
		print("Basis: {}".format(map(lambda x: x+1, basis)))

if __name__ == '__main__':
	model = Model()

	model.A = np.array([[10,7],[1,1]])
	model.b = np.array([40,5])
	model.c = np.array([17,12])
	model.constraints = np.array(["<=", "<="])

	model = branch_and_bound(model, verbose=False)

	print(model)
