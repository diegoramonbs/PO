#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego

import numpy as np 


class Model(object):

	def __init__(self, A=None, b=None, c=None, constraints=None, direction="max", M=10**10):
		self._A = A 
		self._b = b 
		self._c = c 
		self._constraints = constraints 
		self._direction = direction 
		self._M = M
		self._primal_objetive = None 
		self._dual_objetive = None 
		self._primal_variables = None 
		self._dual_variables = None
		self._is_primal_infeasible = None 
		self._is_dual_infeasible = None 
		self._basis = None 
		self._iterations = None 
		self._time_execution = None 
		self._history = []
		self._solved = False

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
	def primal_objetive(self):
		return self._primal_objetive 

	@primal_objetive.setter 
	def primal_objetive(self, x):
		self._primal_objetive = x

	@property
	def dual_objetive(self):
		return self._dual_objetive 

	@dual_objetive.setter 
	def dual_objetive(self, x):
		self._dual_objetive = x

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
		return self._is_dual_feasible 

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

	@property 
	def solved(self, x):
		self._solved = x


	def write(self, filename):
		pass 

	def read(self, filename):
		pass

	def __str__(self):
		pass 

	def __repr__(self):
		pass


if __name__ == '__main__':
	A = np.array([[1,0],[0,2],[3,2]])
	b = np.array([4,12,18])
	c = np.array([3,5])
	k = np.array(["<=", "<=", "<="])
	m = Model(A, b, c, k)







