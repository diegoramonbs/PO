#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego

import numpy as np
import time

np.set_printoptions(precision=5)
np.set_printoptions(suppress=True)


class Simplex(object):
	def __init__(self, A, b, c, s=None, type="max", M=10**10):

		self.m, self.n = A.shape
		self.nn = 0
		self.A = A
		self.b = b 
		self.c = c
		self.s = s
		self.type = type
		self.M = M
		self.tableau = self.build_table(A, b, c, s=s, type=type)
		self.basis = [ self.n+i for i in range(self.m)] 
		self.history = {}
		self.tableau_start = self.tableau[:]
		self.is_solved = False
		self.log = {}


	def get_primal_objetive(self):
		if self.is_solved:
			return self.log['primal objective']

	def set_primal_objetive(self, value):
		self.log['primal objective'] = value

	def is_infeasible(self):
		if self.is_solved:
			return self.log['primal infeasible']

	def get_primal_solution(self):
		if self.is_solved:
			return self.log['primal solution']

	def get_model(self):
		return self.A, self.b, self.c, self.s, self.type 

	def get_log(self):
		return self.log

	def sensitivity_analysis(self, Am, bm, cm):
		pass

	def build_table(self, A, b, c, s, type):
		tableau  = np.zeros(shape=(self.m+1, self.m+self.n+1), dtype=np.float64)
		tableau[0, 0:self.n] = -c if type == "max" else -c * (-1)
		tableau[1:self.m+1,0:self.n] = A
		tableau[1:self.m+1, self.n:self.m+self.n] = np.eye(self.m)
		tableau[1:self.m+1, -1] = b

		if s is None:
			return tableau 

		for i, k in enumerate(s):
			if k == "<=" or k == "=<":
				continue
			elif k == "=":
				tableau[0, i+self.n] = self.M
				tableau[0,:] = tableau[0,:] - self.M*tableau[i+1]
			elif k == ">=" or k == "=>":
				tableau[0, i+self.n] = self.M
				tableau[0,:] = tableau[0,:] - self.M*tableau[i+1]
				tableau[i+1, i+self.n] = -1
				z = np.zeros((self.m+1, 1), dtype=np.float64)
				z[-1] = 1
				tableau  = np.column_stack((tableau , z))
				tableau[:,-2],  tableau [:,-1] = tableau[:,-1],tableau[:,-2].copy()
				tableau[0, i+self.n] = self.M
				self.n += 1
				self.nn += 1
			else:
				raise Exception("Invalid Symbol!")
		return tableau 

	def bland(self):
		q, l = -1, -1
		for i in range(0, self.m+self.n):
			z = self.tableau[0][i]
			if z > 0:
				continue
			z = abs(z)
			if z > l:
			    l = z
			    q = i
		return q

	def min_ratio(self, q):
		p = -1
		for i in range(1, self.m+1):
			# Avoid zero division
			if self.tableau [i][q] <= 0:
			    continue
			elif p == -1:
			    p = i
			elif self.tableau [i][-1]/ self.tableau[i][q] < self.tableau[p][-1] / self.tableau[p][q]:
			    p = i
		return p

	def pivot(self, p, q):
		for i in range(0, self.m+1):
			for j in range(0, self.m+self.n+1 ):
				if i != p and j != q:
					self.tableau[i][j] -= self.tableau[p][j] * self.tableau[i][q] / self.tableau[p][q]
		for i in range(0, self.m+1):
			if i != p:
				self.tableau [i][q] = 0.0
		for j in range(0, self.m+self.n+1 ):
			if j != q:
				self.tableau[p][j] /= self.tableau[p][q]
		self.tableau[p][q] = 1.0

	def is_not_optimal(self):
		return np.any(self.tableau[0, 0:-1] < 0)

	def solve(self):
		start = time.time()
		k = 0
		self.history[k] =(self.tableau.copy(), self.basis[:])

		while self.is_not_optimal():

			q = self.bland()
			p = self.min_ratio(q)

			if p == -1:
				raise Exception("Unbounded problem.")

			self.pivot(p, q)
			k += 1
			self.basis[p-1] = q
			self.history[k] = (self.tableau.copy(), self.basis[:])

		# Primal solution
		x = [None] * (self.n-self.nn)
		for i in range(self.m):
			if self.basis[i] < self.n-self.nn:
				x[self.basis[i]] = self.tableau[i+1][-1]

		# Ax <=> b and x > 0
		is_primal_feasible = False
		for ax, b, s in zip(np.dot(self.A, x), self.b, self.s):
			if s == "=":
				if ax != b: 
					is_primal_feasible = True
				break
			elif s == "<=" or s == "=<":
				if ax > b:
					is_primal_feasible = True
					break
			elif s == "=>" or s == ">=":
				if ax < b:
					is_primal_feasible = True
					break

		is_primal_feasible = np.all(x > 0) and not is_primal_feasible

		# Dual solution 
		y = self.tableau[0,self.n:-1]
		# Ay <=> c and y > 0
		is_dual_feasible = False
		for ax, b, s in zip(np.dot(self.A.T, y), self.c, self.s):
			if s == "=":
				if ax != b: 
					is_dual_feasible = True
				break
			elif s == "<=" or s == "=<":
				if ax > b:
					is_dual_feasible = True
					break
			elif s == "=>" or s == ">=":
				if ax < b:
					is_dual_feasible = True
					break

		is_dual_feasible = np.all(y > 0) and not is_dual_feasible
		total_time = time.time() - start


		self.is_solved = True 

		self.log = { 
			"primal objective" : self.tableau[0][self.m+self.n], 
			"dual objective" : self.tableau[0][self.m+self.n], 
			"primal infeasible": not is_primal_feasible, 
			"dual infeasible": 	not is_dual_feasible,
			"iterations" : k,
			"history": self.history,
			"basis": self.basis,
			"primal solution": x,
			"dual solution": y,
			"time execution": total_time,

		}
		return self.log


def print_matrix(g):
	print('\n'.join(['\t'.join(['{:+7.3e}'.format(item) for item in row])
        for row in g]))

def format_log(log):
	print("Primal infeasible: \t\t {}".format(str(bool(log["primal infeasible"]))))
	print("Dual infeasible: \t\t {}".format(str(bool(log["primal infeasible"]))))
	print("Primal objective: \t\t {}".format(log["primal objective"]))
	print("Primal solution: \t\t {}".format(log["primal solution"]))
	print("Dual objective: \t\t {}".format(log["dual objective"]))
	print("Dual solution:  \t\t {}".format(log["dual solution"]))
	print("Basis:          \t\t {}".format(map(lambda x: x+1, log["basis"])))
	print("Time execution:  \t\t {:.10f} seconds".format(log["time execution"]))


	for iteration, (tableau, basis) in log["history"].items():
		print("\nIteration: {}\n".format(iteration))
		print_matrix(tableau)
		print("Basis: {}".format(map(lambda x: x+1, basis)))



if __name__ == '__main__':

	# Exemplo mais geral, sendo um problema de minimização e com restrições de todos os tipos 
	A = np.array([[1,0],[0,2],[3,2]])
	b = np.array([4,12,18])
	c = np.array([3,5])
	s = np.array(["<=", "<=", "<="])
	log = Simplex(A, b, c, s, "max").solve()
	format_log(log)
