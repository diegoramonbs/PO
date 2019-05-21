#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego

from docplex.mp.model import Model
from docplex.mp.environment import Environment
from docplex.util.environment import get_environment
from collections import namedtuple
from operator import itemgetter
from readB2P import read_input_b2p
import os



def build_model(data, width, height, timelimit=2*60*60):
	
	# Model
	mdl = Model(name='BinPacking2D')

	number_of_items = len(data)

	# Space of indexes for each variable
	X = [(i, j) for i in range(1, number_of_items+1) for j in range(1,number_of_items+1)]
	Y = [i for i in range(1, number_of_items+1)]
	Z = [(i, j) for i in range(1, number_of_items+1) for j in range(1,number_of_items+1)]
	Q = [i for i in range(1, number_of_items+1)]

	x = mdl.binary_var_dict(X, name="x")
	y = mdl.binary_var_dict(Y, name="y")
	z = mdl.binary_var_dict(Z, name="z")
	q = mdl.binary_var_dict(Q, name="q")

	for j in range(1, number_of_items+1):
		mdl.add_constraint(
			mdl.sum(
				x[(i,j)] for i in range(1, j)
				) + y[j] == 1
			)

	for i in range(1, number_of_items+1):
		mdl.add_constraint(
			mdl.sum(
				z[(k,i)] for k in range(1, i) 
			) + q[i] == y[i]
		)
    
	for i in range(1, number_of_items):
		mdl.add_constraint(
			mdl.sum(
				data[j-1][0] * x[(i,j)] for j in range(i+1, number_of_items+1)
			) <= (width - data[i-1][0])*y[i]
		)

	for k in range(1, number_of_items):
		mdl.add_constraint(
			mdl.sum(
				data[i-1][1] * z[(k,i)] for i in range(k+1, number_of_items+1)
			) <= (height - data[k-1][1])*q[k]
		)


	# Minimize cost  q1 + q2 + .... + qn
	mdl.minimize(
		mdl.sum(
			q[i] for i in Q
		)
	)


	mdl.parameters.timelimit = timelimit
	return mdl

def main(instances, verbose=False):
	env = Environment()
	env.print_information()

	dir_path = os.path.dirname(os.path.realpath(__file__))

	solutions = []

	for index, instance in enumerate(instances):

		data, h, w, problem_class = instance 

		data = sorted(data, key=itemgetter(0), reverse=True)

		try:
			mdl = build_model(data, w, h)

			if verbose:
				print(mdl.export_to_string())
				mdl.print_information()

			mdl.export_as_lp(dir_path + "/" + mdl.name + "_" + str(index) + '.lp')
			if mdl.solve(log_output=verbose):
				if verbose:
					mdl.float_precision = 3
					print("* model solved as function:")
					mdl.print_solution()
					mdl.report_kpis()
					
				# Save the CPLEX solution as "solution.json" program output
				with get_environment().get_output_stream(dir_path + "/" + mdl.name + "_" + str(index) + "_solution.json") as fp:
					mdl.solution.export(fp, "json")
			else:
				print("* model has no solution")

		except Exception as e:
			# Problem with more than 1000 variables
			print(e)



if __name__ == '__main__':


	import argparse

	parser = argparse.ArgumentParser(description='BinPacking2D')
	parser.add_argument('-i', '--input', help='Input file name')
	parser.add_argument('-v', '--verbose', help='Enable verbose', action='store_true')


  	args = parser.parse_args()

  	if args.input:
  		instances = read_input_b2p(args.input)
  		main(instances, verbose=args.verbose)
  		


	

