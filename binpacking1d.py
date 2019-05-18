from docplex.mp.model import Model
from docplex.util.environment import get_environment
from collections import namedtuple



DATA = (2, 5, 4, 7, 1, 3, 8)
number_of_items = 7
capacity = 10


def build_model(items, timelimit=2*60*60):
	mdl = Model('BinPacking1D')

	N = [i for i in range(1, number_of_items+1)]
	M = [(i, j) for i in range(1, number_of_items+1) for j in range(1, number_of_items+1)]

	x = mdl.binary_var_dict(M, name="x")
	y = mdl.binary_var_dict(N, name="y")


	# Objetive function
	mdl.minimize(
		mdl.sum(y[i] for i in N)
		)

	# Constraints
	for j in range(1, number_of_items+1):
		mdl.add_constraint(
			mdl.sum(
				x[(i, j)] for i in range(1, number_of_items+1)
			) == 1
		)

	for i in range(1, number_of_items+1):
		mdl.add_constraint(
			mdl.sum( 
				items[j-1] * x[(i,j)] for j in range(1, number_of_items+1)
			) <= capacity * y[i])
	
	return mdl

if __name__ == '__main__':
	mdl = build_model(DATA)
	print(mdl.export_to_string())
	mdl.print_information()
	mdl.export_as_lp()
	solution = mdl.solve(log_output=True)
	if solution:
		mdl.float_precision = 3
		print("* model solved as function:")
		mdl.print_solution()
		mdl.report_kpis()

		# Save the CPLEX solution as "solution.json" program output
		with get_environment().get_output_stream("solution.json") as fp:
			mdl.solution.export(fp, "json")
	else:
		print("* model has no solution")
