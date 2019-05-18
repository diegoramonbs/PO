from docplex.mp.model import Model
from docplex.util.environment import get_environment
from collections import namedtuple

W = 15
H = 10


DATA = ( (0, 0),(10, 7), (9, 5), (7, 5), (7, 5), (3, 5), (3, 5), (10, 4), (5, 4), (5, 4), (8, 3) )

number_of_items = 10

W = 15
H = 10


def build_model(data, timelimit=2*60*60):
	mdl = Model('BinPacking2D')

	X = [(i, j) for i in range(1, number_of_items+1) for j in range(1,number_of_items+1)]
	Y = [i for i in range(1, number_of_items+1)]
	Z = [(i, j) for i in range(1, number_of_items+1) for j in range(1,number_of_items+1)]
	Q = [i for i in range(1, number_of_items+1)]

	x = mdl.binary_var_dict(X, name="x")
	y = mdl.binary_var_dict(Y, name="y")
	z = mdl.binary_var_dict(Z, name="z")
	q = mdl.binary_var_dict(Q, name="q")

	mdl.minimize(mdl.sum(q[i] for i in Q))

	for j in range(1, number_of_items+1):
		if j == 1:
			mdl.add_constraint(mdl.sum(y[1]) == 1)
		else:
			mdl.add_constraint(mdl.sum(x[(i,j)] for i in range(1, j))+ y[j] == 1)

	for i in range(1, number_of_items+1):
		if i == 1:
			mdl.add_constraint(mdl.sum(q[1]) == y[i])
		else:
			mdl.add_constraint(mdl.sum(z[(k,i)] for k in range(1, i) ) + q[i] == y[i])
    
	for i in range(1, number_of_items):
		mdl.add_constraint(mdl.sum(data[j][0] * x[(i,j)] for j in range(i+1, number_of_items+1)) <= (W - data[i][0])*y[i])

	for k in range(1, number_of_items):
		mdl.add_constraint(mdl.sum(data[i][1] * z[(k,i)] for i in range(k+1, number_of_items+1)) <= (H - data[k][1])*q[k])

	
	mdl.parameters.timelimit = timelimit
	return mdl

if __name__ == '__main__':
	mdl = build_model(DATA)
	print(mdl.export_to_string())
	mdl.print_information()
	mdl.export_as_lp()
	if mdl.solve(log_output=True):
		mdl.float_precision = 3
		print("* model solved as function:")
		mdl.print_solution()
		mdl.report_kpis()

		# Save the CPLEX solution as "solution.json" program output
		with get_environment().get_output_stream("solution.json") as fp:
			mdl.solution.export(fp, "json")
	else:
		print("* model has no solution")
