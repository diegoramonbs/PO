from docplex.mp.model import Model
from docplex.util.environment import get_environment
from collections import namedtuple
import readCVRP


(nb_customers, truck_capacity, distance_matrix, 
	distance_warehouses, demands) = readCVRP.read_input_cvrp("instancias/P-VRP/P-n16-k8.vrp")


n_vertices = len(distance_matrix[0])


def build_model():
	mdl = Model('VehicleRoutingProblem')

	arcs = [(i, j) for i in range(n_vertices) for j in range(n_vertices) ]

	x = mdl.binary_var_dict(arcs, name="x")
	y = mdl.binary_var_dict(arcs, name="y")


	mdl.minimize(mdl.sumsq(distance_matrix[i][j] * x[(i, j)] for i, j in arcs ))

	for i in range(1, n_vertices):
		mdl.add_constraint(mdl.sum(x[(i,j)] for j in range(0,  n_vertices)) == 1)

	for i in range(0, n_vertices):
		mdl.add_constraint(mdl.sum(x[(i,j)] for j in range(0,  n_vertices)) - mdl.sum(x[(j,i)] for j in range(0,  n_vertices)) == 0)

	for i in range(1, n_vertices):
		mdl.add_constraint(mdl.sum(y[(i,j)] for j in range(0,  n_vertices) if j!=i) -
		 mdl.sum(y[(j,i)] for j in range(0,  n_vertices) if j!=i) == demands[i])

	mdl.add_constraints(y[(i,j)] <= truck_capacity*x[i,j] for i, j in arcs)


	return mdl


if __name__ == '__main__':
	mdl = build_model()
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