import maximumCut
import networkx as nx

problem = maximumCut.MaximumCutProblem.get_problem_from_input()
problem.prepare()
response = problem.sample_hybrid()
problem.print_result(response)
problem.print_min_energy(response)

'''
graph = nx.fast_gnp_random_graph(10,0.5)
problem = maximumCut.MaximumCutProblem(None,None,graph)
problem.prepare()
response = problem.sample_advantage(100)
problem.print_result(response)
'''

#maximumCut.MaximumCutProblem.test_advantage(100)

#graph = nx.fast_gnp_random_graph(10,0.5)
'''
graph = nx.Graph()
nodes = range(7)
graph.add_nodes_from(nodes)
edges = [(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (3, 4), (5, 6), (0, 6)]
graph.add_edges_from(edges)
maximumCut.MaximumCutProblem.solve_classically(graph)
'''