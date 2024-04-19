import maximumCut
import networkx as nx

'''
problem = maximumCut.MaximumCutProblem.get_problem_from_input()
problem.prepare()
response = problem.sample_hybrid()
problem.print_result(response)
problem.print_min_energy(response)


graph = nx.fast_gnp_random_graph(10,0.5)
problem = maximumCut.MaximumCutProblem(graph)
problem.prepare()
response = problem.sample_advantage(100)
problem.print_result(response)
'''

maximumCut.MaximumCutProblem.test_advantage(90)