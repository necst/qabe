import graphColoring

import networkx as nx
import matplotlib.pyplot as plt





# Create a graph object
G = nx.Graph()

# Add vertices
num_vertices = 5
G.add_nodes_from(range(1, num_vertices + 1))

# Add edges
edges = [
    (2, 1),
    (1, 5),
    (2, 4),
    (2, 5),
    (5, 4),
    (4, 3),
    (3, 2)
]

G.add_edges_from(edges)



problem = graphColoring.GraphColoringProblem(3,None,None,G)
#problem.prepare()
problem.solve_classically()
#response = problem.sample_advantage(100)
#problem.print_result(response)



#problem = graphColoring.GraphColoringProblem.get_problem_from_input()
#problem.prepare()



#graphColoring.GraphColoringProblem.test_advantage(5,3)