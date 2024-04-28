import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 
import networkx as nx
import dimod
import matplotlib.pyplot as plt


def find_all_partitions(node_list):
    if len(node_list) == 0:
        return [[]]
    partitions = []
    node = node_list.pop()
    smaller_partitions = find_all_partitions(node_list)
    for subset in smaller_partitions:
        partitions.append(subset)
        partitions.append(subset + [node])
    return partitions


class MaximumCutProblem:
    """
    This class represents a Maximum Cut problem.
    It takes a graph G = {V,E} as input and partitions the set
    of vertices V in two complementary subsets S and T such
    that the number of edges between S and T is as large as possible
    """
    def __init__(self, edges, vertices_num):
        """
        Constructor of the MaximumCutProblem class.
        Params:
        edges: edges of the graph
        vertices_num: number of vertices in the graph
        """
        self.e = edges
        self.v = vertices_num
        self.q = defaultdict(int)

    def __init__(self, graph):

        if not isinstance(graph,  nx.classes.graph.Graph):
            raise TypeError("A networkx graph is required")

        self.v = graph.number_of_nodes()
        edges = []
        for u,v in graph.edges(data=False):
            edges.append([int(u),int(v)])
        self.e = edges
        self.q = defaultdict(int)

    def prepare(self):
        for i in range(len(self.e)):
            self.q[(self.e[i][0],self.e[i][0])] += -1
            self.q[(self.e[i][1],self.e[i][1])] += -1
            self.q[(self.e[i][0],self.e[i][1])] += 2 
    
    def sample_advantage(self, num_of_reads, chain_strength = None):

        sampler = EmbeddingComposite(DWaveSampler())

        if chain_strength == None:
            chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(self.q, offset = 0.0), sampler)

        print("Computing results on advantage...")

        sample_set = sampler.sample_qubo(self.q,
                               chain_strength=chain_strength,
                               num_reads=num_of_reads,
                               label='Maximum Cut')
        return sample_set
    
    def sample_hybrid(self):

        sampler = LeapHybridSampler()

        sample_set = sampler.sample_qubo(self.q, label="Maximum_Cut")

        return sample_set
    
    def sample_2000Q(self, num_of_reads, chain_strength = None):

        sampler = EmbeddingComposite(DWaveSampler(solver={'topology__type': 'chimera'}))

        if chain_strength == None:
            chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(self.q, offset = 0.0), sampler)

        print("Computing results on advantage...")

        sample_set = sampler.sample_qubo(self.q,
                               chain_strength=chain_strength,
                               num_reads=num_of_reads,
                               label='Maximum Cut')
        return sample_set

    def get_problem_from_input():
            print("Insert the number of vertices: ")
            n_vertices = input()
            if( int(n_vertices) <= 0):
                raise ValueError("The number of vertices must be positive")
            print("Insert the number of edges: ")
            n_edges = input()
            edges = []
            for i in range(int(n_edges)):
                 flag = 0
                 while(flag == 0):
                    print("Insert the first vertex of the edge: ")
                    first_vertex = int(input())
                    print("Insert the second vertex of the edge: ")
                    second_vertex = int(input())
                    if [first_vertex,second_vertex] in edges:
                        print("Edge already inserted, try again")
                    else:
                        edges.append([first_vertex,second_vertex])
                        flag = 1

            return MaximumCutProblem(edges,n_vertices)
    
    def test_advantage(number_of_nodes):
        graph = nx.fast_gnp_random_graph(number_of_nodes,0.5)
        problem = MaximumCutProblem(graph)
        problem.prepare()
        response = problem.sample_advantage(100)
        problem.print_result(response)

    def test_hybrid(number_of_nodes):
        graph = nx.fast_gnp_random_graph(number_of_nodes,0.5)
        problem = MaximumCutProblem(graph)
        problem.prepare()
        response = problem.sample_hybrid()
        problem.print_result(response)

    def solve_classically(graph):
        if not isinstance(graph,  nx.classes.graph.Graph):
            raise TypeError("A networkx graph is required")
        nx.draw(graph)
        plt.show()
        max_cut_size = 0
        max_partition = None
        partitions = find_all_partitions(sorted(graph.nodes))
        for partition in partitions:
            cut_size = nx.algorithms.cuts.cut_size(graph, partition)
            if cut_size > max_cut_size:
                max_cut_size = cut_size
                max_partition = partition

        print("The solution is: ")

        print("{", end = "")
        sorted_nodes = sorted(graph.nodes)
        for node in sorted_nodes:
            if node in max_partition and node+1 in sorted_nodes:
                print("%s: %s" % (node, 1), end = ", ")
            elif node in max_partition and node+1 not in sorted_nodes:
                print("%s: %s" % (node, 1), end = "")
            elif node not in max_partition and node+1 in sorted_nodes:
                print("%s: %s" % (node, 0), end = ", ")
            elif node not in max_partition and node+1 not in sorted_nodes:
                print("%s: %s" % (node, 0), end = "")
        print("}")         

    def print_result(self,response):
        lut = response.first.sample
        for i in range(1, int(self.v) + 1, 1):
            if i not in lut:
                lut[i] = 0

        print("The solution is: ")

        print("{", end = "")
        for key in sorted(lut):
            if key + 1 in lut:
                print("%s: %s" % (key, lut[key]), end = ", ")
            else:
                print("%s: %s" % (key, lut[key]), end = "")  
        print("}")

    def print_min_energy(self,response):
        max_energy = response.first.energy
        print("minimum energy: " + str(max_energy))

