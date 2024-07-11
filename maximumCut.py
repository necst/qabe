import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 
import networkx as nx
import dimod
import time


def find_all_partitions(node_list):
    """
    Given a list of nodes returns all the subsets of the node's set
    """
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
    def __init__(self, edges=None, vertices_num=None, graph = None):
        """
        Constructor of the MaximumCutProblem class.

        Params:

        edges: edges of the graph in the form of an array of arrays
                e.g. [[i,j],[j,k],...] where i-j and j-k are two edges
                connecting the nodes i,j and the nodes j,k
        vertices_num: number of vertices in the graph

        or 

        graph: a NetworkX graph
        """

        if edges == None and vertices_num == None and graph != None:
            if not isinstance(graph,  nx.classes.graph.Graph):
                raise TypeError("A Networkx graph is required")
            
            self.v = graph.number_of_nodes()
            edges = []
            for u,v in graph.edges(data=False):
                edges.append([int(u),int(v)])
            self.e = edges
            self.q = defaultdict(int)

        elif edges != None and vertices_num != None and graph == None:
            self.e = edges
            self.v = vertices_num
            self.q = defaultdict(int)
        else:
            raise TypeError("Specify the arguments in one of the two formats: (None,None,graph) or (edges,vertices_number,None)")
        

    def prepare(self):
        """
        Builds the Q matrix manipulating the problem
        optimization function and constraints. It returns
        the time needed to perform the matrix construction 
        in microseconds.
        """
        start_time = time.perf_counter()

        for i in range(len(self.e)):
            self.q[(self.e[i][0],self.e[i][0])] += -1
            self.q[(self.e[i][1],self.e[i][1])] += -1
            self.q[(self.e[i][0],self.e[i][1])] += 2

        end_time = time.perf_counter()

        return (end_time - start_time)*1000000
    
    def sample_advantage(self, num_of_reads, chain_strength = None):
        """
        Performs the sampling using the D-Wave
        using the given number of reads and chian strength
        Advantage QPU and returns its response
        """
        sampler = EmbeddingComposite(DWaveSampler())

        if chain_strength == None:
            chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(self.q, offset = 0.0), sampler)

        print("Computing results on Advantage...")

        sample_set = sampler.sample_qubo(self.q,
                               chain_strength=chain_strength,
                               num_reads=num_of_reads,
                               label='Maximum Cut')
        return sample_set
    
    def sample_hybrid(self):
        """
        Performs the sampling using the D-Wave
        hybrid solver and returns its response
        """
        sampler = LeapHybridSampler(solver={'category': 'hybrid'})

        print("Computing results on Hybrid...")

        sample_set = sampler.sample_qubo(self.q, label="Maximum_Cut")

        return sample_set


    def get_problem_from_input():
            """
            Builds (and returns) a Maximum Cut problem
            instance using the parameters obtained from 
            the input.
            """
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

            return MaximumCutProblem(edges,n_vertices,None)
    
    def test_advantage(number_of_nodes):
        """
        Takes in input the number of nodes to build a random
        instance of the Maximum Cut problem.
        Then it proceeds to solve the problem
        using the D-Wave Advantage QPU and prints the results.
        """
        graph = nx.fast_gnp_random_graph(number_of_nodes,0.5)
        problem = MaximumCutProblem(graph)
        problem.prepare()
        response = problem.sample_advantage(100)
        problem.print_result(response)

    def test_hybrid(number_of_nodes):
        """
        Takes in input the number of nodes to build a random
        instance of the Maximum Cut problem.
        Then it proceeds to solve the problem
        using the D-Wave hybrid solver and prints the results.
        """
        graph = nx.fast_gnp_random_graph(number_of_nodes,0.5)
        problem = MaximumCutProblem(graph)
        problem.prepare()
        response = problem.sample_hybrid()
        problem.print_result(response)

    def solve_classically(self):
        """
        Takes in input an instance of a Maximum Cut Problem
        and solves it using a classical brute force algorithm.
        Then it prints the results and returns the time in microseconds
        needed to obtain the solution classically 
        """
        start_time = time.perf_counter()

        graph = nx.Graph()
        nodes = range(self.v)
        graph.add_nodes_from(nodes)
        edges=[]
        for edge in self.e:
            edges.append(edge)
        graph.add_edges_from(edges)

        max_cut_size = 0
        max_partition = None
        partitions = find_all_partitions(sorted(graph.nodes))
        
        for partition in partitions:
            cut_size = nx.algorithms.cuts.cut_size(graph, partition)
            if cut_size > max_cut_size:
                max_cut_size = cut_size
                max_partition = partition

        end_time = time.perf_counter()

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

        return (end_time-start_time)*1000000

    def print_result(self,response):
        """
        Prints the solution of the Maximum Cut problem
        instance after the sampling
        """
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
