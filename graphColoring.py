import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 
import networkx as nx
import dimod
import re
from itertools import product


def valid_coloring(edges, coloring):

    for (u, v) in edges:
        if coloring[u] == coloring[v]:
            return False
    return True


def squared_pol(coeff):

  n = len(coeff) - 1
  squared_coeff = {}

  #squared terms(a_i^2)
  for i in range(n):
    squared_coeff[f"x{i+1}@2"] = coeff[i] * coeff[i]

  #cross terms(a_i * a_j)
  for i in range(n):
    for j in range(i + 1, n):
      squared_coeff[f"x{i+1}x{j+1}"] = 2 * coeff[i] * coeff[j]

  #costant multiplying terms (a_i * c)
  for i in range(n):
    squared_coeff[f"x{i+1}c"] = 2 * coeff[i] * coeff[-1]

  #constant squared (c^2)
  squared_coeff["c^2"] = coeff[-1] ** 2

  return squared_coeff


class GraphColoringProblem:
    """
    This class represents a Graph Coloring problem.
    It takes a graph G = {V,E} and an integer K as
    input and returns (ATTEMPTS to return) an assignment
    of K colors in such a way that adjacent nodes receive
    different colors.
    """
    def __init__(self, colors_num , edges=None, vertices_num=None, graph = None):
        """
        Constructor of the graph Coloring class.
        Params:
        colors_num: number of colors
        edges: edges of the graph in the form of an array of arrays
                e.g. [[i,j],[j,k],...] where i-j and j-k are two edges
                connecting the nodes i,j and the nodes j,k
        vertices_num: number of vertices in the graph

        or 
        
        colors_num: number of colors
        graph: a NetworkX graph (the others are None)
        """

        if edges == None and vertices_num == None and graph != None and colors_num != None:
            if not isinstance(graph,  nx.classes.graph.Graph):
                raise TypeError("A Networkx graph is required")
            
            self.v = graph.number_of_nodes()
            edges = []
            add = 0
            for u,v in graph.edges(data=False):
                if int(u)==0 or int(v)==0:
                    add = 1
            for u,v in graph.edges(data=False):
                edges.append([int(u)+add,int(v)+add])
            self.e = edges
            self.q = defaultdict(int)
            self.k = colors_num

        elif edges != None and vertices_num != None and graph == None and colors_num != None:
            self.e = edges
            self.v = int(vertices_num)
            self.q = defaultdict(int)
            self.k = int(colors_num)
        else:
            raise TypeError("Specify the arguments in one of the two formats: (k,None,None,graph,penalty) or (k,edges,num_vertices,None,penalty)")
        
    def prepare(self):

        #any positive value for the penalty will do since we do not have an objective function
        penalty = 4

        variables_number = self.v * self.k

        pol = []

        for i in range(self.k):
            pol.append(1)

        pol.append(-1)

        w_dict = squared_pol(pol)

        for key in w_dict.keys():
           w_dict[key] = 4 * w_dict[key]

        for i in range(0,variables_number,self.k):

            for term,coefficient in w_dict.items():
                if (re.match(r'x(\d+)x(\d+)',term)):
                    match = re.match(r'x(\d+)x(\d+)',term)
                    row = int(match.group(1)) - 1 + i
                    col = int(match.group(2)) - 1 + i

                    #adding the corresponding term to the Q matrix
                    if (col >= row):
                        self.q[(row,col)] = self.q[(row,col)] + coefficient
                        self.q[(col,row)] = self.q[(col,row)] + coefficient

                elif (re.match(r'x(\d+)@2', term)):
                    match = re.match(r'x(\d+)@2',term)
                    index = int(match.group(1)) - 1 + i
                        
                    #adding the corresponding term to the Q matrix
                    self.q[(index,index)] = self.q[(index,index)] + coefficient

                elif (re.match(r'x(\d+)c', term)):      
                    match = re.match(r'x(\d+)c',term)
                    index = int(match.group(1)) - 1 + i
                        
                    #adding the corresponding term to the Q matrix
                    self.q[(index,index)] = self.q[(index,index)] + coefficient


        for edge in self.e:
            v1 = int(edge[0])
            v2 = int(edge[1])

            if (v2 < v1):
                v1, v2 = v2, v1

            for k in range (1, self.k + 1):
                #the last subtraction (-1) is due to the fact that the matrix starts from (0,0)
                v1_index = (self.k * (v1-1) + k) - 1
                v2_index = (self.k * (v2-1) + k) - 1

                self.q[(v1_index,v2_index)] += penalty 


    def sample_advantage(self, num_of_reads, chain_strength = None):

        sampler = EmbeddingComposite(DWaveSampler())

        if chain_strength == None:
            chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(self.q, offset = 0.0), sampler)

        print("Computing results on Advantage...")

        sample_set = sampler.sample_qubo(self.q,
                               chain_strength=chain_strength,
                               num_reads=num_of_reads,
                               label='Graph Coloring')
        return sample_set
    
    def sample_hybrid(self):

        sampler = LeapHybridSampler(solver={'category': 'hybrid'})

        print("Computing results on Hybrid...")

        sample_set = sampler.sample_qubo(self.q, label="Graph Coloring")

        return sample_set
    
    def test_advantage(number_of_nodes,colors):
        graph = nx.fast_gnp_random_graph(number_of_nodes,0.5)
        problem = GraphColoringProblem(colors,None,None,graph)
        problem.prepare()
        response = problem.sample_advantage(100)
        problem.print_result(response)

    def test_hybrid(number_of_nodes,colors):
        graph = nx.fast_gnp_random_graph(number_of_nodes,0.5)
        problem = GraphColoringProblem(colors,None,None,graph)
        problem.prepare()
        response = problem.sample_hybrid()
        problem.print_result(response)
    

    def get_problem_from_input():
            print("Insert the number of colors: ")
            n_colors = input()
            if( int(n_colors) <= 0):
                raise ValueError("The number of colors must be positive")
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

            return GraphColoringProblem(n_colors,edges,n_vertices,None)
    
    def solve_classically(self):

        vertices = set()
        for edge in self.e:
            vertices.update(edge)
        vertices = list(vertices)

        for coloring in product(range(self.k), repeat=self.v):
            vertex_color = {vertex: coloring[i] for i, vertex in enumerate(vertices)}
        
            if valid_coloring(self.e, vertex_color):
                print("The solution is:")
                for key in vertex_color:
                    print("The vertex %s has color: %s" % (key, vertex_color[key]+1))
                return
        
        print("No accettable solution has been found for the problem")

    def print_result(self,response):
        lut = response.first.sample

        print("The solution is: ")

        for key in lut:
            if lut[key] == 1:
                print("The vertex %s has color: %s" % ((key//self.k)+1, self.k-((key+1)%self.k)))