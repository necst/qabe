import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 
import networkx as nx
import dimod
import re




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
    def __init__(self, colors_num , edges=None, vertices_num=None, graph = None, penalty = 10):
        """
        Constructor of the graph Coloring class.
        Params:
        colors_num: number of colors
        edges: edges of the graph in the form of an array of arrays
                e.g. [[i,j],[j,k],...] where i-j and j-k are two edges
                connecting the nodes i,j and the nodes j,k
        vertices_num: number of vertices in the graph
        penalty: the chosen penalty for the constraints


        or 
        
        colors_num: number of colors
        graph: a NetworkX graph (the others are None except for )
        penalty: the chosen penalty for the constraints

        
        """

        if edges == None and vertices_num == None and graph != None and colors_num != None:
            if not isinstance(graph,  nx.classes.graph.Graph):
                raise TypeError("A Networkx graph is required")
            
            self.v = graph.number_of_nodes()
            edges = []
            for u,v in graph.edges(data=False):
                edges.append([int(u),int(v)])
            self.e = edges
            self.q = defaultdict(int)
            self.k = colors_num
            self.pen = penalty

        elif edges != None and vertices_num != None and graph == None and colors_num != None:
            self.e = edges
            self.v = int(vertices_num)
            self.q = defaultdict(int)
            self.k = int(colors_num)
            self.pen = penalty
        else:
            raise TypeError("Specify the arguments in one of the two formats: (k,None,None,graph,penalty) or (k,edges,num_vertices,None,penalty)")
        
    def prepare(self):

        variables_number = self.v * self.k

        pol = []

        for i in range(self.k):
            pol.append(1)

        pol.append(-1)

        w_dict = squared_pol(pol)

        for key in w_dict.keys():
           w_dict[key] = -self.pen * w_dict[key]

        for i in range(0,self.v*self.k,self.k):

            for j in range(i,i + self.k):

                for term,coefficient in w_dict.items():
                    if (re.match(r'x(\d+)x(\d+)',term)):
                        match = re.match(r'x(\d+)x(\d+)',term)
                        row = int(match.group(1)) - 1 + j
                        col = int(match.group(2)) - 1 + j

                        print("adding %s to %s,%s", coefficient, row, col)

                        #adding the corresponding term to the Q matrix
                        if (col >= row):
                            self.q[(row,col)] = self.q[(row,col)] + coefficient
                            self.q[(col,row)] = self.q[(col,row)] + coefficient

                    elif (re.match(r'x(\d+)@2', term)):
                        match = re.match(r'x(\d+)@2',term)
                        index = int(match.group(1)) - 1 + j
                        
                        #adding the corresponding term to the Q matrix
                        self.q[(index,index)] = self.q[(index,index)] + coefficient

                    elif (re.match(r'x(\d+)c', term)):      
                        match = re.match(r'x(\d+)c',term)
                        index = int(match.group(1)) - 1 + j
                        
                        #subtracting the corresponding term to the Q matrix
                        #because we didnt account for the minus before
                        self.q[(index,index)] = self.q[(index,index)] - coefficient

        print(self.q)
        


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

            print("Insert the penalty: ")
            penalty = int(input())

            return GraphColoringProblem(n_colors,edges,n_vertices,None,penalty)
    