import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 
import networkx as nx
import dimod



class MaximumCutProblem:
    """
    This class represents a Graph Coloring problem.
    It takes a graph G = {V,E} and an integer K as
    input and returns (ATTEMPTS to return) an assignment
    of K colors in such a way that adjacent nodes receive
    different colors.
    """
    def __init__(self, k , edges=None, vertices_num=None, graph = None):
        """
        Constructor of the MaximumCutProblem class.
        Params:
        edges: edges of the graph
        vertices_num: number of vertices in the graph

        or 


        """

        if edges == None and vertices_num == None and graph != None and k != None:
            if not isinstance(graph,  nx.classes.graph.Graph):
                raise TypeError("A Networkx graph is required")
            
            self.v = graph.number_of_nodes()
            edges = []
            for u,v in graph.edges(data=False):
                edges.append([int(u),int(v)])
            self.e = edges
            self.q = defaultdict(int)

        elif edges != None and vertices_num != None and graph == None and k != None:
            self.e = edges
            self.v = vertices_num
            self.q = defaultdict(int)
        else:
            raise TypeError("Specify the arguments in the format: (k,None,None,graph)")
        
