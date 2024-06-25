import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 


class QuadraticAssignmentProblem:
    """
    This class represents a Quadratic Assignment problem.
    It takes two matrices of values representing respectively
    the flow of material between facilities and the distances
    between the sites as input and returns an assignment
    of the facilities to locations to minimize the weighted
    flow across the system.
    """
    def __init__(self, facilities, locations, penalty = 10):
        """
        Constructor of the graph Coloring class.
        Params:
        colors_num: number of colors
        facilities: matrix indicating the flow of material between
                facilities in the form of a list of lists
                where facilities[i][j] represents the flow
                between facility i and facility j
        locations: matrix indicating the distances between
                locations in the form of a list of lists
                where locations[i][j] represents the distance
                between location i and location j
        """

        self.f = facilities
        self.l = locations
        self.pen = penalty
        self.q = defaultdict(int)


    def prepare(self):

        var_number = len(self.f)

        print(var_number)