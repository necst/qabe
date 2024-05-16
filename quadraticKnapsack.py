import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 

class quadraticKnapsackProblem:
    """
    This class represents a Quadratic Knapsack problem.
    It takes a sum of profits given as a function to maximize.
    The function is expressed using linear terms to describe
    the profit of selecting a particular item and quadratic 
    cross terms to describe the profit of selecting both of them.
    The choice is subject to a constraint formalizing the
    capacity of the knapsack, also given as input.
    It produces the best choice for the items in the knapsack.
    """
    def __init__(self, profits, weights):
        
        """
        Constructor of the quadraticKnapsackProblem class.
        Params:
        profits: coefficients of the function to maximize
        weights: coefficients of the capacity constraint expressed as a polynomial, e.g.
                    a1x1 + a2x2 + a3x3 + ... + aNxN - C <= 0 becomes
                    [a1 a2 a3 ... an C] where C is the capacity limit
        """
        self.p = profits
        self.w = weights
        self.q = defaultdict(int)