import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict
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

  # Calcola il quadrato della costante (c^2)
  squared_coeff["c^2"] = coeff[-1] ** 2

  return squared_coeff


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
    def __init__(self, profits, weights, penalty = 10):
        
        """
        Constructor of the quadraticKnapsackProblem class.
        Params:
        profits: list of list of profits (polynomial) in the form:
                profits[0][0] = a1 + a11 where a1 is the coefficient of x1 and
                                a11 is the coefficient of (x1)^2 in the polynomial
                profits[0][1]= a12 where a12 is the coefficient of 
                            x1x2 (=x2x1) in the polynomial 
        weights: coefficients of the capacity constraint expressed as a polynomial, e.g.
                    a1x1 + a2x2 + a3x3 + ... + aNxN - C <= 0 becomes
                    [a1 a2 a3 ... an C] where C is the capacity limit
        """
        self.p = profits
        self.w = weights
        self.pen = penalty
        self.q = defaultdict(int)


    def prepare(self):
       
       #adding the profits
       for i in range(len(self.p)):
            for j in range(len(self.p[i])):
                self.q[(i,j)] += self.p[i][j]
                if i!=j:
                  self.q[(j,i)] += self.p[i][j]
       
       #computing weights coefficient
       weights = self.w

       #adding slack variables
       boundary = weights.pop()
       weights.append(1)
       weights.append(2)
       weights.append(boundary)

       #square the weights obtaining a dictionary
       w_dict = squared_pol(weights)

       for key in w_dict.keys():
           w_dict[key] = -self.pen * w_dict[key];

       #here the constant is already negative
       constant = w_dict["c^2"]

       for term,coefficient in w_dict.items():
          if (re.match(r'x(\d)x(\d)',term)):
             match = re.match(r'x(\d)x(\d)',term)
             row = int(match.group(1)) - 1
             col = int(match.group(2)) - 1

             #adding the corresponding term to the Q matrix
             self.q[(row,col)] = self.q[(row,col)] + coefficient
             self.q[(col,row)] = self.q[(col,row)] + coefficient

          elif (re.match(r'x(\d)@2', term)):
             match = re.match(r'x(\d)@2',term)
             index = int(match.group(1)) - 1
            
             #adding the corresponding term to the Q matrix
             self.q[(index,index)] = self.q[(index,index)] + coefficient

          elif (re.match(r'x(\d)c', term)):      
             match = re.match(r'x(\d)c',term)
             index = int(match.group(1)) - 1
             
             #subtracting the corresponding term to the Q matrix
             #because we didnt account for the minus before
             self.q[(index,index)] = self.q[(index,index)] - coefficient
    
       for key in self.q.keys():
           self.q[key] += constant

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
    
    def print_result(self,response):
        lut = response.first.sample
        for i in range(1, len(self.w) - 1, 1):
            if i not in lut:
                lut[i] = 0
            elif ((i == len(self.w) - 2) or (i == len(self.w) - 3)) and (i in lut):
                lut[i] = -99

        print("The solution is: ")

        print("{", end = "")
        for key in sorted(lut):
            if key + 1 in lut:
                print("%s: %s" % (key, lut[key]), end = ", ")
            else:
                print("%s: %s" % (key, lut[key]), end = "")  
        print("}")
