import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict 
import math
import re
import itertools
import random
import time


def generate_matrix(n):
    """
    The function generates a matrix of random integers with a
    null diagonal given the size received in input.
    """
    matrix = []

    for i in range(n):

        row = []

        for j in range(n):
            if i == j:
                row.append(0)
            else:
                row.append(random.randint(1, 15))

        matrix.append(row)

    return matrix

def modified_squared_pol(coeff,offset):
  """
  Returns a dictionary containing all the terms and
  their coefficients (as key value pairs) of the squared 
  polynomial represented by the coefficients' list in input.
  To the terms considered is applied an offset to shift their
  indexes accordingly.
  """
  n = len(coeff) - 1
  squared_coeff = {}

  #squared terms(a_i^2)
  for i in range(n):
    index = i*offset
    squared_coeff[f"x{index+1}@2"] = coeff[i] * coeff[i]

  #cross terms(a_i * a_j)
  for i in range(n):
    for j in range(i + 1, n):
      index_i = i*offset
      index_j = j*offset
      squared_coeff[f"x{1+index_i}x{1+index_j}"] = 2 * coeff[i] * coeff[j]

  #costant multiplying terms (a_i * c)
  for i in range(n):
    index = i*offset
    squared_coeff[f"x{index+1}c"] = 2 * coeff[i] * coeff[-1]

  #constant squared (c^2)
  squared_coeff["c^2"] = coeff[-1] ** 2

  return squared_coeff


def squared_pol(coeff):
  """
  Returns a dictionary containing all the terms and
  their coefficients (as key value pairs) of the squared 
  polynomial represented by the coefficients' list in input.
  """
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


def compute_cost(permutation, flow, distance, n):
    """
    Computes and returns the cost of the single facilities-locations
    assignment returning the sum of the flow times distance for each
    pair of facilities
    """
    cost = 0

    for i in range(n):
        for j in range(n):
            cost += flow[i][j] * distance[permutation[i]][permutation[j]]

    return cost


class QuadraticAssignmentProblem:
    """
    This class represents a Quadratic Assignment problem.
    It takes two matrices of values representing respectively
    the flow of material between facilities and the distances
    between the sites as input and returns an assignment
    of the facilities to locations to minimize the weighted
    flow across the system.
    """
    def __init__(self, facilities, distances, penalty = 10):
        """
        Constructor of the graph Coloring class.
        Params:
        colors_num: number of colors
        facilities: matrix indicating the flow of material between
                facilities in the form of a list of lists
                where facilities[i][j] represents the flow
                between facility i and facility j
        distances: matrix indicating the distances between
                locations in the form of a list of lists
                where locations[i][j] represents the distance
                between location i and location j
        """

        self.f = facilities
        self.d = distances
        self.pen = penalty
        self.q = defaultdict(int)


    def prepare(self):
        """
        Builds the Q matrix manipulating the problem
        optimization function and constraints. It returns
        the time needed to perform the matrix construction 
        in microseconds.
        """
        start_time = time.perf_counter()

        var_number = len(self.f)

        squared_var_number = var_number*var_number

        for k in range(1,squared_var_number):

            for s in range(k+1,squared_var_number+1):

                #from 1 to n
                i = math.floor((k-1) / var_number) + 1
                j = (k-1) % var_number + 1

                #from 1 to n
                m = math.floor((s-1)/var_number) + 1
                n = (s-1)%var_number + 1

    
                if (i != m and j != n):
                    # Using x-1 to access matrices from 0 to n
                    self.q[(k-1,s-1)] = self.q[(k-1,s-1)] + 2*self.f[i-1][m-1]*self.d[j-1][n-1]

        pol = []

        for i in range(var_number):
            pol.append(1)

        pol.append(-1)

        w_dict = squared_pol(pol)

        for key in w_dict.keys():
           w_dict[key] = self.pen * w_dict[key]

        for i in range(0,squared_var_number,var_number):

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


        pol = []

        for i in range(var_number):
            pol.append(1)

        pol.append(-1)

        w_dict = modified_squared_pol(pol,var_number)

        for key in w_dict.keys():
           w_dict[key] = self.pen * w_dict[key]

        for i in range(0,var_number):
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
                               label='Graph Coloring')
        return sample_set
    
    def sample_hybrid(self):
        """
        Performs the sampling using the D-Wave
        hybrid solver and returns its response
        """
        sampler = LeapHybridSampler(solver={'category': 'hybrid'})

        print("Computing results on Hybrid...")

        sample_set = sampler.sample_qubo(self.q, label="Graph Coloring")

        return sample_set
    
    def test_advantage(number_of_facilities, penalty = 200):
        """
        Takes in input the number of facilities (the
        number of locations is the same) and the penalty values
        to build a random instance of the Quadratic Assignment problem.
        Then it proceeds to solve the problem
        using the D-Wave Advantage QPU and prints the results.
        """
        flow = generate_matrix(number_of_facilities)
        distance = generate_matrix(number_of_facilities)
        problem = QuadraticAssignmentProblem(flow,distance,penalty)
        print(flow)
        print(distance)
        problem.prepare()
        response = problem.sample_advantage(100)
        problem.print_result(response)
    
    def test_hybrid(number_of_facilities, penalty = 200):
        """
        Takes in input the number of facilities (the
        number of locations is the same) and the penalty values
        to build a random instance of the Quadratic Assignment problem.
        Then it proceeds to solve the problem
        using the D-Wave hybrid solver and prints the results.
        """
        flow = generate_matrix(number_of_facilities)
        distance = generate_matrix(number_of_facilities)
        problem = QuadraticAssignmentProblem(flow,distance,penalty)
        problem.prepare()
        response = problem.sample_hybrid()
        problem.print_result(response)

    def solve_classically(self):
        """
        Takes in input an instance of a Quadratic Assignment Problem
        and solves it using a classical brute force algorithm.
        Then it prints the results and returns the time in microseconds
        needed to obtain the solution classically 
        """
        start_time = time.perf_counter()
       
        var_number = len(self.f)

        permutations = list(itertools.permutations(range(var_number)))
        best_cost = float('inf')
        best_permutation = None

        for perm in permutations:
          current_cost = compute_cost(perm, self.f, self.d, var_number)
          if current_cost < best_cost:
              best_cost = current_cost
              best_permutation = perm
        
        end_time = time.perf_counter()
        
        print("The solution is: ")

        for i in range (len(best_permutation)):
            print("The facility %s is assigned to location: %s" % (i+1, best_permutation[i]+1))

        return (end_time - start_time)*1000000

        

    def print_result(self,response):
        """
        Prints the solution of the Quadratic Assignment problem
        instance after the sampling
        """
        lut = response.first.sample
        print(lut)

        print("The solution is: ")

        for key in lut:
            if lut[key] == 1:
                print("The facility %s is assigned to location: %s" % (math.floor((key) / len(self.d)) + 1, (key) % len(self.d) + 1))
    
        

