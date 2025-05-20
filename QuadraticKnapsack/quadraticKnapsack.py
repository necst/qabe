import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from collections import defaultdict
import re
import itertools
import random
import time
import math


def min_slack_variables(coefs, A):
    """
    Calculates the minimum number of slack variables required to transform the inequality
    into an equality by adding slack variables using binary expansion.
    """
    max_sum = sum(coefs)
    diff = max_sum - A
    
    if diff < 0:
        return 0,0
    
    # Calculate the number of slack variables required (binary representation of the difference)
    num_slack_vars = math.ceil(math.log2(diff + 1)) if diff > 0 else 0
    
    #binary coefficients
    slack_var_coeffs = [2**i for i in range(num_slack_vars)]
    
    return num_slack_vars, slack_var_coeffs

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
        profits: list of lists of profits (polynomial) in the form:
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
       """
        Builds the Q matrix manipulating the problem
        optimization function and constraints. It returns
        the time needed to perform the matrix construction 
        in microseconds.
       """
       start_time = time.perf_counter()
       
       #adding the profits
       for i in range(len(self.p)):
            for j in range(len(self.p[i])):
                if j >= i:
                    self.q[(i,j)] += self.p[i][j]
                    if i!=j:
                        self.q[(j,i)] += self.p[i][j]
       
       #computing weights coefficient
       weights = self.w

       #adding slack variables
       boundary = weights.pop()

       print(weights)
       print(boundary)

       slack_vars_number, slack_var_coeffs = min_slack_variables(weights, boundary)

       print(slack_var_coeffs)

       for i in range(slack_vars_number):
           weights.append(slack_var_coeffs[i])

       weights.append(boundary)
       #square the weights obtaining a dictionary
       w_dict = squared_pol(weights)

       for key in w_dict.keys():
           w_dict[key] = -self.pen * w_dict[key]

       #here the constant is already negative
       constant = w_dict["c^2"]

       for term,coefficient in w_dict.items():
          if (re.match(r'x(\d+)x(\d+)',term)):
             match = re.match(r'x(\d+)x(\d+)',term)
             row = int(match.group(1)) - 1
             col = int(match.group(2)) - 1

             #adding the corresponding term to the Q matrix
             if (col >= row):
                self.q[(row,col)] = self.q[(row,col)] + coefficient
                self.q[(col,row)] = self.q[(col,row)] + coefficient

          elif (re.match(r'x(\d+)@2', term)):
             match = re.match(r'x(\d+)@2',term)
             index = int(match.group(1)) - 1
            
             #adding the corresponding term to the Q matrix
             self.q[(index,index)] = self.q[(index,index)] + coefficient

          elif (re.match(r'x(\d+)c', term)):      
             match = re.match(r'x(\d+)c',term)
             index = int(match.group(1)) - 1
            
             #subtracting the corresponding term to the Q matrix
             #because we didnt account for the minus before
             self.q[(index,index)] = self.q[(index,index)] - coefficient

       #making Q an upper diagonal matrix (and) for a minimaztion problem
       for (i, j) in self.q.keys():
            if j < i:
                self.q[(i, j)] = 0
            else:
                self.q[(i, j)] *= -1

       end_time = time.perf_counter()

       # Define matrix size (15x15 based on the dictionary)
       size = 7

    # Create a matrix and initialize all values to 0
       matrix = [[0] * size for _ in range(size)]

    # Populate the matrix based on the dictionary
       for (i, j), value in self.q.items():
           matrix[i][j] = value

# Print the matrix in rows and columns
       for row in matrix:
           print(row)

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
                               label='Quadratic Knapsack')
        return sample_set
    
    
    def sample_hybrid(self):
        """
        Performs the sampling using the D-Wave
        hybrid solver and returns its response
        """
        sampler = LeapHybridSampler(solver={'category': 'hybrid'})

        print("Computing results on Hybrid...")

        sample_set = sampler.sample_qubo(self.q, label="Quadratic Knapsack")

        return sample_set
    

    def test_advantage(number_of_variables):
        """
        Takes in input the number of variables of the problem
        to build a random instance of the Quadratic Assignment problem.
        Then it proceeds to solve the problem
        using the D-Wave Advantage QPU and prints the results.
        """
        if (number_of_variables <= 1):
            raise TypeError("The number of variables must be at least 2")
        penalty = 10
        profits = [[random.randint(0, 10) for i in range(number_of_variables)] for j in range(number_of_variables)]
        #Set to have at least weight one for each 
        weights = [random.randint(1, 10) for i in range(number_of_variables)]
        capacity = random.randint((number_of_variables - 1) * 3, number_of_variables * 3)
        weights.append(capacity)
        problem = quadraticKnapsackProblem(profits,weights,penalty)
        problem.prepare()
        response = problem.sample_advantage(100)
        problem.print_result(response)

    def test_hybrid(number_of_variables):
        """
        Takes in input the number of variables of the problem
        to build a random instance of the Quadratic Assignment problem.
        Then it proceeds to solve the problem
        using the D-Wave hybrid solver and prints the results.
        """
        if (number_of_variables <= 1):
            raise TypeError("The number of variables must be at least 2")
        penalty = 10
        profits = [[random.randint(0, 10) for i in range(number_of_variables)] for j in range(number_of_variables)]
        #Set to have at least weight one for each 
        weights = [random.randint(1, 10) for i in range(number_of_variables)]
        capacity = random.randint((number_of_variables - 1) * 3, number_of_variables * 3)
        weights.append(capacity)
        problem = quadraticKnapsackProblem(profits,weights,penalty)
        problem.prepare()
        response = problem.sample_hybrid()
        problem.print_result(response)



    def print_result(self,response):
        """
        Prints the solution of the Quadratic Knapsack problem
        instance after the sampling
        """
        lut = response.first.sample
        for i in range(0, len(self.w) - 1, 1):
            if i not in lut:
                lut[i] = 0
            elif ((i == len(self.w) - 2) or (i == len(self.w) - 3)) and (i in lut):
                del lut[i]

        print("The solution is: ")

        print("{", end = "")
        for key in sorted(lut):
            if key + 1 in lut:
                print("%s: %s" % (key + 1, lut[key]), end = ", ")
            else:
                print("%s: %s" % (key + 1, lut[key]), end = "")  
        print("}")


    def solve_classically(self):
        """
        Takes in input an instance of a Quadratic Knapsack Problem
        and solves it using a classical brute force algorithm.
        Then it prints the results and returns the time in microseconds
        needed to obtain the solution classically 
        """
        start_time = time.perf_counter()

        weights = self.w
        profits = self.p
        
        capacity = weights.pop()
        n = len(profits[0])

        best_profit = 0
        best_combination = None

        #check all possible combinations
        for combination in itertools.product([0, 1], repeat=n):
            
            total_weight = sum(weights[i] for i, selected in enumerate(combination) if selected == 1)

            #iterate only if the capacity constraint is respected
            if total_weight <= capacity:

                total_profit = 0
                for i in range(n):
                    if combination[i] == 1:
                        total_profit += profits[i][i]
                        for j in range(i+1, n):
                            if combination[j] == 1:
                                total_profit += profits[i][j]

                #check if this combination has a better profit
                if total_profit > best_profit:
                    best_profit = total_profit
                    best_combination = combination

        end_time = time.perf_counter()

        print('The solution is:')
        print("{", end = "")
        for i in range(len(best_combination)):
            if i + 1  < len(best_combination):
                print("%s: %s" % (i + 1, best_combination[i]), end = ", ")
            else:
                print("%s: %s" % (i + 1, best_combination[i]), end = "")
        print("}")

        return (end_time - start_time)*1000000