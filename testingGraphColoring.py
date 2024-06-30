import graphColoring
import dimod.binary_quadratic_model
from dwave.system import EmbeddingComposite
from dwave.embedding.chain_strength import *
from dwave.system import DWaveSampler
from dwave.embedding import *
from dimod import BinaryQuadraticModel
from dwave.system import DWaveSampler
from dwave.system import LeapHybridSampler
from dwave.system import EmbeddingComposite
from dwave.inspector import *
from dwave.embedding.chimera import *
from dwave.embedding.chain_strength import *
import networkx as nx
import math

def get_parameters(target_variables, p):
    number_of_nodes = 1
    while True:
        number_of_colors = calculate_colors(number_of_nodes, p)
        total_variables = number_of_nodes * number_of_colors
        if total_variables >= target_variables:
            return number_of_nodes, number_of_colors, total_variables
        number_of_nodes += 1
    

def calculate_colors(number_of_nodes, p):
    avg_degree = p * (number_of_nodes - 1)
    return math.ceil(avg_degree / 2)

f = open("graphColoring.csv", "a")
f.write("TESTING ADVANTAGE\n")
f.write("numvar, numqubit, minenergy, maxchainlength, chainstrength, qpusamplingtime, qpuaccesstime, qpuprogrammingtime, preparetime, classicaltime\n")



for i in range(3,11,1):
    var_number = i*i
    edge_probability = 0.3
    nodes, colors, variables = get_parameters(var_number, edge_probability)
    graph = nx.fast_gnp_random_graph(nodes, edge_probability)
    problem = graphColoring.GraphColoringProblem(colors,None,None,graph)
    classical_time = 0.0#problem.solve_classically()
    prepare_time = problem.prepare()
    sampler = EmbeddingComposite(DWaveSampler())
    chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(problem.q, offset = 0.0), sampler)
    sample_set = sampler.sample_qubo(problem.q,
                               chain_strength=chain_strength,
                               num_reads=100,
                               label='Graph Coloring')
    embedding = sample_set.info['embedding_context']['embedding']
    lengths = [len(chain) for chain in embedding.values()]
    num_qubit = sum(lengths)
    max_chain_length = max(lengths)
    f.write("%d, %d, %f, %d, %f, %s, %s, %s, %f, %f\n" % \
                ( variables, num_qubit, sample_set.first.energy, max_chain_length, chain_strength, \
                str(sample_set.info['timing']['qpu_sampling_time']),\
                str(sample_set.info['timing']['qpu_access_time']), str(sample_set.info['timing']['qpu_programming_time']), \
                      prepare_time, classical_time))


f.write("TESTING HYBRID\n")
f.write("numvar, minenergy, qpuaccesstime, classicaltime\n")
for i in range(3,11,1):
    var_number = i*i
    edge_probability = 0.3
    nodes, colors, variables = get_parameters(var_number, edge_probability)
    graph = nx.fast_gnp_random_graph(nodes, edge_probability)
    problem = graphColoring.GraphColoringProblem(colors,None,None,graph)
    classical_time = 0.0
    prepare_time = problem.prepare()
    sampler = LeapHybridSampler(solver={'category': 'hybrid'})
    sample_set = sampler.sample_qubo(problem.q, label="Quadratic Knapsack")
    f.write("%d, %f, %s, %f\n" % \
                (var_number, sample_set.first.energy, \
                str(sample_set.info['qpu_access_time']), classical_time))