import quadraticKnapsack
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
from dwave.embedding import uniform_torque_compensation
import random


f = open("quadraticKnapsack.csv", "a")
f.write("TESTING ADVANTAGE1\n")
f.write("numvar, numqubit, minenergy, maxchainlength, chainstrength, qpusamplingtime, qpuaccesstime, qpuprogrammingtime, preparetime, classicaltime\n")


for i in range(3,11,1):
    var_number = i*i
    profits = [[random.randint(0, 10) for i in range(var_number)] for j in range(var_number)]
    weights = [random.randint(1, 10) for i in range(var_number)]
    capacity = random.randint((var_number - 1) * 3, var_number * 3)
    weights.append(capacity)
    problem = quadraticKnapsack.quadraticKnapsackProblem(profits,weights,10)
    classical_time = 0.0#problem.solve_classically()
    prepare_time = problem.prepare()
    sampler = EmbeddingComposite(DWaveSampler())
    chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(problem.q, offset = 0.0), sampler)
    sample_set = sampler.sample_qubo(problem.q,
                               chain_strength=chain_strength,
                               num_reads=100,
                               label='Quadratic Knapsack')
    embedding = sample_set.info['embedding_context']['embedding']
    lengths = [len(chain) for chain in embedding.values()]
    num_qubit = sum(lengths)
    max_chain_length = max(lengths)
    f.write("%d, %d, %f, %d, %f, %s, %s, %s, %f, %f\n" % \
                (var_number, num_qubit, sample_set.first.energy, max_chain_length, chain_strength, \
                str(sample_set.info['timing']['qpu_sampling_time']),\
                str(sample_set.info['timing']['qpu_access_time']), str(sample_set.info['timing']['qpu_programming_time']), \
                      prepare_time, classical_time))


'''
f.write("TESTING HYBRID1\n")
f.write("numvar, chainstrength, minenergy, qpuaccesstime, classicaltime\n")
for i in range(3,11,1):
    var_number = i*i
    profits = [[random.randint(0, 10) for i in range(var_number)] for j in range(var_number)]
    weights = [random.randint(1, 10) for i in range(var_number)]
    capacity = random.randint((var_number - 1) * 3, var_number * 3)
    weights.append(capacity)
    problem = quadraticKnapsack.quadraticKnapsackProblem(profits,weights,10)
    classical_time = 0.0#problem.solve_classically()
    prepare_time = problem.prepare()
    sampler = LeapHybridSampler(solver={'category': 'hybrid'})
    chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(problem.q, offset = 0.0), sampler)
    sample_set = sampler.sample_qubo(problem.q, label="Quadratic Knapsack")
    f.write("%d, %f, %f, %s, %f\n" % \
                (var_number, chain_strength, sample_set.first.energy, \
                str(sample_set.info['qpu_access_time']), classical_time))
'''
'''
profits = [[0] * 4 for _ in range(4)]
profits[0] = [2,8,6,10]
profits[1] = [0,5,2,6]
profits[2] = [0,0,2,4]
profits[3] = [0,0,0,4]

weights = [8,6,5,3,16]

problem = quadraticKnapsack.quadraticKnapsackProblem(profits,weights,10)
prepare_time = problem.prepare()
'''