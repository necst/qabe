import quadraticAssignment
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
import random

# Copyright 2020 D-Wave Systems Inc.
#
#    Licensed under the Apache License, Version 2.0 (the "License");
#    you may not use this file except in compliance with the License.
#    You may obtain a copy of the License at
#
#        http://www.apache.org/licenses/LICENSE-2.0
#
#    Unless required by applicable law or agreed to in writing, software
#    distributed under the License is distributed on an "AS IS" BASIS,
#    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#    See the License for the specific language governing permissions and
#    limitations under the License.

"""Utility functions for calculating chain strength.

Examples:
    This example uses :func:`uniform_torque_compensation`, given a prefactor of 2, 
    to calculate a chain strength that :class:`EmbeddingComposite` then uses.

    >>> from functools import partial
    >>> from dwave.system import EmbeddingComposite, DWaveSampler
    >>> from dwave.embedding.chain_strength import uniform_torque_compensation
    ...
    >>> Q = {(0,0): 1, (1,1): 1, (2,3): 2, (1,2): -2, (0,3): -2}
    >>> sampler = EmbeddingComposite(DWaveSampler())
    >>> # partial() can be used when the BQM or embedding is not accessible
    >>> chain_strength = partial(uniform_torque_compensation, prefactor=2)
    >>> sampleset = sampler.sample_qubo(Q, chain_strength=chain_strength, return_embedding=True)
    >>> sampleset.info['embedding_context']['chain_strength']
    1.224744871391589

"""
import math

__all__ = ['uniform_torque_compensation', 'scaled']

def uniform_torque_compensation(bqm, embedding=None, prefactor=1.414):
    """Chain strength that attempts to compensate for torque that would break
    the chain.

    The RMS of the problem's quadratic biases is used for calculation.

    Args: 
        bqm (:obj:`.BinaryQuadraticModel`):
            A binary quadratic model.

        embedding (dict/:class:`.EmbeddedStructure`, default=None):
            Included to satisfy the `chain_strength` callable specifications 
            for `embed_bqm`. 

        prefactor (float, optional, default=1.414):
            Prefactor used for scaling. For non-pathological problems, the recommended 
            range of prefactors to try is [0.5, 2].

    Returns:
        float: The chain strength, or 1 if chain strength is not applicable.

    """
    if bqm.num_interactions > 0:
        squared_j = (j ** 2 for j in bqm.quadratic.values())
        rms = math.sqrt(sum(squared_j)/bqm.num_interactions)
        avg_degree = bqm.degrees(array=True).mean()

        return prefactor * rms * math.sqrt(avg_degree)
    else:
        return 1    # won't matter (chain strength isn't needed to embed this problem)


def scaled(bqm, embedding=None, prefactor=1.0):
    """Chain strength that is scaled to the problem bias range.

    Args:
        bqm (:obj:`.BinaryQuadraticModel`):
            A binary quadratic model.

        embedding (dict/:class:`.EmbeddedStructure`, default=None):
            Included to satisfy the `chain_strength` callable specifications 
            for `embed_bqm`. 

        prefactor (float, optional, default=1.0):
            Prefactor used for scaling. 

    Returns:
        float: The chain strength, or 1 if chain strength is not applicable.

    """  
    if bqm.num_interactions > 0:
        max_bias = max(max(bqm.linear.max(), -bqm.linear.min()), 
                       max(bqm.quadratic.max(), -bqm.quadratic.min()))
        return prefactor * max_bias
    else:
        return 1    # won't matter (chain strength isn't needed to embed this problem)


def generate_matrix(n):

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


f = open("quadraticAssignment.csv", "a")
f.write("TESTING ADVANTAGE1\n")
f.write("numvar, numqubits, minenergy, maxchainlength, chainstrength, qpusamplingtime, qpuaccesstime, qpuprogrammingtime, preparetime, classicaltime\n")

for i in range(3,11,1):
    var_number = i*i
    flow = generate_matrix(i)
    distance = generate_matrix(i)
    problem = quadraticAssignment.QuadraticAssignmentProblem(flow,distance)
    classical_time = 0.0#problem.solve_classically()
    prepare_time = problem.prepare()
    sampler = EmbeddingComposite(DWaveSampler())
    chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(problem.q, offset = 0.0), sampler)
    sample_set = sampler.sample_qubo(problem.q,
                                     chain_strength=chain_strength,
                                    num_reads=100,
                                    label="Quadratic Assignment")
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
    flow = generate_matrix(i)
    distance = generate_matrix(i)
    problem = quadraticAssignment.QuadraticAssignmentProblem(flow,distance)
    classical_time = 0.0#problem.solve_classically()
    prepare_time = problem.prepare()
    sampler = LeapHybridSampler(solver={'category': 'hybrid'})
    chain_strength = uniform_torque_compensation(dimod.BinaryQuadraticModel.from_qubo(problem.q, offset = 0.0), sampler)
    sample_set = sampler.sample_qubo(problem.q, label="Quadratic Knapsack")
    f.write("%d, %f, %f, %s, %f\n" % \
                (var_number, chain_strength, sample_set.first.energy, \
                str(sample_set.info['qpu_access_time']), classical_time))

'''