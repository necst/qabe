import quadraticKnapsack
import array as ar
#f = x1 + x1x3 +x2 + x3
'''
list = [[1,0,1],[0,1,0],[1,0,1]]
b = ar.array('i',(0 for i in range(0,4)))
b[0] = 1
b[1] = 2
b[2] = 3
b[3] = 1
'''
list = [[2,8,6,10],[0,5,2,6],[0,0,2,4],[0,0,0,4]]
b = ar.array('i',(0 for i in range(0,5)))
b[0] = 8
b[1] = 6
b[2] = 5
b[3] = 3
b[4] = 16
#problem = quadraticKnapsack.quadraticKnapsackProblem(list,b,10)
#problem.prepare()
#response = problem.sample_advantage(100)
#problem.print_result(response)
#quadraticKnapsack.quadraticKnapsackProblem.solve_classically(list,b)

quadraticKnapsack.quadraticKnapsackProblem.test_hybrid(9)