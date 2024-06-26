import quadraticAssignment

list = [[0,5,2],[5,0,3],[2,3,0]]
list2 = [[0,8,15],[8,0,13],[15,13,0]]
problem = quadraticAssignment.QuadraticAssignmentProblem(list,list2,200)
problem.prepare()
response = problem.sample_advantage(100)
problem.print_result(response)
