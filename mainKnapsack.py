import quadraticKnapsack
import array as ar
#f = x1 + x1x3 +x2 + x3
list = [[1,0,1],[0,1,0],[0,0,1]]
print(list[0][0])
b = ar.array('i',(0 for i in range(0,4)))
b[0] = 1
b[1] = 2
b[2] = 3
b[3] = 1
quadraticKnapsack.quadraticKnapsackProblem(list,b,1).prepare()
'''
coeff_polinomio = b
coeff_quadrato = quadraticKnapsack.squared_pol(coeff_polinomio)

print("Coefficienti del polinomio:", coeff_polinomio)
print("Coefficienti del polinomio al quadrato:")
for termine, coefficiente in coeff_quadrato.items():
  print(f"{termine}: {coefficiente}")'''