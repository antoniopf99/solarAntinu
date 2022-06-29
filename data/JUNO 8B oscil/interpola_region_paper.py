import numpy as np
from scipy import interpolate
from scipy.interpolate import fitpack

nome_entrada = "total.txt"#input("Entre com o nome do arquivo para interpolar: ")
data = np.loadtxt(nome_entrada)

nome_saida = nome_entrada.replace(".txt", "_int.txt")
saida = open(nome_saida, "w")

X = data[:,(0,1)]
Z = data[:,2]

interpola = interpolate.NearestNDInterpolator(X, Z, rescale=True)

N = 200

X_axis = data[:,0]
Y_axis = data[:,1]

min_x = min(X_axis)
max_x = max(X_axis)
min_y = min(Y_axis)
max_y = max(Y_axis)

h_x = (max_x-min_x)/N
h_y = (max_y-min_y)/N

for x in np.arange(min_x, max_x, h_x):
	for y in np.arange(min_y, max_y, h_y):
		z = float(interpola(x, y))
		saida.write(str(x)+" "+str(y)+" "+str(z)+"\n")
	saida.write("\n")