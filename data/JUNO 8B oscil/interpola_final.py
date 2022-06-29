import numpy as np
from scipy import interpolate
from scipy.interpolate import fitpack

nome_entrada = "total_int.txt"#input("Entre com o nome do arquivo para interpolar: ")
data = np.loadtxt(nome_entrada)

nome_saida = nome_entrada.replace(".txt", "_int.txt")
saida = open(nome_saida, "w")

X_pre = data[:,0]
Y_pre = data[:,1]
Z_pre = data[:,2]

X = np.array([])
Y = np.array([])

X = np.append(X, X_pre[0])
for i in range(np.size(X_pre)):
	if X[-1] < X_pre[i]:
		X = np.append(X, X_pre[i])

i = 0
Y = np.append(Y, Y_pre[0])
foi = False
while not foi:
	if Y[-1] < Y_pre[i]:
		Y = np.append(Y, Y_pre[i])
	if Y[-1] > Y_pre[i]:
		foi = True
	i += 1

X_dim = np.size(X)
Y_dim = np.size(Y)

Z = np.empty((X_dim,Y_dim))

print(X_dim, Y_dim)

for i in range(X_dim):
	for j in range(Y_dim):
		Z[i][j] = Z_pre[j+i*Y_dim]

print(Z.shape)

interpola = interpolate.RectBivariateSpline(X, Y, Z, kx=1, ky=1)

N = 1_000

min_x = 0.2#min(X)
max_x = 0.45#max(X)
min_y = 2e-5#min(Y)
max_y = 1e-4#max(Y)

h_x = (max_x-min_x)/N
h_y = (max_y-min_y)/N

for x in np.arange(min_x, max_x, h_x):
	for y in np.arange(min_y, max_y, h_y):
		z = float(interpola(x, y))
		saida.write(str(x)+" "+str(y)+" "+str(z)+"\n")
	saida.write("\n")