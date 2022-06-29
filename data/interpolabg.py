import numpy as np
from scipy import interpolate

nome_entrada = input("Entre com o nome do arquivo: ")

data = np.loadtxt(nome_entrada)

x = data[:,0]
y = data[:,1]

interpola = interpolate.interp1d(x,y,kind="nearest",fill_value="extrapolate")

nome_saida = nome_entrada.replace(".txt", "_int.txt")
saida = open(nome_saida, "w")

E = 2
binsize = 0.1

for i in range(140):
	if(E <= x[-1]):
		saida.write("{:.1f} {:f}\n".format(E, interpola(E+binsize/2)))
		saida.write("{:.1f} {:f}\n".format(E+binsize, interpola(E+binsize/2)))
	else:
		saida.write("{:.1f} 0\n".format(E))
		saida.write("{:.1f} 0\n".format(E+binsize))
	E += binsize

saida.close()
