import numpy as np

mais = 1

bg_total = np.zeros(140)

nome_padrao = "JUNO_BG_FLAG_int.txt"

while mais == 1:
	nome_entrada = input("Entre com o tipo de BG: ")
	nome_entrada = nome_padrao.replace("FLAG", nome_entrada)
	data = np.loadtxt(nome_entrada)
	bg = data[:,1]
	for i in range(140):
		bg_total[i] += bg[2*i]
	mais = int(input("Mais arquivos?: "))

saida = open("JUNO_BG_total.txt", "w")

for i in range(140):
	saida.write("{:f}\n".format(bg_total[i]))
