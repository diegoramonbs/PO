#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego


import numpy as np

np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)

def simplex(A, b, c, s, type="max", M=100):
	m, n = A.shape

	t = np.zeros(shape=(m+1, m+n+1), dtype=np.float64)

	t[0, 0:n] = -c if type == "max" else -c * (-1)
	t[1:m+1,0:n] = A
	t[1:m+1, n:m+n] = np.eye(m)
	t[1:m+1, -1] = b

	for i, k in enumerate(s):
		if k == "<=" or k == "=<":
			continue
		elif k == "=":
			# Aplica a regra do "Big-M": introduz uma variável artificial de valor muito grande 
			# para ser possível fechar as dimensões das variáveis básicas, ou seja, faz aparecer 
			# a matriz identidade. 
			t[0, i+n] = M

			# Deixa a função objetivo em função das variáveis originais 
			t[0,:] = t[0,:] - M*t[i+1]
		elif k == ">=" or k == "=>":
			t[0, i+n] = M
			t[0,:] = t[0,:] - M*t[i+1]
			t[i+1, i+n] = -1
			z = np.zeros((m+1, 1), dtype=np.float64)
			z[-1] = 1
			t = np.column_stack((t, z))
			t[:,-2],  t[:,-1] = t[:,-1],t[:,-2].copy()
			t[0, i+n] = M
			n += 1
		else:
			raise Exception("Símbolo Inválido!")

	print("-------------------------------------------------")
	print("Tabela inicial")
	print("-------------------------------------------------\n")
	print(t)

	# Interação atual
	k = 0

	# Lista de variáveis para visualização do processo de entrada e saída da base
	v = ["x" + str(i) for i in range(1, m+n+1)]

	# Se cada coeficiente da eq(0) for maior do que zero, então a solução é ótima
	while np.any(t[0,0:-1] < 0):
	    
	    # Seleciona o coeficiente de maior valor absoluto da eq(0), ou seja, selecionado a 
	    # variável que entrará na base
	    q = -1
	    l = -1
	    for i in range(0, m+n):
	    	z = t[0][i]
	    	if z > 0:
	    		continue
	        z = abs(z)
	        if z > l:
	            l = z
	            q = i
         
	    print("\nVariável que entrará na base: {}, valor: {}".format(v[q], t[0][q]))
	    # Determina a coluna pivô pelo teste da razão mínima       
	    p = -1
	    for i in range(1, m+1):
	        # Avoid zero division
	        if t[i][q] <= 0:
	            continue
	        elif p == -1:
	            p = i
	        elif t[i][-1]/t[i][q] < t[p][-1] / t[p][q]:
	            p = i
	    
	    print("Variável que sairá na base: {}, valor: {}".format(v[p+1], t[p][q]))
	    # Aplica eliminação de Gauss-Jordan no elemento pivô (p, q)
	    for i in range(0, m+1):
	        for j in range(0, m+n+1 ):
	            if i != p and j != q:
	                t[i][j] -= t[p][j] * t[i][q] / t[p][q]
	    for i in range(0, m+1):
	        if i != p:
	            t[i][q] = 0.0
	    for j in range(0, m+n+1 ):
	        if j != q:
	            t[p][j] /= t[p][q]
	    t[p][q] = 1.0

	    v[p+1], v[q] = v[q],  v[p+1]
	    print(v)
	    
	    print("\n-------------------------------------------------")
	    print("Tabela Interação - {}".format(k))
	    print("-------------------------------------------------\n")
	    print(t)
	    k += 1
	print("\n-------------------------------------------------")
	print("Solução ótima: {}".format(t[0][m+n]))
