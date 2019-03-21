#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego


import numpy as np

np.set_printoptions(precision=3)
np.set_printoptions(suppress=True)


def simplex(A, b, c):
    m, n = A.shape
    
    # Cria a tabela inicial
    t = np.zeros(shape=(m+1,m+n+1), dtype=np.float64)
    t[0, 0:n] = -c
    t[1:m+1,0:n] = A
    t[1:m+1, n:m+n] = np.eye(m)
    t[1:m+1, -1] = b
    
    #[x1 x2 x3 x4 x5]
    v = ["x" + str(i) for i in range(1, m+n+1)]
    
    print("-------------------------------------------------")
    print("Tabela inicial")
    print("-------------------------------------------------\n")
    print(t)
    
    k = 0
    
    # Se cada coeficiente da eq(0) for maior do que zero, então a solução é ótima
    while np.any(t[0,] < 0):
        
        # Seleciona o coeficiente de maior valor absoluto da eq(0), ou seja, selecionado a 
        # variável que entrará na base
        q = -1
        l = -1
        for i in range(0, m+n):
            z = abs(t[0][i])
            if z > l:
                l = z
                q = i
                
        print("\nVariável que entrará na base: {}".format(v[q]))
                
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
        
        print("Variável que sairá na base: {}".format(v[p+1]))

        # Aplica eliminação de Gauss-Jordan no elemento pivô (p, q)
        for i in range(0, m+1):
            for j in range(0, m+n+1):
                if i != p and j != q:
                    t[i][j] -= t[p][j] * t[i][q] / t[p][q]
        for i in range(0, m+1):
            if i != p:
                t[i][q] = 0.0
        for j in range(0, m+n+1):
            if j != q:
                t[p][j] /= t[p][q]
        t[p][q] = 1.0
        
        print("\n-------------------------------------------------")
        print("Tabela Interação - {}".format(k))
        print("-------------------------------------------------\n")
        print(t)
        k += 1
    print("\n-------------------------------------------------")
    print("Solução ótima: {}".format(t[0][m+n]))


A = np.array([[1,0],[0,2],[3,2]])
b = np.array([4,12,18])
c = np.array([3,5])

simplex(A, b,c)
