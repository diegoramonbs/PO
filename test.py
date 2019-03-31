#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego


from simplex import simplex


# Exemplo com restrinções do tipo igualdade 
A = np.array([[1,0],[0,2],[3,2]])
b = np.array([4,12,18])
c = np.array([3,5])
s = np.array(["<=", "<=", "="])
simplex(A, b, c, s, "min", verbose=False)


# Exemplo de problema ilimitado
A = np.array([[1,-1],[-1,1]])
b = np.array([4, 4])
c = np.array([-1, -1])
s = np.array(["<=", "<="])
simplex(A, b, c, s, "min", verbose=False)



# Exemplo mais geral, sendo um problema de minimização e com restrições de todos os tipos 
A = np.array([[0.3,0.1],[0.5,0.5],[0.6,0.4]])
b = np.array([2.7, 6, 6])
c = np.array([0.4, 0.5])
s = np.array(["<=", "=", ">="])
simplex(A, b, c, s, "min", verbose=False)


