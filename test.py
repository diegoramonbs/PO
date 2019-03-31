#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Author: Diego


from simplex import simplex

# Exemplo padrão do slide  (Valor ótimo: 36.0)
A = np.array([[1,0],[0,2],[3,2]])
b = np.array([4,12,18])
c = np.array([3, 5])
s = np.array(["<=", "<=", "<="])

# Exemplo com restrições do tipo igualdade  (Valor ótimo: 36.0)
A = np.array([[1,0],[0,2],[3,2]])
b = np.array([4,12,18])
c = np.array([3,5])
s = np.array(["<=", "<=", "="])
simplex(A, b, c, s, "min", verbose=False)

# Exemplo de problema ilimitado (Exception)
A = np.array([[1,-1],[-1,1]])
b = np.array([4, 4])
c = np.array([-1, -1])
s = np.array(["<=", "<="])
simplex(A, b, c, s, "min", verbose=False)


# Exemplo mais geral, sendo um problema de minimização e com restrições de todos os tipos (Valor ótimo: -5.25)
A = np.array([[0.3,0.1],[0.5,0.5],[0.6,0.4]])
b = np.array([2.7, 6, 6])
c = np.array([0.4, 0.5])
s = np.array(["<=", "=", ">="])
simplex(A, b, c, s, "min", verbose=False)


