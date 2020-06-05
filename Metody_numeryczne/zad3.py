import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la

E12 = 25
E23 = 50
E34 = 50
E35 = 25
Qa = 200
Ca = 2
Qb = 300
Cb = 2
Qc = 150
Qd = 350
Ws = 1500
Wg = 2500

A = np.array([[E12 + Qa, -E12, 0, 0, 0],
              [-E12 - Qa, E12 + E23 + Qa + Qb, -E23, 0, 0],
              [0, -E23 - Qa - Qb, E23 + E34 + E35 + Qa + Qb, -E34, -E35],
              [0, 0, -E34 - Qa - Qb + Qd, E34 + Qc, 0],
              [0, 0, -E35 - Qa - Qb + Qc, 0, E35 + Qd]])

B = np.array([Ws + Qa*Ca, Qb*Cb, 0, 0, Wg])

# Wyznaczanie macierzy L i U dla macierzy A
_, L, U = la.lu(A)

# Obliczanie wektora stężeń CO na początku
z = la.solve_triangular(L, B, lower=True)
x = la.solve_triangular(U, z)

# Obliczanie wektora stężeń CO po zmniejszeniu Ws i Wg
Ws1 = 800
Wg1 = 1200
B1 = np.array([Ws1 + Qa*Ca, Qb*Cb, 0, 0, Wg1])

z1 = la.solve_triangular(L, B1, lower=True)
x1 = la.solve_triangular(U, z1)

# Obliczanie macierzy odwrotnej metodą LU
temp = la.solve_triangular(L, [1, 0, 0, 0, 0], lower=True)
row1 = la.solve_triangular(U, temp)
temp = la.solve_triangular(L, [0, 1, 0, 0, 0], lower=True)
row2 = la.solve_triangular(U, temp)
temp = la.solve_triangular(L, [0, 0, 1, 0, 0], lower=True)
row3 = la.solve_triangular(U, temp)
temp = la.solve_triangular(L, [0, 0, 0, 1, 0], lower=True)
row4 = la.solve_triangular(U, temp)
temp = la.solve_triangular(L, [0, 0, 0, 0, 1], lower=True)
row5 = la.solve_triangular(U, temp)

A_inv = np.array([row1, row2, row3, row4, row5])
A_inv = np.transpose(A_inv)
A_inv_np = np.linalg.inv(A)

# Obliczanie procentowago udziału CO w pokoju dla dzieci na początku
proc_grill = (A_inv[3, 4] * Wg) / x[3]
proc_palacz = (A_inv[3, 0] * Ws) / x[3]
proc_ulica = (A_inv[3, 0] * Qa * Ca + A_inv[3, 1] * Qb * Cb) / x[3]

# Obliczanie procentowago udziału CO w pokoju dla dzieci po zmniejszeniu Wg i Ws
proc_grill1 = (A_inv[3, 4] * Wg1) / x1[3]
proc_palacz1 = (A_inv[3, 0] * Ws1) / x1[3]
proc_ulica1 = (A_inv[3, 0] * Qa * Ca + A_inv[3, 1] * Qb * Cb) / x1[3]

print('Na początku:')
print('Wektor stężenia CO w stanie ustalonym:')
print('|  P1  |  P1  |  P1  |  P1  |  P1  |')
print('|{:6.2f}|{:6.2f}|{:6.2f}|{:6.2f}|{:6.2f}|'.format(x[0], x[1], x[2], x[3], x[4]))
print('Udział grilla: {:.2f}%'.format(proc_grill*100))
print('Udział papierosow: {:.2f}%'.format(proc_palacz*100))
print('Udział ulicy: {:.2f}%'.format(proc_ulica*100))
print('\n')

print('Po obniżeniu Wg i Ws:')
print('Wektor stężenia CO w stanie ustalonym:')
print('|  P1  |  P1  |  P1  |  P1  |  P1  |')
print('|{:6.2f}|{:6.2f}|{:6.2f}|{:6.2f}|{:6.2f}|'.format(x1[0], x1[1], x1[2], x1[3], x1[4]))
print('Udział grilla: {:.2f}%'.format(proc_grill1*100))
print('Udział papierosow: {:.2f}%'.format(proc_palacz1*100))
print('Udział ulicy: {:.2f}%'.format(proc_ulica1*100))
print('\n')

np.set_printoptions(precision=3)
print('A-1 metoda LA:')
print(A_inv)
print('A-1 obliczone w numpy:')
print(A_inv_np)
