import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('dane/measurements23.txt')
px = data[:, 0]
py = data[:, 1]

T = 1
F = np.array([[1, 0, T, 0],
             [0, 1, 0, T],
             [0, 0, 1, 0],
             [0, 0, 0, 1]])

G = np.array([[0, 0],
              [0, 0],
              [1, 0],
              [0, 1]])

H = np.array([[1, 0, 0, 0],
              [0, 1, 1, 1]])

P = 5 * np.identity(4)
Q = 0.25 * np.identity(2)
R = 2 * np.identity(2)
s = np.array([[px[0]], [py[0]], [0], [0]])

zx = [0]
zy = [0]
# Implementacja filtru Kalmana
for i in range(1, len(px)):
    s = np.dot(F, s)
    P = np.dot(np.dot(F, P), np.transpose(F)) + np.dot(np.dot(G, Q), np.transpose(G))
    z = np.dot(H, s)
    e = np.array([[px[i]], [py[i]]]) - z
    S = np.dot(np.dot(H, P), np.transpose(H)) + R
    K = np.dot(np.dot(P, np.transpose(H)), np.linalg.inv(S))
    s = s + np.dot(K, e)
    P = np.dot((np.identity(4) - np.dot(K, H)), P)
    zx.append(z[0])
    zy.append(z[1])


zx_p = [zx[len(zx) - 1]]
zy_p = [zy[len(zy) - 1]]
# Predykcja przyszłych 5 próbek za pomocą wyznaczonego filtru Kalmana
for i in range(5):
    s = np.dot(F, s)
    P = np.dot(np.dot(F, P), np.transpose(F)) + np.dot(np.dot(G, Q), np.transpose(G))
    z = np.dot(H, s)
    zx_p.append(z[0])
    zy_p.append(z[1])

zx = np.asarray(zx)
zy = np.asarray(zy)
zx_p = np.asarray(zx_p)
zy_p = np.asarray(zy_p)
line1, = plt.plot(px, py, 'rx')
line2, = plt.plot(zx, zy, 'b-')
line3, = plt.plot(zx_p, zy_p, 'b--')
plt.plot(zx_p[-1], zy_p[-1], 'bo')
plt.legend([line1, line2, line3], ['Punkty pomiarowe', 'Wygładzona trajektoria', 'Przewidywana przyszła trajektoria'])
plt.show()
