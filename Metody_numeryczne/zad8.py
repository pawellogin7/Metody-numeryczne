import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('dane/data23.txt')
x = data[:, 0]
y = data[:, 1]

# Z powodu, że w otrzymanych danych znajduje się 11 punktów należy wyznaczyć wielomian interpolacyjny Newtona 10 rzędu.
# Istnieje tylko jeden taki wielomian 10 rzędu, który będzie przechodził przez wszystkie 11 punktów
def interpolacja_newtona(points_x, points_y):
    x = np.asarray(points_x)
    y = np.asarray(points_y)
    y_new = np.copy(y)
    x_new = np.copy(x)
    b = []
    b.append(y[0])
    for leng in range(len(x) - 1, 0, -1):
        y_cp = np.copy(y_new)
        y_new = np.zeros(leng)
        for i in range(leng):
            leng1 = len(x) - leng
            y_new[i] = (y_cp[i + 1] - y_cp[i]) / (x[i + leng1] - x_new[i])
        b.append(y_new[0])
    return b


# Używamy obliczonych wyżej wag wielomianu do wyznaczenia przebiegu funkcji interpolacyjnej
def wielomian(x, b):
    x1 = np.arange(np.amin(x), np.amax(x), 0.001)
    y1 = np.zeros_like(x1)
    for i in range(len(b)):
        iloczyn = np.ones_like(y1)
        for j in range(i):
            iloczyn *= x1 - x[j]
        y1 += b[i]*iloczyn
    return x1, y1


# Wyznaczamy wagi funkcji sklejanych 3 rzędu na każdym odcinku pomiędzy dwoma punktami
def funkcje_sklejane(points_x, points_y):
    x = np.asarray(points_x)
    y = np.asarray(points_y)
    A = np.zeros((len(x), len(x)))
    B = np.zeros(len(x))
    A[0, 0] = 1
    A[len(x) - 1, len(x) - 1] = 1
    for i in range(1, len(x) - 1):
        A[i, i - 1] = x[i] - x[i - 1]
        A[i, i] = 2*((x[i] - x[i - 1]) + (x[i + 1] - x[i]))
        A[i, i + 1] = x[i + 1] - x[i]
        B[i] = 3*((y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]))

    ci = np.linalg.solve(A, B)
    index = np.arange(len(x) - 1, dtype=np.uint8)
    hi = np.zeros(len(index))
    ai = np.zeros(len(index))
    bi = np.zeros(len(index))
    di = np.zeros(len(index))
    hi[index] = x[index + 1] - x[index]
    ai[index] = y[index]
    bi[index] = ((y[index + 1] - y[index]) / hi[index]) - (hi[index] / 3) * (2*ci[index] + ci[index + 1])
    di[index] = (ci[index + 1] - ci[index]) / (3*hi[index])
    return ai, bi, ci, di


# Używamy obliczonych wyżej wag funkcji sklejanych do wyznaczenia przebiegu funkcji interpolacyjnej
# Dla każdego punktu pi, x[i] <= pi < x[i+1] używamy wag a[i], b[i], c[i] oraz d[i]
def wielomian_sklejane(x, a, b, c, d):
    x = np.asarray(x)
    a = np.asarray(a)
    b = np.asarray(b)
    c = np.asarray(c)
    d = np.asarray(d)
    x1 = np.arange(np.amin(x), np.amax(x), 0.001)
    index = np.zeros_like(x1, dtype=np.uint8)
    for i in range(len(x)):
        index[x1 >= x[i]] = i
    y1 = a[index] + b[index]*(x1 - x[index]) + c[index]*np.power((x1 - x[index]), 2) + d[index]*np.power((x1 - x[index]), 3)
    return x1, y1


b = interpolacja_newtona(x, y)
x1, y1 = wielomian(x, b)

a1, b1, c1, d1 = funkcje_sklejane(x, y)
x2, y2 = wielomian_sklejane(x, a1, b1, c1, d1)

line1, = plt.plot(x, y, 'ro')
line2, = plt.plot(x1, y1, 'g-')
line3, = plt.plot(x2, y2, 'b-')
plt.legend([line1, line2, line3],
           ['Punkty z danych', 'Interpolacja Newtona wielomianem 10 rzędu', 'Interpolacja funkcjami sklejanymi 3 rzędu'])
plt.show()
