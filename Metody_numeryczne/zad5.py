import numpy as np
import matplotlib.pyplot as plt


# Funkcja zwracająca transpozycję macierzy. Używana dla wektorów, jako że musimy najpierw nadać im 2 wymiary.
def transpose(mat):
    return np.transpose(np.atleast_2d(mat))


A = np.array([[4.4, -4.65, 1.35],
              [1, 0, 0],
              [0, 1, 0]])

B = np.array([1, 0, 0])
C = np.array([1, 1, 1])
D = 0

# Rysowanie odpowiedzi układu na sygnał skokowy
liczba_probek = 50
x0 = 0
n = np.arange(0, liczba_probek, 1)
x = np.zeros((liczba_probek, 3))
u = np.ones(liczba_probek)
y = np.zeros(liczba_probek)
x[0, :] = 0
for i in range(liczba_probek):
    if i < liczba_probek - 1:
        x[i + 1] = np.dot(A, x[i]) + B * u[i]
    y[i] = np.dot(C, x[i]) + D * u[i]

fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 8))
ax1.plot(n, y)
ax1.set_title("Odpowiedź skokowa układu bez sterownika LQR - układ jest niestabilny")


# Obliczanie parametrów sterownika LQR
# Zwiększanie stałej c1 zwiększa szybkość układu ale go destabilizuje i powoduje większe przeregulowania
# Dodatkowo wzrost stałej c1 powoduje zmniejszenie się amplitudy sygnału wyjściowego
# Zwiększanie stałej c2 zwiększa zapas stabilności układu, ale spowalnia układ
# Dodatkowo wzrost stałej c1 powoduje zwiększenie się amplitudy sygnału wyjściowego
max_iter = 10
c1 = 1
c2 = 1
Q = c1 * np.identity(3)
R = c2
P = np.zeros((3, 3))

# Obliczanie macierzy P w sterowniku
# Poprawność wyliczeń macierzy P zgodnie ze wzorem z wykładu 6 została sprawdzona w programie Octave
B = transpose(B)
for i in range(max_iter):
    mat1 = np.power(R + np.dot(np.dot(transpose(B), P), B), -1)
    mat2 = np.dot(P, np.dot(B, mat1))
    mat3 = P - np.dot(mat2, np.dot(transpose(B), P))
    P = Q + np.dot(np.dot(transpose(A), mat3), A)
np.set_printoptions(precision=3)
print('Macierz P sterownika obliczona dla {} iteracji oraz wartości stałych c1 = {} i c2 = {}:'.format(max_iter, c1, c2))
print(P)
print('\n')

# Implementacja sterownika do układu ze sprzężeniem od stanu
# Obliczanie nowej macierzy A
# Poprawność wyliczeń macierzy F zgodnie ze wzorem z wykładu 6 została sprawdzona w programie Octave
F = np.power(R + np.dot(transpose(B), np.dot(P, B)), -1)
F = np.dot(F, np.dot(transpose(B), np.dot(P, A)))
A1 = A - np.dot(B, F)
np.set_printoptions(formatter={'float': lambda l: "{0:0.6f}".format(l)})
print('Nowa macierz A układu po implementacji sterownika LQR i sprzężenia od stanu:'.format(max_iter, c1, c2))
print(A1)

# Rysowanie odpowiedzi skokowej układu ze sterownikiem LQR
B = transpose(B)
liczba_probek = 50
n1 = np.arange(0, liczba_probek, 1)
u1 = np.ones(liczba_probek)
x1 = np.zeros((liczba_probek, 3))
y1 = np.zeros(liczba_probek)
x1[0, :] = 0
for i in range(liczba_probek):
    if i < liczba_probek - 1:
        x1[i + 1] = np.dot(A1, x1[i]) + B * u1[i]
    y1[i] = np.dot(C, x1[i]) + D * u1[i]

wyj, = ax2.plot(n1, y1)
wej, = ax2.plot(n1, u1)
ax2.legend((wej, wyj), ('Sygnał sterujący', 'Odpowiedź ukł1du'))
ax2.set_title('Odpowiedź skokowa układu ze sterownikiem LQR dla c1 = {} i c2 = {}'.format(c1, c2))

plt.show()
