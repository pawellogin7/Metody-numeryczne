import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize


# Obliczanie odp. skokowej układu oscylacyjnego II rzędu dla danej chwili czasu i określonych parametrach
def odp_skok(params):
    t, tau, zeta = params
    odp = 1 - np.exp(-zeta*t/tau) / np.sqrt(1 - np.square(zeta)) * \
          np.sin((np.sqrt(1 - np.square(zeta)) * t / tau + np.arctan(np.sqrt(1 - np.square(zeta)) / zeta)))
    return odp


# Obliczanie odp. impulsowej układu oscylacyjnego II rzędu dla danej chwili czasu i określonych parametrach
def odp_imp(params):
    t, tau, zeta = params
    odp = np.exp(-zeta*t/tau) / tau * np.sqrt(1 - np.square(zeta)) * np.sin(np.sqrt(1 - np.square(zeta)) * t / tau)
    return odp


# Obliczanie odp. skokowej układu danego w zadaniu dla danej chwili czasu i określonych parametrach
def odp_skokowa_uklad(params):
    t, k, tau_z, tau, zeta = params
    odp = k * (tau_z * odp_imp([t, tau, zeta]) + odp_skok([t, tau, zeta]))
    return odp


# Regresja nieliniwa z kryterium najmniejszych kwadratów dla odpowiedzi skokowej układu
def regresja(params):
    k, tau_z, tau, zeta = params
    sum = 0
    for t, A in zip(time, amp):
        # Dla każdej chwili czasu od odp. skokowej z danych odejmujemy odp. skokową wyznaczoną z wzoru
        sum += np.square(A - odp_skokowa_uklad([t, k, tau_z, tau, zeta]))
    return sum


# Otwieranie pliku i wczytywanie danych
with open('dane/data16.txt', 'r') as file:
    data = []
    for line in file:
        data_temp = []
        for x in line.split():
            data_temp.append(float(x))
        data_temp = np.asarray(data_temp)
        data.append(data_temp)
    data = np.asarray(data, dtype=np.float32)
    time = data[:, 0]
    amp = data[:, 1]


# Poszukiwanie parametrów układu
min = optimize.fmin(regresja, [1.0, -1.0, 1, 0.5])
print('Parametry układu:')
print('k = {}'.format(min[0]))
print('tau z = {}'.format(min[1]))
print('tau = {}'.format(min[2]))
print('zeta = {}'.format(min[3]))

# Rysowanie teoretycznej odp. skokowej układu
odp_s = []
for t in time:
    odp_s.append(odp_skokowa_uklad([t, min[0], min[1], min[2], min[3]]))
odp_s = np.asarray(odp_s)

line1, = plt.plot(time, amp)
line2, = plt.plot(time, odp_s)
plt.legend([line1, line2], ['Odp. skokowa układu z danych', 'Odp. skokowa układu wyznaczonego'])
plt.show()
