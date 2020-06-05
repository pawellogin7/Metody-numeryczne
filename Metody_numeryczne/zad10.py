import numpy as np
import matplotlib.pyplot as plt


def funkcja(y, t):
    return (3 * t - 4 * np.power(t, 2)) * np.sqrt(y)
    return y*(1 - y)*(2*t - 3)


def pochodna_analityczna(t):
    return np.power(((-1) * (4 / 3) * np.power(t, 3) + (3 / 4) * np.power(t, 2) + 1), 2)
    return np.power(np.e, np.square(t)) / (2*np.power(np.e, 3*t) + np.power(np.e, np.square(t)))


def pochodna_eulera(y0, t0, tk, n):
    h = (tk - t0) / n
    tn = np.arange(t0, tk + h, h)
    yn = np.zeros(len(tn))
    yn[0] = y0
    ad = np.arange(1, n)
    yn[ad] = yn[ad - 1] + funkcja(yn[ad - 1], tn[ad - 1]) * h
    for i in range(1, len(yn), 1):
        yn[i] = yn[i - 1] + funkcja(yn[i - 1], tn[i - 1]) * h
    return yn


def pochodna_hena(y0, t0, tk, n):
    h = (tk - t0) / n
    tn = np.arange(t0, tk + h, h)
    yn = np.zeros(len(tn))
    yn[0] = y0
    for id in range(1, len(yn), 1):
        y_temp = yn[id - 1] + funkcja(yn[id - 1], tn[id - 1]) * h
        yn[id] = yn[id - 1] + (funkcja(yn[id - 1], tn[id - 1]) + funkcja(y_temp, tn[id])) * h / 2
    return yn


def pochodna_srodkowy(y0, t0, tk, n):
    h = (tk - t0) / n
    tn = np.arange(t0, tk + h, h)
    yn = np.zeros(len(tn))
    yn[0] = y0
    for i in range(1, len(yn), 1):
        y_srod = yn[i - 1] + funkcja(yn[i - 1], tn[i - 1]) * (h / 2)
        dy_srod = funkcja(y_srod, (tn[i - 1] + tn[i])/2)
        yn[i] = yn[i - 1] + dy_srod*h
    return yn


t0 = 0
y0 = 1
n = 100
print('Podaj wartość tk, dla której liczona bedzie pochodna: ')
while True:
    try:
        tk = float(input())
        if 0 < tk < 14:
            break
        else:
            print('Podana wartość wychodzi poza przedział. Podaj nową wartość:')
    except ValueError:
        print('Podana wartość nie jest liczbą! Podaj prawidłową wartość:')

y_a = pochodna_analityczna(tk)
y_e = pochodna_eulera(y0, t0, tk, n)
y_h = pochodna_hena(y0, t0, tk, n)
y_s = pochodna_srodkowy(y0, t0, tk, n)

# print('Pochodna funkcji dla t = {} obliczona analitycznie: {}'.format(tk, y_a))
# print('Pochodna funkcji dla t = {} obliczona metodą Eulera: {}'.format(tk, y_e))
# print('Pochodna funkcji dla t = {} obliczona metodą Hena: {}'.format(tk, y_h))
# print('Pochodna funkcji dla t = {} obliczona metodą punktu środkowego: {}'.format(tk, y_s))
t = np.arange(t0, tk + (tk - t0)/n, (tk - t0)/n )
y_a = pochodna_analityczna(t)

plt.plot(t, y_a)
plt.plot(t, y_e, 'g--')
plt.plot(t, y_h, 'ro')
plt.plot(t, y_s, 'b.')
plt.show()
