import numpy as np


# Funkcja oblicza wartość wielomianu max 5 stopnia w punkcie x
def wielomian(a5, a4, a3, a2, a1, a0, x):
    return a5*np.power(x, 5) + a4*np.power(x, 4) + a3*np.power(x, 3) + a2*np.power(x, 2) + a1*x + a0


# Funkcja oblicza całkę metodą trapezów dla wielomianu max 5 stopnia, dla przedziału <a, b> oraz n segmentów
def całka_trapezy(a5, a4, a3, a2, a1, a0, a, b, n):
    h = (b - a) / n
    x_tab = np.arange(a, b + h, h)
    y1 = wielomian(a5, a4, a3, a2, a1, a0, x_tab[:-1])
    y2 = wielomian(a5, a4, a3, a2, a1, a0, x_tab[1:])
    return np.sum((y1 + y2)*h / 2)


# Dane
a4 = -0.06543
a3 = 0.6341
a2 = -0.7652
a1 = -3.65
a0 = 1.543
x1 = -2
x2 = 7

# Całka obliczona analitycznie
aw_5 = a4 / 5
aw_4 = a3 / 4
aw_3 = a2 / 3
aw_2 = a1 / 2
aw_1 = a0 / 1
calka_wynik = wielomian(aw_5, aw_4, aw_3, aw_2, aw_1, 0, x2) - wielomian(aw_5, aw_4, aw_3, aw_2, aw_1, 0, x1)
print("Całka obliczona analitycznie: {}".format(np.round(calka_wynik, 6)))
print("Wspolczynniki wielomianu wynikowego całki:")
aw = np.round([aw_5, aw_4, aw_3, aw_2, aw_1], 3)
print("||a5 = {}||a4 = {}||a3 = {}||a2 = {}||a1 = {}||\n".format(aw[0], aw[1], aw[2], aw[3], aw[4]))

# Całka obliczona metodą Romberga
print("Tabela przedstawiająca wyniki całek obliczonych numerycznie(metoda Romberga i  trzypunktowekwadratury Gaussa):")
print("||Metoda ||Iter||  Wynik  || Błąd[%] ||")
print("||-------++----++---------++---------||")
blad = 1
iter = 1
while blad > 0.2:
    calka_t = całka_trapezy(0, a4, a3, a2, a1, a0, x1, x2, np.power(2, iter,))
    calka_t2 = całka_trapezy(0, a4, a3, a2, a1, a0, x1, x2, np.power(2, iter + 1))
    calka_romberg = 4/3*calka_t2 - 1/3*calka_t
    blad = np.abs((calka_wynik - calka_romberg) / calka_wynik) * 100
    print("||Romberg|| {}  ||{:0.6f}||{:<9}||".format(iter, calka_romberg, np.round(blad, 5)))
    iter += 1

# Całka obliczona metodą trzypunktowej kwadratury Gaussa
c = np.array([5/9, 8/9, 5/9])
t = np.array([-1*np.sqrt(3/5), 0, np.sqrt(3/5)])
x_gauss = ((x2 + x1) + (x2 - x1)*t) / 2
dx_gauss = (x2 - x1) / 2
calka_gauss = np.sum(c*wielomian(0, a4, a3, a2, a1, a0, x_gauss)*dx_gauss)
blad = np.abs((calka_wynik - calka_gauss) / calka_wynik) * 100
# print("Całka obliczona metodą trzypunktowej kwadratury Gaussa:")
# print("||  Wynik  || Błąd[%] ||")
print("|| Gauss || -- ||{:0.6f}||{:<9}||".format(calka_gauss, np.round(blad, 14)))

