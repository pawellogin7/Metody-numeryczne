import numpy as np
import matplotlib.pyplot as plt


min_x = -3
max_x = 3

# Pierwiastki: (x1, y1) = (-1.45, -0.31); (x2, y2) = (-0.33, 0); (x3, y3) = (0.46, 1.65)
# Rozwiązanie graficzne funkcji
x1 = np.arange(min_x, max_x, 0.001)
y1 = np.square(x1) + 2*x1 + 0.5

x2a = np.arange(min_x, 3/7, 0.001)
y2a = 2*np.square(x2a) / (7*x2a - 3)
x2b = np.arange(3/7 + 0.001, max_x, 0.001)
y2b = 2*np.square(x2b) / (7*x2b - 3)

fun1, = plt.plot(x1, y1, 'b-')
fun2, = plt.plot(x2a, y2a, 'g-')
plt.plot(x2b, y2b, 'g-')
plt.ylim(-3, 3)

plt.legend((fun1, fun2), ('y = x^2 + 2x + 0.5', 'Y = 2x^2 / (7x - 3)'))
plt.title('Rozwiązanie graficzne układu równań')

print('Rozwiązania układu równań to: (x1, y1) = (-1.45, -0.31); (x2, y2) = (-0.33, 0); (x3, y3) = (0.46, 1.65)')

# Metoda iteracyjnego podstawiania
# Wyznaczamy rozwiązanie metody dla pierwiastka (x3, y3) = (0.46, 1.65)
max_iter = 10
x0 = 1
y0 = 0.5
xi = x0
yi = y0
for i in range(max_iter):
    yi = xi*xi + 2*xi + 0.5
    xi = (3*yi + 2*xi*xi) / (7*yi)

print('Rozwiązanie układu równań dla metody iteracyjnego podstawiania dla punktu startowego ({}, {}) '
      'po {} iteracjach wynosi ({:.3f}, {:.3f})'.format(x0, y0, max_iter, xi, yi))


# Metoda Newtona-Raphsona
# Wyznaczamy rozwiązanie metody dla 3 pierwiastków:
# (x1, y1) = (-1.45, -0.31)
# (x2, y2) = (-0.33, 0
# (x3, y3) = (0.46, 1.65)
max_iter = 10
max_pierw = 3
x0 = [2, 0, -2]
y0 = [2, 0, -1]
for i in range(max_pierw):
    xi = x0[i]
    yi = y0[i]
    for j in range(max_iter):
        f1i_x = -2*xi - 2
        f1i_y = 1
        f2i_x = 4*xi - 7*yi
        f2i_y = 3 - 7*xi
        det = f1i_x * f2i_y - f1i_y * f2i_x
        f1i = yi - xi*xi - 2*xi - 0.5
        f2i = 3*yi + 2*xi*xi - 7*xi*yi
        xi = xi - (f1i*f2i_y - f2i*f1i_y) / det
        yi = yi - (f2i*f1i_x - f1i*f2i_x) / det
    print('Rozwiązanie układu równań dla metody N-R dla punktu startowego ({}, {}) '
          'po {} iteracjach wynosi ({:.3f}, {:.3f})'.format(x0[i], y0[i], max_iter, xi, yi))

plt.show()
