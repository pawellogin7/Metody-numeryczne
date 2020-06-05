import numpy as np
import matplotlib.pyplot as plt


def dx_dt(x, y, z):
    return -10*x + 10*y


def dy_dt(x, y, z):
    return 28*x - y - x*z


def dz_dt(x, y, z):
    return -8/3*z + x*y


def metoda_RK(t, h, x0, y0, z0):
    x = np.zeros_like(t)
    y = np.zeros_like(t)
    z = np.zeros_like(t)
    x[0] = x0
    y[0] = y0
    z[0] = z0
    # Obliczanie wartości zmiennych x, y i z dla całego przedziału t
    for i in range(len(t) - 1):
        # k1
        k1_x = dx_dt(x[i], y[i], z[i])
        k1_y = dy_dt(x[i], y[i], z[i])
        k1_z = dz_dt(x[i], y[i], z[i])
        # k2
        x_temp = x[i] + k1_x * h / 2
        y_temp = y[i] + k1_y * h / 2
        z_temp = z[i] + k1_z * h / 2
        k2_x = dx_dt(x_temp, y_temp, z_temp)
        k2_y = dy_dt(x_temp, y_temp, z_temp)
        k2_z = dz_dt(x_temp, y_temp, z_temp)
        # k3
        x_temp = x[i] + k2_x * h / 2
        y_temp = y[i] + k2_y * h / 2
        z_temp = z[i] + k2_z * h / 2
        k3_x = dx_dt(x_temp, y_temp, z_temp)
        k3_y = dy_dt(x_temp, y_temp, z_temp)
        k3_z = dz_dt(x_temp, y_temp, z_temp)
        # k4
        x[i+1] = x[i] + k3_x * h
        y[i+1] = y[i] + k3_y * h
        z[i+1] = z[i] + k3_z * h
        k4_x = dx_dt(x[i+1], y[i+1], z[i+1])
        k4_y = dy_dt(x[i+1], y[i+1], z[i+1])
        k4_z = dz_dt(x[i+1], y[i+1], z[i+1])
        # x[i+1], y[i+1], z[i+1]
        x[i+1] = x[i] + 1/6*(k1_x + 2*k2_x + 2*k3_x + k4_x)*h
        y[i+1] = y[i] + 1 / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y) * h
        z[i+1] = z[i] + 1 / 6 * (k1_z + 2 * k2_z + 2 * k3_z + k4_z) * h
    return x, y, z


# inicjalizacja zmiennych
t0 = 0
tk = 25
h = 0.03125
x0 = y0 = z0 = 5
t = np.arange(t0, tk + h, h)

# obliczanie x, y, z za pomocą metody RK-4
x_rk, y_rk, z_rk = metoda_RK(t, h, x0, y0, z0)

# Rysowanie przebiegów x, y i z w czasie
fig, axes = plt.subplots(3, 1)
fig.set_size_inches(8, 6)
line1, = axes[0].plot(t, x_rk, 'b-')
line2, = axes[1].plot(t, y_rk, 'g-')
line3, = axes[2].plot(t, z_rk, 'r-')
plt.Axes.set_title(axes[0], 'Przebieg zmiennej x w czasie')
plt.Axes.set_title(axes[1], 'Przebieg zmiennej y w czasie')
plt.Axes.set_title(axes[2], 'Przebieg zmiennej z w czasie')

# Rysowanie trajektorii fazowych w przestrzeni 3D
fig3d = plt.figure()
ax = fig3d.add_subplot(111, projection='3d')
ax.plot(x_rk, y_rk, z_rk)
plt.Axes.set_title(ax, 'Trajektoria fazowa w przestrzeni trójwymiarowej')

plt.show()
