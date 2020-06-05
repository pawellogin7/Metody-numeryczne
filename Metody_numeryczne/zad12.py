import numpy as np
import matplotlib.pyplot as plt


Ta = 100
h = 0.05
dx = 2
dy = 2

border_temp_up = np.arange(400, 290, -10)
border_temp_left = np.arange(400, 290, -10)
border_temp_down = np.arange(300, 190, -10)
border_temp_right = np.arange(300, 190, -10)


A = np.zeros((81, 81))
for i in range(len(A)):
    A[i, i] -= 2 / np.square(dx) + 2 / np.square(dx) + h
    if i % 9 != 0:
        A[i, i - 1] += 1 / np.square(dx)
    if i % 9 != 8:
        A[i, i + 1] += 1 / np.square(dx)
    if i // 9 != 0:
        A[i, i - 9] += 1 / np.square(dy)
    if i // 9 != 8:
        A[i, i + 9] += 1 / np.square(dy)

B = np.zeros(81)
for i in range(len(B)):
    B[i] -= h * Ta
    if i % 9 == 0:
        B[i] -= border_temp_left[i % 9 + 1] / np.square(dx)
    if i % 9 == 8:
        B[i] -= border_temp_right[i % 9 + 1] / np.square(dx)
    if i // 9 == 0:
        B[i] -= border_temp_up[i % 9 + 1] / np.square(dy)
    if i // 9 == 8:
        B[i] -= border_temp_down[i % 9 + 1] / np.square(dy)

T = np.linalg.solve(A, B).reshape(9, 9)

dim = np.mgrid[2:22:2, 2:22:2]
dim_x = dim[1]
dim_y = dim[0]

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
plot = ax.pcolor(dim_x, dim_y, T, cmap='inferno', vmin=150, vmax=350)
plt.Axes.set_xlabel(ax, 'X')
plt.Axes.set_ylabel(ax, 'Y')
plt.Axes.set_title(ax, 'Temperatura na płytce w stanie ustalonym')
cbar = fig.colorbar(plot)
cbar.set_label('Temperatura w stopniach')

plt.show()
