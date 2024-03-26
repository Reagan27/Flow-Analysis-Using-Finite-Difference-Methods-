import numpy as np
import matplotlib.pyplot as plt


Lx = 1.0
Ly = 1.0
Nx = 50
Ny = 50
U = 1.0
R = 1.0


x = np.linspace(0, Lx, Nx)
y = np.linspace(0, Ly, Ny)
X, Y = np.meshgrid(x, y)
dx = Lx / (Nx - 1)
dy = Ly / (Ny - 1)


psi = np.sin(2 * np.pi * X) * np.cos(np.pi * Y)
omega = np.cos(2 * np.pi * X) * np.sin(np.pi * Y)
def psi_x(psi, i, j):
    return (psi[i, j+1] - psi[i, j]) / dx

def psi_y(psi, i, j):
    return (psi[i+1, j] - psi[i, j]) / dy

def omega_x(omega, i, j):
    return (omega[i, j+1] - omega[i, j]) / dx

def omega_y(omega, i, j):
    return (omega[i+1, j] - omega[i, j]) / dy


def update_psi(psi, omega, i, j):
    return psi[i, j] + dx * dy * (omega_x(omega, i, j) - psi_x(psi, i, j))

def update_omega(psi, omega, i, j):
    return omega[i, j] + dx * dy * (R * (omega_x(omega, i, j) / dx + omega_y(omega, i, j) / dy)
                                    - psi_x(psi, i, j) * psi_y(psi, i, j))


iterations = 100
for _ in range(iterations):
    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            psi[i, j] = update_psi(psi, omega, i, j)


for _ in range(iterations):
    for i in range(1, Ny - 1):
        for j in range(1, Nx - 1):
            omega[i, j] = update_omega(psi, omega, i, j)


plt.figure(figsize=(12, 5))


plt.subplot(1, 2, 1)
plt.contourf(X, Y, psi, cmap='viridis')
plt.colorbar(label='Stream Function')
plt.title('Stream Function')


plt.subplot(1, 2, 2)
plt.contourf(X, Y, omega, cmap='viridis')
plt.colorbar(label='Vorticity')
plt.title('Vorticity')

plt.tight_layout()
plt.show()
