import numpy as np
import matplotlib.pyplot as plt

# Constants
d = 0.001          # m  (string diameter)
E = 2.0e11         # Pa
rho = 7850.0       # kg/m^3
f = 440.0          # Hz (target frequency)


A = 0.25 * np.pi * d**2         # cross-sectional area (m^2)
I = np.pi * d**4 / 64.0         # second moment of area (m^4)
mu = rho * A                    # linear mass density (kg/m)

# Frequency equation (Hz) with inharmonicity
def f_stiff(T, L):
    """Return frequency (Hz) including inharmonicity."""
    B = (np.pi**2 * E * I) / (T * L**2)
    return (1/ (2 * L)) * np.sqrt(T / mu) * np.sqrt(1 + B)

# Meshgrid: L on x-axis, T on y-axis (figure 3)
L_vals = np.linspace(0.2, 1.2, 100)      # m
T_vals = np.linspace(100, 1500, 100)   # N
LL, TT = np.meshgrid(L_vals, T_vals, indexing='xy')

# Compute frequency surface
FF = f_stiff(TT, LL)
FF_plane = 440 * np.ones_like(FF)       # reference plane at 440 Hz
"""
np.ones_like() returns an array with the
same shape and type as the given array
"""
# 3D Plot
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(111, projection='3d')

# Frequency surface
surf = ax.plot_surface(LL, TT, FF, cmap='viridis', alpha=0.9, linewidth=0)

# 440 Hz plane
ax.plot_surface(LL, TT, FF_plane, color='r', alpha=0.25, linewidth=0)

# Axis labels
ax.set_xlabel('String Length L [m]')
ax.set_ylabel('Tension T [N]')
ax.set_zlabel('Frequency f [Hz]')
ax.set_title('f(L, T) Surface with Inharmonicity (d = 1 mm)')
ax.view_init(elev=25, azim=-70)

plt.tight_layout()
plt.show()

# Stiff-string equation solved for T(L)
# f = (1/(2L)) * sqrt(T/mu) * sqrt(1 + (pi^2 * E * I)/(T * L^2))
# Rearranged closed form:
# T(L) = 4*mu*f^2*L^2 - (pi^2 * E * I)/(L^2)
def tension(L):
    return 4.0 * mu * (f**2) * (L**2) - (np.pi**2 * E * I) / (L**2)

# Range of string lengths
L_vals = np.linspace(0.1, 1.5, 400)  # m
T_vals = tension(L_vals)

# Plotting graph (figure 4)
plt.figure(figsize=(7,5))
plt.plot(L_vals, T_vals, color='navy', lw=2)
plt.xlabel('String Length L [m]')
plt.ylabel('Tension T [N]')
plt.title('Tension vs String Length for 440 Hz (with Inharmonicity)')
plt.grid(True)
plt.tight_layout()
plt.show()
        
