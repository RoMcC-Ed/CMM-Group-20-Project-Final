import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# Material constants for the model
E = 2.0e11        # Young's modulus [Pa]
rho = 7850.0      # density [kg/m^3]
sigma_y = 2.1e9   # yield stress of music wire [Pa]

# Defining length, working stress, target frequency, and diameter
def g_L(L, sigma_w, f_target, d):
    """
    Nonlinear equation for L:
    4 rho f^2 L^4 - sigma_w L^2 - (pi^2 E d^2)/16 = 0
    """
    return 4.0 * rho * f_target**2 * L**4 - sigma_w * L**2 - (np.pi**2 * E * d**2) / 16.0

def solve_L(f_target, SF, d):
    """
    For a given frequency f_target (Hz), safety factor SF and diameter d (m),
    solve g(L, sigma_w) = 0 for the speaking length L using Brent's method.
    """
    A = 0.25 * np.pi * d**2       # cross-sectional area
    sigma_w = sigma_y / SF        # working stress

    # Initial bracket for L [m]
    L_min, L_max = 0.1, 2.0

    g_min = g_L(L_min, sigma_w, f_target, d)
    g_max = g_L(L_max, sigma_w, f_target, d)

    # Expand the bracket until a sign change is found
    while g_min * g_max > 0:
        L_max += 0.5
        if L_max > 5.0:
            raise RuntimeError("Could not bracket a root for L - check inputs.")
        g_max = g_L(L_max, sigma_w, f_target, d)

    # Brent's root finding for L*. Attempting inverse quadratic interpolation, secant, or bisection.
    L_star = brentq(g_L, L_min, L_max, args=(sigma_w, f_target, d))
    return L_star

# Real data table (Figure 10) 
# Data array for frequency
f = np.array([
    880.00, 830.61, 784.00, 740.00,
    698.46, 659.26, 622.25, 587.33,
    554.37, 523.25, 493.88, 466.16,
    440.00, 415.30, 392.00, 369.99,
    349.23, 329.63, 311.13, 293.66,
    277.18, 261.63, 246.94, 233.08,
    220.00
])

# Data array for length in cm
L_cm = np.array([
    22.22, 23.52, 24.95, 26.50,
    27.20, 28.90, 30.56, 32.50,
    34.35, 36.40, 37.50, 39.80,
    42.20, 44.70, 47.30, 50.10,
    51.80, 54.80, 58.20, 61.60,
    63.70, 67.30, 69.70, 73.80,
    76.30
])

# Convert length from cm to m
L = L_cm / 100.0


# Linear regression: f vs 1/L  (Figure 8)
x_lin = 1.0 / L        # independent variable x = 1/L
y = f                  # dependent variable y = f

# Least-squares linear regression: y ≈ a1 * x + a0
coef_lin = np.polyfit(x_lin, y, 1)
a1, a0 = coef_lin

# Fitted line for x values and graph plotting
x_fit_lin = np.linspace(x_lin.min(), x_lin.max(), 300)
f_fit_lin = a1 * x_fit_lin + a0

plt.figure(figsize=(6, 4))
plt.plot(x_lin, y, 'ko', label='Scale data')              # scatter points
plt.plot(x_fit_lin, f_fit_lin, 'b-', label='Linear fit')  # fitted straight line
plt.xlabel('1 / L (1/m)')
plt.ylabel('Frequency f (Hz)')
plt.grid(True)
plt.legend()
plt.tight_layout()

# Quadratic regression: f vs L  (Figure 9)
# Least-squares quadratic regression: f(L) ≈ a * L^2 + b * L + c
coef_quad = np.polyfit(L, f, 2)
a, b, c = coef_quad

# Fitting curve
L_fit = np.linspace(L.min(), L.max(), 300)
f_fit_quad = a * L_fit**2 + b * L_fit + c

# Modelling curve: f–L relation from root-finding model
SF_model = 2.5          # choosing a safety factor
d_model = 0.00095       # 0.95 mm in metres (typical diameter in this range: b51-f#46)

# Uses frequency range similar to the data
f_model = np.linspace(f.min(), f.max(), 200)
L_model = []

for f_i in f_model:
    L_star = solve_L(f_i, SF_model, d_model)
    L_model.append(L_star)

L_model = np.array(L_model)

# Plotting data, quadratic regression, and analytical model (Figure 9)
plt.figure(figsize=(6, 4))
plt.plot(L, f, 'ko', label='Scale data')                    # scatter points
plt.plot(L_fit, f_fit_quad, 'r--', label='Quadratic fit')   # quadratic fit to data
plt.plot(L_model, f_model, 'g-', label='Analytical model')  # model curve

plt.xlabel('Speaking length L (m)')
plt.ylabel('Frequency f (Hz)')
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.show()

print("Linear fit (f vs 1/L): a1, a0 =", coef_lin)
print("Quadratic fit f(L): a, b, c =", coef_quad)
