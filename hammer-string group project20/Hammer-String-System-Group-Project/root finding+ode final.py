
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq


# Material constants

E = 2.0e11        # Pa
rho = 7850.0      # kg/m^3
sigma_y = 2.1e9   # Pa, yield stress of music wire


# Nonlinear equation g(L) = 0

# 4 rho f^2 L^4 - sigma_w L^2 - (pi^2 E d^2)/16 = 0

def g_L(L, sigma_w, f_target, d):
    return 4.0 * rho * f_target**2 * L**4 - sigma_w * L**2 - (np.pi**2 * E * d**2) / 16.0

def solve_L_T(f_target, SF, d):
    """
    Given target frequency f_target (Hz), safety factor SF, and diameter d (m),
    return (L*, T*) for the string using Brent's method.
    """
    # cross-sectional area & working stress
    A = 0.25 * np.pi * d**2
    sigma_w = sigma_y / SF

    # initial bracket for L [m]
    L_min, L_max = 0.1, 2.0

    # ensure there is a sign change in the bracket
    g_min = g_L(L_min, sigma_w, f_target, d)
    g_max = g_L(L_max, sigma_w, f_target, d)

    while g_min * g_max > 0:
        L_max += 0.5
        if L_max > 5.0:
            raise RuntimeError("Could not bracket a root for L - check inputs.")
        g_max = g_L(L_max, sigma_w, f_target, d)

    # Brent's root finding
    L_star = brentq(g_L, L_min, L_max, args=(sigma_w, f_target, d))

    # Corresponding tension from T = sigma_w * A
    T_star = sigma_w * A
    return L_star, T_star


# RK4 ODE solver: Using L, d, f, T60 to simulate decay

def rk4_decay(f, L, d, T60):
    """
    Free decay of the string for given design (f, L, d) and sustain time T60.
    Uses the damped oscillator x¨ + 2ζω0 x˙ + ω0^2 x = 0 with RK4.
    """

    # mass m = ρ A L
    A = 0.25 * np.pi * d**2
    m = rho * A * L
    
    #natural frequency
    omega0 = 2.0 * np.pi * f

    # damping ratio ζ from chosen T60:  x(T60)/x(0) = 1/1000
    zeta = np.log(1000.0) / (omega0 * T60)

    # Degree one： y = [x, v]
    def rhs(t, y):
        x, v = y
        return np.array([v, -2.0*zeta*omega0*v - (omega0**2)*x])

    # time
    t0   = 0.0
    t_end = T60 + 1             # simulate to T60+1
    dt   = 1.0 / (200.0 * f)      # ~200 steps per period
    N    = int((t_end - t0)/dt) + 1

    t = np.linspace(t0, t_end, N)
    y = np.zeros((N, 2))

    # initial condition：x(0)=0, v(0)=1/m（Unit impulsive force）
    x0 = 0.0
    v0 = 1 / m
    y[0, :] = [x0, v0]

    # RK4 main loop
    for k in range(N-1):
        tk = t[k]
        yk = y[k]

        k1 = rhs(tk,          yk)
        k2 = rhs(tk+dt/2.0,   yk + dt*k1/2.0)
        k3 = rhs(tk+dt/2.0,   yk + dt*k2/2.0)
        k4 = rhs(tk+dt,       yk + dt*k3)

        y[k+1] = yk + dt * (k1 + 2*k2 + 2*k3 + k4) / 6.0

    x = y[:, 0]
    return t, x


# Main script: design helper + damping exploration
if __name__ == "__main__":

    # Graph for Safety Factor against L*
    f_A4 = 440.0
    d_A4 = 0.001

    SF_vals = np.linspace(2.0, 3.5, 20)
    L_vals, T_vals = [], []

    # Calculating values to create an array
    for SF in SF_vals:
        L_star, T_star = solve_L_T(f_A4, SF, d_A4)
        L_vals.append(L_star)
        T_vals.append(T_star)

    L_vals = np.array(L_vals)
    
    # Plotting graph (figure 5)
    plt.figure()
    plt.plot(SF_vals, L_vals, "o-")
    plt.xlabel("Safety factor SF")
    plt.ylabel("Speaking length L*(SF) [m]")
    plt.title("Speaking length vs safety factor (f = 440 Hz, d = 1 mm)")
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # User input (f, SF, d, with multiple T60)
    print("Piano string design helper (enter 'q' to quit)\n")

    while True:
        #root finding first to find L* & d_user
        s = input("Target frequency f (Hz): ")
        if s.lower().startswith("q"):
            break
        try:
            f_user = float(s)
        except ValueError:
            print(" Invalid frequency, try again.\n")
            continue

        s = input("Safety factor SF (e.g. 2.5): ")
        if s.lower().startswith("q"):
            break
        try:
            SF_user = float(s)
        except ValueError:
            print(" Invalid SF, try again.\n")
            continue

        s = input("String diameter d (mm, e.g. 1.0): ")
        if s.lower().startswith("q"):
            break
        try:
            d_mm   = float(s)
            d_user = d_mm / 1000.0   # mm → m
        except ValueError:
            print(" Invalid diameter, try again.\n")
            continue

        try:
            L_star, T_star = solve_L_T(f_user, SF_user, d_user)
        except Exception as e:
            print(" Error solving for L,T:", e, "\n")
            continue

        A_user  = 0.25 * np.pi * d_user**2
        sigma_w = sigma_y / SF_user

        print(f"\nFor f = {f_user:.1f} Hz, SF = {SF_user:.2f}, d = {d_mm:.3f} mm:")
        print(f" Working stress σ_w = {sigma_w/1e6:.1f} MPa")
        print(f" Speaking length L* = {L_star:.3f} m")
        print(f" Tension T*        = {T_star:.1f} N")
        print(f" Cross-section A   = {A_user*1e6:.3f} mm^2\n")

       # Different T60 values
        print("Now explore damping for this design (press ENTER to go back to a new design).")
        while True:
            s = input("  Sustain time T60 (s), or just ENTER for new design: ")
            if s == "":        # Press enter to go back to last program
                print("")
                break
            if s.lower().startswith("q"):   # Input q causes program to terminate
                running = False
                print("\nExiting program.\n")
                break

            try:
                T60_user = float(s)
            except ValueError:
                print("   Invalid T60, try again.\n")
                continue

           # Using the current design from user inputs (f_user, L_star, d_user) with this T60, graph plotted
            t, x = rk4_decay(f_user, L_star, d_user, T60_user)

            plt.figure()
            plt.plot(t, x)
            plt.xlabel("t (s)")
            plt.ylabel("x (m)")
            plt.title(
                f"Free decay: f = {f_user:.1f} Hz, "
                f"L* = {L_star:.3f} m, d = {d_mm:.2f} mm, T60 = {T60_user:.1f} s"
            )
            plt.grid(True)
            plt.tight_layout()
            plt.show()


