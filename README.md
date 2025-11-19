# CMM-Group-20-Project-Final
Modern Grand Piano hammer string mechanism


## General model


No input requied from user, simply run.

This code produces figure 3 and 4 modelling the behaviour of a piano string and visualises how string length and tension (with inharmonictiy) produce a target frequency of 440Hz.

numpy and matplotlib.pyplot is imported

After defining the constants and frequency equation with inharmonicity 
(1/ (2 * L)) * np.sqrt(T / mu) * np.sqrt(1 + B), 
A 3D graph is generated via the use of arrays.

Along with this, our reference plane of 440Hz intercepts the function.
From this, a 2D graph is created with the derived equation
T(L) = 4*mu*f^2*L^2 - (pi^2 * E * I)/(L^2),
as seen by figure 4.


## Root finding+ode final


With the use of root finding(Brent's method), this file initially finds the optimum length and tension with input variables of safety factor and diameter.
[Brent's method lines 30-48]

It forms oscillating graph by solving the ODE using RK4 method.
[RK4 lines 53-101]

User inputs frequency, safe factor and diameter which outputs an optimum tension and speaking length.
Then asks user to input a sustain time which then outputs an oscillating graph.
To terminate this loop, user is required to enter 'q' into sustain time to exit damping exploration. To exit entire program, 'q' again must be entered into target frequency.

numpy, matpotlib.pyplot and scipy.optimize(brentq) are imported.

After defining constants and non linear equation, 
4 rho f^2 L^4 - sigma_w L^2 - (pi^2 E d^2)/16 = 0, 
an initial bracket is created for length and ensuring there is a sign change within it to enable root finding.
it then runs Brent's method: firstly attempting inverse quadratic interpolation, if this fails it attempts secant methods. If both faster steps are unreliable, it defaults to a bisection method.

After determining the optimum length and tension, RK4 ODE solver is used with length, diameter, frequency and T60 to simulate decay.
RK4 first uses K1 to find the slope at the beginning of the interval. K2 and K3 are midpoint methods of the slope using its previous points. K4 is the end points of the interval and act as the predictor corrector.

To create a graph exploring damping, arrays of safety values and length are created to produce a graph (Figure 5).
The user is then asked for inputs to create a oscillating graph(Figure 6 and 7).


## Linear regression and intepolation


This file explores real word data(figure 10) from Pohlmann's scale design notes ranging from a61 to a37 to our model and produces figure 8 and 9.

No input requied from user, simply run.

numpy, matpotlib.pyplot and scipy.optimize(brentq) are imported.

It defines material characteristics along with the equations for length, working stress, target frequency and diameter as seen before.
After an initial bracket is created, Brent's method is again used for finding optimum length(L*).
[Brent's method line 26-41]
Data from figure 10 is used in creating arrays for frequency and length.
Linear regression(least-squares method) is used to plot this data: f vs 1/L (figure 8).
[Linear regression line 70-80]

To compare our model with data obtained, least-squares quadratic regression is used(Figure 9).
[Quadratic regression 91-112]







