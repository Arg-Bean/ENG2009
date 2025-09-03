import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# Define the system of differential equations
def model(y, t):
    x1, x1_prime, x2, x2_prime = y
    dx1dt = x1_prime
    dx1dt_prime = -x1_prime - 2 - 4 * x1 + x2_prime
    dx2dt = x2_prime
    dx2dt_prime = -x2_prime + f(t)
    return [dx1dt, dx1dt_prime, dx2dt, dx2dt_prime]

# Define the impulse force function f(t)
def f(t):
    # Unit impulse of duration 0.1 seconds starting from t = 0
    return 10.0 if 0 <= t < 0.1 else 0.0

# Euler method for numerical integration
def euler_method(model, y0, t):
    dt = t[1] - t[0]
    y = np.zeros((len(t), len(y0)))
    y[0] = y0
    for i in range(1, len(t)):
        y[i] = y[i - 1] + np.array(model(y[i - 1], t[i - 1])) * dt
    return y

# Initial conditions
y0 = [0, 0, 0, 0]

# Time points for simulation
t = np.linspace(0, 10, 1000)

# Perform Euler simulation
euler_solution = euler_method(model, y0, t)

# Extract the solution for each variable
x1_euler, x1_prime_euler, x2_euler, x2_prime_euler = euler_solution.T

# Use odeint for comparison
odeint_solution = odeint(model, y0, t)
x1_odeint, x1_prime_odeint, x2_odeint, x2_prime_odeint = odeint_solution.T

# Plot the results
plt.figure(figsize=(12, 8))

# Euler method plot
plt.subplot(2, 1, 1)
plt.plot(t, x1_euler, label='Euler: x1(t)')
plt.plot(t, x2_euler, label='Euler: x2(t)')
plt.title('Euler Method Simulation')
plt.legend()
plt.grid(True)

# odeint plot
#plt.subplot(2, 1, 2)
#plt.plot(t, x1_odeint, label='odeint: x1(t)')
#plt.plot(t, x2_odeint, label='odeint: x2(t)')
#plt.title('odeint Simulation')
#plt.legend()
#plt.grid(True)

plt.tight_layout()
plt.show()
