import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Lorenz system parameters
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0

# Lorenz system differential equations.
def lorenz(t, state):
    x, y, z = state
    return [sigma * (y - x), x * (rho - z) - y, x * y - beta * z]

# Initial state (near the origin but not exactly at it).
initial_state = [0.1, 0.0, 0.0]

# Time points to solve for (from 0 to 40 in 10000 steps).
t_points = np.linspace(0, 40, 10000)

# Solve the Lorenz system.
solution = solve_ivp(lorenz, [0, 40], initial_state, t_eval=t_points)

# Extract the solution for plotting.
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

# Normalize the time points to use for color mapping.
norm = plt.Normalize(t_points.min(), t_points.max())
colors = plt.cm.jet(norm(t_points))

# Set up the figure and axis with black background.
fig = plt.figure(figsize=(12, 10))
fig.patch.set_facecolor('black')
ax = fig.add_subplot(111, projection='3d', facecolor='black')

# Plot the Lorenz attractor with the color gradient.
ax.plot(x, y, z, color='r', alpha=0.4, lw=0.6)

# Plot projections as lines on the planes.
ax.plot(x, y, min(z)-0.5*np.ones_like(z), color='w', alpha=0.3, lw=0.5)  # XY plane projection
ax.plot(x, max(y)+0.5*np.ones_like(y), z, color='w', alpha=0.3, lw=0.5)  # YZ plane projection
ax.plot(min(x)-0.5*np.ones_like(x), y, z, color='w', alpha=0.3, lw=0.5)  # XZ plane projection

# Set labels with white color.
ax.set_xlabel("X Axis", color='white')
ax.set_ylabel("Y Axis", color='white')
ax.set_zlabel("Z Axis", color='white')

# Set the tick colors to white.
ax.tick_params(axis='x', colors='white')
ax.tick_params(axis='y', colors='white')
ax.tick_params(axis='z', colors='white')

# Set the color of the axis spines to white.
for spine in ax.spines.values():
    spine.set_color('white')

# Show the plot.
plt.show()
