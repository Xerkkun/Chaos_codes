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

# Initial state and time array
initial_state = [1.0, 1.0, 1.0]
t_span = (0, 40)
t_points = np.linspace(*t_span, 10000)

# Solve the Lorenz system
solution = solve_ivp(lorenz, t_span, initial_state, t_eval=t_points)

# Extract the solution for plotting
x = solution.y[0]
y = solution.y[1]
z = solution.y[2]

# Normalize the time points to use for color mapping
norm = plt.Normalize(t_points.min(), t_points.max())
colors = plt.cm.hsv(norm(t_points))

# Set up the figure and axis with black background
fig = plt.figure(figsize=(15, 5))
ax1 = fig.add_subplot(131, projection='3d', facecolor='black')
ax2 = fig.add_subplot(132, projection='3d', facecolor='black')
ax3 = fig.add_subplot(133, projection='3d', facecolor='black')

# Plot for each subplot
for ax in [ax1, ax2, ax3]:
    # Set the view angle to make the projections visible
    ax.grid(False)  # Turn off the grid
    ax.xaxis.pane.fill = False  # Turn off the pane filling
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('white')  # Set pane edge color to white
    ax.yaxis.pane.set_edgecolor('white')
    ax.zaxis.pane.set_edgecolor('white')
    ax.set_xticks([])  # Remove the ticks
    ax.set_yticks([])
    ax.set_zticks([])
    ax.scatter(x, y, z, c=colors, s=1, alpha=0.6)  # Use scatter plot

# Set the view for each plot
ax1.view_init(elev=20, azim=60)   # Adjust this to get the desired view
ax2.view_init(elev=20, azim=60)
ax3.view_init(elev=20, azim=60)

# Hide the axes
for ax in [ax1, ax2, ax3]:
    ax._axis3don = False

# Tight layout to maximize space
plt.tight_layout()
plt.show()
