import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Lorenz system parameters
sigma = 10.0
rho = 28.0
beta = 8.0 / 3.0

# Lorenz system differential equations
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

# Normalize the x-axis to use for color mapping
norm = plt.Normalize(x.min(), x.max())
colors = plt.cm.hsv(norm(x))

# Set up the figure and axis with black background
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_facecolor('black')

# Plot the Lorenz attractor with a color gradient dependent on the x position
# Here we plot only the XY projection by setting the Z values to 0
sc = ax.scatter(y, z, c=colors, s=0.5, alpha=0.6)

# Hide the axes
ax.set_xticks([])
ax.set_yticks([])
ax.axis('off')  # Turn off the axis completely

# Tight layout and save to PDF with a black background
plt.tight_layout()
plt.savefig("lorenz_attractor_yz_projection.png", format='png', bbox_inches='tight', facecolor='black', edgecolor='none')

# Close the plot
plt.close()
