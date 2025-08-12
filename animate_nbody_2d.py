import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load trajectory data
data = np.loadtxt("trajectory_testing.csv", delimiter=",")
steps, cols = data.shape
N = cols // 3
positions = np.reshape(data, (steps, N, 3))

# Center the system on the Sun (body 0)
positions -= positions[:, 0:1, :]

# Plot setup
fig, ax = plt.subplots()
#colors = ['orange'] + ['blue'] * (N - 1) # to be commented out when not simulating the solar system
#sizes = [80] + [20] * (N - 1)		 # to be commented out when not simulating the solar system
colors = ['orange', 'gray', 'yellow', 'blue', 'red', 'brown', 'gold', 'cyan', 'navy', 'darkred']
sizes = [100, 5, 6, 8, 6, 40, 35, 25, 24, 3]

# Initial positions
x0 = positions[0, :, 0]
y0 = positions[0, :, 1]

# Scatter for bodies
scat = ax.scatter(x0, y0, s=sizes, c=colors)

# Lines for trails (one per body)
trails = [ax.plot([], [], lw=1, color=colors[i])[0] for i in range(N)]

# Axis config
ax.set_xlim(-12, 12)
ax.set_ylim(-12, 12)
ax.set_aspect('equal')
ax.set_title("2D N-Body Simulation with Trails")

# Trail length (number of past frames to show)
trail_length = 100

def update(frame):
    # Update body positions
    x = positions[frame, :, 0]
    y = positions[frame, :, 1]
    scat.set_offsets(np.c_[x, y])

    # Update trails
    for i in range(N):
        trails[i].set_data(positions[:frame+1, i, 0], positions[:frame+1, i, 1])


    return [scat] + trails

ani = FuncAnimation(fig, update, frames=steps, interval=30, blit=True)
ani.save("orbit.mp4", fps=30, dpi=200)
plt.show()
