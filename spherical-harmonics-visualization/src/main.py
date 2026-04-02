import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import plotly.graph_objects as go
import numpy as np
import math

from spherical_harmonics import sph_harm_real
from wigner_d import rotate_sh


n = 100

# Create grid
theta = np.linspace(0, math.pi, num=n)
phi = np.linspace(0, 2 * math.pi, num=n)

theta_grid, phi_grid = np.meshgrid(theta, phi)

# Flatten ONLY for computation
theta_flat = theta_grid.flatten()
phi_flat = phi_grid.flatten()

coords = np.stack((theta_flat, phi_flat), axis=-1)

a_0 = np.array([0, 0, 0, 0, math.sqrt(7 / 12), 0, 0, 0, math.sqrt(5 / 12)]) # reference frame

# rotation in euler angles
alpha = math.pi / 6
beta = math.pi / 2
gamma = math.pi / 4

R = rotate_sh(alpha, beta, gamma)

print(R)

a = R @ a_0

print(a)

l = (len(a_0) - 1) // 2

# Compute results (vectorized instead of loop)
results = np.array([
    sum(a[j + l] * sph_harm_real(l, j, phi_flat[i], theta_flat[i])
        for j in range(-l, l + 1))
    for i in range(coords.shape[0])
])

# Reshape back to grid
results = results.reshape((n, n))

# Compute cartesian coordinates
offset = 0.5
r = offset + results

x = r * np.sin(theta_grid) * np.cos(phi_grid)
y = r * np.sin(theta_grid) * np.sin(phi_grid)
z = r * np.cos(theta_grid)

# Plot
fig = go.Figure(data=[go.Surface(x=x, y=y, z=z, surfacecolor=results)])

camera = dict(
    up=dict(x=0, y=0, z=1),
    center=dict(x=0, y=0, z=0),
    eye=dict(x=1.25, y=1.25, z=1.25)
)

fig.update_layout(
    title="Spherical harmonics plot",
    scene_camera=camera
)

fig.show()