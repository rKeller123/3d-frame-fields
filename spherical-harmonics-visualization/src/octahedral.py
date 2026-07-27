import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import plotly.graph_objects as go
import numpy as np
import math
from dash import Dash, html, dcc, callback, Output, Input

from spherical_harmonics import sph_harm_real
from rotations import rotate_octahedral

def canonical_octahedral():
    return np.array([0, 0, 0, 0, math.sqrt(7 / 12), 0, 0, 0, math.sqrt(5 / 12)])

def generate_sh_values_from_coordinates_octa(coordinates):
    n = 100

    # Create grid
    theta = np.linspace(0, math.pi, num=n)
    phi = np.linspace(0, 2 * math.pi, num=n)

    theta_grid, phi_grid = np.meshgrid(theta, phi)

    # Flatten ONLY for computation
    theta_flat = theta_grid.flatten()
    phi_flat = phi_grid.flatten()

    l = 4

    # Compute results (vectorized instead of loop)
    sh_values = np.array([
        sum(coordinates[j + l] * sph_harm_real(l, j, theta_flat[i], phi_flat[i])
            for j in range(-l, l + 1))
        for i in range(theta_flat.size)
    ])

    # Reshape back to grid
    sh_values = sh_values.reshape((n, n))

    # Compute cartesian coordinates
    offset = 0.8
    r = offset + sh_values

    x = r * np.sin(theta_grid) * np.cos(phi_grid)
    y = r * np.sin(theta_grid) * np.sin(phi_grid)
    z = r * np.cos(theta_grid)

    return sh_values, x, y, z

def generate_octa_coordinates(alpha, beta, gamma):
    a_0 = canonical_octahedral()
    a = rotate_octahedral(a_0, alpha, beta, gamma)
    return a