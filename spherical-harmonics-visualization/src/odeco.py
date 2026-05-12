import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import plotly.graph_objects as go
import numpy as np
import math
from dash import Dash, html, dcc, callback, Output, Input

from spherical_harmonics import sph_harm_real
from rotations import rotate_odeco


def canonical_odeco(l_x, l_y, l_z):
    # assert l_x >= 0, "stretching ratio x has to be positive"
    # assert l_y >= 0, "stretching ratio y has to be positive"
    # assert l_z >= 0, "stretching ratio z has to be positive"

    # reference frame expressed in SH basis
    a_0 = [2/5 * math.sqrt(math.pi) * (l_x + l_y + l_z), # band 0
           0, 0, 4 / 7 * math.sqrt(math.pi / 5) * (2 * l_z - (l_x + l_y)), 0, 4 / 7 * math.sqrt(3 / 5 * math.pi) * (l_x - l_y), # band 2
           0, 0, 0, 0, 2 / 35 * math.sqrt(math.pi) * (l_x + l_y + 8 / 3 * l_z), 0, -4 / 21 * math.sqrt(math.pi / 5) * (l_x - l_y), 0, 2 / 3 * math.sqrt(math.pi / 35) * (l_x + l_y) # band 4
           ]

    assert len(a_0) == 15, "odeco tensor should have dimension 15"

    return a_0

def generate_sh_values_from_coordinates(coordinates):
    n = 100

    # Create grid
    theta = np.linspace(0, math.pi, num=n)
    phi = np.linspace(0, 2 * math.pi, num=n)

    theta_grid, phi_grid = np.meshgrid(theta, phi)

    theta_flat = theta_grid.flatten()
    phi_flat = phi_grid.flatten()

    bands = [0, 2, 4]
    sh_values = []

    for k in range(n * n):
        sh_value = 0
        i = 0
        for l in bands:
            for m in range(-l, l + 1):
                sh_value += coordinates[i] * sph_harm_real(l, m, theta_flat[k], phi_flat[k])
                i += 1
        sh_values.append(sh_value)

    sh_values = np.array(sh_values)

    sh_values = sh_values.reshape((n, n))
    x = sh_values * np.sin(theta_grid) * np.cos(phi_grid)
    y = sh_values * np.sin(theta_grid) * np.sin(phi_grid)
    z = sh_values * np.cos(theta_grid)
    return sh_values, x, y, z

def generate_coordinates(l_x, l_y, l_z, alpha, beta, gamma):
    a_0 = canonical_odeco(l_x, l_y, l_z)
    a = rotate_odeco(a_0, alpha, beta, gamma)
    return a

def generate_sh_values(l_x, l_y, l_z, alpha, beta, gamma):
    a = generate_coordinates(l_x, l_y, l_z, alpha, beta, gamma)
    sh_values, x, y, z = generate_sh_values_from_coordinates(a)
    return sh_values, x, y, z