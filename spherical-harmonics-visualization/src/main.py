import sys
import os

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import plotly.graph_objects as go
import numpy as np
import math
from dash import Dash, html, dcc, callback, Output, Input

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

def generate_sh_values(alpha, beta, gamma):
    a_0 = np.array([0, 0, 0, 0, math.sqrt(7 / 12), 0, 0, 0, math.sqrt(5 / 12)]) # reference frame

    # rotation in euler angles
    R = rotate_sh(alpha, beta, gamma)
    a = R @ a_0

    l = (len(a_0) - 1) // 2

    # Compute results (vectorized instead of loop)
    sh_values = np.array([
        sum(a[j + l] * sph_harm_real(l, j, phi_flat[i], theta_flat[i])
            for j in range(-l, l + 1))
        for i in range(coords.shape[0])
    ])

    # Reshape back to grid
    sh_values = sh_values.reshape((n, n))

    # Compute cartesian coordinates
    offset = 0.5
    r = offset + sh_values

    x = r * np.sin(theta_grid) * np.cos(phi_grid)
    y = r * np.sin(theta_grid) * np.sin(phi_grid)
    z = r * np.cos(theta_grid)

    return sh_values, x, y, z

# Plot
app = Dash()

app.layout = html.Div([
    html.H1("Octahedral Frame in SH representation"),

    html.Div([
        html.Div([
            html.Label("Alpha", style={"fontWeight": "bold"}),
            dcc.Slider(id="alpha", min=0, max=0.5*math.pi, step=0.01, value=0,
                       marks={0: '0', math.pi/4: 'π/4', math.pi/2: 'π/2'}),
        ], style={"flex": 1}),

        html.Div([
            html.Label("Beta", style={"fontWeight": "bold"}),
            dcc.Slider(id="beta", min=0, max=0.5*math.pi, step=0.01, value=0,
                       marks={0: '0', math.pi/4: 'π/4', math.pi/2: 'π/2'}),
        ], style={"flex": 1}),

        html.Div([
            html.Label("Gamma", style={"fontWeight": "bold"}),
            dcc.Slider(id="gamma", min=0, max=0.5*math.pi, step=0.01, value=0,
                       marks={0: '0', math.pi/4: 'π/4', math.pi/2: 'π/2'}),
        ], style={"flex": 1}),

    ], style={"display": "flex", "gap": "40px", "padding": "16px 24px"}),

    dcc.Graph(id="3d-surface-plot", style={"height": "90vh"})
])

@callback(
    Output("3d-surface-plot", "figure"),
    Input("alpha", "value"),
    Input("beta",  "value"),
    Input("gamma", "value")
)
def update_graph(alpha, beta, gamma):
    sh_values, x, y, z = generate_sh_values(alpha, beta, gamma)

    figure = go.Figure(
        data=[go.Surface(x=x, y=y, z=z, surfacecolor=sh_values)]
    )

    figure.update_layout(
        uirevision="fixed",
        scene_camera=dict(
            up=dict(x=0, y=0, z=1),
            center=dict(x=0, y=0, z=0),
            eye=dict(x=1.25, y=1.25, z=1.25)
        )
    )

    return figure

if __name__ == "__main__":
    app.run(debug=True)