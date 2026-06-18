import subprocess
import numpy as np
import os
import platform
import math
from time import time

PLATFORM = platform.system()
executable_path = os.path.abspath("./3d-odeco-frame-fields/build/cli")
if (PLATFORM == "Windows"):
    executable_path = os.path.abspath("./3d-odeco-frame-fields/out/build/x64-debug/cli.exe")


import plotly.graph_objects as go
from dash import Dash, html, dcc, callback, Output, Input, exceptions, ALL

from odeco import generate_sh_values_from_coordinates, generate_coordinates

def run_cli(mode: int, values: np.ndarray, d: np.ndarray):
    values = np.asarray(values).flatten()

    if values.size != 15:
        raise ValueError(f"Expected 15 values, got {values.size}")

    parts = [str(mode)]

    if mode == 1:
        parts.extend(map(str, d))

    parts.extend(map(str, values))

    arg = ";".join(parts)

    result = subprocess.run(
        [executable_path, arg],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        raise RuntimeError(result.stderr)

    out = result.stdout.strip().split(";")
    num_iter = int(out[0])
    values = np.array([float(x) for x in out[1:]])
    return num_iter, values

app = Dash()

app.layout = html.Div([
    html.H1("Odeco Frame Projection"),

    dcc.Store(id="coordinates-store"),

    html.Div([
        html.Div([
            html.H3("Controls"),

            html.Div([
                dcc.Checklist(id="z-aligned",
                              options=[
                                  {"label": "aligned projection", "value": "z_proj"}
                              ],
                              value=[]),
                dcc.Input(id="d_x_input", type="number", value=0),
                dcc.Input(id="d_y_input", type="number", value=0),
                dcc.Input(id="d_z_input", type="number", value=1),
                html.Br(), html.Br(),
            ]),

            html.Div([
                html.Label("Scale x"),
                dcc.Slider(id="l_x", min=0, max=10, step=0.01, value=1),
                dcc.Input(id="l_x_input", type="number", value=1),
            ]),

            html.Div([
                html.Label("Scale y"),
                dcc.Slider(id="l_y", min=0, max=10, step=0.01, value=1),
                dcc.Input(id="l_y_input", type="number", value=1),
            ]),

            html.Div([
                html.Label("Scale z"),
                dcc.Slider(id="l_z", min=0, max=10, step=0.01, value=1),
                dcc.Input(id="l_z_input", type="number", value=1),
            ]),

            html.Hr(),

            html.Div([
                html.Label("Rotate x"),
                dcc.Slider(id="alpha", min=-0.5 * math.pi, max=0.5 * math.pi, step=0.01, value=0),
                dcc.Input(id="alpha_input", type="number", value=0),
            ]),

            html.Div([
                html.Label("Rotate y"),
                dcc.Slider(id="beta", min=-0.5 * math.pi, max=0.5 * math.pi, step=0.01, value=0),
                dcc.Input(id="beta_input", type="number", value=0),
            ]),

            html.Div([
                html.Label("Rotate z"),
                dcc.Slider(id="gamma", min=-0.5 * math.pi, max=0.5 * math.pi, step=0.01, value=0),
                dcc.Input(id="gamma_input", type="number", value=0),
            ]),

            html.Hr(),

            html.H3("Override Coordinates"),

            html.Div([
                dcc.Input(
                    id={"type": "override-coord", "index": i},
                    type="number",
                    step=0.01,
                    value=None,
                    style={"width": "100%"}
                )
                for i in range(15)
            ], style={
                "display": "flex",
                "flexDirection": "column",
                "gap": "6px"
            })
        ], style={"width": "20%", "padding": "20px", "borderRight": "1px solid #ddd", "display": "flex", "flexDirection": "column"}),

        html.Div([
            html.H3("Target"),
            dcc.Graph(id="input-plot", style={"height": "80vh"})
        ], style={"width": "40%"}),

        html.Div([
            html.H3("Projection"),
            html.P(id="info"),
            dcc.Graph(id="output-plot", style={"height": "80vh"})
        ], style={"width": "40%"}),

    ], style={"display": "flex"})
])


for slider_id in ["l_x", "l_y", "l_z", "alpha", "beta", "gamma"]:
    @app.callback(
        Output(f"{slider_id}_input", "value"),
        Input(slider_id, "value"),
        prevent_initial_call=True
    )
    def sync_slider(v):
        return v

@callback(
    Output("coordinates-store", "data"),
    Input("l_x_input", "value"),
    Input("l_y_input", "value"),
    Input("l_z_input", "value"),
    Input("alpha_input", "value"),
    Input("beta_input", "value"),
    Input("gamma_input", "value"),
)
def compute_coordinates(l_x, l_y, l_z, alpha, beta, gamma):
    coordinates = generate_coordinates(l_x, l_y, l_z, alpha, beta, gamma)
    return coordinates.tolist()

@callback(
    Output({"type": "override-coord", "index": ALL}, "value"),
    Input("coordinates-store", "data"),
)
def sync_override_boxes(coords):
    if coords is None:
        raise exceptions.PreventUpdate
    return coords

@callback(
    Output("input-plot", "figure"),
    Input("coordinates-store", "data"),
    Input({"type": "override-coord", "index": ALL}, "value"),
)
def update_input_plot(coords, overrides):

    if coords is None:
        raise exceptions.PreventUpdate

    final = np.array(coords)

    if overrides and any(v is not None for v in overrides):
        for i, v in enumerate(overrides):
            if v is not None:
                final[i] = v

    sh_values, x, y, z = generate_sh_values_from_coordinates(final)

    fig = go.Figure(data=[go.Surface(x=x, y=y, z=z, surfacecolor=sh_values)])
    return fig

@callback(
    Output("output-plot", "figure"),
    Output("info", "children"),
    Input("input-plot", "figure"),
    Input("coordinates-store", "data"),
    Input({"type": "override-coord", "index": ALL}, "value"),
    Input("z-aligned", "value"),
    Input("d_x_input", "value"),
    Input("d_y_input", "value"),
    Input("d_z_input", "value"),
)
def update_output(_, coords, overrides, z_aligned, d_x, d_y, d_z):

    if coords is None or z_aligned is None:
        raise exceptions.PreventUpdate

    final = np.array(coords)
    mode = int('z_proj' in z_aligned)
    d = np.array([d_x, d_y, d_z])

    if overrides and any(v is not None for v in overrides):
        for i, v in enumerate(overrides):
            if v is not None:
                final[i] = v

    start = time()
    num_iter, projection = run_cli(mode, final, d)
    end = time()
    sh_values, x, y, z = generate_sh_values_from_coordinates(projection)

    fig = go.Figure(data=[go.Surface(x=x, y=y, z=z, surfacecolor=sh_values)])
    return fig, f"{num_iter} iterations | {end - start:.2f} s"


app.run(debug=True)