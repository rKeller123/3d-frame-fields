import subprocess
import numpy as np
import os
import math

import plotly.graph_objects as go
from dash import Dash, html, dcc, callback, Output, Input, exceptions

from odeco import generate_sh_values, generate_sh_values_from_coordinates, generate_coordinates

executable_path = os.path.abspath("./3d-odeco-frame-fields/out/build/x64-debug/cli.exe")

def run_cli(values: np.ndarray):
    values = np.asarray(values)
    values = values.flatten()

    if values.size != 15:
        raise ValueError(f"Expected 15 values, got {values.size}")

    arg = ";".join(map(str, values))

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

    html.H1("Odeco Frame in SH representation"),

    dcc.Store(id="coordinates-store"),

    html.Div([
        html.Div([

            html.H3("Controls"),

            html.Div([
                html.Label("Stretch x", style={"fontWeight": "bold"}),
                dcc.Slider(id="l_x", min=0, max=10, step=0.01, value=1),
                dcc.Input(id="l_x_input", type="number", value=1, step=0.01),
            ], style={"marginBottom": "15px"}),

            html.Div([
                html.Label("Stretch y", style={"fontWeight": "bold"}),
                dcc.Slider(id="l_y", min=0, max=10, step=0.01, value=1),
                dcc.Input(id="l_y_input", type="number", value=1, step=0.01),
            ], style={"marginBottom": "15px"}),

            html.Div([
                html.Label("Stretch z", style={"fontWeight": "bold"}),
                dcc.Slider(id="l_z", min=0, max=10, step=0.01, value=1),
                dcc.Input(id="l_z_input", type="number", value=1, step=0.01),
            ], style={"marginBottom": "25px"}),

            html.Hr(),

            html.Div([
                html.Label("Rotate x"),
                dcc.Slider(id="alpha", min=-0.5 * math.pi, max=0.5 * math.pi,
                           step=0.01, value=0),
                dcc.Input(id="alpha_input", type="number", value=0, step=0.01),
            ], style={"marginBottom": "15px"}),

            html.Div([
                html.Label("Rotate y"),
                dcc.Slider(id="beta", min=-0.5 * math.pi, max=0.5 * math.pi,
                           step=0.01, value=0),
                dcc.Input(id="beta_input", type="number", value=0, step=0.01),
            ], style={"marginBottom": "15px"}),

            html.Div([
                html.Label("Rotate z"),
                dcc.Slider(id="gamma", min=-0.5 * math.pi, max=0.5 * math.pi,
                           step=0.01, value=0),
                dcc.Input(id="gamma_input", type="number", value=0, step=0.01),
            ]),

        ], style={
            "width": "20%",
            "padding": "20px",
            "borderRight": "1px solid #ddd",
            "overflowY": "auto",
            "height": "85vh"
        }),
        html.Div([
            html.H3("Input"),
            dcc.Graph(id="input-plot", style={"height": "80vh"})
        ], style={
            "width": "40%",
            "padding": "10px"
        }),
        html.Div([
            html.H3("Output"),
            dcc.Graph(id="output-plot", style={"height": "80vh"})
        ], style={
            "width": "40%",
            "padding": "10px"
        }),

    ], style={
        "display": "flex",
        "width": "100%"
    })

])


for slider_id in ["l_x", "l_y", "l_z", "alpha", "beta", "gamma"]:
    @app.callback(
        Output(f"{slider_id}_input", "value"),
        Input(slider_id, "value"),
        prevent_initial_call=True,
    )
    def sync_slider_to_box(val):
        return val

@callback(
    Output("input-plot", "figure"),
    Output("coordinates-store", "data"),
    Input("l_x_input", "value"),
    Input("l_y_input", "value"),
    Input("l_z_input", "value"),
    Input("alpha_input", "value"),
    Input("beta_input", "value"),
    Input("gamma_input", "value"),
)
def update_graph(l_x, l_y, l_z, alpha, beta, gamma):
    sh_values, x, y, z = generate_sh_values(l_x, l_y, l_z, alpha, beta, gamma)
    coordinates = generate_coordinates(l_x, l_y, l_z, alpha, beta, gamma)

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
    return figure, coordinates.tolist()

@callback(
    Output("output-plot", "figure"),
    Input("coordinates-store", "data"),
    config_prevent_initial_callbacks=True
)
def update_output_graph(coordinates):
    if coordinates is None:
        raise exceptions.PreventUpdate
    num_iter, projection = run_cli(np.array(coordinates))
    print(num_iter)
    sh_values, x, y, z = generate_sh_values_from_coordinates(projection)

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

app.run(debug=True)