import plotly.graph_objects as go
from dash import Dash, html, dcc, callback, Output, Input, exceptions, ALL
import numpy as np
from odeco import generate_coordinates, canonical_odeco
from plotutils import plot_odeco
from rotations import rotate_z, rotate_x

# Define the bases you want to compare, by name
bases = {
    "octa": rotate_x(np.pi / 2) @ np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0, 0, 0, 0]
    ),
    "odeco": canonical_odeco(1, 1, 1),
    "odeco spinner": rotate_x(np.pi / 2) @ np.array(
        [1, 0, 0, 0.2, 0, 0, 0, 0, 0, 0, 0.6, 0, 0, 0, 0]
    ),
    "odeco only x": canonical_odeco(1, 0, 0)
}

app = Dash()

app.layout = html.Div(
    [
        html.Div(
            [
                html.Label("n: "),
                dcc.Input(id="n-input", type="number", value=3, min=1, step=1),
            ],
            style={"padding": "10px"},
        ),
        html.Div(
            [
                html.Div(
                    [
                        html.H3(name),
                        html.Div(
                            [
                                html.Div(
                                    [
                                        html.H4("Individual odecos"),
                                        dcc.Graph(
                                            id={"type": "odecos-plot", "index": name},
                                            style={"height": "65vh"},
                                        ),
                                    ],
                                    style={"flex": "1 1 0", "minWidth": "350px"},
                                ),
                                html.Div(
                                    [
                                        html.H4("Average"),
                                        dcc.Graph(
                                            id={"type": "average-plot", "index": name},
                                            style={"height": "65vh"},
                                        ),
                                    ],
                                    style={"flex": "1 1 0", "minWidth": "350px"},
                                ),
                            ],
                            style={"display": "flex", "gap": "10px"},
                        ),
                    ],
                    style={"marginBottom": "30px"},
                )
                for name in bases
            ]
        ),
    ]
)


def compute_odeco_set(base, n):
    angles = np.linspace(0, np.pi, n, endpoint=False)
    average = np.zeros(shape=(15,), dtype=np.float64)
    odeco_coords = []
    for a in angles:
        R = rotate_z(a)
        coords = R @ base
        odeco_coords.append(coords)
        average += (1 / n) * coords
    return odeco_coords, average


@callback(
    Output({"type": "odecos-plot", "index": ALL}, "figure"),
    Output({"type": "average-plot", "index": ALL}, "figure"),
    Input("n-input", "value"),
)
def update_plots(n):
    if not n or n < 1:
        raise exceptions.PreventUpdate

    odecos_figures = []
    average_figures = []

    for name, base in bases.items():
        odeco_coords, average = compute_odeco_set(base, n)

        odecos_figures.append(
            plot_odeco(odeco_coords, opacities=[0.3] * len(odeco_coords))
        )
        average_figures.append(plot_odeco([average], opacities=[1.0]))

    return odecos_figures, average_figures


app.run(debug=True)