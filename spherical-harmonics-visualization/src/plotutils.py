
import plotly.graph_objects as go
from dash import Dash, html, dcc, callback, Output, Input, exceptions, ALL

from odeco import generate_sh_values_from_coordinates

def plot_odeco(coords_list, opacities=None, colorscale="Viridis", showscale=False):
    if not isinstance(coords_list, list):
        coords_list = [coords_list]
    if opacities is None:
        opacities = [1.0] * len(coords_list)

    fig = go.Figure()
    for coords, opacity in zip(coords_list, opacities):
        sh_values, x, y, z = generate_sh_values_from_coordinates(coords)
        fig.add_trace(
            go.Surface(
                x=x,
                y=y,
                z=z,
                surfacecolor=sh_values,
                opacity=opacity,
                colorscale=colorscale,
                showscale=showscale,
            )
        )
    return fig