import numpy as np
import pyvista as pv

from odeco import canonical_odeco, generate_sh_values_from_coordinates
from rotations import rotate_z, rotate_x

n = 5

bases = {
    "octa": rotate_x(np.pi / 2) @ np.array(
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.6, 0, 0, 0, 0]
    ),
    # "odeco": canonical_odeco(1, 1, 1),
    "odeco spinner": rotate_x(np.pi / 2) @ np.array(
        [1, 0, 0, 0.2, 0, 0, 0, 0, 0, 0, 0.6, 0, 0, 0, 0]
    ),
    "odeco only x": canonical_odeco(1, 0, 0),
}


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


def add_odeco_surface(plotter, coords, opacity=1.0):
    sh_values, x, y, z = generate_sh_values_from_coordinates(coords)
    grid = pv.StructuredGrid(x, y, z)
    # StructuredGrid flattens points in Fortran order; match that for the scalars.
    grid["sh_values"] = np.asarray(sh_values).ravel(order="F")
    plotter.add_mesh(
        grid,
        scalars="sh_values",
        opacity=opacity,
        cmap="viridis",
        show_scalar_bar=False,
        smooth_shading=True,
    )


plotter = pv.Plotter(shape=(len(bases), 2), window_size=(1200, 600 * len(bases)))

for row, (name, base) in enumerate(bases.items()):
    odeco_coords, average = compute_odeco_set(base, n)

    plotter.subplot(row, 0)
    plotter.add_text(f"{name} - individual odecos", font_size=10)
    for coords in odeco_coords:
        add_odeco_surface(plotter, coords, opacity=0.3)

    plotter.subplot(row, 1)
    plotter.add_text(f"{name} - average", font_size=10)
    add_odeco_surface(plotter, average, opacity=1.0)

plotter.link_views()
plotter.show(screenshot="./average.png")