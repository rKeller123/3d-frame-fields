import re
import subprocess

import numpy as np
import pyvista as pv
from tqdm import tqdm

from odeco import generate_sh_values_from_coordinates

CLI = "../snail_field_cxx/build/examples/averaging"  # the unmodified C++ program
NUM = re.compile(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?")

LINK_VIEWS = True
ZOOM = 1.1


def floats(line):
    return np.array([float(x) for x in NUM.findall(line)])


def run_cpp():
    out = subprocess.run([CLI], capture_output=True, text=True, check=True).stdout
    experiments, current = [], None
    for line in out.splitlines():
        line = line.strip()
        if not line:
            continue
        key, sep, rest = line.partition(":")
        if not sep:
            current = {"name": line}
            experiments.append(current)
        elif current is not None and key in ("base", "average", "projection", "scales"):
            vals = floats(rest)
            if vals.size:
                current[key] = vals
    return experiments


def add_odeco_surface(plotter, coords, opacity=1.0):
    sh_values, x, y, z = generate_sh_values_from_coordinates(coords)
    grid = pv.StructuredGrid(x, y, z)
    grid["sh_values"] = np.asarray(sh_values).ravel(order="F")
    plotter.add_mesh(
        grid,
        scalars="sh_values",
        opacity=opacity,
        cmap="viridis",
        show_scalar_bar=False,
        smooth_shading=True,
    )


def fmt(v):
    return "(" + ", ".join(f"{x:.3f}" for x in v) + ")"


experiments = run_cpp()

# Columns: base, average (only if the C++ printed it), projection
has_average = any("average" in e and e["average"].size == 15 for e in experiments)
columns = ["base"] + (["average"] if has_average else []) + ["projection"]

plotter = pv.Plotter(
    shape=(len(experiments), len(columns)),
    window_size=(600 * len(columns), 600 * len(experiments)),
)

with tqdm(total=len(experiments) * len(columns), desc="Plotting") as pbar:
    for row, exp in enumerate(experiments):
        for col, key in enumerate(columns):
            pbar.set_postfix(experiment=exp["name"], column=key)
            plotter.subplot(row, col)

            title = f"{exp['name']} - {key}"
            if key == "projection" and "scales" in exp:
                title += f"\nscales: {fmt(exp['scales'])}"
            plotter.add_text(title, font_size=10)

            if key in exp and exp[key].size == 15:
                add_odeco_surface(plotter, exp[key])

            plotter.add_axes()      # coordinate system widget for this subplot
            plotter.reset_camera()  # fit this subplot's own geometry
            plotter.camera.zoom(ZOOM)
            pbar.update(1)

if LINK_VIEWS:
    plotter.link_views()

plotter.show(screenshot="./average.png")