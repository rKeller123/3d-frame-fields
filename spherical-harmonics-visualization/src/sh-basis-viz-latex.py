"""Plot functions given by coefficients in an even-order real SH basis (l = 0, 2, 4),
as several frames next to each other.

Basis ordering (15 coefficients), the usual one for symmetric SH bases:
    index j = l (l + 1) / 2 + m,  for l = 0, 2, 4 and m = -l..l
    j = 0          -> l=0, m=0
    j = 1..5       -> l=2, m=-2..2
    j = 6..14      -> l=4, m=-4..4
So j=10 is (l=4, m=0) and j=14 is (l=4, m=4).

f(theta, phi) = sum_j c_j * Y_{l_j m_j}(theta, phi)

Surface radius = |f|, color = signed value of f (viridis; see COLOR_RANGE).
All frames share ONE spatial scale and ONE color scale, so they are directly
comparable. A single colorbar is drawn below the frames.

Output: sh_combination.pdf (vector text, rasterized surfaces) and a 400 dpi PNG.
Requires: numpy, matplotlib. Expects spherical_harmonics.py and odeco.py in the same folder.
"""
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource, Normalize

from spherical_harmonics import sph_harm_real
from odeco import generate_coordinates

# ----------------------------------------------------------------- input
# One entry per frame. Each entry is passed to generate_coordinates(*entry).
COORDINATES = [
    (1, 1, 1, 0, 0, 0),
    (2, 1, 1.5, 0, 0, 0),
    (0, 1, 2, 0, 0, 0),
]

# Optional titles above the frames (same length as COORDINATES), or None.
# Math is fine, e.g. [r"$\lambda=(1,2,1)$", ...]
LABELS = None

N_COLS = None   # frames per row; None = all frames in one row

# Alternatively, give the 15 coefficients directly instead of using odeco, e.g.
#   octahedral canonical: [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, math.sqrt(7 / 12), 0, 0, 0, math.sqrt(5 / 12)]
# BASES = [ [..15 numbers..], [..15 numbers..] ]
BASES = [generate_coordinates(*c) for c in COORDINATES]
# BASES = [[2, 0, 0, 0, 0, 0, 0, 0, 0, 0, math.sqrt(7 / 12), 0, 0, 0, math.sqrt(5 / 12)]]

N = 100  # subdivisions of theta and phi

# ----------------------------------------------------------------- settings
CMAP_NAME = "plasma"
# Color range shared by all frames:
#   "auto"      : if f >= 0 everywhere use [0, max]; otherwise symmetric [-max, +max]
#   "symmetric" : always [-max, +max] (zero at the center of the colormap)
#   "data"      : [min f, max f]
COLOR_RANGE = "auto"
PANEL = 2.3                 # width and height of one frame, inches
CBAR_BLOCK = 0           # space reserved for the colorbar at the bottom, inches
TITLE_STRIP = 0.14          # fraction of a frame's height reserved for its title
VIEW_ELEV, VIEW_AZIM = 22, -55
ZOOM = 1.3
FS_TITLE, FS_TICK = 13, 12

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.unicode_minus": True,
})

_BASIS_LM = [(l, m) for l in (0, 2, 4) for m in range(-l, l + 1)]
_Y_CACHE = {}


def basis_indices(l_values=(0, 2, 4)):
    """(l, m) pair for each coefficient index j = l(l+1)/2 + m."""
    return [(l, m) for l in l_values for m in range(-l, l + 1)]


def basis_function(l, m, theta, phi):
    """Y_lm on the (theta, phi) grid; cached so frames can reuse it."""
    key = (l, m)
    if key not in _Y_CACHE:
        Y = np.empty((len(theta), len(phi)))
        for i, t in enumerate(theta):
            for j, p in enumerate(phi):
                Y[i, j] = sph_harm_real(l, m, t, p)
        _Y_CACHE[key] = Y
    return _Y_CACHE[key]


def evaluate(coeffs, theta, phi):
    """f(theta, phi) on the grid, using only the nonzero coefficients."""
    lm = basis_indices()
    assert len(coeffs) == len(lm), f"expected {len(lm)} coefficients, got {len(coeffs)}"
    f = np.zeros((len(theta), len(phi)))
    for c, (l, m) in zip(coeffs, lm):
        if c == 0:
            continue
        f += c * basis_function(l, m, theta, phi)
    return f


def main():
    n = len(BASES)
    if LABELS is not None:
        assert len(LABELS) == n, "LABELS must have one entry per frame"
    n_cols = n if N_COLS is None else min(N_COLS, n)
    n_rows = math.ceil(n / n_cols)

    theta = np.linspace(0, np.pi, N)       # polar angle
    phi = np.linspace(0, 2 * np.pi, N)     # azimuth
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")

    fs = [evaluate(b, theta, phi) for b in BASES]
    fmax = max(np.abs(f).max() for f in fs)       # shared by all frames
    fmin_data = min(f.min() for f in fs)
    fmax_data = max(f.max() for f in fs)
    for k, f in enumerate(fs):
        print(f"frame {k}: f in [{f.min():.4f}, {f.max():.4f}]")

    if COLOR_RANGE == "symmetric" or (COLOR_RANGE == "auto" and fmin_data < -1e-9):
        norm = Normalize(vmin=-fmax, vmax=fmax)
    elif COLOR_RANGE == "auto":                   # nonnegative everywhere
        norm = Normalize(vmin=0.0, vmax=fmax)
    elif COLOR_RANGE == "data":
        norm = Normalize(vmin=fmin_data, vmax=fmax_data)
    else:
        raise ValueError(f"unknown COLOR_RANGE: {COLOR_RANGE!r}")
    cmap = plt.get_cmap(CMAP_NAME)
    light = LightSource(azdeg=300, altdeg=55)

    fig_w = n_cols * PANEL
    fig_h = n_rows * PANEL + CBAR_BLOCK
    fig = plt.figure(figsize=(fig_w, fig_h))

    for k, f in enumerate(fs):
        row, col = divmod(k, n_cols)
        r = np.abs(f)
        x = r * np.sin(THETA) * np.cos(PHI)
        y = r * np.sin(THETA) * np.sin(PHI)
        z = r * np.cos(THETA)

        # cell rectangle in figure fractions
        x0 = col * PANEL / fig_w
        y0 = (CBAR_BLOCK + (n_rows - 1 - row) * PANEL) / fig_h
        w, h = PANEL / fig_w, PANEL / fig_h
        title_h = TITLE_STRIP if LABELS is not None else 0.0

        ax = fig.add_axes([x0, y0, w, h * (1 - title_h)], projection="3d")
        ax.plot_surface(x, y, z, facecolors=cmap(norm(f)), shade=True, lightsource=light,
                        rstride=1, cstride=1, linewidth=0, antialiased=False,
                        rasterized=True)
        ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
        ax.set_xlim(-fmax, fmax); ax.set_ylim(-fmax, fmax); ax.set_zlim(-fmax, fmax)
        ax.set_box_aspect((1, 1, 1), zoom=ZOOM)
        ax.set_axis_off()

        if LABELS is not None:
            fig.text(x0 + w / 2, y0 + h * (1 - title_h / 2), LABELS[k],
                     fontsize=FS_TITLE, ha="center", va="center")

    # one shared colorbar, centered under the frames
    # cb_w = min(0.5, 3.0 / fig_w)   # fraction of the figure width (max 3 in wide)
    # cax = fig.add_axes([0.5 - cb_w / 2, 0.36 / fig_h, cb_w, 0.14 / fig_h])
    # cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
    #                     orientation="horizontal")
    # cbar.ax.tick_params(labelsize=FS_TICK, length=4, width=0.8)
    # cbar.outline.set_linewidth(0.8)

    fig.savefig("sh_combination.pdf", dpi=400)
    fig.savefig("sh_combination.png", dpi=400)
    print("Saved sh_combination.pdf and sh_combination.png")
    plt.show()


if __name__ == "__main__":
    main()