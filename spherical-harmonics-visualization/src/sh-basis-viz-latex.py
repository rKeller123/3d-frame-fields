"""Plot a function given by coefficients in an even-order real SH basis (l = 0, 2, 4).

Basis ordering (15 coefficients), the usual one for symmetric SH bases:
    index j = l (l + 1) / 2 + m,  for l = 0, 2, 4 and m = -l..l
    j = 0          -> l=0, m=0
    j = 1..5       -> l=2, m=-2..2
    j = 6..14      -> l=4, m=-4..4
So j=10 is (l=4, m=0) and j=14 is (l=4, m=4).

f(theta, phi) = sum_j c_j * Y_{l_j m_j}(theta, phi)

Surface radius = |f|, color = signed value of f (viridis, symmetric scale).

Output: sh_combination.pdf (vector text, rasterized surface) and a 400 dpi PNG.
Requires: numpy, matplotlib. Expects spherical_harmonics.py in the same folder.
"""
import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource, Normalize

from spherical_harmonics import sph_harm_real

# --------------------------------------------------------------- coefficients
# use sh basis in orders 0, 2 and 4. 15 coefficients.
basis = [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, math.sqrt(7 / 12), 0, 0, 0, math.sqrt(5 / 12)] # octahedral canonical

N = 100  # subdivisions of theta and phi

# ----------------------------------------------------------------- settings
CMAP_NAME = "viridis"
FIG_SIZE = (5.0, 5.6)
VIEW_ELEV, VIEW_AZIM = 22, -55
FS_CBAR, FS_TICK = 14, 12

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.unicode_minus": True,
})


def basis_indices(l_values=(0, 2, 4)):
    """(l, m) pair for each coefficient index j = l(l+1)/2 + m."""
    return [(l, m) for l in l_values for m in range(-l, l + 1)]


def evaluate(coeffs, theta, phi):
    """f(theta, phi) on the grid, using only the nonzero coefficients."""
    lm = basis_indices()
    assert len(coeffs) == len(lm), f"expected {len(lm)} coefficients, got {len(coeffs)}"
    f = np.zeros((len(theta), len(phi)))
    for c, (l, m) in zip(coeffs, lm):
        if c == 0:
            continue
        for i, t in enumerate(theta):
            for j, p in enumerate(phi):
                f[i, j] += c * sph_harm_real(l, m, t, p)
    return f


def main():
    theta = np.linspace(0, np.pi, N)       # polar angle
    phi = np.linspace(0, 2 * np.pi, N)     # azimuth
    THETA, PHI = np.meshgrid(theta, phi, indexing="ij")

    f = evaluate(basis, theta, phi)
    r = np.abs(f)
    x = r * np.sin(THETA) * np.cos(PHI)
    y = r * np.sin(THETA) * np.sin(PHI)
    z = r * np.cos(THETA)

    fmax = np.abs(f).max()
    norm = Normalize(vmin=-fmax, vmax=fmax)
    cmap = plt.get_cmap(CMAP_NAME)
    light = LightSource(azdeg=300, altdeg=55)

    fig = plt.figure(figsize=FIG_SIZE)
    ax = fig.add_axes([0.0, 0.19, 1.0, 0.81], projection="3d")
    ax.plot_surface(x, y, z, facecolors=cmap(norm(f)), shade=True, lightsource=light,
                    rstride=1, cstride=1, linewidth=0, antialiased=False,
                    rasterized=True)
    ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
    ax.set_xlim(-fmax, fmax); ax.set_ylim(-fmax, fmax); ax.set_zlim(-fmax, fmax)
    ax.set_box_aspect((1, 1, 1), zoom=1.35)
    ax.set_axis_off()

    cax = fig.add_axes([0.2, 0.125, 0.6, 0.028])
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                        orientation="horizontal")
    cbar.ax.tick_params(labelsize=FS_TICK, length=4, width=0.8)
    cbar.outline.set_linewidth(0.8)

    fig.savefig("sh_combination.pdf", dpi=400)
    fig.savefig("sh_combination.png", dpi=400)
    print("Saved sh_combination.pdf and sh_combination.png")
    plt.show()


if __name__ == "__main__":
    main()