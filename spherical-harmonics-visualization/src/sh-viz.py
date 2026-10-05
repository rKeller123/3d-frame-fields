"""Publication-quality visualization of real spherical harmonics, l = 0..3.

Surface radius = |Y_lm(theta, phi)|, surface color = signed value of Y_lm
(viridis, one symmetric color scale shared by all panels). All panels also
share the same spatial scale, so relative lobe sizes are directly comparable.

Outputs: spherical_harmonics.pdf (vector text, rasterized surfaces) and
         spherical_harmonics.png (400 dpi).

Requires: numpy, matplotlib. Expects spherical_harmonics.py in the same folder.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LightSource, Normalize
from matplotlib.gridspec import GridSpec

from spherical_harmonics import sph_harm_real

# ----------------------------------------------------------------- settings
L_MAX = 3
N_THETA, N_PHI = 120, 240          # mesh resolution
CMAP_NAME = "viridis"
FIG_WIDTH, FIG_HEIGHT = 7.5, 6.3   # inches; design at final printed size
VIEW_ELEV, VIEW_AZIM = 20, -55
ZOOM = 1.5                         # enlarges each 3D panel within its cell

FS_ROW = 15       # row labels (l = ...)
FS_COL = 13       # panel titles (m = ...)
FS_CBAR = 14      # colorbar label
FS_TICK = 12      # colorbar ticks

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "axes.unicode_minus": True,
})

# ------------------------------------------------------------------- grid
theta = np.linspace(0, np.pi, N_THETA)      # polar angle
phi = np.linspace(0, 2 * np.pi, N_PHI)      # azimuth
THETA, PHI = np.meshgrid(theta, phi, indexing="ij")


def evaluate(l: int, m: int) -> np.ndarray:
    """Evaluate Y_lm on the whole (theta, phi) grid using the scalar function."""
    Y = np.empty_like(THETA)
    for i in range(N_THETA):
        for j in range(N_PHI):
            Y[i, j] = sph_harm_real(l, m, theta[i], phi[j])
    return Y


def to_cartesian(r):
    return (r * np.sin(THETA) * np.cos(PHI),
            r * np.sin(THETA) * np.sin(PHI),
            r * np.cos(THETA))


def main():
    data = {(l, m): evaluate(l, m)
            for l in range(L_MAX + 1) for m in range(-l, l + 1)}
    vmax = max(np.abs(Y).max() for Y in data.values())
    norm = Normalize(vmin=-vmax, vmax=vmax)
    cmap = plt.get_cmap(CMAP_NAME)
    light = LightSource(azdeg=300, altdeg=55)

    n_rows, n_cols = L_MAX + 1, 2 * L_MAX + 1
    fig = plt.figure(figsize=(FIG_WIDTH, FIG_HEIGHT))
    gs = GridSpec(n_rows, n_cols, figure=fig, left=0.075, right=0.995,
                  top=0.99, bottom=0.17, hspace=0.0, wspace=0.0)

    for (l, m), Y in data.items():
        x, y, z = to_cartesian(np.abs(Y))
        # manual placement: lower part of the cell for the plot, top strip for the title
        cell = gs[l, m + L_MAX].get_position(fig)
        ax = fig.add_axes([cell.x0, cell.y0, cell.width, cell.height * 0.84],
                          projection="3d")
        fig.text(cell.x0 + cell.width / 2, cell.y0 + cell.height * 0.90,
                 f"$m={m}$", fontsize=FS_COL, ha="center", va="center")
        ax.plot_surface(
            x, y, z, facecolors=cmap(norm(Y)), shade=True, lightsource=light,
            rstride=1, cstride=1, linewidth=0, antialiased=False,
            rasterized=True,
        )
        ax.view_init(elev=VIEW_ELEV, azim=VIEW_AZIM)
        # common limits across all panels -> sizes are comparable
        ax.set_xlim(-vmax, vmax); ax.set_ylim(-vmax, vmax); ax.set_zlim(-vmax, vmax)
        ax.set_box_aspect((1, 1, 1), zoom=ZOOM)
        ax.set_axis_off()

    # row labels (l = ...) to the left of each row
    for l in range(n_rows):
        pos = gs[l, 0].get_position(fig)
        fig.text(0.012, pos.y0 + pos.height * 0.42, rf"$\ell={l}$",
                 fontsize=FS_ROW, ha="left", va="center")

    # horizontal colorbar below the figure
    cax = fig.add_axes([0.22, 0.095, 0.56, 0.022])
    cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                        orientation="horizontal")
    cbar.set_label(r"$Y_{\ell}^m(\theta,\phi)$", fontsize=FS_CBAR, labelpad=4)
    cbar.ax.tick_params(labelsize=FS_TICK, length=4, width=0.8)
    cbar.outline.set_linewidth(0.8)

    fig.savefig("spherical_harmonics.pdf", dpi=400)
    fig.savefig("spherical_harmonics.png", dpi=400)
    print("Saved spherical_harmonics.pdf and spherical_harmonics.png")
    plt.show()


if __name__ == "__main__":
    main()