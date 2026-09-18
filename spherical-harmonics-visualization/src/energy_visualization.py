"""
Parse `energy_current=<value>` lines out of the four snail_field_cxx
energy log files and plot all four series on one chart.

Usage:
    python plot_energies.py
"""

import os
import re
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# File locations (same paths as the original snippet)
# ---------------------------------------------------------------------------
odeco_sphere_path = os.path.abspath("../snail_field_cxx/build/sphere_energy.txt")
odeco_torus_path = os.path.abspath("../snail_field_cxx/build/torus_energy.txt")
octa_sphere_path = os.path.abspath("../snail_field_cxx_main/build/sphere_energy.txt")
octa_torus_path = os.path.abspath("../snail_field_cxx_main/build/torus_energy.txt")

FILES = {
    "odeco sphere": odeco_sphere_path,
    "odeco torus": odeco_torus_path,
    "octa sphere": octa_sphere_path,
    "octa torus": octa_torus_path,
}

# Matches lines like: energy_current=2.16134
ENERGY_RE = re.compile(r"^energy_current=([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*$")


def parse_energy_blocks(filepath: str) -> list[list[float]]:
    """Extract energy_current values from a log file, grouped into blocks.

    Each run of consecutive `energy_current=<num>` lines is one block
    (these correspond to one refinement/tree-depth stage, e.g. everything
    between two "Depth N tree has ..." lines). Any other line (mesh info,
    tree stats, timing, blank lines, etc.) ends the current block. Blocks
    are returned in file order; each block is the list of energy values
    for that stage.
    """
    blocks: list[list[float]] = []
    current: list[float] = []
    with open(filepath, "r") as f:
        for line in f:
            match = ENERGY_RE.match(line.strip())
            if match:
                current.append(float(match.group(1)))
            elif current:
                blocks.append(current)
                current = []
    if current:
        blocks.append(current)
    return blocks


def blocks_to_xy(blocks: list[list[float]]) -> tuple[list[float], list[float]]:
    """Flatten blocks into (x, y) with each block normalized to unit width.

    Block i occupies the x-range [i, i+1), with its samples evenly spaced
    across that range regardless of how many samples the block has. This
    way, block/refinement-stage boundaries line up at the same x position
    (0, 1, 2, ...) across files even when files have different numbers of
    samples per stage.
    """
    xs: list[float] = []
    ys: list[float] = []
    for i, block in enumerate(blocks):
        n = len(block)
        if n == 1:
            positions = [0.0]
        else:
            positions = [j / (n - 1) for j in range(n)]
        xs.extend(i + p for p in positions)
        ys.extend(block)
    return xs, ys


def main() -> None:
    plt.figure(figsize=(11, 6))

    max_blocks = 0
    for label, path in FILES.items():
        if not os.path.exists(path):
            print(f"Warning: file not found, skipping: {path}")
            continue

        blocks = parse_energy_blocks(path)
        if not blocks:
            print(f"Warning: no energy_current values found in: {path}")
            continue

        max_blocks = max(max_blocks, len(blocks))
        xs, ys = blocks_to_xy(blocks)
        total_n = sum(len(b) for b in blocks)
        plt.plot(xs, ys, label=f"{label} ({len(blocks)} stages, n={total_n})", linewidth=1.2)

    # Vertical guides + integer ticks at each refinement-stage boundary
    for i in range(max_blocks + 1):
        plt.axvline(i, color="gray", linewidth=0.5, alpha=0.3)

    plt.xlabel("Refinement stage (depth), each stage normalized to unit width")
    plt.ylabel("energy_current")
    plt.title("Energy convergence: odeco vs octa, sphere vs torus (aligned by refinement stage)")
    plt.legend()
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()

    out_path = "energy_plot.png"
    plt.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")
    plt.show()


if __name__ == "__main__":
    main()