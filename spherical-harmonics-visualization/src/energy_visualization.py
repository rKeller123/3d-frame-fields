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


def parse_energy(filepath: str) -> list[float]:
    """Extract energy_current values from a log file, flattened to one list
    """

    energy: list[float] = []
    with open(filepath, "r") as f:
        for line in f:
            match = ENERGY_RE.match(line.strip())
            if match:
                energy.append(float(match.group(1)))
    return energy


def main() -> None:
    plt.figure(figsize=(11, 6))

    for label, path in FILES.items():
        if not os.path.exists(path):
            print(f"Warning: file not found, skipping: {path}")
            continue

        energy = parse_energy(path)
        if not energy:
            print(f"Warning: no energy_current values found in: {path}")
            continue

        total_n = len(energy)
        plt.plot(energy, label=f"{label} (n={total_n})", linewidth=1.2)

    plt.xlabel("Iteration")
    plt.ylabel("energy_current")
    plt.title("Energy convergence: odeco vs octa, sphere vs torus")
    plt.legend()
    plt.grid(True, axis="y", alpha=0.3)
    plt.tight_layout()

    out_path = "energy_plot.png"
    plt.savefig(out_path, dpi=150)
    print(f"Saved plot to {out_path}")
    plt.show()


if __name__ == "__main__":
    main()