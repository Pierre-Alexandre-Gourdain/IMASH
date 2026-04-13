import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="Directory containing IMASH output files")
args = parser.parse_args()

base = Path(args.dir).resolve()

# ============================================================
# Helpers
# ============================================================

def reshape_on_grid(x, y, f):
    xvals = np.unique(x)
    yvals = np.unique(y)

    nx = len(xvals)
    ny = len(yvals)

    if f.size != nx * ny:
        raise RuntimeError(
            f"reshape_on_grid: field has size {f.size}, but nx*ny = {nx*ny}"
        )

    F = f.reshape(nx, ny).T
    XX, YY = np.meshgrid(xvals, yvals, indexing="xy")
    return xvals, yvals, XX, YY, F


def plot_rho1(ax):
    ax.contour(RR, ZZ, rho2, levels=[1.0], colors="white", linewidths=2.5)
    ax.contour(RR, ZZ, rho2, levels=[1.0], colors="black", linewidths=1.0)


# ============================================================
# Load core / equilibrium-grid data
# ============================================================

core = np.loadtxt(base / "core_on_equilibrium_grid.dat", comments="#")
Rc   = core[:, 0]
Zc   = core[:, 1]
rho  = core[:, 2]

Rvals, Zvals, RR, ZZ, rho2 = reshape_on_grid(Rc, Zc, rho)

# ============================================================
# Load raw edge data
# ============================================================

edge_raw = np.loadtxt(base / "edge_cells_2d.dat", comments="#")
Re_raw   = edge_raw[:, 1]
Ze_raw   = edge_raw[:, 2]
ne_raw   = edge_raw[:, 3]
Te_raw   = edge_raw[:, 4]

tri_raw = mtri.Triangulation(Re_raw, Ze_raw)

# ============================================================
# Load interpolated edge data (C++ output)
# ============================================================

edge_int = np.loadtxt(base / "edge_cells_2d_interpolated.dat", comments="#")

Ri_int     = edge_int[:, 0]
Zi_int     = edge_int[:, 1]
ne_int     = edge_int[:, 2]
Te_int     = edge_int[:, 3]
rho_int    = edge_int[:, 4]
theta_int  = edge_int[:, 5]

Rvals_i, Zvals_i, RRi, ZZi, ne_i2 = reshape_on_grid(Ri_int, Zi_int, ne_int)
_,       _,       _,   _,   Te_i2 = reshape_on_grid(Ri_int, Zi_int, Te_int)
_,       _,       _,   _,   rho_i2 = reshape_on_grid(Ri_int, Zi_int, rho_int)
_,       _,       _,   _,   theta_i2 = reshape_on_grid(Ri_int, Zi_int, theta_int)

ne_i2 = np.maximum(np.nan_to_num(ne_i2, nan=0.0), 0.0)
Te_i2 = np.maximum(np.nan_to_num(Te_i2, nan=0.0), 0.0)

# ============================================================
# Diagnostics
# ============================================================

print("Raw edge ne min/max =", np.nanmin(ne_raw), np.nanmax(ne_raw))
print("Raw edge Te min/max =", np.nanmin(Te_raw), np.nanmax(Te_raw))
print("Interpolated ne min/max =", np.nanmin(ne_i2), np.nanmax(ne_i2))
print("Interpolated Te min/max =", np.nanmin(Te_i2), np.nanmax(Te_i2))

# ============================================================
# Plots
# ============================================================

fig, axs = plt.subplots(2, 3, figsize=(12, 10), constrained_layout=True)

# ------------------------------------------------------------
# Column 1: rho / theta
# ------------------------------------------------------------

# rho
pcm0 = axs[0, 0].pcolormesh(RRi, ZZi, rho_i2, shading="auto")
plot_rho1(axs[0, 0])
fig.colorbar(pcm0, ax=axs[0, 0], label="$\\rho$")
axs[0, 0].set_title("Normalized $\\rho$")

# theta
pcm1 = axs[1, 0].pcolormesh(RRi, ZZi, theta_i2, shading="auto")
plot_rho1(axs[1, 0])
fig.colorbar(pcm1, ax=axs[1, 0], label="$\\theta$ (rd)")
axs[1, 0].set_title("$\\theta$")

# ------------------------------------------------------------
# Column 2: density
# ------------------------------------------------------------

# Raw edge ne
tpc0 = axs[0, 1].tricontourf(tri_raw, np.log10(ne_raw + 1.0e14), levels=50)
plot_rho1(axs[0, 1])
fig.colorbar(tpc0, ax=axs[0, 1], label="$\\log_{10}n_e\\,(m^{-3})$")
axs[0, 1].set_title("Raw edge density\n(ignore core values)")

# Interpolated edge ne
pcm2 = axs[1, 1].pcolormesh(RRi, ZZi, np.log10(ne_i2 + 1.0e14), shading="auto")
plot_rho1(axs[1, 1])
fig.colorbar(pcm2, ax=axs[1, 1], label="$\\log_{10}n_e\\,(m^{-3})$")
axs[1, 1].set_title("Interpolated edge density")

# ------------------------------------------------------------
# Column 3: temperature
# ------------------------------------------------------------

# Raw edge Te
tpc1 = axs[0, 2].tricontourf(tri_raw, Te_raw, levels=50)
plot_rho1(axs[0, 2])
fig.colorbar(tpc1, ax=axs[0, 2], label="T (eV)")
axs[0, 2].set_title("Raw edge temperature\n(ignore core values)")

# Interpolated edge Te
pcm3 = axs[1, 2].pcolormesh(RRi, ZZi, Te_i2, shading="auto")
plot_rho1(axs[1, 2])
fig.colorbar(pcm3, ax=axs[1, 2], label="T (eV)")
axs[1, 2].set_title("Interpolated edge temperature")

# cosmetics
for ax in axs.flat:
    ax.set_xlabel("R")
    ax.set_ylabel("Z")
    ax.set_aspect("equal")

plt.show()