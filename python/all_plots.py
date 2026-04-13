import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("dir", help="Directory containing IMASH output files")
args = parser.parse_args()

base = Path(args.dir).resolve()

# ============================================================
# Helpers
# ============================================================

# Separatrix
def plot_rho1(ax):
    ax.contour(RR, ZZ, rho2, levels=[1.0],
               colors='white', linewidths=2.5)
    ax.contour(RR, ZZ, rho2, levels=[1.0],
               colors='black', linewidths=1)

# ============================================================
# Load core data
# ============================================================
core = np.loadtxt(base / "core_on_equilibrium_grid.dat", comments="#")
Rc   = core[:, 0]
Zc   = core[:, 1]
rho  = core[:, 2]
ne_c = core[:, 3]
Te_c = core[:, 4]

Rvals = np.unique(Rc)
Zvals = np.unique(Zc)

nR = len(Rvals)
nZ = len(Zvals)

rho2  = rho.reshape(nR, nZ).T
ne_c2 = ne_c.reshape(nR, nZ).T
Te_c2 = Te_c.reshape(nR, nZ).T

RR, ZZ = np.meshgrid(Rvals, Zvals, indexing="xy")

# ============================================================
# Load interpolated edge data
# ============================================================
edge = np.loadtxt(base / "edge_cells_2d_interpolated.dat", comments="#")

# Expected columns:
# R Z ne Te rho theta
Re    = edge[:, 0]
Ze    = edge[:, 1]
ne_e  = edge[:, 2]
Te_e  = edge[:, 3]

Rvals_e = np.unique(Re)
Zvals_e = np.unique(Ze)

nRe = len(Rvals_e)
nZe = len(Zvals_e)

ne_e2 = ne_e.reshape(nRe, nZe).T
Te_e2 = Te_e.reshape(nRe, nZe).T

RRe, ZZe = np.meshgrid(Rvals_e, Zvals_e, indexing="xy")

ne_e2 = np.nan_to_num(ne_e2, nan=0.0)
Te_e2 = np.nan_to_num(Te_e2, nan=0.0)

ne_e2 = np.maximum(ne_e2, 0.0)
Te_e2 = np.maximum(Te_e2, 0.0)

# Total = core inside + edge outside
ne_T = ne_c2.copy()
Te_T = Te_c2.copy()

outside_sep = rho2 > 1.0
ne_T[outside_sep] = ne_e2[outside_sep]
Te_T[outside_sep] = Te_e2[outside_sep]

# ============================================================
# Load magnetic field data
# ============================================================
bf   = np.loadtxt(base / "bfield.dat", comments="#")
Rb   = bf[:, 0]
Zb   = bf[:, 1]
BR   = bf[:, 2]
BZ   = bf[:, 3]
Bphi = bf[:, 4]
Bpol = bf[:, 5]
Bmag = bf[:, 6]

Rvals_b = np.unique(Rb)
Zvals_b = np.unique(Zb)

nRb = len(Rvals_b)
nZb = len(Zvals_b)

BR2   = BR.reshape(nRb, nZb).T
BZ2   = BZ.reshape(nRb, nZb).T
Bpol2 = Bpol.reshape(nRb, nZb).T
Bmag2 = Bmag.reshape(nRb, nZb).T
Bphi2 = Bphi.reshape(nRb, nZb).T

RRb, ZZb = np.meshgrid(Rvals_b, Zvals_b, indexing="xy")

# ============================================================
# Cold / hot electron masses
# ============================================================
eps0 = 8.8541878128e-12
echarge = 1.602176634e-19
m0e = 9.1093837015e-31
c = 299792458.0

m_cold = m0e * np.ones_like(Te_T)
m_hot  = m0e * np.sqrt(1.0 + 5.0 * echarge * Te_T / (m0e * c**2))

# ============================================================
# Frequencies
# ============================================================
omega_pe_cold = np.sqrt(echarge**2 * ne_T / (eps0 * m_cold))
omega_pe_hot  = np.sqrt(echarge**2 * ne_T / (eps0 * m_hot))

Omega_e_cold = echarge * Bmag2 / m_cold
Omega_e_hot  = echarge * Bmag2 / m_hot

omega_R_cold = 0.5 * Omega_e_cold + np.sqrt(0.25 * Omega_e_cold**2 + omega_pe_cold**2)
omega_R_hot  = 0.5 * Omega_e_hot  + np.sqrt(0.25 * Omega_e_hot**2  + omega_pe_hot**2)

# Convert to GHz
GHz = 1.0e9
f_O_cold_GHz = omega_pe_cold / (2.0 * np.pi * GHz)
f_O_hot_GHz  = omega_pe_hot  / (2.0 * np.pi * GHz)

f_R_cold_GHz = omega_R_cold / (2.0 * np.pi * GHz)
f_R_hot_GHz  = omega_R_hot  / (2.0 * np.pi * GHz)

# Common contour levels in GHz
def build_levels(a, b, nlev=8):
    vals = np.concatenate([
        a[np.isfinite(a) & (a > 0.0)],
        b[np.isfinite(b) & (b > 0.0)]
    ])
    vmin = np.percentile(vals, 15)
    vmax = np.percentile(vals, 85)
    return np.linspace(vmin, vmax, nlev)

levels_O = build_levels(f_O_cold_GHz, f_O_hot_GHz, nlev=8)
levels_R = build_levels(f_R_cold_GHz, f_R_hot_GHz, nlev=8)

# ============================================================
# Make combined figure
# ============================================================
fig, axs = plt.subplots(2, 4, figsize=(16, 10), constrained_layout=True)

# ------------------------------------------------------------
# Core ne
# ------------------------------------------------------------
pcm0 = axs[0, 0].pcolormesh(RR, ZZ, ne_c2, shading="auto")
plot_rho1(axs[0, 0])
fig.colorbar(pcm0, ax=axs[0, 0], label="$n_e\\,(m^{-3})$")
axs[0, 0].set_title("Core density")
axs[0, 0].set_xlabel("R")
axs[0, 0].set_ylabel("Z")
axs[0, 0].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Core Te
# ------------------------------------------------------------
pcm1 = axs[0, 1].pcolormesh(RR, ZZ, Te_c2 * 1e-3, shading="auto")
plot_rho1(axs[0, 1])
fig.colorbar(pcm1, ax=axs[0, 1], label="$T_e\\,(keV)$")
axs[0, 1].set_title("Core temperature")
axs[0, 1].set_xlabel("R")
axs[0, 1].set_ylabel("Z")
axs[0, 1].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Poloidal magnetic field magnitude
# ------------------------------------------------------------
pcm2 = axs[0, 2].pcolormesh(RRb, ZZb, Bpol2, shading="auto")
fig.colorbar(pcm2, ax=axs[0, 2], label="$B_{pol}\\,(T)$")
axs[0, 2].set_title("Poloidal magnetic field magnitude")
axs[0, 2].set_xlabel("R")
axs[0, 2].set_ylabel("Z")
axs[0, 2].set_aspect("equal", adjustable="box")
axs[0, 2].streamplot(Rvals_b, Zvals_b, BR2, BZ2, density=1.2, color="red")
plot_rho1(axs[0, 2])

# ------------------------------------------------------------
# Interpolated edge ne
# ------------------------------------------------------------
pcm_edge_ne = axs[1, 0].pcolormesh(RRe, ZZe, ne_e2, shading="auto")
plot_rho1(axs[1, 0])
fig.colorbar(pcm_edge_ne, ax=axs[1, 0], label="$n_e\\,(m^{-3})$")
axs[1, 0].set_title("Edge density")
axs[1, 0].set_xlabel("R")
axs[1, 0].set_ylabel("Z")
axs[1, 0].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Interpolated edge Te
# ------------------------------------------------------------
pcm_edge_te = axs[1, 1].pcolormesh(RRe, ZZe, Te_e2, shading="auto")
plot_rho1(axs[1, 1])
fig.colorbar(pcm_edge_te, ax=axs[1, 1], label="$T_e\\,(eV)$")
axs[1, 1].set_title("Edge temperature")
axs[1, 1].set_xlabel("R")
axs[1, 1].set_ylabel("Z")
axs[1, 1].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Total density on equilibrium grid
# ------------------------------------------------------------
pcm3 = axs[1, 2].pcolormesh(RR, ZZ, np.log10(np.abs(ne_T)+np.average(ne_c2)/20.), shading="auto")
plot_rho1(axs[1, 2])
fig.colorbar(pcm3, ax=axs[1, 2], label="$\\log_{10}n_e\\,(m^{-3})$")
axs[1, 2].set_title("Total electron density")
axs[1, 2].set_xlabel("R")
axs[1, 2].set_ylabel("Z")
axs[1, 2].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Toroidal magnetic field with omega_R overlay
# ------------------------------------------------------------
pcm4 = axs[0, 3].pcolormesh(RR, ZZ, Bphi2, shading="auto")
fig.colorbar(pcm4, ax=axs[0, 3], label="$B_\\phi\\,(T)$")

cR_cold = axs[0, 3].contour(
    RR, ZZ, f_R_cold_GHz,
    levels=levels_R, colors="blue", linewidths=1.5
)
cR_hot = axs[0, 3].contour(
    RR, ZZ, f_R_hot_GHz,
    levels=levels_R, colors="red", linewidths=1.2
)

axs[0, 3].clabel(cR_cold, fmt="%.1f GHz", fontsize=7, colors="blue")
axs[0, 3].clabel(cR_hot,  fmt="%.1f GHz", fontsize=7, colors="red")

plot_rho1(axs[0, 3])
axs[0, 3].set_title(r"$B_\phi$ with $\omega_R$ contours")
axs[0, 3].set_xlabel("R")
axs[0, 3].set_ylabel("Z")
axs[0, 3].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Total density with omega_O overlay
# Split contours into core-only and edge-only
# ------------------------------------------------------------

pcm5 = axs[1, 3].pcolormesh(
    RR, ZZ,
    np.log10(np.abs(ne_T) + np.average(ne_c2)/20.0),
    shading="auto"
)
fig.colorbar(pcm5, ax=axs[1, 3], label="$\\log_{10}n_e\\,(m^{-3})$")

mask_core = rho2 <= 1.0
mask_edge = rho2 > 1.0

# Masked fields
f_O_cold_core = np.where(mask_core, f_O_cold_GHz, np.nan)
f_O_hot_core  = np.where(mask_core, f_O_hot_GHz,  np.nan)

f_O_cold_edge = np.where(mask_edge, f_O_cold_GHz, np.nan)
f_O_hot_edge  = np.where(mask_edge, f_O_hot_GHz,  np.nan)

def build_levels_masked(a, b, mask, nlev=8, pmin=10, pmax=90):
    vals = np.concatenate([
        a[mask & np.isfinite(a) & (a > 0.0)],
        b[mask & np.isfinite(b) & (b > 0.0)]
    ])
    if vals.size == 0:
        return None
    vmin = np.percentile(vals, pmin)
    vmax = np.percentile(vals, pmax)
    if vmax <= vmin:
        vmax = vmin + 1.0e-6
    return np.linspace(vmin, vmax, nlev)

# Separate contour levels
levels_O_core = build_levels_masked(f_O_cold_GHz, f_O_hot_GHz, mask_core, nlev=8, pmin=10, pmax=90)
levels_O_edge = build_levels_masked(f_O_cold_GHz, f_O_hot_GHz, mask_edge, nlev=8, pmin=10, pmax=90)

# Core contours
if levels_O_core is not None:
    cO_cold_core = axs[1, 3].contour(
        RR, ZZ, f_O_cold_core,
        levels=levels_O_core, colors="blue", linewidths=1.5
    )
    cO_hot_core = axs[1, 3].contour(
        RR, ZZ, f_O_hot_core,
        levels=levels_O_core, colors="red", linewidths=1.2
    )

    axs[1, 3].clabel(cO_cold_core, fmt="%.1f GHz", fontsize=7, colors="blue")
    axs[1, 3].clabel(cO_hot_core,  fmt="%.1f GHz", fontsize=7, colors="red")

# Edge contours
if levels_O_edge is not None:
    cO_cold_edge = axs[1, 3].contour(
        RR, ZZ, f_O_cold_edge,
        levels=levels_O_edge, colors="cyan", linewidths=1.2, linestyles="--"
    )
    cO_hot_edge = axs[1, 3].contour(
        RR, ZZ, f_O_hot_edge,
        levels=levels_O_edge, colors="magenta", linewidths=1.0, linestyles="--"
    )

    axs[1, 3].clabel(cO_cold_edge, fmt="%.1f GHz", fontsize=7, colors="cyan")
    axs[1, 3].clabel(cO_hot_edge,  fmt="%.1f GHz", fontsize=7, colors="magenta")

plot_rho1(axs[1, 3])
axs[1, 3].set_title(r"Total density with split $\omega_O=\omega_{pe}$ contours")
axs[1, 3].set_xlabel("R")
axs[1, 3].set_ylabel("Z")
axs[1, 3].set_aspect("equal", adjustable="box")

plt.show()