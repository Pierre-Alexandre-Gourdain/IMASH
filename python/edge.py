import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

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

def wrap_angle_diff(a, b):
    d = a - b
    return (d + np.pi) % (2.0*np.pi) - np.pi

# ============================================================
# Load core / equilibrium-grid data
# ============================================================

core = np.loadtxt("core_on_equilibrium_grid.dat", comments="#")
Rc   = core[:, 0]
Zc   = core[:, 1]
rho  = core[:, 2]
ne_c = core[:, 3]
Te_c = core[:, 4]

Rvals, Zvals, RR, ZZ, rho2  = reshape_on_grid(Rc, Zc, rho)
_,     _,     _,  _,  ne_c2 = reshape_on_grid(Rc, Zc, ne_c)
_,     _,     _,  _,  Te_c2 = reshape_on_grid(Rc, Zc, Te_c)

# ============================================================
# Load magnetic axis
# ============================================================

axis = np.loadtxt("magnetic_axis.dat")
R_axis = axis[0]
Z_axis = axis[1]

print("Magnetic axis =", R_axis, Z_axis)

## ============================================================
# Use magnetic axis as angular center
# ============================================================

R0 = R_axis
Z0 = Z_axis

theta2 = np.arctan2(ZZ - Z0, RR - R0)
# ============================================================
# Load raw edge points
# ============================================================

edge = np.loadtxt("edge_cells_2d.dat", comments="#")
Re   = edge[:, 1]
Ze   = edge[:, 2]
ne_e = edge[:, 3]
Te_e = edge[:, 4]

# ============================================================
# Interpolate rho from equilibrium grid to edge points
# ============================================================

eq_points = np.column_stack((Rc, Zc))
rho_e = griddata(eq_points, rho, (Re, Ze), method="linear")

# Remove points that failed rho interpolation
good = np.isfinite(rho_e)
Re   = Re[good]
Ze   = Ze[good]
ne_e = ne_e[good]
Te_e = Te_e[good]
rho_e = rho_e[good]

# Edge angle
theta_e = np.arctan2(Ze - Z0, Re - R0)

# ============================================================
# Duplicate in theta to handle periodicity
# ============================================================

rho_ext   = np.concatenate([rho_e, rho_e, rho_e])
theta_ext = np.concatenate([theta_e - 2.0*np.pi, theta_e, theta_e + 2.0*np.pi])
ne_ext    = np.concatenate([ne_e, ne_e, ne_e])
Te_ext    = np.concatenate([Te_e, Te_e, Te_e])

edge_param_points = np.column_stack((rho_ext, theta_ext))

# ============================================================
# Evaluate edge interpolation on equilibrium grid in (rho,theta)
# ============================================================

rho_flat   = rho2.ravel()
theta_flat = theta2.ravel()

target_points = np.column_stack((rho_flat, theta_flat))

ne_e_flat = griddata(edge_param_points, ne_ext, target_points, method="linear")
Te_e_flat = griddata(edge_param_points, Te_ext, target_points, method="linear")

ne_e2 = ne_e_flat.reshape(rho2.shape)
Te_e2 = Te_e_flat.reshape(rho2.shape)

# Clean up
ne_e2 = np.nan_to_num(ne_e2, nan=0.0)
Te_e2 = np.nan_to_num(Te_e2, nan=0.0)

ne_e2 = np.maximum(ne_e2, 0.0)
Te_e2 = np.maximum(Te_e2, 0.0)

# ============================================================
# Merge: core inside, edge outside
# ============================================================

ne_T = ne_c2.copy()
Te_T = Te_c2.copy()

outside_sep = rho2 > 1.0
ne_T[outside_sep] = ne_e2[outside_sep]
Te_T[outside_sep] = Te_e2[outside_sep]

# ============================================================
# Diagnostics
# ============================================================

print("R0, Z0 =", R0, Z0)
print("edge rho min/max =", np.nanmin(rho_e), np.nanmax(rho_e))
print("edge theta min/max =", np.nanmin(theta_e), np.nanmax(theta_e))
print("edge ne min/max =", np.nanmin(ne_e2), np.nanmax(ne_e2))
print("edge Te min/max =", np.nanmin(Te_e2), np.nanmax(Te_e2))
print("Max edge density inside core:",np.max(ne_e2[rho2 <= 1.0]))
print("Max edge density outisde core:",np.max(ne_e2[rho2 > 1.0]))
fig, axs = plt.subplots(1, 4, figsize=(18, 5), constrained_layout=True)

# ------------------------------------------------------------
# Edge interpolation
# ------------------------------------------------------------
pcm0 = axs[0].pcolormesh(RR, ZZ, np.log10(ne_e2 + 1.0e14), shading="auto")
plot_rho1(axs[0])
fig.colorbar(pcm0, ax=axs[0], label=r"$\log_{10} n_{e,\mathrm{edge}}$")
axs[0].set_title("Edge density $(\\rho,\\theta)$")
axs[0].set_xlabel("R")
axs[0].set_ylabel("Z")
axs[0].set_aspect("equal")

# ------------------------------------------------------------
# Total field
# ------------------------------------------------------------
pcm1 = axs[1].pcolormesh(RR, ZZ, np.log10(ne_T + 1.0e14), shading="auto")
plot_rho1(axs[1])
fig.colorbar(pcm1, ax=axs[1], label=r"$\log_{10} n_e$")
axs[1].set_title("Total density")
axs[1].set_xlabel("R")
axs[1].set_ylabel("Z")
axs[1].set_aspect("equal")

# ------------------------------------------------------------
# Poloidal angle
# ------------------------------------------------------------
pcm2 = axs[2].pcolormesh(RR, ZZ, theta2, shading="auto")
plot_rho1(axs[2])
fig.colorbar(pcm2, ax=axs[2], label=r"$\theta$")
axs[2].set_title("Poloidal angle")
axs[2].set_xlabel("R")
axs[2].set_ylabel("Z")
axs[2].set_aspect("equal")

# ------------------------------------------------------------
# NEW: rho plot
# ------------------------------------------------------------
pcm3 = axs[3].pcolormesh(RR, ZZ, rho2, shading="auto")
plot_rho1(axs[3])
fig.colorbar(pcm3, ax=axs[3], label=r"$\rho$")
axs[3].set_title("Normalized flux $\\rho$")
axs[3].set_xlabel("R")
axs[3].set_ylabel("Z")
axs[3].set_aspect("equal")

plt.show()