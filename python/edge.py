import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
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
    return (d + np.pi) % (2.0 * np.pi) - np.pi


def local_interp_rho_theta(rho_q, theta_q,
                           rho_s, theta_s, val_s,
                           sigma_rho, sigma_theta,
                           max_drho=None, max_dtheta=None,
                           min_weight_sum=1.0e-12):
    out = np.zeros_like(rho_q)

    for j in range(rho_q.shape[0]):
        for i in range(rho_q.shape[1]):
            rq = rho_q[j, i]
            tq = theta_q[j, i]

            drho = rho_s - rq
            dtheta = wrap_angle_diff(theta_s, tq)

            if max_drho is None:
                mask_rho = np.ones_like(drho, dtype=bool)
            else:
                mask_rho = np.abs(drho) <= max_drho

            if max_dtheta is None:
                mask_theta = np.ones_like(dtheta, dtype=bool)
            else:
                mask_theta = np.abs(dtheta) <= max_dtheta

            mask = mask_rho & mask_theta

            if not np.any(mask):
                out[j, i] = 0.0
                continue

            dr = drho[mask] / sigma_rho
            dt = dtheta[mask] / sigma_theta

            d2 = dr * dr + dt * dt
            w = np.exp(-0.5 * d2)

            wsum = np.sum(w)
            if wsum < min_weight_sum:
                out[j, i] = 0.0
            else:
                out[j, i] = np.sum(w * val_s[mask]) / wsum

    return out


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
R0 = axis[0]
Z0 = axis[1]

print("Magnetic axis =", R0, Z0)

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

good = np.isfinite(rho_e)
Re   = Re[good]
Ze   = Ze[good]
ne_e = ne_e[good]
Te_e = Te_e[good]
rho_e = rho_e[good]

theta_e = np.arctan2(Ze - Z0, Re - R0)

tri = mtri.Triangulation(Re, Ze)

# ============================================================
# Local interpolation in (rho, theta)
# ============================================================

sigma_rho   = 0.0025
sigma_theta = 0.1

max_drho   = 3.0 * sigma_rho
max_dtheta = 3.0 * sigma_theta

ne_e2 = local_interp_rho_theta(
    rho2, theta2,
    rho_e, theta_e, ne_e,
    sigma_rho=sigma_rho,
    sigma_theta=sigma_theta,
    max_drho=max_drho,
    max_dtheta=max_dtheta
)

Te_e2 = local_interp_rho_theta(
    rho2, theta2,
    rho_e, theta_e, Te_e,
    sigma_rho=sigma_rho,
    sigma_theta=sigma_theta,
    max_drho=max_drho,
    max_dtheta=max_dtheta
)

outside_sep = rho2 > 1.0
ne_e2[~outside_sep] = 0.0
Te_e2[~outside_sep] = 0.0

ne_e2 = np.nan_to_num(ne_e2, nan=0.0)
Te_e2 = np.nan_to_num(Te_e2, nan=0.0)

ne_e2 = np.maximum(ne_e2, 0.0)
Te_e2 = np.maximum(Te_e2, 0.0)

# ============================================================
# Merge: core inside, edge outside
# ============================================================

ne_T = ne_c2.copy()
Te_T = Te_c2.copy()

ne_T[outside_sep] = ne_e2[outside_sep]
Te_T[outside_sep] = Te_e2[outside_sep]

# ============================================================
# Diagnostics
# ============================================================

print("edge rho min/max =", np.nanmin(rho_e), np.nanmax(rho_e))
print("edge theta min/max =", np.nanmin(theta_e), np.nanmax(theta_e))
print("edge ne min/max =", np.nanmin(ne_e2), np.nanmax(ne_e2))
print("edge Te min/max =", np.nanmin(Te_e2), np.nanmax(Te_e2))
print("Max edge density inside core:", np.max(ne_e2[rho2 <= 1.0]))
print("Max edge density outside core:", np.max(ne_e2[rho2 > 1.0]))
print("Max edge temperature inside core:", np.max(Te_e2[rho2 <= 1.0]))
print("Max edge temperature outside core:", np.max(Te_e2[rho2 > 1.0]))

# ============================================================
# Plots
# ============================================================

fig, axs = plt.subplots(2, 4, figsize=(15, 10), constrained_layout=True)

# ------------------------------------------------------------
# Row 1: density
# ------------------------------------------------------------

# Raw edge ne
tpc0 = axs[0, 0].tricontourf(tri, np.log10(ne_e + 1.0e14), levels=50)
plot_rho1(axs[0, 0])
fig.colorbar(tpc0, ax=axs[0, 0], label=r"$\log_{10} n_{e,\mathrm{raw}}$ (m$^{-3}$)")
axs[0, 0].set_title("Raw edge density")
axs[0, 0].set_xlabel("R")
axs[0, 0].set_ylabel("Z")
axs[0, 0].set_aspect("equal")

# Interpolated edge ne
pcm0 = axs[0, 1].pcolormesh(RR, ZZ, np.log10(ne_e2 + 1.0e14), shading="auto")
plot_rho1(axs[0, 1])
fig.colorbar(pcm0, ax=axs[0, 1], label=r"$\log_{10} n_{e,\mathrm{interp}}$ (m$^{-3}$)")
axs[0, 1].set_title("Interpolated edge density")
axs[0, 1].set_xlabel("R")
axs[0, 1].set_ylabel("Z")
axs[0, 1].set_aspect("equal")

# Total ne
pcm1 = axs[0, 2].pcolormesh(RR, ZZ, np.log10(ne_T + 1.0e14), shading="auto")
plot_rho1(axs[0, 2])
fig.colorbar(pcm1, ax=axs[0, 2], label=r"$\log_{10} n_e$ (m$^{-3}$)")
axs[0, 2].set_title("Total density")
axs[0, 2].set_xlabel("R")
axs[0, 2].set_ylabel("Z")
axs[0, 2].set_aspect("equal")

# theta
pcm2 = axs[0, 3].pcolormesh(RR, ZZ, theta2, shading="auto")
plot_rho1(axs[0, 3])
fig.colorbar(pcm2, ax=axs[0, 3], label=r"$\theta$")
axs[0, 3].set_title("Poloidal angle")
axs[0, 3].set_xlabel("R")
axs[0, 3].set_ylabel("Z")
axs[0, 3].set_aspect("equal")

# ------------------------------------------------------------
# Row 2: temperature
# ------------------------------------------------------------

# Raw edge Te
tpc1 = axs[1, 0].tricontourf(tri, Te_e, levels=50)
plot_rho1(axs[1, 0])
fig.colorbar(tpc1, ax=axs[1, 0], label=r"$T_{e,\mathrm{raw}}\;(\mathrm{eV})$")
axs[1, 0].set_title("Raw edge temperature")
axs[1, 0].set_xlabel("R")
axs[1, 0].set_ylabel("Z")
axs[1, 0].set_aspect("equal")

# Interpolated edge Te
pcm4 = axs[1, 1].pcolormesh(RR, ZZ, Te_e2, shading="auto")
plot_rho1(axs[1, 1])
fig.colorbar(pcm4, ax=axs[1, 1], label=r"$T_{e,\mathrm{interp}}\;(\mathrm{eV})$")
axs[1, 1].set_title("Interpolated edge temperature")
axs[1, 1].set_xlabel("R")
axs[1, 1].set_ylabel("Z")
axs[1, 1].set_aspect("equal")

# Total Te
pcm5 = axs[1, 2].pcolormesh(RR, ZZ, Te_T, shading="auto")
plot_rho1(axs[1, 2])
fig.colorbar(pcm5, ax=axs[1, 2], label=r"$T_e\;(\mathrm{eV})$")
axs[1, 2].set_title("Total temperature")
axs[1, 2].set_xlabel("R")
axs[1, 2].set_ylabel("Z")
axs[1, 2].set_aspect("equal")

# rho again
pcm7 = axs[1, 3].pcolormesh(RR, ZZ, rho2, shading="auto")
plot_rho1(axs[1, 3])
fig.colorbar(pcm7, ax=axs[1, 3], label=r"$\rho$")
axs[1, 3].set_title("Normalized flux $\\rho$")
axs[1, 3].set_xlabel("R")
axs[1, 3].set_ylabel("Z")
axs[1, 3].set_aspect("equal")

plt.show()