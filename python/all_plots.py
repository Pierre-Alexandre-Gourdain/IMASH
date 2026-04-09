import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

# ============================================================
# Load core data
# ============================================================
core = np.loadtxt("core_on_equilibrium_grid.dat", comments="#")
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
# Load edge data
# ============================================================
edge = np.loadtxt("edge_cells_2d.dat", comments="#")
Re   = edge[:, 1]
Ze   = edge[:, 2]
ne_e = edge[:, 3]
Te_e = edge[:, 4]

tri = Triangulation(Re, Ze)

# ============================================================
# Load magnetic field data
# ============================================================
bf   = np.loadtxt("bfield.dat", comments="#")
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

RRb, ZZb = np.meshgrid(Rvals_b, Zvals_b, indexing="xy")

# ============================================================
# Plot the separatrix
# ============================================================

def plot_rho1(ax):
    ax.contour(RR, ZZ, rho2, levels=[1.0],
               colors='white', linewidths=2.5)
    ax.contour(RR, ZZ, rho2, levels=[1.0],
               colors='black', linewidths=1)

# ============================================================
# Make combined figure
# ============================================================
fig, axs = plt.subplots(2, 3, figsize=(16, 10), constrained_layout=True)

# ------------------------------------------------------------
# Core ne
# ------------------------------------------------------------
pcm0 = axs[0, 0].pcolormesh(RR, ZZ, ne_c2, shading="auto")
plot_rho1(axs[0, 0])
fig.colorbar(pcm0, ax=axs[0, 0], label="$n_e ( m^{-3} )$")
axs[0, 0].set_title("Core density")
axs[0, 0].set_xlabel("R")
axs[0, 0].set_ylabel("Z")
axs[0, 0].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Core Te
# ------------------------------------------------------------
pcm1 = axs[0, 1].pcolormesh(RR, ZZ, Te_c2*1e-3, shading="auto")
plot_rho1(axs[0, 1])
fig.colorbar(pcm1, ax=axs[0, 1], label="$T_e (keV)$")
axs[0, 1].set_title("Core temperature")
axs[0, 1].set_xlabel("R")
axs[0, 1].set_ylabel("Z")
axs[0, 1].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Magnetic field magnitude
# ------------------------------------------------------------
pcm2 = axs[0, 2].pcolormesh(RRb, ZZb, Bpol2, shading="auto")
plot_rho1(axs[0, 2])
fig.colorbar(pcm2, ax=axs[0, 2], label="$B_{pol} (T)$")
axs[0, 2].set_title("Poloidal magnetic field magnitude")
axs[0, 2].set_xlabel("R")
axs[0, 2].set_ylabel("Z")
axs[0, 2].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Edge ne
# ------------------------------------------------------------
tpc0 = axs[1, 0].tricontourf(tri, ne_e, levels=50)
plot_rho1(axs[1, 0])
fig.colorbar(tpc0, ax=axs[1, 0], label="$n_e ( m^{-3} )$")
axs[1, 0].set_title("Edge density")
axs[1, 0].set_xlabel("R")
axs[1, 0].set_ylabel("Z")
axs[1, 0].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Edge Te
# ------------------------------------------------------------
tpc1 = axs[1, 1].tricontourf(tri, Te_e, levels=50)
plot_rho1(axs[1, 1])
fig.colorbar(tpc1, ax=axs[1, 1], label="$T_e (eV)$")
axs[1, 1].set_title("Edge temperature")
axs[1, 1].set_xlabel("R")
axs[1, 1].set_ylabel("Z")
axs[1, 1].set_aspect("equal", adjustable="box")

# ------------------------------------------------------------
# Poloidal field lines
# ------------------------------------------------------------

# Bottom-right: poloidal field lines
axs[1, 2].streamplot(Rvals_b, Zvals_b, BR2, BZ2, density=1.2)
plot_rho1(axs[1, 2])
axs[1, 2].set_title("Poloidal magnetic field lines")
axs[1, 2].set_xlabel("R")
axs[1, 2].set_ylabel("Z")
axs[1, 2].set_aspect("equal", adjustable="box")

# Add a dummy colorbar only to keep panel widths identical
sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0, vmax=1), cmap="Greys")
sm.set_array([])
cbar = fig.colorbar(sm, ax=axs[1, 2])
cbar.set_label("")
cbar.ax.set_yticklabels([])

# Force same plotting window on all panels
Rmin, Rmax = Rvals_b.min(), Rvals_b.max()
Zmin, Zmax = Zvals_b.min(), Zvals_b.max()

for ax in axs.flat:
    ax.set_xlim(Rmin, Rmax)
    ax.set_ylim(Zmin, Zmax)
    ax.set_aspect("equal", adjustable="box")

plt.show()