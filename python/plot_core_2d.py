import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) < 2:
    print("Usage: python plot_core_on_eq.py <datafile>")
    sys.exit(1)

filename = sys.argv[1]

data = np.loadtxt(filename, comments="#")

R   = data[:, 0]
Z   = data[:, 1]
rho = data[:, 2]
ne  = data[:, 3]
Te  = data[:, 4]

# infer structured grid sizes from repeated R/Z pattern
R_unique = np.unique(R)
Z_unique = np.unique(Z)

nR = len(R_unique)
nZ = len(Z_unique)

if nR * nZ != len(R):
    raise RuntimeError("Data do not appear to form a structured R-Z grid")

R2   = R.reshape(nR, nZ)
Z2   = Z.reshape(nR, nZ)
rho2 = rho.reshape(nR, nZ)
ne2  = ne.reshape(nR, nZ)
Te2  = Te.reshape(nR, nZ)

plt.figure()
plt.pcolormesh(R2, Z2, ne2, shading="auto")
plt.colorbar(label="ne")
plt.xlabel("R")
plt.ylabel("Z")
plt.title(f"Core electron density on equilibrium grid\n{filename}")
plt.gca().set_aspect("equal")
plt.tight_layout()

plt.figure()
plt.pcolormesh(R2, Z2, Te2, shading="auto")
plt.colorbar(label="Te")
plt.xlabel("R")
plt.ylabel("Z")
plt.title(f"Core electron temperature on equilibrium grid\n{filename}")
plt.gca().set_aspect("equal")
plt.tight_layout()

plt.figure()
plt.pcolormesh(R2, Z2, rho2, shading="auto")
plt.colorbar(label="rho")
plt.xlabel("R")
plt.ylabel("Z")
plt.title(f"Mapped rho on equilibrium grid\n{filename}")
plt.gca().set_aspect("equal")
plt.tight_layout()

plt.show()