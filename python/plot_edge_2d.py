import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

if len(sys.argv) < 2:
    print("Usage: python plot_edge_2d.py <datafile>")
    sys.exit(1)

filename = sys.argv[1]

data = np.loadtxt(filename, comments="#")

k  = data[:, 0]
R  = data[:, 1]
Z  = data[:, 2]
ne = data[:, 3]
Te = data[:, 4]

triang = tri.Triangulation(R, Z)

plt.figure()
plt.tricontourf(triang, ne, levels=50)
plt.colorbar(label="ne")
plt.xlabel("R")
plt.ylabel("Z")
plt.title(f"2D edge electron density ({filename})")
plt.gca().set_aspect("equal")
plt.tight_layout()

plt.figure()
plt.tricontourf(triang, Te, levels=50)
plt.colorbar(label="Te")
plt.xlabel("R")
plt.ylabel("Z")
plt.title(f"2D edge electron temperature ({filename})")
plt.gca().set_aspect("equal")
plt.tight_layout()

plt.show()