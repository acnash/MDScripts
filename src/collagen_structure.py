import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# Load trajectory
u = mda.Universe(r"C:\Users\Anthony\Downloads\collagen.gro")
selection = u.select_atoms("all")
rama = Ramachandran(selection).run()

# Extract phi/psi data
phi_psi = rama.angles.reshape(-1, 2)
phi = phi_psi[:, 0]
psi = phi_psi[:, 1]

# Plot with Seaborn KDE
plt.figure(figsize=(7, 7))
sns.kdeplot(x=phi, y=psi, cmap="viridis", fill=True, thresh=0.01, levels=100)

# Scatter overlay
plt.plot(phi, psi, 'k.', markersize=1, alpha=0.5)

# Labels and formatting
plt.xlim(-180, 180)
plt.ylim(-180, 180)
plt.xlabel("Φ (phi)")
plt.ylabel("Ψ (psi)")
plt.title("Ramachandran Plot with KDE Contour")
plt.grid(True)
plt.show()
