import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis as HBA
import matplotlib.pyplot as plt

import numpy as np

# mda.analysis.base.AnalysisBase.parallelizable = False

# if you're using a GRO trajectory you still need topology information. Make a subset of the original
# tpr file to match the trajectory:
# gmx convert-tpr -s topol.tpr -n index.ndx -o collagen.tpr
u = mda.Universe(r"../trajectories/collagen.tpr", r"../trajectories/collagen.trr")

hbonds = HBA(universe=u)
hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
hbonds.run()

# Extract H-bond data: shape = (n_bonds, 6), where 6 columns include frame, donor/acceptor/resids etc.
hbonds_array = hbonds.results.hbonds

# Count number of H-bonds per frame
unique_frames, counts = np.unique(hbonds_array[:, 0], return_counts=True)

# Convert frame numbers to actual time using trajectory time step
frame_times = []
for ts in u.trajectory:
    if ts.frame in unique_frames:
        frame_times.append(ts.time)  # In ps

# Plot time vs number of H-bonds
plt.figure(figsize=(10, 6))
plt.plot(frame_times, counts, marker='o', linestyle='-', color='darkblue')
plt.xlabel("Time (ps)", fontsize=12)
plt.ylabel("Number of Hydrogen Bonds", fontsize=12)
#plt.title("Hydrogen Bonds vs Time", fontsize=14)
plt.grid(True)
plt.tight_layout()
plt.show()