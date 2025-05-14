import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Load the trajectory
u = mda.Universe("../trajectories/collagen.tpr", "../trajectories/collagen.trr")

# Prepare arrays to store data
n_frames = len(u.trajectory)
box_lengths = np.zeros((n_frames, 3))  # columns: lx, ly, lz
times = np.zeros(n_frames)

# Loop over trajectory
for i, ts in enumerate(u.trajectory):
    box_lengths[i] = ts.dimensions[:3]  # lx, ly, lz in Ångström
    times[i] = ts.time  # time in ps

# Plot box lengths over time
plt.figure(figsize=(10, 6))
plt.plot(times, box_lengths[:, 0], label='Lx')
plt.plot(times, box_lengths[:, 1], label='Ly')
plt.plot(times, box_lengths[:, 2], label='Lz')
plt.xlabel("Time (ps)")
plt.ylabel("Box Length (Å)")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("unit_cell_dimensions.jpg", dpi=600)
plt.show()
