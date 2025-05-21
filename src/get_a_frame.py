import MDAnalysis as mda

# Load the topology and trajectory
u = mda.Universe("../trajectories/collagen.tpr", "../trajectories/collagen.trr")

# Jump to the first frame explicitly (optional, but good practice)
u.trajectory[0]

# Write the first frame to a PDB file
with mda.Writer("../temp/first_frame.pdb") as pdb_writer:
    pdb_writer.write(u.atoms)
