import numpy as np
import MDAnalysis as mda
from MDAnalysis.coordinates.base import Timestep

# Step 1: Load the two structures
uA = mda.Universe(r"C:\Users\Anthony\sensaas\MMP1_COLLAGEN\pocket4_atm.pdb")
uB = mda.Universe(r"C:\Users\Anthony\sensaas\MMP1_COLLAGEN\pocket2_atm.pdb")

# Step 2: Read the 4x4 transformation matrix
transform_matrix = np.loadtxt(r"C:\Users\Anthony\sensaas\tran_4_2.txt").reshape(4, 4)

# Extract the 3x3 rotation and 3x1 translation components
rotation = transform_matrix[:3, :3]
translation = transform_matrix[:3, 3]

# Step 3: Apply transformation to FileA.pdb coordinates
atomsA = uA.atoms
coords = atomsA.positions  # shape: (n_atoms, 3)

# Apply rotation and translation
transformed_coords = np.dot(coords, rotation.T) + translation

# Store transformed coordinates in the universe
atomsA.positions = transformed_coords

# Optional: Write out transformed structure
uA.atoms.write(r"..\temp\FileA_transformed.pdb")

# Now FileA_transformed.pdb is superimposed onto FileB.pdb using the provided matrix