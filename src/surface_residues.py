import MDAnalysis as mda
import freesasa

# Load the structure
u = mda.Universe("../trajectories/4auo_single_no_collagen.pdb")

# Select the protein atoms
protein = u.select_atoms("protein")
protein.write("temp_protein.pdb")

# Use freesasa to load and calculate SASA
structure = freesasa.Structure("temp_protein.pdb")
result = freesasa.calc(structure)

# Get per-residue SASA
residue_areas = result.residueAreas()

# Threshold for what we consider "surface"
surface_threshold = 30.0  # Å²

# Collect surface residues
surface_residues = []
for chain in residue_areas:
    for resnum, resinfo in residue_areas[chain].items():
        total_area = resinfo.total  # attribute
        if total_area > surface_threshold:
            surface_residues.append((chain, resinfo.residueType, resnum))

# Print the surface residues
for chain, resname, resnum in surface_residues:
    print(f"{resname} {resnum} {chain}")