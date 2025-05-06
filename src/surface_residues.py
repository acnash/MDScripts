import MDAnalysis as mda
import freesasa

# Load the structure
u = mda.Universe("../trajectories/cluster_40360_free.pdbqt")

# Select the protein atoms
protein = u.select_atoms("protein")
protein.write("temp_protein.pdb")

# Use freesasa to load and calculate SASA
structure = freesasa.Structure("temp_protein.pdb")
result = freesasa.calc(structure)

# Get per-residue SASA
residue_areas = result.residueAreas()

# Threshold for what we consider "surface"
# < 5 Å²	Very buried (ignore)
# 5–30 Å²	Partially exposed residues
# >30–50+ Å²	Fully exposed residues
surface_threshold = 30.0  # Å²

# Collect surface residues
surface_residues = []
for chain in residue_areas:
    for resnum, resinfo in residue_areas[chain].items():
        total_area = resinfo.total  # attribute
        if total_area > surface_threshold:
            surface_residues.append((chain, resinfo.residueType, resnum))

selection = []
for chain, restype, resnum in surface_residues:
    if chain == "X":
        selection.append(f"resid {resnum}")
    else:
        selection.append(f"segid {chain} and resid {resnum}")

surface_sel = u.select_atoms(" or ".join(selection))

surface_sel.write("surface_residues.pdb")

# Print the surface residues
with open("MMP1_free_residues.txt", "w") as file:
    for chain, resname, resnum in surface_residues:
        ca_atom = u.select_atoms(f"resid {resnum} and name CA")
        ca_pos = ca_atom.positions[0]
        file.write(f"{resname} {resnum} {chain} {ca_pos[0]:.3f} {ca_pos[1]:.3f} {ca_pos[2]:.3f}\n")