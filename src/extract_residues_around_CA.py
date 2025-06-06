import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

# === Input ===
pdb_file = r'..\temp\similarity\with_collagen\with_collagen.pdb'            # Your PDB file
residue_number = 386               # Target residue number (as in the PDB, 1-based)
distance_cutoff = 17.6            # Radius in Ångstroms
output_file = r'..\temp\similarity\with_collagen\R386_nearby_residues.pdb'  # Output PDB file

# === Load the structure ===
u = mda.Universe(pdb_file)

# === Select the Cα atom of the specified residue ===
target_res = u.select_atoms(f"resid {residue_number} and name CA")
if len(target_res) != 1:
    raise ValueError(f"Could not uniquely identify CA atom of residue {residue_number}")

target_ca = target_res.positions[0]

# === Compute distances from the Cα to all atoms in the protein ===
protein = u.select_atoms("protein")
dists = distance_array(target_ca.reshape(1, 3), protein.positions)[0]

# === Identify residues with at least one atom within the cutoff ===
within_cutoff = protein[dists < distance_cutoff]
unique_residues = within_cutoff.residues

# === Write the nearby residues to a new PDB file ===
with mda.Writer(output_file, multiframe=False) as w:
    w.write(unique_residues.atoms)

print(f"{len(unique_residues)} residues written to {output_file}")