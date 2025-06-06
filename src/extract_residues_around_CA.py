import argparse

import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array


def main():
    parser = argparse.ArgumentParser(description="Silent Wealth")

    parser.add_argument("--input", type=str, required=True,
                        help="Input PDB file (include file path). Make sure it's clean..")
    parser.add_argument("--output", type=str, required=True, help="Output PDB file location (include file path).")
    parser.add_argument("--residue_number", type=int,
                        help="Target residue number to extract residues around (extraction includes the target "
                             "residue). Must be unique in the file. Can't have multiple chains.")
    parser.add_argument("--distance_cutoff", type=float, help="Inclusion metric from the carbon alpha of the residue.")

    args = parser.parse_args()
    pdb_file = args.input
    output_file = args.output
    residue_number = args.residue_number
    distance_cutoff = args.distance_cutoff

    # === Input ===
    #pdb_file = r'..\temp\similarity\with_collagen\with_collagen.pdb'            # Your PDB file
    #residue_number = 386               # Target residue number (as in the PDB, 1-based)
    #distance_cutoff = 17.6            # Radius in Ångstroms
    #output_file = r'..\temp\similarity\with_collagen\R386_nearby_residues.pdb'  # Output PDB file

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


if __name__ == "__main__":
    main()
