import argparse

import MDAnalysis as mda
import freesasa

def main():
    parser = argparse.ArgumentParser(description="Silent Wealth")

    parser.add_argument("--input", type=str, required=True,
                        help="Input PDB file (include file path). Make sure it's clean..")
    parser.add_argument("--output", type=str, required=True, help="Output PDB file location (include file path).")
    parser.add_argument("--surface_threshold", type=int, required=True, help="Threshold in Å² that defines when a residue becomes a surface residue. Set to 30.")
    parser.add_argument("--residue_file", type=str, required=True, help="Output text file location (include file path) of stored residue IDs that make up the surface.")

    args = parser.parse_args()
    pdb_file = args.input
    output_file = args.output
    residue_file = args.residue_file
    surface_threshold = args.surface_threshold

    # Load the structure
    u = mda.Universe(pdb_file)

    # Select the protein atoms
    protein = u.select_atoms("protein")
    protein.write("../temp/temp_protein.pdb")

    # Use freesasa to load and calculate SASA
    structure = freesasa.Structure("../temp/temp_protein.pdb")
    result = freesasa.calc(structure)

    # Get per-residue SASA
    residue_areas = result.residueAreas()

    # Threshold for what we consider "surface"
    # < 5 Å²	Very buried (ignore)
    # 5–30 Å²	Partially exposed residues
    # >30–50+ Å²	Fully exposed residues
    #surface_threshold = 30.0  # Å²

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

    #surface_sel = u.select_atoms(" or ".join(selection))
    surface_sel = u.atoms[[]]
    for sel in selection:
        surface_sel += u.select_atoms(sel)

    surface_sel.write(output_file)

    # Print the surface residues
    with open(residue_file, "w") as file:
        for chain, resname, resnum in surface_residues:
            ca_atom = u.select_atoms(f"resid {resnum} and name CA")
            if ca_atom:
                ca_pos = ca_atom.positions[0]
                file.write(f"{resname} {resnum} {chain} {ca_pos[0]:.3f} {ca_pos[1]:.3f} {ca_pos[2]:.3f}\n")

if __name__ == "__main__":
    main()
