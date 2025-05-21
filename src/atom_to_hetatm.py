import os

# --- User Inputs ---
input_file = "../temp/zinc14_pmf_1000_window00.pdb"  # Replace with your actual filename
target_resname = "LIG"        # Residue name to target

# --- Generate Output File Path ---
base, ext = os.path.splitext(input_file)
output_file = f"{base}_replaced.pdb"

# --- Process the PDB ---
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        if line.startswith("ATOM  ") and target_resname in line:
            line = line.replace("ATOM  ", "HETATM", 1)
        outfile.write(line)

print(f"Finished. Modified file saved as: {output_file}")
