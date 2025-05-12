import subprocess
import os

# Path to the executable (adjust the path as needed)
vina_exe = "../external/vina_1.2.5_win.exe"

receptor_file = "../trajectories/MMP1.pdbqt"
#ligand_file = "../trajectories/ZINC000001149132.pdbqt"
ligand_directory = "../trajectories/top_50_ARG405_collagen/"
output_directory = "../temp/with_collagen/"
data_file = "../temp/MMP1_collagen_residues.txt"

# go through each file in the ligand_directory
# create a directory in temp/with_collagen

for pdbqt_filename in os.listdir(ligand_directory):
    pdbqt_full_path = os.path.join(ligand_directory, pdbqt_filename)
    zinc_name = os.path.splitext(pdbqt_filename)[0]
    zinc_output_folder = os.path.join(output_directory, zinc_name)
    os.makedirs(zinc_output_folder, exist_ok=True)

    with open(data_file, "r") as file:

        for line in file:
            if line.strip():
                #PHE 81 A 16.868 132.433 28.508
                parts = line.strip().split()
                resname = parts[0]
                resid = parts[1]
                chain = parts[2]
                x_pos = parts[3]
                y_pos = parts[4]
                z_pos = parts[5]

                args = ["--receptor", f"{receptor_file}",
                        "--ligand", f"{pdbqt_full_path}",
                        "--center_x", f"{x_pos}",
                        "--center_y", f"{y_pos}",
                        "--center_z", f"{z_pos}",
                        "--size_x", "25",
                        "--size_y", "25",
                        "--size_z", "25",
                        "--exhaustiveness", "36",
                        "--out", f"{zinc_output_folder}/{resid}_{resname}.pdbqt"]
                print(" ".join(args))
                result = subprocess.run([vina_exe] + args, capture_output=True, text=True)

                # Print stdout and stderr
                print("STDOUT:\n", result.stdout)
                print("STDERR:\n", result.stderr)
