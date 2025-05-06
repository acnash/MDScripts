import subprocess

# Path to the executable (adjust the path as needed)
vina_exe = "../external/vina_1.2.5_win.exe"

receptor_file = "../trajectories/MMP1.pdbqt"
ligand_file = "../trajectories/ZINC000001149132.pdbqt"

data_file = "../temp/MMP1_collagen_residues.txt"

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

            args = ["--receptor", "{receptor_file}",
                    "--ligand", "{ligand_file}",
                    "--center_x", "{x_pos}",
                    "--center_y", "{y_pos}",
                    "--center_z", "{z_pos}",
                    "--size_x", "25",
                    "--size_y", "25",
                    "--size_z", "25",
                    "--exhaustiveness", "36",
                    "--out", f"../temp/with_collagen/ZINC11/{resname}_{resid}.pdbqt"]
            result = subprocess.run([vina_exe] + args, capture_output=True, text=True)

            # Print stdout and stderr
            print("STDOUT:\n", result.stdout)
            print("STDERR:\n", result.stderr)
