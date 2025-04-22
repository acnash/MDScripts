import MDAnalysis as mda

u = mda.Universe(r"C:\Users\Anthony\Downloads\collagen.tpr", r"C:\Users\Anthony\Downloads\collagen.trr")

for ts in u.trajectory:
    box = ts.dimensions  # returns [lx, ly, lz, alpha, beta, gamma]
    print(f"Time {ts.time} ps - Box: {box[:3]} Å")


