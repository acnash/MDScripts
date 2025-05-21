import os
import glob
from collections import defaultdict
import numpy as np

# Define the folder containing your files
folder_path = '../trajectories/ZINC_14_1000kj_pullx'
target_value = 0.3998

# Get all text files
file_paths = glob.glob(os.path.join(folder_path, '*.xvg'))

# Dictionary to group second-column values by identifier
data_by_id = defaultdict(list)

# Extract identifier and collect second column values
for file_path in file_paths:
    filename = os.path.basename(file_path)

    # Extract identifier from filename, e.g., "_window27_"
    parts = filename.split('_')
    identifier = '_'.join(part for part in parts if 'window' in part)

    # Read second column, skipping metadata lines
    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('@') or line.startswith('#'):
                continue
            cols = line.split()
            if len(cols) >= 2:
                try:
                    value = float(cols[1])
                    data_by_id[identifier].append(value)
                except ValueError:
                    continue  # Skip non-numeric values

# Calculate average per identifier
averages = {k: np.mean(v) for k, v in data_by_id.items()}

# Find the identifier with the closest average to target_value
closest_id = min(averages, key=lambda k: abs(averages[k] - target_value))
closest_value = averages[closest_id]

# Output all averages
for identifier, avg in sorted(averages.items()):
    print(f'{identifier}: {avg:.4f}')

# Output the closest match
print(f'\nClosest to target ({target_value}): {closest_id} with average {closest_value:.4f}')