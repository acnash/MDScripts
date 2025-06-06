# MDScripts

Scripts for analysing Molecular Dynamics trajectories and molecular models.

## Spherical residue extraction
All residues within a specified distance (in angstroms) from a specified residue are extracted from an input PDB file and saved to an output residue file. 

### Options
- `input_file` - Path and file name of a clean PDB input file.
- `output_file` - Path and file name of where to save residues.
- `residue_number` - The residue number of the residue in which the carbon alpha is used to derive an inclusion sphere. The residue number can only be present once in the input file, otherwise multiple sites of inclusion will be saved.
- `distance_cutoff` - A radius (angstrom) from the alpha carbon of the specified residue.

## Extract SASA residues
Uses FreeSASA to isolate the solvent accessible surface areas of a protein (from a PDB file) and saves those residues to a PDB file. This will not discriminate chains. Every molecule in the input file is parsed.

A dependency on having FreeSASA installed. This works on my Mac laptop, but not on my Windows machine.

### Options
- `input_file` - Path and file name of a clean PDB input file.
- `output_file` - Path and file name of where to save residues.