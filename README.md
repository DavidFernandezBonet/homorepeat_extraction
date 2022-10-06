# Homorepeat Analysis

Requirements:
- Ubuntu OS (because of directories "/" as opposed to "\\"). TODO: fix this, should be universal
- DSSP

Input: folders containing .pdb files
For each folder, expected result is 
- files_containing_homorepeats: list of pdb files that present homorepeats (cut-off at n=4, or 4 aminoacids in a row)
- homorepeat information: if an homorepeat is found, annotate aminoacid, homorepeat length, position and secondary structure percentage ([%alpha, %beta, %coil]) --> e.g. GLU 4 714-717 [0.0, 75.0, 25.0]
- homorepeat analysis: count % of homorepeats out of the total chains. Also count % of homorepeat aminoacid type
- dihedral data: compute dihedral data for each input folder and plot results 

Procedure:
1. Generate data with "dihedral and homo.py"  --> homorepeat information, dihedral data
2. Analyze percentage of homorepeats with "cound homorepeats kind.py"
3. Plot dihedral data with "Read dihedral data and plot.py"

