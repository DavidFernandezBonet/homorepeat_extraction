# Homorepeat Extraction

<img src="homorepeats.jpeg" width="200" height="200" />


## Purpose:
Extract aminoacid homorepeats (several aminoacids in a row) from protein sequences using large databases. Furthermore, gather miscellaneous data such as the ratio of homorepeats given the aminoacid type or the secondary structure of the homorepeats. 

## Requirements:
- biopython==1.78
- contourpy==1.0.5
- cycler==0.11.0
- fonttools==4.37.4
- kiwisolver==1.4.4
- matplotlib==3.6.0
- numpy==1.23.3
- packaging==21.3
- Pillow==9.2.0
- pyparsing==3.0.9
- python-dateutil==2.8.2
- six==1.16.0
- DSSP (for e.g. Ubuntu follow these [steps](https://zoomadmin.com/HowToInstall/UbuntuPackage/dssp)

## Usage: 

It is designed to be simple: open the terminal and run a Python script with two arguments:
`python main.py [homorepeat_threshold] [path_to_pdb]`
- **Homorepeat threshold**: The minimum number of aminoacids in a row to be considered a homorepeat. For example, AAAAT contains a "A" homorepeat if the threshold is 4, but not if it is 5.
- **path_to_pdb**: Path to the directory containing .pdb files.


## Output:

- `files_with_homorepeats.txt`: list of pdb files that present homorepeats 
- `data` folder containing:
	- homorepeat information: if an homorepeat is found, annotate aminoacid, homorepeat length, position and secondary structure percentage ([%alpha, %beta, %coil]) --> e.g. GLU 4 714-717 [0.0, 75.0, 25.0] means that a homorepeat of consisting of 4 GLU aminoacids was found at positions 714-715-716-717 with a 75% beta structure and a 25% alpha structure.
	
	- homorepeat analysis:  annotate the ratio of found homorepeats out of the total number of chains. 
	- dihedral data: annotate dihedral angles for each folder and draw a psi-phi plot




