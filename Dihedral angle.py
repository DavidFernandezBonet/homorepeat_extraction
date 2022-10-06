from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import DSSP
from Bio.PDB.vectors import Vector, rotmat
from Bio.PDB.Atom import Atom
import numpy as np

def find_dihedral(r1,r2,r3,r4):
    """" Using 4 points calculate the dihedral angle in radians """""

    #Transform into numpy format for easy operations  --> They already are
    r1 = np.array(r1)
    r2 = np.array(r2)
    r3 = np.array(r3)
    r4 = np.array(r4)

    v1 = r2-r1
    v2 = r3-r2
    v3 = r4-r3

    v1v2= np.cross(v1,v2)
    v2v3= np.cross(v2,v3)

    cos_theta = np.dot(v1v2,v2v3)/(np.linalg.norm(v1v2)*np.linalg.norm(v2v3))
    theta = np.arccos(cos_theta)

    # Change to degrees
    #theta = (theta/np.pi)*180


    return theta

r1= [-14.333,-13.53,1.956]
r2= [-14.913,-12.838,0.828]
r3= [ -14.297,-13.087,-0.542]
r4= [-14.469, -14.359,  -0.959]
theta= find_dihedral(r1,r2,r3,r4)

print("la mega teta és:", theta)

path = 'deca_ALA.pdb'
p = PDBParser(PERMISSIVE=True, QUIET=True)
structure = p.get_structure("tmp", path)

# Problem: identify the atoms needed.
# Phi: C-N-CA-C   (first C is i-1)
# Psi: N-CA-C-N   (last N is i+1)

# With "element" we can know simply the atom
# With "get_parent" we get the amino (but as a class)
# With "get_id" we have the atom name with no bullshit

for model in structure:
    for chain in model:
        for residue in chain:
           for atom in residue:

               # Way to know the amino acid position?
               serial = Atom.get_serial_number(atom)
               coord = Atom.get_coord(atom)
               #ALL = Atom.get_full_id(atom)
               altloc = Atom.get_parent(atom)
               print(altloc)
               print(atom,coord,serial,altloc)