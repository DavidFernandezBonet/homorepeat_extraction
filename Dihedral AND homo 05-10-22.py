import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from itertools import groupby
import csv
from Bio.PDB import DSSP
import os
import time

### BIOPYTHON V 1.78
### ATTENTION!!! DSSP NEEDS TO BE INSTALLED https://zoomadmin.com/HowToInstall/UbuntuPackage/dssp
### MUST BE RUN IN UBUNTU FOR DIRECTORIES : "\", TODO: solve for windows with "//"
# FIND HOMOREPEATS
# ----------------------------------------
def annotate_dictionary(my_list):


    """"" 
    Creates a dictionary from a list. It indicates the element of the list (key)
    and how many times it repeats in a row.
    """""
    # TO DO: add also the place where it starts.
    # https://stackoverflow.com/questions/773/how-do-i-use-itertools-groupby maybe unsterstand it better with this?

    #No sé si és possible perquè el bucle passa 99 vegades en una list de llargada 111.
    #Diria que fa skip a coses repetides
    #Counter is not good enough

    d = dict()

    #Check if the list is empty
    #If it is not, the function begins
    if not my_list:
        print("Empty sequence")

    else:


        #Creates a counter_list, that is, the times a certain aminoacid repeats
        #It is useful to later determine the position of the aminoacids
        counter_list=[]
        i=0
        initial_position= my_list[0][0] #first element of the list with index 1 is position
        counter= initial_position

        # determine final position (as there are some lists that are continuously repeating)
        positions=[]
        for o in range(len(my_list)):
            positions.append(my_list[o][0])
        final_position= max(positions)

        #print("Initial Position", initial_position, ", Final Position", final_position, ", Length:", final_position-initial_position+1)


        #Create a list with only the aminoacid sequence
        list_sequence= []
        for j in range(len(my_list)):
            list_sequence.append(my_list[j][1])  #append aminoacid



        for k,v in groupby(list_sequence):
            counter_list.append(len(list(v)))


        #Creates a dictionary with 'AMINOACID': (number of repeats, starting position of the repeat)
        for k, v in groupby(list_sequence):

            #setdefault returns the value of a certain key
            #The value is converted in a list

            d.setdefault(k, []).append( (len(list(v)), counter ))
            counter += counter_list[i]
            if counter >= final_position:
                counter = initial_position
            i += 1

    return d
    #It creates a dictionary: key--> the element of the list provided
    #                         value--> list with the times it repeats in a row
def unique(list1):
    # intilize a null list
    unique_list = []

    # traverse for all elements
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)

    return unique_list
def aminoacid_abreviation(letter):
    """""
    Translates aminoacids abreviation (one letter) into the usual three leters
    For example "G" --> "GLY"
    """""


    if letter == 'A':
        letter = 'ALA'

    elif letter == 'R':
        letter = 'ARG'

    elif letter == 'N':
        letter = 'ASN'

    elif letter == 'D':
        letter = 'ASP'

    elif letter == 'C':
        letter = 'CYS'

    elif letter == 'Q':
        letter = 'GLN'

    elif letter == 'E':
        letter = 'GLU'

    elif letter == 'G':
        letter = 'GLY'

    elif letter == 'H':
        letter = 'HIS'

    elif letter == 'I':
        letter = 'ILE'

    elif letter == 'L':
        letter = 'LEU'

    elif letter == 'K':
        letter = 'LYS'

    elif letter == 'M':
        letter = 'MET'

    elif letter == 'F':
        letter = 'PHE'

    elif letter == 'P':
        letter = 'PRO'

    elif letter == 'O':
        letter = 'PYL'

    elif letter == 'S':
        letter = 'SER'

    elif letter == 'U':
        letter = 'SEC'

    elif letter == 'T':
        letter = 'THR'

    elif letter == 'W':
        letter = 'TRP'

    elif letter == 'Y':
        letter = 'TYR'

    elif letter == 'V':
        letter = 'VAL'

    elif letter == 'B':
        letter = 'ASX'

    elif letter == 'Z':
        letter = 'GLX'

    elif letter == 'X':
        letter = 'XAA'

    elif letter == 'J':
        letter = 'XLE'



    else:
        print("caca de la vaca, alguna no s'ha detectat")

    return letter
def structure_abbreviation(letter):
    if (letter == 'H') or (letter == 'H') or (letter == 'H'):
        letter= "HELIX"
    elif (letter== 'E'):
        letter= "SHEET"
    else:
        letter= "COIL"
    return letter
def are_there_homorepeats(d,homorepeat_min,pos_struc_X,final_result):
    #Tells if there are homorepeats, which aminoacid and the length and the position
    #The convention is the first aminoacid has the postiion number 1 (not 0)
    homorepeats = False

    #pos_struc_X: [(757, 'COIL'), ..., ]
    #the new one has (757, aminoacid, structure) !!!

    h_positions = []
    for item in d.items():

        for i in range(len(item[1])): #length of te values
            #item[0] --> key
            #item[1] --> tuples of the form: (number of repeats, position)
            #item[1][i][0] --> number of repeats
            if item[1][i][0] >= homorepeat_min:
                aminoacid= item[0]
                length= item[1][i][0]
                start_pos= item[1][i][1]
                final_pos = start_pos + length -1

                # Fill the h_positions list
                for iteration in range(length):
                    h_positions.append((aminoacid,start_pos+iteration)) # (append aminoacid, position)

                helix_count=0
                sheet_count=0
                coil_count=0

                # ------------------------------
                #Structure determination
                for p in range(length):
                    #locating the index of the starting amino acid in the pos_struc_X list
                    index= (start_pos+p)-pos_struc_X[0][0]
                    structure = pos_struc_X[index][2]
                    if structure == 'HELIX':
                        helix_count += 1
                    if structure == 'SHEET':
                        sheet_count += 1
                    if structure == 'COIL':
                        coil_count +=1

                percentage_helix= (helix_count/length)*100
                percentage_sheet= (sheet_count/length)*100
                percentage_coil= (coil_count/length)*100
                # ------------------------------

                # print("Homorepeat", aminoacid, ". Length:", length, ". Position:",
                #       start_pos, "to", final_pos, ". Its structure is:", percentage_helix, "% HELIX, ",
                #       percentage_sheet, "% SHEET, ", percentage_coil, "% COIL.")

                with open(final_result, 'a', newline='') as f:
                    writer = csv.writer(f, delimiter='\t')
                    writer.writerow(["","","",aminoacid, length, str(start_pos) + "-"+ str(final_pos), [percentage_helix,
                                     percentage_sheet, percentage_coil]])


                homorepeats= True




    if homorepeats == False:
        with open(final_result, 'a', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["","","","There are no homorepeats"])


    return h_positions, homorepeats
def find_structure(title, directory):
    # Returns a list with (position, aminoacid, structure)
    path = directory + "/" + title
    multiple_chains = True  # ACTIVAR SI NOMÉS VULL LA PRIMERA CADENA
    all = []
    with open(path, "r") as f:
        # with warnings.catch_warnings():
        #     warnings.simplefilter("error")  # warnings become exceptions
        #     try:
        #         p = PDBParser(PERMISSIVE=True, QUIET=True)
        #         structure = p.get_structure("tmp", path)
        #     except:
        #     #except PDBConstructionWarning:
        #         pass
                # The Model in the PDB file consists of multiple Chains.
                # Calling DSSP with the current Model will concatenate
                # all Chains, and this is unavoidable because it's
                # external to Python.  The DSSP result will need to be
                # truncated after the fact.
                #warnings.simplefilter("ignore")
                #structure = PDBParser().get_structure("tmp", path)
                #multiple_chains = True
        #model = next(structure.get_models())  # I always work with Model 1

        try:
            p = PDBParser(PERMISSIVE=True, QUIET=True)
            structure = p.get_structure("tmp", path)
            model = next(structure.get_models())  # I always work with Model 1
            output = DSSP(model, path)  # DSSP expects a Model, why?


            if multiple_chains:
                # A DSSP object is dictionary-like, but with ordered keys.
                # The keys are nested tuples which are structured like PDB
                # Residue identifiers.  The first element of the tuple is
                # the Chain name (typically, "A").

                # Appending all unique labels in the label_list
                k0_list = []
                start_pos_list = []  # List in which we store the first start pos of a chain
                for k in output.keys():
                    k0_list.append(k[0])
                label_list = unique(k0_list)


                for i in range(len(label_list)):
                    start_pos_logic = False

                    for j in range(len(output.keys())):
                        if output.keys()[j][0] == label_list[i]:

                            if start_pos_logic == False:
                                start_pos = output.keys()[j][1][1]
                                start_pos_list.append(start_pos)
                                start_pos_logic = True



                # subset_keys = [k for k in output.keys() if k[0] == chain_a]     #Chain B happens if k[0]="B"

                subset_keys_list = []  # list to store the different sequences. It will support "A", "B", "C", "D"

                for i in range(len(label_list)):
                    subset_keys = [k for k in output.keys() if k[0] == label_list[i]]
                    output_append = [output[k] for k in subset_keys]
                    subset_keys_list.append(output_append)

            all = []  # num, amino, dssp, sequence A,B,C,D

            for i in range(len(subset_keys_list)):  # same length as label list (length = number of chains)
                seq_label = label_list[i]
                start_pos = start_pos_list[i]
                for (n, a, d, *unused) in subset_keys_list[i]:
                    all.append([start_pos, a, d, seq_label])
                    start_pos += 1
        except:
            print("ADEU ADEU")
            pass


    return all
def correct_sequence(sequence):
    #From a list that has (position,aminoacid,structure)
    # correct the position, the aminoacid 3 letter, the structure helix/sheet/coil
    for i in range(len(sequence)):

        sequence[i][1]= aminoacid_abreviation(sequence[i][1])    #aminoacid 3 letters
        sequence[i][2]= structure_abbreviation(sequence[i][2])   #correct structure (only 3 kinds)

    return sequence
def final_message_structure(homorepeat_min,dct, final_result, directory):
    h_positions_list = []  # list for the homo positions for EACH chain
    homorepeats_in_chain = []
    for key in dct.keys():         # The keys are the labels of the chain!
        with open(final_result, 'a',newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["CHAIN", str(key)])

        sequence = dct.get(key)
        d_A = annotate_dictionary(sequence)
        h_positions, homorepeats = are_there_homorepeats(d_A, homorepeat_min, sequence,final_result)  #list-tupple [(amino, position),...,(amino,position)]
        h_positions_list.append(h_positions)    # list of lists of tuples
        homorepeats_in_chain.append(homorepeats)
    # WE GET THE PHI_PSI
    get_homo_phipsi(title,h_positions_list,dihedral_t, directory)
    return homorepeats_in_chain

def read_file_find_homorepeats(title,homorepeat_min,directory,final_result):


    # Retrieve info about sequences (position, amino, structure, seq_label)
    all_sequences= find_structure(title, directory)
    all_sequences = correct_sequence(all_sequences)

    # Separate them by label (max 4 sequences)

    k0_list = []
    for i in range(len(all_sequences)):
        k0_list.append(all_sequences[i][3])
    label_list = unique(k0_list)

    dct = {}
    for i in label_list:
        dct[i] = []               # Dictionary with the labels

    # And now we append to the dictionary
    for u in range(len(label_list)):
        for i in range(len(all_sequences)):
            position= all_sequences[i][0]
            amino = all_sequences[i][1]
            dssp = all_sequences[i][2]
            seq_label = all_sequences[i][3]
            if all_sequences[i][3] == label_list[u]:
                dct[label_list[u]].append([position,amino,dssp])
    homorepeats_in_chain = final_message_structure(homorepeat_min,dct, final_result, directory)
    return homorepeats_in_chain
def read_all_files(homorepeat_min,directory,final_result):
    with open(final_result, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Homorepeat", "Length", "Position", "%HELIX", "%SHEET", "%COIL"])
    # directory = os.path.join(mydir, myfile)
    counter = 0
    for filename in os.listdir(directory):


        if filename.endswith(".pdb"):    #end with .pdb or .ent
            counter += 1
            print(counter)
            title = filename
            #print(filename)
            with open(final_result, 'a', newline="") as f:
                writer = csv.writer(f, delimiter='\t')
                writer.writerow([]) #empty row
                writer.writerow(["-------------------------------------------------------------------------"])
                writer.writerow([counter, filename])

            homorepeats_chain = read_file_find_homorepeats(title, homorepeat_min, directory,final_result)



            continue
        else:
            continue

#-----------------------------------------

# FIND PHI_PSI
def phi_psi_degrees(phi_psi_index):
    """
    Converts (phi,psi) from radians to degrees)
    """
    phi_psi_index = list(phi_psi_index)
    for i in range (len(phi_psi_index)):
        if phi_psi_index[i] == None:
            continue
        else:

            phi_psi_index[i] = phi_psi_index[i]*(180/np.pi)

    return phi_psi_index



def get_homo_phipsi(title, h_positions_list, dihedral_t, directory):
    p = PDBParser(PERMISSIVE=True, QUIET=True)
    path = directory + "/" + title
    structure = p.get_structure("tmp", path)
    counter_chain = 0
    print(title)
    for model in structure:
        for chain in model:
            if not h_positions_list:   # if there is no homorepeat continue
                continue
            elif counter_chain >= len(h_positions_list):   # avoid getting out of index
                continue
            else:

                h_positions = h_positions_list[counter_chain]   # list containing the tuples
                counter_chain += 1
                poly = Bio.PDB.Polypeptide.Polypeptide(chain)
                #print ("Model %s Chain %s" % (str(model.id), str(chain.id)))
                print("Chain %s" % (str(chain.id)))
                phi_psi = poly.get_phi_psi_list()
                for res_index, residue in enumerate(poly):
                    res_name = "%s, %i," % (residue.resname, residue.id[1])
                    name_only = residue.resname   # CAREFUL, I THINK IT DOES NOT GET THE XAA and it classifies it as unknown.


                    for iter in range(len(h_positions)):
                        if residue.id[1] == h_positions[iter][1]:   # the positions are in index 1
                            phi_psi_index = phi_psi_degrees(
                                phi_psi[res_index])  # get the phi,psi of this particular index
                            with open(dihedral_t, 'a', newline='') as f:
                                writer = csv.writer(f, delimiter='\t')
                                writer.writerow([h_positions[iter][0], phi_psi_index])
                            print(h_positions[iter][0], phi_psi_index)

                    # if residue.id[1] in h_positions:  # if the position is in the h_positions
                    #
                    #     phi_psi_index = phi_psi_degrees(phi_psi[res_index])    # get the phi,psi of this particular index
                    #     with open(dihedral_t, 'a', newline='') as f:
                    #         writer = csv.writer(f, delimiter='\t')
                    #         writer.writerow([name_only, phi_psi_index])
                    #     print(res_name, phi_psi_index)
    return


def create_folder_if_not_exist(parent_dir, folder_name):
    path = parent_dir + "/" + folder_name
    isExist = os.path.exists(path)
    # Check whether the specified path exists or not
    if not isExist:
        # Create a new directory because it does not exist
        os.makedirs(path)
        print(f"The new directory {path} is created!")


rootdir = r"/home/david/Desktop/swissprot_pdb_v2_grouped"  # Directory with PDB files
code_dir = os.getcwd()



# #run your code
homorepeat_min= 4 #Threshold for the homorepeat


#directory=  r"C:\Users\David\PycharmProjects\Internships\Internship Dilraj\AF-test new database. 08-2022\AF-test"
#title = 'deca_ALA.pdb'



counter = 0
# 1. Gets structure 2. Orders it in comprehensive list. 3. Summons final_message_structure
# final_message_structure: 1. Writes chain A,B... 2. Summons annotate_dict and are_there_homorepeats


# # Don't remember how I obtained homo_files... I think this is done for efficiency, only reading the files that really contain homorepeats
# homo_files_list = []        # list of each pdb file with homorepeat
# with open("Homo_files.txt", 'r') as f:
#     reader = csv.reader(f, delimiter='\t')
#
#     for row in reader:
#         homo_files_list.append(row[1])
#
#
# with open(dihedral_t, 'w', newline='') as f:
#     writer = csv.writer(f, delimiter='\t')
#     writer.writerow(["Amino Acid, [Phi,Psi]"])


# # Read all files
# for filename in os.listdir(directory):
#
#     #if filename.endswith(".pdb"):  # end with .pdb or .ent
#     #if filename in homo_files_list:                   # DEACTIVATED: UPDATE 2022-08-08. IDK why there was a homo_files_list
#
#     counter += 1
#     print(counter)
#     title = filename
#     #print(filename)
#     with open(final_result, 'a', newline="") as f:
#         writer = csv.writer(f, delimiter='\t')
#         read_file_find_homorepeats(title, homorepeat_min, directory, final_result)


homorepeat_in_pdb_file = f'{code_dir}/data/homorepeats'

with open(homorepeat_in_pdb_file, 'a', newline="") as h:
    writer2 = csv.writer(h, delimiter='\t')
    writer2.writerow(["File Name", "Chain position of homorepeat"])

# Read all files
time_start = time.perf_counter()
for subdir, dirs, files in os.walk(rootdir):
    folder_name = subdir[-4:]
    # create_folder_if_not_exist(code_dir, folder_name)  # Do not need to create folders, for each data folder there is one file
    dihedral_t = f'{code_dir}/data/dihedral_data_{folder_name}.txt'
    final_result = f"{code_dir}/data/{folder_name}.txt"

    with open(dihedral_t, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["Amino Acid, [Phi,Psi]"])

    for file in files:
        counter += 1
        print(counter)


        with open(final_result, 'a', newline="") as f:
            title = str(file)
            writer = csv.writer(f, delimiter='\t')

            homorepeats_chain = read_file_find_homorepeats(title, homorepeat_min, subdir, final_result)
        if True in homorepeats_chain:
            with open(homorepeat_in_pdb_file, 'a', newline="") as h:
                writer2 = csv.writer(h, delimiter='\t')
                writer2.writerow([file, [i for i, x in enumerate(homorepeats_chain) if x]])

        print(os.path.join(subdir, file))

    time_elapsed = (time.perf_counter() - time_start)
    print("time_elapsed", time_elapsed, "current_folder:", folder_name)


time_elapsed = (time.perf_counter() - time_start)
print("time_elapsed", time_elapsed)

 # CAREFUL, I THINK IT DOES NOT GET THE XAA and it classifies it as unknown.
