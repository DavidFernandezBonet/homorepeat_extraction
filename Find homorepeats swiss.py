from itertools import groupby
from itertools import chain
import csv

import re
import ast  #for literal eval
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import DSSP
from Bio.PDB.PDBExceptions import PDBConstructionWarning
import warnings
import os
import time


time_start = time.perf_counter()

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

        i = 0
        counter= 0


        #print("Initial Position", initial_position, ", Final Position", final_position, ", Length:", final_position-initial_position+1)





        for k,v in groupby(my_list):
            counter_list.append(len(list(v)))


        #Creates a dictionary with 'AMINOACID': (number of repeats, starting position of the repeat)
        for k, v in groupby(my_list):

            #setdefault returns the value of a certain key
            #The value is converted in a list

            d.setdefault(k, []).append( (len(list(v)), counter ))
            counter += counter_list[i]

            i += 1

    return d
    #It creates a dictionary: key--> the element of the list provided
    #                         value--> list with the times it repeats in a row


# The other part of the problem is knowing how to read the pbd files
# The keyword "SEQRES" starts the sequence of aminoacids
# Then the keyword "HELIX" or "SHEET" determine the structure
# This is for reading multiple files: https://towardsdatascience.com/the-best-practice-of-reading-text-files-in-python-509b1d4f5a4


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





def are_there_homorepeats(d,homorepeat_min,final_result):
    #Tells if there are homorepeats, which aminoacid and the length and the position
    #The convention is the first aminoacid has the postiion number 1 (not 0)

    homorepeats = False

    #pos_struc_X: [(757, 'COIL'), ..., ]
    #the new one has (757, aminoacid, structure) !!!

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

                with open(final_result, 'a', newline='') as file:
                    writer = csv.writer(file, delimiter='\t')
                    writer.writerow([aminoacid, length, str(start_pos) + "-"+ str(final_pos)])


                homorepeats= True




    if homorepeats == False:
        with open(final_result, 'a', newline='') as file:
            writer = csv.writer(file, delimiter='\t')
            writer.writerow(["There are no homorepeats"])






    return


def read_fasta_file(homorepeat_min, title, final_result, directory):
    path = directory + "\\" + title


    with open(path, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        counter = 0
        for row in reader:
            #seguent = next(reader)[0]

            if counter ==0:    # First row
                if row[0][0] == ">":
                    counter += 1
                    amino_list = []
                    cut= row[0].find("|", 4)+1   # Cuts the proper protein name
                    protein_name = row[0][0:cut]


                    with open(final_result, 'a', newline="") as file:

                        writer = csv.writer(file, delimiter='\t')
                        writer.writerow([protein_name])  # maybe it needs row[0]

                else:

                    amino_list.append(list(row[0]))

                    if len(row[0]) != 60:     # final row    (sometimes if the final row has exactly 60 it does not work)

                        amino_list = list(chain.from_iterable(amino_list))# unnest list [[1,2]] --> [1,2]

                        for i in range(len(amino_list)):
                            amino_list[i] = aminoacid_abreviation(amino_list[i])  # Transform abreviations to normal amino

                        d = annotate_dictionary(amino_list)
                        are_there_homorepeats(d,homorepeat_min,final_result)


            if counter > 0:
                if row[0][0] == ">":
                    amino_list = list(chain.from_iterable(amino_list))  # unnest list [[1,2]] --> [1,2]

                    for i in range(len(amino_list)):
                        amino_list[i] = aminoacid_abreviation(amino_list[i])  # Transform abreviations to normal amino

                    d = annotate_dictionary(amino_list)
                    are_there_homorepeats(d, homorepeat_min, final_result)

                    counter += 1
                    amino_list = []
                    cut= row[0].find("|", 4)+1   # Cuts the proper protein name
                    protein_name = row[0][0:cut]


                    with open(final_result, 'a', newline="") as file:

                        writer = csv.writer(file, delimiter='\t')
                        writer.writerow([protein_name])  # maybe it needs row[0]

                else:

                    amino_list.append(list(row[0]))







            if counter % 1000 == 0:
                print(counter)


    return




#run your code

title = "1629465504.fasta"
directory = r"C:\Users\David\Desktop\FastaSP2\Swissprot_1"
final_result = "uniprot_sprot_2_analysis.txt"
homorepeat_min= 5 #Threshold for the homorepeat
#read_all_files(homorepeat_min,directory,final_result)

with open(final_result, 'w', newline="") as file:
    writer = csv.writer(file, delimiter='\t')

    writer.writerow(["FASTA ANALYSIS"])
read_fasta_file(homorepeat_min, title, final_result,directory)

time_elapsed = (time.perf_counter() - time_start)


print("time_elapsed", time_elapsed)


# Limits: only handles up to 4 chains and the last three have to be named ("B", "C", "D")