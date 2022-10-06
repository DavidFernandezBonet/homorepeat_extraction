import csv
import os

amino_list=["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
            "SER","THR","TRP","TYR","VAL","XAA"]







def initialize_dictionary(amino_list):
    dct_homo_counter= {}
    for i in amino_list:
        dct_homo_counter[i] = 0      # Establishing keys and setting value as 0
    return dct_homo_counter




# fer diccionari per contar homorepeats (key: amino, value: counter)
def count_homorepeats(data_dir, file_name, dct_homo_counter):


    file_dir = f"{data_dir}/Homorepeat information/{file_name}"
    chain_counter = 0
    homo_counter = 0
    with open(file_dir, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            #print(len(row))
            if "CHAIN" in row:
                chain_counter += 1
            if len(row) == 7:   # data row (length = 7)           # FOR FASTA FILES IT IS MORE THAN 7 (PROB 10)
                homo_counter += 1
                if row[3] in dct_homo_counter:                       # FOR FASTA FILES IT IS 6
                    dct_homo_counter[row[3]] += 1    # normally 3    # FOR FASTA FILES IT IS 6
                else:
                    dct_homo_counter["OTHERS"] += 1


    print(chain_counter)
    dct_homo_average = dct_homo_counter.copy()

    for key in dct_homo_average.keys():
        dct_homo_average[key] = (dct_homo_average[key] / homo_counter)*100  # returns the percentage
    fraction_homo= (homo_counter/chain_counter)*100

    print("There are", fraction_homo, "% of homorepeats out of a total of", chain_counter, "chains")
    print("The % of homorepeats for each kind are:")
    print(dct_homo_average)


    doc_name= f"{data_dir}/Homorepeat analysis/homorepeat_analysis_{file}.csv"
    with open(doc_name, 'w', newline='') as f:

        writer = csv.writer(f, delimiter=' ')

        writer.writerow([f'There are {fraction_homo} % of homorepeats out of a total of {chain_counter} chains'])
        writer.writerow(['The % of homorepeats of each kind are:'])
        writer.writerow([])
        for key,value in dct_homo_average.items():
            writer.writerow([key,value])

# Directory:
code_dir = os.getcwd()
data_dir = str(code_dir) + "/data"
# Run this in a loop for all files, to run single file just call functions inside loop
files = [file for file in os.listdir(f"{data_dir}/Homorepeat information/")]  # Get files that are not the dihedral angle data

for file in files:
    dct_homo_counter = initialize_dictionary(amino_list)
    count_homorepeats(data_dir, file, dct_homo_counter)



