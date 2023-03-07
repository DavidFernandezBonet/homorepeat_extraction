import csv
import ast
import matplotlib.pyplot as plt
import numpy as np
import os


def read_dihedral_and_plot_main():
    amino_list=["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
                "SER","THR","TRP","TYR","VAL","XAA"]



    def initialize_dictionary(amino_list):
        dct_homo_counter= {}
        for i in amino_list:
            dct_homo_counter[i] = []      # Establishing keys and setting value as empty list. For the dihedral data we will have:
                                         # AMINO : [[phi1,psi1], ..., [phi_n,psi_n]]
        return dct_homo_counter



    # fer diccionari per contar homorepeats (key: amino, value: counter)
    final_result = "dihedral_combined.txt"

    def create_folder_if_not_exist(parent_dir, folder_name):
        path = parent_dir + "/" + folder_name
        isExist = os.path.exists(path)
        # Check whether the specified path exists or not
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(path)
            print(f"The new directory {path} is created!")

    def dihedral_plots(data_dir, file_name, dct_homo_counter):
        file_dir = f"{data_dir}/Dihedral data/{file_name}"
        with open(file_dir, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                try:
                    if row[0] in dct_homo_counter:
                        amino_acid = row[0]
                        phi_psi = ast.literal_eval(row[1])

                        if (phi_psi[0] or phi_psi[1]) == None:   # sometimes the value is None because its the first element or something
                            continue
                        else:
                            dct_homo_counter[amino_acid].append(phi_psi)   # dct[label_list[u]].append([position,amino,dssp])
                except:
                    print("additional space")
                    continue
                else:
                    continue


        #print(dct_homo_counter)


        # Unzip the points to plot:
        create_folder_if_not_exist(f"{data_dir}/Dihedral data/", f"{file_name}_Plots")
        for key in dct_homo_counter.keys():
            zipped_points = dct_homo_counter[key]

            x_val = [x[0] for x in zipped_points]
            y_val = [y[1] for y in zipped_points]

            plt.title(key)
            plt.plot(x_val, y_val, ".")
            plt.xlabel("$\phi$")
            plt.ylabel("$\psi$")
            plt.xticks(np.arange(-180,181, 60))
            plt.yticks(np.arange(-180,181, 60))
            plt.xlim([-180, 181])
            plt.ylim([-180, 181])

            plt.savefig(f'{data_dir}/Dihedral data/{file_name}_Plots/phi_psi_{key}.pdf')
            plt.close()
            #plt.show()







    # Directory:
    code_dir = os.getcwd()
    data_dir = str(code_dir) + "/data"

    # Run this in a loop
    files = [file for file in os.listdir(f"{data_dir}/Dihedral data/") if file[-4:] == f".txt"]  # Get files that are not the dihedral angle data

    print(files)
    for file in files:
        dct_homo_counter = initialize_dictionary(amino_list)
        dihedral_plots(data_dir, file, dct_homo_counter)






