import csv
import os

amino_list=["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO",
            "SER","THR","TRP","TYR","VAL","XAA"]

dct_homo_counter= {}
for i in amino_list:
    dct_homo_counter[i] = 0      # Establishing keys and setting value as 0



# fer diccionari per contar homorepeats (key: amino, value: counter)
final_result = "uniprot_sprot_2_analysis.txt"

chain_counter = 0
homo_counter = 0
with open(final_result, 'r') as f:
    reader = csv.reader(f, delimiter='\t')

    for row in reader:
        #print(len(row))
        if ">" in row[0]:
            chain_counter += 1


        if len(row) == 3:   # data row (length = 3). It found a homorepeat
            homo_counter += 1
            if row[0] in dct_homo_counter:
                dct_homo_counter[row[0]] += 1
            else:
                print("not found")


print(chain_counter)
dct_homo_average = dct_homo_counter.copy()

for key in dct_homo_average.keys():

    dct_homo_average[key] = (dct_homo_average[key] / homo_counter)*100  # returns the percentage

fraction_homo= (homo_counter/chain_counter)*100

print("There are", fraction_homo, "% of homorepeats out of a total of", chain_counter, "chains")
print("The % of homorepeats of each kind are:")
print(dct_homo_average)

doc_name= "homo_ratio_uniprot_sprot_2.csv"
with open(doc_name, 'w', newline='') as f:

    writer = csv.writer(f, delimiter=' ')

    writer.writerow(['There are',fraction_homo,'% of homorepeats out of a total of', chain_counter,'chains'])
    writer.writerow(['The % of homorepeats of each kind are:'])
    writer.writerow([])
    for key,value in dct_homo_average.items():
        writer.writerow([key,value])





