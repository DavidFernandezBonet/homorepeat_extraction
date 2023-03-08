from homorepeat_and_dihedral_data_extraction_051022 import *
from count_homorepeats_kind import *
from read_dihedral_data_and_plot import *
import sys

if __name__ == "__main__":


    homorepeat_min = int(sys.argv[1])
    rootdir = sys.argv[2]

    print(homorepeat_min)
    print(rootdir)
    print(type(homorepeat_min))
    if type(homorepeat_min) != int:
        raise ValueError("Homorepeat threshold should be an integer")

    if type(rootdir) != str:
        raise ValueError("Path containing the .pdb files should be a string")
    # # User parameters
    # rootdir = r"/home/david/Desktop/swissprot_pdb_v2_grouped/test"  # Directory with PDB files
    # homorepeat_min = 4  # Threshold for the homorepeat

    code_dir = os.getcwd()
    homorepeat_in_pdb_file = f'{code_dir}/data/files_with_homorepeats.txt'




    create_folder_if_not_exist(code_dir, "data")
    create_folder_if_not_exist(f"{code_dir}/data/", "Homorepeat information")
    create_folder_if_not_exist(f"{code_dir}/data/", "Homorepeat analysis")
    create_folder_if_not_exist(f"{code_dir}/data/", "Dihedral data")
    create_folder_if_not_exist(f"{code_dir}/data/", "test")

    with open(homorepeat_in_pdb_file, 'w', newline="") as h:
        writ = csv.writer(h, delimiter='\t')
        writ.writerow(["File Name", "Chain position of homorepeat"])

    # Read all files
    time_start = time.perf_counter()
    counter = 0
    for subdir, dirs, files in os.walk(rootdir):
        folder_name = subdir[-4:]


        dihedral_t = f'{code_dir}/data/Dihedral data/dihedral_data_{folder_name}.txt'
        final_result = f"{code_dir}/data/Homorepeat information/{folder_name}.txt"
        with open(dihedral_t, 'w', newline='') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(["Amino Acid, [Phi,Psi]"])

        for file in files:
            counter += 1
            print(counter)

            with open(final_result, 'a', newline="") as f:
                title = str(file)
                writer = csv.writer(f, delimiter='\t')
                homorepeats_chain = read_file_find_homorepeats(title, homorepeat_min, subdir, final_result, dihedral_t)
            if True in homorepeats_chain:
                with open(homorepeat_in_pdb_file, 'a', newline="") as h:
                    writ = csv.writer(h, delimiter='\t')
                    writ.writerow([file, [i for i, x in enumerate(homorepeats_chain) if x]])
            print(os.path.join(subdir, file))

        time_elapsed = (time.perf_counter() - time_start)
        print("time_elapsed", time_elapsed, "current_folder:", folder_name)

    time_elapsed = (time.perf_counter() - time_start)
    print("time elapsed (seconds):", time_elapsed)

    count_homorepeats_main()
    read_dihedral_and_plot_main()

    time_elapsed = (time.perf_counter() - time_start)
    print("time elapsed (seconds):", time_elapsed)