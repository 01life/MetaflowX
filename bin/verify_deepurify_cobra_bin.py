import sys
import os
import shutil
from pathlib import Path


deepurify_path = Path(sys.argv[1]).resolve()
deepurify_files = {}
cobra_path = Path(sys.argv[2]).resolve()
cobra_files = {}
outpath = Path(sys.argv[3]).resolve()
os.makedirs(outpath)

refineBinF = open(sys.argv[4],'w')


# Step 1: Retrieve the .fa files from the deepurify directory.

for root, _, files in os.walk(deepurify_path):
    for file in files:
        if file.endswith('.fa'):
            file_path = os.path.join(root, file)
            file_name = os.path.splitext(file)[0].lstrip("Deepurify_")  # Get the file name without the extension.
            deepurify_files[file_name] = file_path

# Step 2: Retrieve the .fa files from the cobra directory.

for root, _, files in os.walk(cobra_path):
    for file in files:
        if file.endswith('.fa'):
            file_path = os.path.join(root, file)
            file_name = os.path.splitext(file)[0].lstrip("COBRA_")  # Get the file name.
            cobra_files[file_name] = file_path

# Step 3: Loop through all files in the cobra directory and check if the file size is 0; if so, replace it.

with open(os.path.join(outpath,"Refine_bin_choose_info.txt"),'w') as outFile:
    for file_name, cobra_file_path in cobra_files.items():
        # Check if the file size is 0.
        if os.path.getsize(cobra_file_path) == 0:
            if file_name in deepurify_files:
                deepurify_file_path = deepurify_files[file_name]
                # Replace with the same named file from the deepurify directory.
                newPath =  os.path.join(outpath, "Refine_"+file_name+".fa")
                shutil.copy(deepurify_file_path, newPath)
                refineBinF.write(f"{file_name}\t{str(newPath)}\n")
                outFile.write(f"Refine_{file_name}.fa\tDeepurify\n")
                print(f"bin {cobra_file_path} is empty, repleace as {deepurify_file_path}")
            else:
                print(f"bin {cobra_file_path} is empty, but deepurify not exit")
        else:
            newPath =  os.path.join(outpath, "Refine_"+file_name,".fa")
            shutil.copy(cobra_file_path, newPath)
            refineBinF.write(f"{file_name}\t{str(newPath)}\n")
            outFile.write(f"Refine_{file_name}.fa\tCOBRA\n")
            print(f"bin {cobra_file_path} not empty, will ussing it.")

refineBinF.close()