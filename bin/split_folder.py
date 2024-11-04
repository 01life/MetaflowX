#!/usr/bin/env python

import os
import sys
import shutil

def move_files_to_subfolders(base_dir, num):
    # Ensure the directory exists
    if not os.path.exists(base_dir):
        print(f"Directory '{base_dir}' does not exist.")
        return

    # Get all files in the base directory (excluding subfolders)
    file_list = [f for f in os.listdir(base_dir) if os.path.isfile(os.path.join(base_dir, f))]

    filenum = len(file_list)
    touchFile = open(base_dir +"/"+str(filenum) +".txt",'w')
    touchFile.close()

    # If the file count is less than or equal to 'num', move all to 'folder1'
    if len(file_list) <= num:
        new_folder = os.path.join(base_dir, 'folder1')

        # If the target folder already exists, alert and exit
        if os.path.exists(new_folder):
            print(f"'{new_folder}' already exists, exiting to avoid conflicts.")
            sys.exit(1)

        os.mkdir(new_folder)

        for file in file_list:
            shutil.move(os.path.join(base_dir, file), os.path.join(new_folder, file))

        print(f"Moved {len(file_list)} files to 'folder1'.")
        return

    # If the file count is greater than 'num', split into subfolders
    num_files = (len(file_list) // num) + (1 if len(file_list) % num != 0 else 0)

    print(f"Total files: {len(file_list)}")
    print(f"Creating {num_files} subfolders.")

    # Create and distribute files among subfolders
    for i in range(num_files):
        new_folder = os.path.join(base_dir, f'floder{i + 1}')

        if os.path.exists(new_folder):
            print(f"'{new_folder}' already exists, exiting to avoid conflicts.")
            sys.exit(1)

        os.mkdir(new_folder)

        files_to_move = file_list[i * num : (i + 1) * num]

        for file in files_to_move:
            shutil.move(os.path.join(base_dir, file), os.path.join(new_folder, file))

    print("Task completed.")

    

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <directory> <files_per_folder>")
        sys.exit(1)

    base_dir = sys.argv[1]
    num = int(sys.argv[2])

    move_files_to_subfolders(base_dir, num)


