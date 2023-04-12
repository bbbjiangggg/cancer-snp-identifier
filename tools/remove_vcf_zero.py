import os

# get the current directory
root_dir = os.getcwd()

# loop over all subdirectories in the root directory
for dirpath, dirnames, filenames in os.walk(root_dir):
    
    # check if the directory name contains "RR"
    if "RR" in dirpath:
        
        # loop over all files in the directory
        for filename in filenames:
            
            # check if the file ends with "-final.vcf"
            if filename.endswith("-final.vcf"):
                
                # construct the full file path
                filepath = os.path.join(dirpath, filename)
                
                # check if the file is at zero bytes
                if os.path.getsize(filepath) == 0:
                    # delete the file
                    os.remove(filepath)
                    print(f"{filepath} has been deleted.")

