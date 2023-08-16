import os

# Iterate over the items in the current directory
for item in os.listdir():
    # Check if the item is a directory and if it starts with "SRR"
    if os.path.isdir(item) and item.startswith('SRR'):
        # Initialize a variable to check if a file with the specific extension is found
        bam_file_found = False
        
        # Iterate over the files in the directory
        for filename in os.listdir(item):
            # Check if the file ends with "_mapped.sorted.bam"
            if filename.endswith("_mapped.sorted.bam"):
                bam_file_found = True
            
            # Check if the file ends with ".fastq" and delete it if true
            if filename.endswith(".fastq"):
                file_to_delete = os.path.join(item, filename)
                os.remove(file_to_delete)
                print(f"Deleted {file_to_delete}")

        # If no file with "_mapped.sorted.bam" was found, print the directory name
        if not bam_file_found:
            print(f"The directory {item} does not contain a file ending with '_mapped.sorted.bam'")

