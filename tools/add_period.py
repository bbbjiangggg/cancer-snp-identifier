import os

# Define the starting directory to search
starting_dir = '/media/aguilarr@campus.wra.net/WD_GREEN/ovcasc2'

# Loop through all subdirectories of the starting directory
for root, dirs, files in os.walk(starting_dir):
    # Check if the directory name starts with SRR or ERR
    if os.path.basename(root).startswith(('SRR', 'ERR')):
        # Loop through all files in the directory
        for filename in files:
            # Check if the file name ends with _mapped.var-final.vcf
            if filename.endswith('_mapped.var-final.vcf'):
                # Rename the file to end with _mapped.var.-final.vcf instead
                new_filename = filename.replace('_mapped.var-final.vcf', '_mapped.var.-final.vcf')
                os.rename(os.path.join(root, filename), os.path.join(root, new_filename))

