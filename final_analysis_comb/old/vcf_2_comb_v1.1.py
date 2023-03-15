import os
from os import listdir
import shutil
import re
import requests
import io

#This program will copy all vcf files from all SRA directories to a new directory
chr = input('\033[1;45m Enter the chrosome number of this analysis: \033[0;0;0m ')
print('\n')
cancer = input('\033[1;45m Enter an abbreviation for the cancer type (e.g. pnca): \033[0;0;0m ')
print('\n')

directory = 'ch' + chr + '_' + cancer + '_vcf'
print("This is your directory's name: " + directory)
print('\n')
homedir = input('\033[1;45m Enter the path where the SRA directories are currently stored: \033[0;0;0m ')
print('\n')

isecdir = os.path.join(homedir,directory)

#make new directory
if directory in os.listdir():
    print('Directory already exists, removing directory.')
    os.system('rm -r ' + directory)
    os.mkdir(directory)
else:
    os.mkdir(directory)
    print('New directory has been created: ' + directory)
    print('\n')
input('\033[1;45mPress enter to continue...\033[0;0;0m ')
print('\n')
os.chdir(isecdir)
#copy all vcf files to new directory
print('Copying all vcf files to: ' + directory + '\n')

for item in listdir(homedir):
    if 'RR' in item:
        subdir = os.path.join(homedir, item)
        for content in listdir(subdir):
            if '.vcf' in content:
                shutil.move(os.path.join(subdir, content),
                         os.path.join(isecdir, content))
            else: pass
    else: pass

print('All vcf files have been copied to: ' + directory)
print('\n')
print('As backup, a copy of ' + directory + ' is being created...')
print('\n')

copyName = 'copy_'+directory
if copyName in os.listdir():
    print('Directory already exists, removing directory.')
    os.system('rm -r ' + copyName)
    os.mkdir(copyName)
else:
    #copydir = os.path.join(homedir, copyName)
    shutil.copytree(isecdir, copyName)

print('A copy of ' + directory + ' has been created.')
print('\n')
input('\033[1;35;40m Press enter to continue...\033[0m')
print('\n')

#Combining vcf reports
print('\033[1;45m Combining all vcf files in: \033[0;0;0m ' + directory + '\n')
#this command will bgzip all .vcf files, this program must be located in the same directory
#if all files are already bgzipped then it will throw back:
#"ls: cannot access '*.vcf': No such file or directory"

print('Compressing all vcf files to gzip...')


os.system('ls *.vcf | xargs -n1 -P0 bgzip')
print('\n')
print('Organizing all compressed vcfgz files...')
print('\n')

#move all .vcf files to a directory called isec_vcfgz_files

os.system('mkdir isec_vcfgz_files')
os.system('mv *.gz* isec_vcfgz_files')

#change directory to isec_vcfgz_files
os.chdir('isec_vcfgz_files')



#this command will Write all file names into a txt file, same line, 
#one space, named "isec_tools_commands.txt"
os.system('ls -1 | paste -sd " " ->> isec_tools_commands.txt')

#deleting string that results from previous command
infile = "isec_tools_commands.txt"
outfile = "isec_tools_commands2.txt"

delete_string = ["isec_tools_commands.txt"]
fin = open(infile)
fout = open(outfile, "w+")
for line in fin:
    for word in delete_string:
        line = line.replace(word, "")
    fout.write(line)
fin.close()
fout.close()



#read the existing text from file in READ mode 
with open('isec_tools_commands2.txt','r') as src:
    fline='bcftools isec -n +2 '
    #Prepending string
    oline=src.readlines()
    #prepend the string we want to, on first line 
    oline.insert(0,fline) 
#open the file in WRITE mode  
with open('isec_tools_commands2.txt','w') as src: 
   src.writelines(oline) 

#vcf_file = input('\033[1;45m Enter the name you wish to give the combined vcf files (e.g. ch22_pnca_comb): \033[0m')

vcf_file = 'ch' + chr + '_' + cancer + '_comb'
print('\033[1;45m This is your combined vcf file name: \033[0;0;0m' + vcf_file)
print('\n')
print('Combining all vcf files...')
print('\n')

#open a file with access mode 'a+'web: https://stackabuse.com/file-handling-in-python/
with open('isec_tools_commands2.txt', 'a+') as file_object:
    # Append "| bgzip -c >bgzip -c > vcf_file.vcf.gz" at the end of file
    file_object.write('| bgzip -c >bgzip -c > ' + vcf_file + '.vcf.gz')
    
#remove all white spaces > 1 and save it as "isec_tools_commands2.txt"
with open('isec_tools_commands2.txt', 'r') as file_object, open ('isec_tools_commands3.txt', 'w') as file_object2:
    for line in file_object:
        file_object2.write(re.sub('\s+',' ',line))

print('\n')
print('Indexing vcf files using tabix...')


#tabix all vcfgz files
os.system('for f in ./*.vcf.gz; do tabix -p vcf -f $f;done')
print('\n')
print('Creating intersections, unions and complements (isec)...')

#run isec on isec_tools_commands2.txt files
os.system('cat isec_tools_commands3.txt | bash')

#unzip gzip file
os.system('gunzip '+ vcf_file + '.vcf.gz')

#open unzipped file
print('\n')
print('\033[1;45m Done! Your final combined vcf file is located here:\033[0;0;0m' + directory + '/isec_vcfgz_files/' + vcf_file + '.vcf ') 
print('\n')
print('Open the file using the command: cat ' + directory + '/isec_vcfgz_files/' + vcf_file + '.vcf | less')
print('\n')


os.system('rm bgzip isec_tools_commands.txt isec_tools_commands2.txt isec_tools_commands3.txt')

#move comb.vcf file to parent directory
parent_path = os.path.abspath(os.path.join(os.getcwd(), os.pardir))

for filename in os.listdir("."):
    if filename.endswith("comb.vcf"):
        # Move the file to the parent directory
        shutil.move(filename, parent_path)

# Create a dictionary of Google Drive links
link_list = {
    "1": "https://drive.google.com/file/d/1-09an3LHEJYnf-0s4S70ATis81dZPV1G/view?usp=sharing",
    "2": "https://drive.google.com/file/d/1xP8ELyW5Kr7fay0NulbgiPQOAx5Eq9E3/view?usp=sharing",
    "3": "https://drive.google.com/file/d/1-JNXwA4VvwcSOoulcmRmkLq8qYueltyQ/view?usp=sharing",
    "4": "https://drive.google.com/file/d/1n2T5xzMOqOobm5w97uthlvW5d7e3nV9M/view?usp=sharing",
    "5": "https://drive.google.com/file/d/1s8OPif3dKAesCgoOTIQGdAr9CdcUhe2q/view?usp=sharing",
    "6": "https://drive.google.com/file/d/1IzA-rUe8ELPyFtCePiCQgQeWDnOEPJxN/view?usp=sharing",
    "7": "https://drive.google.com/file/d/16eGfwwe1yjvyY6TO-hWj1Rx7YMvXuSbO/view?usp=sharing",
    "8": "https://drive.google.com/file/d/1uSqwMSvQpyBxPwiHhgUcNVkL7YdkcOHg/view?usp=sharing",
    "9": "https://drive.google.com/file/d/1chh47ez6Vqe1W0bdoKqxPHO8iOgMMiJe/view?usp=sharing",
    "10": "https://drive.google.com/file/d/12NXqBfLFJZwifp7_LLAzZFsfpE60czgN/view?usp=sharing",
    "11": "https://drive.google.com/file/d/1nsicoTeVg3t4AL592QGLa0QCiOmrhZz4/view?usp=sharing",
    "12": "https://drive.google.com/file/d/1hhX-NaqjzkpsA9cvFMj4F-1cEPCl6JmK/view?usp=sharing",
    "13": "https://drive.google.com/file/d/112aOwZiP4kOJhrBsBockjM0q2i1a1lx-/view?usp=sharing",
    "14": "https://",
    "15": "https://drive.google.com/file/d/1-umaSF-rAFqqSwVeQI6djyW7QKxBm5A7/view?usp=sharing",
    "16": "https://drive.google.com/file/d/1cOFSi6ujbJV_I144je6PnBUmiAZL1-J0/view?usp=share_link",
    "17": "https://drive.google.com/file/d/1yhKZrhhw2YlBMySWa4L-DYaJ4UVYH-y8/view?usp=share_link",
    "18": "https://drive.google.com/file/d/1yvj7_SQE93sJ2kwceFr5ikYsFstOZ36I/view?usp=sharing",
    "19": "https://drive.google.com/file/d/1hRKEmwtw7uBkx8n_5wBWvdY2WqQjuOZz/view?usp=sharing",
    "20": "http",
    "21": "http",
    "22": "https://drive.google.com/file/d/1zqOFpGV8HEVBo_uDIqVi0_YVQOCL6LZZ/view?usp=share_link",
    "X": "https",
    "Y": "https",
    # and so on
}

while True:
    normal_file = input('\033[1;45mWould you like to download a normal combined VCF file (yes or no)?\033[0;0;0m ')
    print('\n')
    if normal_file.lower() == 'yes':
        # Get chromosome number from user input
        chromosome = input('\033[1;45mEnter the chromosome number of interest (1-22, X or Y): \033[0;0;0m ')

        # Check if chromosome number is in link list
        if chromosome in link_list:
            # Extract the file ID from the Google Drive link
            file_id = link_list[chromosome].split("/")[-2]

            # Download the VCF file as bytes
            print("Downloading chromosome " + chromosome + " normal combined VCF file...")
            url = f"https://drive.google.com/uc?id={file_id}"
            response = requests.get(url)
            vcf_bytes = io.BytesIO(response.content).read()

            # Write the VCF bytes to a file
            filename = f"ch{chromosome}_norm_comb.vcf"
            with open(filename, "wb") as f:
                f.write(vcf_bytes)
            print(f"{filename} downloaded successfully!")
            shutil.move(filename, parent_path)
            break
        else:
            print("Chromosome number not found in link list.")
    elif normal_file.lower() == 'no':
        break
    else:
        print("Invalid input. Please enter 'yes' or 'no'.")