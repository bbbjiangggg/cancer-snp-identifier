#!/usr/bin/env python3
import os



#THIS PROGRAM IS FOR TRIMMED FILES ONLY

#must have sendemail installed on terminal
#for ubuntu use: $ sudo apt-get install libio-socket-ssl-perl libnet-ssleay-perl sendemail
# for Mac, use: brew install sendemail

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

# Getting the current working directory
src_dir = os.getcwd()

# Printing current directory
print(f'{MAGENTA}\nCurrent working directory:{RESET} ' + src_dir)

# Add the email to be notified when the process is done
user = input('Enter the email address to be notified once the analysis is complete: ')

# Add the job title
job = input('Enter a job name: ')

def replace_in_file(file_path, old_string, new_string):
    with open(file_path, 'r') as file:
        filedata = file.read()

    new_data = filedata.replace(old_string, new_string)

    with open(file_path, 'w') as file:
        file.write(new_data)

# Add the path to where bowtie files are found (must end in "bowtie/bowtie")
bowtie = input('Copy and paste the complete path to your bowtie files: ')
replace_in_file('trimmed_bash_sra_v1.2.txt', 'bowtie2_path', bowtie)

# Add the path to where reference chromosome is found
ref_chrom = input('Copy and paste the complete path to your reference chromosome: ')
replace_in_file('trimmed_bash_sra_v1.2.txt', 'ref_chrom', ref_chrom)


def find_unanalyzed_files(src_dir):
    files = [f for f in os.listdir(src_dir) if 'RR' in f and '.txt' not in f]

    vcfs = []
    for f in files:
        file_dir = os.path.join(src_dir, f)
        content = os.listdir(file_dir)
        for doc in content:
            if '_mapped.var-final.vcf' in doc:
                name = doc.replace('_mapped.var-final.vcf', '')
                vcfs.append(name)

    return sorted(list(set(files) - set(vcfs)))


def print_unanalyzed_files(unanalyzed_files):
    print(f'{MAGENTA}\nThese are the unanalyzed files in the current directory:{RESET}')
    for f in unanalyzed_files:
        print(f)
    print(f'{MAGENTA}There are{RESET} {len(unanalyzed_files)} {MAGENTA}unanalyzed files.{RESET}')


# Find and print the sorted list of unanalyzed files
new = find_unanalyzed_files(src_dir)
print_unanalyzed_files(new)


# This will store placement numbers into a list
ordinal = lambda n: "%d%s" % (n, "tsnrhtdd"[(n//10%10!=1)*(n%10<4)*n%10::4])
number = int(input('How many SRA sequences do you wish to analyze: '))
placement = [ordinal(n) for n in range(1, number + 1)]

# This will set different variables for different SRR sequences
srr_list = new[:number]

def replace_in_file(file_path, old_string, new_string):
    with open(file_path, 'r') as file:
        filedata = file.read()

    new_data = filedata.replace(old_string, new_string)

    with open(file_path, 'w') as file:
        file.write(new_data)

# These commands will replace each SRR number in the .txt file with 
# each of the accession numbers entered by the user
for index, srr in enumerate(srr_list):
    replace_in_file('trimmed_bash_sra_v1.2.txt', 'number', placement[index])
    replace_in_file('trimmed_bash_sra_v1.2.txt', 'now', srr)

    # Run the commands on the trimmed_bash_sra_v1.2.txt file
    os.system('cat trimmed_bash_sra_v1.2.txt | bash')

    # Replace the changed names back to the original
    replace_in_file('trimmed_bash_sra_v1.2.txt', placement[index], 'number')
    replace_in_file('trimmed_bash_sra_v1.2.txt', srr, 'now')

# Reset the bowtie2_path and ref_chrom path
replace_in_file('trimmed_bash_sra_v1.2.txt', bowtie, 'bowtie2_path')
replace_in_file('trimmed_bash_sra_v1.2.txt', ref_chrom, 'ref_chrom')




#send an email to the user to let them know the analysis is done
email_message = f"Your {job}_name analysis is complete. Please log in to check the results."
os.system(f'sendemail -f sudoroot1775@outlook.com -t {user} -u {job}_name Analysis Complete -m "{email_message}" -s smtp-mail.outlook.com:587 -o tls=yes -xu sudoroot1775@outlook.com -xp ydAEwVVu2s7uENC')
print(f'Sent email to {user}.')


def run_again():
    #runs trimmed_analysis_tools.py again
    print(f'{MAGENTA}\ntrimmed_analysis_tools.py is ready to run.{RESET}')

    while True:
        print(f'{MAGENTA}\nWould you like to run another analysis?{RESET}')
        choice = input('Enter yes or no to continue: ')

        if choice.lower() == 'yes':
            print(f'{MAGENTA}1){RESET} Email address used: ', user)
            print(f'{MAGENTA}2){RESET} Bowtie file path: ', bowtie)
            print(f'{MAGENTA}3){RESET} BWA reference chromosome path: ', ref_chrom)
            
            # execute the analysis script
            os.system('python3 trimmed_analysis_tools_v1.1.py')
        
        elif choice.lower() == 'no':
            # exit the loop and end the program
            print(f'{MAGENTA}\nAnalysis terminated. Goodbye.{RESET}')
            break
        else:
            # handle invalid input by asking the user to enter either yes or no
            print(f'{MAGENTA}\nEnter either yes or no.{RESET}')

run_again()