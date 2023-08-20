#!/usr/bin/env python3

import subprocess
import sys

# Define color codes
RED = '\033[1;31m'
GREEN = '\033[1;32m'
MAGENTA = '\033[1;35m'
BLUE = '\033[1;34m'
RESET = '\033[0m'

def check_install(package):
    print(f"{MAGENTA}\nChecking for the installation of {package}{RESET}")

    if subprocess.call(['which', package]) == 0:
        print(f"{GREEN}{package} is installed{RESET}")
    else:
        print(f"{RED}{package} is not installed{RESET}")
        print(f"{MAGENTA}installing {package}{RESET}")
        subprocess.check_call(['sudo', 'apt', 'install', package, '-y'])
        print(f"{GREEN}{package} is installed{RESET}")

print(f"{RED}\nONLY RUN THIS PROGRAM IN YOUR HOME DIRECTORY{RESET}")
print(f"{MAGENTA}\nTo access your home directory use the command: '$ explorer.exe .' only on Ubuntu WSL{RESET}")

input(f"{BLUE}\nPress enter to continue...{RESET}\n")

print(f"{MAGENTA}\nInstalling Python3 and pip3{RESET}")
subprocess.check_call(['sudo', 'apt-get', '-y', 'install', 'python3', 'python3-pip'])

print(f"{MAGENTA}\nInstalling pyutil and pandas{RESET}")
subprocess.check_call(['pip3', 'install', 'pyutil', 'pandas'])

print(f"{MAGENTA}\nInstalling biopython{RESET}")
subprocess.check_call(['pip3', 'install', 'biopython'])

############################################################################

# Update terminal and install necessaries
print(f"{MAGENTA}\nUpdating terminal{RESET}")
input(f"{BLUE}\nPress enter to continue...{RESET}\n")

#subprocess.check_call(['sudo', 'apt', 'update', '-y'])
#subprocess.check_call(['sudo', 'apt', 'upgrade', '-y'])
print('\n')

# Checking for the installation of wget
check_install('wget')
input(f"{BLUE}\nPress enter to continue...{RESET}\n")
print('\n')

# Check and install unzip
check_install('unzip')
input(f'{BLUE}Press enter to continue...{RESET}\n')

# Check and install sendemail
check_install('sendemail')
input(f'{BLUE}Press enter to continue...{RESET}\n')

# Check and install tabix
check_install('tabix')
input(f'{BLUE}Press enter to continue...{RESET}\n')

########################################################################

print(f'{BLUE}Downloading SRA Toolkit{RESET}')
subprocess.run(['sudo', 'apt', 'install', 'sra-toolkit'])
print(f'{GREEN}SRA Toolkit is downloaded{RESET}')

input(f'{BLUE}Press enter to continue...{RESET}\n')
print('\n')

# Install FASTQC
print(f'{MAGENTA}Installing FASTQC{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')
subprocess.run(['sudo', 'apt-get', 'install', '-y', 'fastqc'])
print(f'{GREEN}FASTQC is installed{RESET}\n')
input(f'{BLUE}Press enter to continue...{RESET}\n')

##############################################################################

# checking for the installation of java
print(f'{MAGENTA}Checking for the installation of java{RESET}')
if subprocess.run(['java', '-version']).returncode == 0:
    print(f'{GREEN}java is installed{RESET}')
else:
    print(f'{RED}java is not installed{RESET}')
    print(f'{MAGENTA}installing java{RESET}')


# check if default-jre is installed
print(f'{MAGENTA}Checking default-jre{RESET}')
result = subprocess.run(['java', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if 'command not found' in result.stderr.decode('utf-8'):
    print(f'{RED}default-jre is not installed{RESET}')
    # install default-jre
    print(f'{MAGENTA}Installing default-jre{RESET}')
    subprocess.run(['sudo', 'apt-get', 'install', 'default-jre', '-y'])
    print(f'{GREEN}default-jre is installed{RESET}')
else:
    print(f'{GREEN}default-jre is already installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')
print()

# check if default-jdk is installed
'''print(f'{MAGENTA}Checking default-jdk{RESET}')
result = subprocess.run(['javac', '-version'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
if 'command not found' in result.stderr.decode('utf-8'):
    print(f'{RED}default-jdk is not installed{RESET}')
    # install default-jdk
    print(f'{MAGENTA}Installing default-jdk{RESET}')
    subprocess.run(['sudo', 'apt-get', 'install', 'default-jdk'])
    print(f'{GREEN}default-jdk is installed{RESET}')
else:
    print(f'{GREEN}default-jdk is already installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')
print()'''


# create directory $HOME/local/bin if it does not exist
subprocess.run(['mkdir', '-p', '$HOME/local/bin'])

# install Trimmomatic
print(f'{MAGENTA}Installing Trimmomatic{RESET}')
subprocess.run(['wget', 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip'])
subprocess.run(['unzip', 'Trimmomatic-0.39.zip'])
subprocess.run(['cp', 'Trimmomatic-0.39/trimmomatic-0.39.jar', '$HOME/local/bin'])
subprocess.run(['chmod', 'a+x', '$HOME/local/bin/trimmomatic-0.39.jar'])
print(f'{GREEN}Trimmomatic is installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')
print()

##############################################################################

# check if bwa is installed
result = subprocess.run(['which', 'bwa'], capture_output=True, text=True)
if not result.stdout.strip():
    print(f'{MAGENTA}Installing BWA{RESET}')
    subprocess.run(['sudo', 'apt-get', 'install', 'bwa', '-y'])
else:
    print(f'{GREEN}BWA is already installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')

# check if bowtie2 is installed
result = subprocess.run(['which', 'bowtie2'], capture_output=True, text=True)
if not result.stdout.strip():
    print(f'{MAGENTA}Installing Bowtie2{RESET}')
    subprocess.run(['sudo', 'apt-get', 'install', 'bowtie2', '-y'])
else:
    print(f'{GREEN}Bowtie2 is already installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')

# check if samtools is installed
result = subprocess.run(['which', 'samtools'], capture_output=True, text=True)
if not result.stdout.strip():
    print(f'{MAGENTA}Installing samtools{RESET}')
    subprocess.run(['sudo', 'apt-get', 'install', 'samtools', '-y'])
else:
    print(f'{GREEN}samtools is already installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')

# check if bcftools is installed
result = subprocess.run(['which', 'bcftools'], capture_output=True, text=True)
if not result.stdout.strip():
    print(f'{MAGENTA}Installing bcftools{RESET}')
    subprocess.run(['sudo', 'apt-get', 'install', 'bcftools', '-y'])
else:
    print(f'{GREEN}bcftools is already installed{RESET}')
input(f'{BLUE}Press enter to continue...{RESET}\n')

print(f'{GREEN}All tools have been installed.{RESET}')
print(f'{MAGENTA}Use the command "vdb-config -i" to configure your SRA Toolkit.{RESET}')