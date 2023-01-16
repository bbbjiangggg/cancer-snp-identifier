#!/usr/bin/env python3

import os
import sys
import subprocess

print("\033[1;45m ONLY RUN THIS PROGRAM IN YOUR HOME DIRECTORY  \033[0;0;0m")
print("\033[1;45m to access your home directory open Finder and then click on Go and Home.  \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')

#checking for the installation of xcode
print("\033[1;45m Checking for the installation of Xcode \033[0;0;0m")
os.system('xcode-select --install')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for installation of homebrew
print("\033[1;45m Checking for installation of Homebrew \033[0;0;0m")
if os.system("brew --version") == 0:
    print("\033[1;45m Homebrew is installed \033[0;0;0m")
else:
    print("\033[1;45m Homebrew is not installed \033[0;0;0m")
    print("\033[1;45m installing Homebrew \033[0;0;0m")
    os.system("/bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)\"")
    print("\033[1;45m Homebrew is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for installation of python3
print("\033[1;45m Checking for installation of Python3 \033[0;0;0m")
if os.system("python3 --version") == 0:
    print("\033[1;45m Python3 is installed \033[0;0;0m")
else:
    print("\033[1;45m Python3 is not installed \033[0;0;0m")
    print("\033[1;45m installing Python3 \033[0;0;0m")
    os.system("brew install python3")
    print("\033[1;45m Python3 is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for installation of pip3
print("\033[1;45m Checking for installation of Pip3 \033[0;0;0m")
if os.system("pip3 --version") == 0:
    print("\033[1;45m Pip3 is installed \033[0;0;0m")
else:
    print("\033[1;45m Pip3 is not installed \033[0;0;0m")
    print("\033[1;45m installing Pip3 \033[0;0;0m")
    os.system("brew install pip3")
    print("\033[1;45m Pip3 is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

import pip 

os.system('pip3 install python-util')
os.system('pip3 install pandas')

#checking for the installation of wget
print("\033[1;45m Checking for installation of Wget \033[0;0;0m")
if os.system("wget --version") == 0:
    print("\033[1;45m Wget is installed \033[0;0;0m")
else:
    print("\033[1;45m Wget is not installed \033[0;0;0m")
    print("\033[1;45m installing Wget \033[0;0;0m")
    os.system("brew install wget")
    print("\033[1;45m Wget is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of sendemail
print("\033[1;45m Checking for installation of Sendemail \033[0;0;0m")
if os.system("sendemail --version") == 0:
    print("\033[1;45m Sendemail is installed \033[0;0;0m")
else:
    print("\033[1;45m Sendemail is not installed \033[0;0;0m")
    print("\033[1;45m installing Sendemail \033[0;0;0m")
    os.system("brew install sendemail")
    print("\033[1;45m Sendemail is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of tabix
print("\033[1;45m Checking for installation of Tabix \033[0;0;0m")
if os.system("tabix --version") == 0:
    print("\033[1;45m Tabix is installed \033[0;0;0m")
else:
    print("\033[1;45m Tabix is not installed \033[0;0;0m")
    print("\033[1;45m installing Tabix \033[0;0;0m")
    os.system("brew install Tabix")
    print("\033[1;45m Tabix is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of gzip
print("\033[1;45m Checking for installation of Gzip \033[0;0;0m")
if os.system("gzip --version") == 0:
    print("\033[1;45m Gzip is installed \033[0;0;0m")
else:
    print("\033[1;45m Gzip is not installed \033[0;0;0m")
    print("\033[1;45m installing Gzip \033[0;0;0m")
    os.system("brew install gzip")
    print("\033[1;45m Gzip is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of coreutils
print("\033[1;45m Checking for installation of Coreutils \033[0;0;0m")
if os.system("gsha256sum --version") == 0:
    print("\033[1;45m Coreutils is installed \033[0;0;0m")
else:
    print("\033[1;45m Coreutils is not installed \033[0;0;0m")
    print("\033[1;45m installing Coreutils \033[0;0;0m")
    os.system("brew install coreutils")
    print("\033[1;45m Coreutils is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#downloading SRA Toolkit
print("\033[1;45m Downloading SRA Toolkit \033[0;0;0m")
os.system('curl -OL --output sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz')
print("\033[1;45m SRA Toolkit downloaded \033[0;0;0m")
print('\n')

#extract SRA toolkit tar file contents
print("\033[1;45m Extracting SRA Toolkit \033[0;0;0m")
os.system('gunzip sratoolkit.current-mac64.tar.gz')
os.system('tar -vxzf sratoolkit.current-mac64.tar')
print("\033[1;45m SRA Toolkit extracted \033[0;0;0m")

#remove SRA Toolkit tar file
os.system('rm sratoolkit.current-mac64.tar')
os.system('export PATH=$PATH:~/sratoolkit/bin')

#install fastqc
print("\033[1;45m Installing FastQC \033[0;0;0m")
os.system("brew install fastqc")

print("\033[1;45m FastQC installed \033[0;0;0m")
print('\n')


#checking for the installation of java
print("\033[1;45m Checking for installation of Java \033[0;0;0m")
if os.system("java -version") == 0:
    print("\033[1;45m Java is installed \033[0;0;0m")
else:
    print("\033[1;45m Java is not installed \033[0;0;0m")
    print("\033[1;45m installing Java \033[0;0;0m")
    os.system("brew cask install java")
    os.system("brew install openjdk")
    print("\033[1;45m Java is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install trimmomatic
print("\033[1;45m Installing Trimmomatic \033[0;0;0m")
os.system('wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip')
os.system('unzip Trimmomatic-0.39.zip')
os.system('rm Trimmomatic-0.39.zip')
os.system('mkdir -p local/bin')
os.system('cp Trimmomatic-0.39/trimmomatic-0.39.jar $HOME/local/bin')
print("\033[1;45m Trimmomatic installed \033[0;0;0m")
input('\033[1;45m Press enter to test Trimmomatic... \033[0;0;0m \n')
print('\033[1;45m Checking if Trimmomatic was properly installed \033[0;0;0m')
os.system('java -jar $HOME/local/bin/trimmomatic-0.39.jar')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')


#install bwa
print("\033[1;45m Installing BWA \033[0;0;0m")
os.system("brew install bwa")
print("\033[1;45m BWA installed \033[0;0;0m")
print('\n')

#install bowtie2
print("\033[1;45m Installing Bowtie2 \033[0;0;0m")
os.system("brew install bowtie2")
print("\033[1;45m Bowtie2 installed \033[0;0;0m")
print('\n')

#install samtools
print("\033[1;45m Installing SAMtools \033[0;0;0m")
os.system("brew install samtools")
print("\033[1;45m SAMtools installed \033[0;0;0m")
print('\n')

#install bcftools
print("\033[1;45m Installing BCFtools \033[0;0;0m")
os.system("brew install bcftools")
print("\033[1;45m BCFtools installed \033[0;0;0m")
print('\n')

print('\033[1;45m All tools have been installed. \033[0;0;0m')
print('\033[1;45m Use the command "vdb-config -i" to configure your SRA Toolkit \033[0;0;0m')
