#!/usr/bin/env python3

import os
import sys
import subprocess

print("\033[1;45m ONLY RUN THIS PROGRAM IN YOUR HOME DIRECTORY  \033[0;0;0m")
print("\033[1;45m to access your home directory open Finder and then click on Go and Home.  \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')

#checking for the installation of xcode
print("\033[1;45m Checking for the installation of xcode \033[0;0;0m")
os.system('/usr/bin/xcodebuild -version')
print("\033[1;45m If you do not have xcode installed, install it from the App Store. \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for installation of homebrew
print("\033[1;45m Checking for installation of homebrew \033[0;0;0m")
if os.system("brew --version") == 0:
    print("\033[1;45m homebrew is installed \033[0;0;0m")
else:
    print("\033[1;45m homebrew is not installed \033[0;0;0m")
    print("\033[1;45m installing homebrew \033[0;0;0m")
    os.system("/bin/bash -c \"$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install.sh)\"")
    print("\033[1;45m homebrew is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for installation of python3
print("\033[1;45m Checking for installation of python3 \033[0;0;0m")
if os.system("python3 --version") == 0:
    print("\033[1;45m python3 is installed \033[0;0;0m")
else:
    print("\033[1;45m python3 is not installed \033[0;0;0m")
    print("\033[1;45m installing python3 \033[0;0;0m")
    os.system("brew install python3")
    print("\033[1;45m python3 is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for installation of pip3
print("\033[1;45m Checking for installation of pip3 \033[0;0;0m")
if os.system("pip3 --version") == 0:
    print("\033[1;45m pip3 is installed \033[0;0;0m")
else:
    print("\033[1;45m pip3 is not installed \033[0;0;0m")
    print("\033[1;45m installing pip3 \033[0;0;0m")
    os.system("brew install pip3")
    print("\033[1;45m pip3 is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

import pip 

os.system('pip3 install python-util')
os.system('pip3 install pandas')

#checking for the installation of wget
print("\033[1;45m Checking for installation of wget \033[0;0;0m")
if os.system("wget --version") == 0:
    print("\033[1;45m wget is installed \033[0;0;0m")
else:
    print("\033[1;45m wget is not installed \033[0;0;0m")
    print("\033[1;45m installing wget \033[0;0;0m")
    os.system("brew install wget")
    print("\033[1;45m wget is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of sendemail
print("\033[1;45m Checking for installation of sendemail \033[0;0;0m")
if os.system("sendemail --version") == 0:
    print("\033[1;45m sendemail is installed \033[0;0;0m")
else:
    print("\033[1;45m sendemail is not installed \033[0;0;0m")
    print("\033[1;45m installing sendemail \033[0;0;0m")
    os.system("brew install sendemail")
    print("\033[1;45m sendemail is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of tabix
print("\033[1;45m Checking for installation of tabix \033[0;0;0m")
if os.system("tabix --version") == 0:
    print("\033[1;45m tabix is installed \033[0;0;0m")
else:
    print("\033[1;45m tabix is not installed \033[0;0;0m")
    print("\033[1;45m installing tabix \033[0;0;0m")
    os.system("brew install tabix")
    print("\033[1;45m tabix is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of gzip
print("\033[1;45m Checking for installation of gzip \033[0;0;0m")
if os.system("gzip --version") == 0:
    print("\033[1;45m gzip is installed \033[0;0;0m")
else:
    print("\033[1;45m gzip is not installed \033[0;0;0m")
    print("\033[1;45m installing gzip \033[0;0;0m")
    os.system("brew install gzip")
    print("\033[1;45m gzip is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

# checking for the installation of coreutils
print("\033[1;45m Checking for installation of coreutils \033[0;0;0m")
if os.system("gsha256sum --version") == 0:
    print("\033[1;45m coreutils is installed \033[0;0;0m")
else:
    print("\033[1;45m coreutils is not installed \033[0;0;0m")
    print("\033[1;45m installing coreutils \033[0;0;0m")
    os.system("brew install coreutils")
    print("\033[1;45m coreutils is installed \033[0;0;0m")

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
print("\033[1;45m Installing fastqc \033[0;0;0m")
os.system("brew install fastqc")

print("\033[1;45m fastqc installed \033[0;0;0m")
print('\n')


#checking for the installation of java
print("\033[1;45m Checking for installation of java \033[0;0;0m")
if os.system("java -version") == 0:
    print("\033[1;45m java is installed \033[0;0;0m")
else:
    print("\033[1;45m java is not installed \033[0;0;0m")
    print("\033[1;45m installing java \033[0;0;0m")
    os.system("brew cask install java")
    print("\033[1;45m java is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for the installation of jdk
print("\033[1;45m Checking for installation of jdk \033[0;0;0m")
if os.system("javac -version") == 0:
    print("\033[1;45m jdk is installed \033[0;0;0m")
else:
    print("\033[1;45m jdk is not installed \033[0;0;0m")
    print("\033[1;45m installing jdk \033[0;0;0m")
    os.system("brew install openjdk")
    print("\033[1;45m jdk is installed \033[0;0;0m")

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install trimmomatic
print("\033[1;45m Installing Trimmomatic \033[0;0;0m")
os.system('wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip')
os.system('unzip Trimmomatic-0.39.zip')
os.system('rm Trimmomatic-0.39.zip')
print("\033[1;45m Trimmomatic installed \033[0;0;0m")
input('\033[1;45m Press enter to test Trimmomatic... \033[0;0;0m \n')
print('\033[1;45m Checking if Trimmomatic was properly installed \033[0;0;0m')
os.system('java -jar $HOME/local/bin/trimmomatic-0.39.jar')
print('\n')

#install bwa
print("\033[1;45m Installing BWA \033[0;0;0m")
os.system("brew install bwa")
print("\033[1;45m BWA installed \033[0;0;0m")
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
print('033[1;45m Use the command "vdb-config -i" to configure your SRA Toolkit \033[0;0;0m')
