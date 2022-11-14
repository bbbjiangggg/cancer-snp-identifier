#!/usr/bin/env python3

import os
import sys
import subprocess

print("\033[1;45m ONLY RUN THIS PROGRAM IN YOUR HOME DIRECTORY  \033[0;0;0m")
print("\033[1;45m to access your home directory open Finder and then click on Go and Home.  \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')

print("\033[1;45m Installing Python3 and pip3 \033[0;0;0m")
subprocess.call(['sudo', 'apt', 'install', 'python3', 'python3-pip'])

import pip 
#Install necessaries 
print("\033[1;45m Install Homebrew, vim, and wget \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"')
os.system('brew install wget')
os.system('brew install vim')
os.system('brew install sendemail')
os.system('brew install gzip')
os.system('brew install python')
os.system('curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py')
os.system('python3 get-pip.py')
os.system('brew install coreutils')

#Check Python version
if sys.version_info.major == 3:
    print('\033[1;45m Python3 is installed. \033[0;0;0m')
else:
    print('\033[0;101m You need to install a current version of Python3 \033[0;0;0m')

#Check Pip version
pip_versn = pip.__version__
print('\033[1;45m Your Pip version is ' + pip_versn + '\033[0;0;0m')
print('\033[1;45m Note, if you do not have Pip installed, install it using the command: curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py \033[0;0;0m')
os.system('pip3 install python-util')


#downloading SRA Toolkit
print('\033[1;45m Downloading SRA Toolkit \033[0;0;0m')
os.system('curl -OL --output sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz')

#extract tar file contents
os.system('gunzip sratoolkit.current-mac64.tar.gz')
os.system('tar -vxzf sratoolkit.current-mac64.tar')
print('\033[1;45m Extracted SRA Toolkit files \033[0;0;0m')

#remove SRA Toolkit tar file
os.system('rm sratoolkit.current-mac64.tar')
print('\n')

#show files in the current directory
print('\033[1;45m These are the files in the current directory: \033[0;0;0m')
os.system('ls')

#add export path to the $path directory
sra_tool = input('\033[1;45m Copy and paste the file name of the sratoolkit: \033[0;0;0m \n')
print('\033[1;45m This is your export path: \033[0;0;0m export PATH=$PATH:$PWD/' + sra_tool + '/bin \n')
print('\033[1;45m The following 6 steps must be followed carefully: \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\033[1;45m 1) Open another terminal window (TW)\033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\033[1;45m 2) Copy and paste the command "sudo vim .zshrc" into the TW and press "Enter" \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\033[1;45m 3) Enter "I" on the TW and press "Enter" \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\033[1;45m 4) Copy and paste the following path below the last line of the TW: \033[0;0;0m export PATH=$PATH:$PWD/' + sra_tool + '/bin')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\033[1;45m 5) Press the following keys on the TW in the correct order: "Esc", ":", "w", "q", "Enter" \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\033[1;45m 6) If Vim closed correctly, you may now close the new terminal window \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')

#install fastqc
print('\033[1;45m Ready to install FASTQC \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('brew install fastqc')
print('\n')

#install java
print('\033[1;45m Ready to install Java (JRE) \033[0;0;0m')
print('\033[1;45m This will open a webpage where you can download the approprite JRE \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m')
os.system('open https://java.com/en/download/')
input('\033[1;45m Press enter after JRE is installed \033[0;0;0m \n')
print('\033[1;45m Ready to install Java SE Development Kit (JDK) \033[0;0;0m')
print('\033[1;45m This will open a webpage where you can download the approprite JDK \033[0;0;0m')
print('\033[1;45m For older Macbooks use the x64 DMG Installer \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('open https://www.oracle.com/java/technologies/downloads/')
input('\033[1;45m Press enter after JDK is installed to test Java... \033[0;0;0m \n')
print('\033[1;45m Testing Java \033[0;0;0m')
os.system('java -version')
print('\n')

#install trimmomatic
print('\033[1;45m Ready to install Trimmomatic \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip')
os.system('unzip Trimmomatic-0.39.zip')
os.system('rm Trimmomatic-0.39.zip')
os.system('mkdir -p local/bin')
os.system('sudo cp Trimmomatic-0.39/trimmomatic-0.39.jar $HOME/local/bin')
input('\033[1;45m Press enter to test Trimmomatic... \033[0;0;0m \n')
print('\033[1;45m Checking if Trimmomatic was properly installed \033[0;0;0m')
os.system('java -jar $HOME/local/bin/trimmomatic-0.39.jar')
print('\n')

#install other programs
print('\033[1;45m Installing BWA \033[0;0;0m')
os.system('brew install bwa')
print('\n')
print('\033[1;45m Installing Bowtie2 \033[0;0;0m')
os.system('brew install bowtie2')
print('\n')
print('\033[1;45m Installing Samtools \033[0;0;0m')
os.system('brew install samtools')
print('\n')
print('\033[1;45m Installing BCFtools \033[0;0;0m')
os.system('brew install bcftools')

print('\n')
print('\033[1;45m All tools have been installed \033[0;0;0m')

