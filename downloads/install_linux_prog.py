#!/usr/bin/env python3

import os
import sys
 
import subprocess

print("\033[1;45m ONLY RUN THIS PROGRAM IN YOUR HOME DIRECTORY  \033[0;0;0m")
print("\033[1;45m to access your home directory use the command: '$ explorer.exe .' only on Ubuntu WSL  \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')

print("\033[1;45m Installing Python3 and pip3 \033[0;0;0m")
subprocess.call(['sudo', 'apt', 'install', 'python3', 'python3-pip'])

import pip

#Update terminal and install necessaries 
print("\033[1;45m Updating terminal  \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('sudo apt-get update -y')
os.system('sudo apt-get upgrade -y')
print('\n')

#checking for the installation of wget
print('\033[1;45m Checking for the installation of wget \033[0;0;0m')
if os.system('wget --version') == 0:
    print('\033[1;45m wget is installed \033[0;0;0m')
else:
    print('\033[1;45m wget is not installed \033[0;0;0m')
    print('\033[1;45m installing wget \033[0;0;0m')
    os.system('sudo apt-get install wget -y')
    print('\033[1;45m wget is installed \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

print('\033[1;45m Checking for the installation of unzip \033[0;0;0m')
os.system('sudo apt-get install unzip')
os.system('sudo apt-get install libio-socket-ssl-perl libnet-ssleay-perl sendemail')
os.system('sudo apt-get install --reinstall xdg-utils')
os.system('sudo apt-get install tabix')


#Check Python version
if sys.version_info.major == 3:
    print('\033[1;45m Python3 is installed. \033[0;0;0m')
else:
    print('\033[0;101m You need to install a current version of Python3 \033[0;0;0m')

#Check Pip version
pip_versn = pip.__version__
print('\033[1;45m Your Pip version is ' + pip_versn + '\033[0;0;0m')
print('\033[1;45m Note, if you do not have Pip installed, install it using the command: sudo apt install python3-pip \033[0;0;0m')
os.system('pip3 install python-util')

#downloading SRA Toolkit, you can also try sudo apt install sra-toolkit (do not needt to enter export path)
print('\033[1;45m Downloading SRA Toolkit \033[0;0;0m')
os.system('wget --output-document sratoolkit.tar.gz http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz')

#extract tar file contents
os.system('tar -vxzf sratoolkit.tar.gz')
print('\033[1;45m Extracted SRA Toolkit files \033[0;0;0m')

#remove SRA Toolkit tar file
os.system('rm sratoolkit.tar.gz')

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
print('\033[1;45m 2) Copy and paste the command "vim ~/.bashrc" into the TW and press "Enter" \033[0;0;0m')
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
os.system('sudo apt-get install -y fastqc')
print('\n')

#install java
print('\033[1;45m Ready to install Java \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('sudo apt install default-jre')
print('\n')
print('\033[1;45m Ready to install Java SE Development Kit (JDK) \033[0;0;0m')
print('\033[1;45m This will open a webpage where you can download the approprite JDK \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('xdg-open https://www.oracle.com/java/technologies/downloads/')
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
os.system('sudo apt-get install -y bwa')
print('\n')
print('\033[1;45m Installing Bowtie2 \033[0;0;0m')
os.system('sudo apt-get install -y bowtie2')
print('\n')
print('\033[1;45m Installing Samtools \033[0;0;0m')
os.system('sudo apt-get install -y samtools')
print('\n')
print('\033[1;45m Installing BCFtools \033[0;0;0m')
os.system('sudo apt-get install -y bcftools')

print('\n')
print('\033[1;45m All tools have been installed \033[0;0;0m')

