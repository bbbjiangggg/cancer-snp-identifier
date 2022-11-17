#!/usr/bin/env python3

import os
import sys
 
import subprocess

print("\033[1;45m ONLY RUN THIS PROGRAM IN YOUR HOME DIRECTORY  \033[0;0;0m")
print("\033[1;45m to access your home directory use the command: '$ explorer.exe .' only on Ubuntu WSL  \033[0;0;0m")
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')

print("\033[1;45m Installing Python3 and pip3 \033[0;0;0m")
subprocess.call(['sudo', 'apt', 'install', 'python3', 'python3-pip'])
print("\033[1;45m Installing pyutil and pandas \033[0;0;0m")
subprocess.call(['pip3', 'install', 'pandas'])
os.system('pip3 install python-util')

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

#checking for the installation of unzip
print('\033[1;45m Checking for the installation of unzip \033[0;0;0m')
if os.system('unzip --version') == 0:
    print('\033[1;45m unzip is installed \033[0;0;0m')
else:
    print('\033[1;45m unzip is not installed \033[0;0;0m')
    print('\033[1;45m installing unzip \033[0;0;0m')
    os.system('sudo apt-get install unzip -y')
    print('\033[1;45m unzip is installed \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for the installation of sendemail
print('\033[1;45m Checking for the installation of sendemail \033[0;0;0m')
if os.system('sendemail --version') == 0:
    print('\033[1;45m sendemail is installed \033[0;0;0m')
else:
    print('\033[1;45m sendemail is not installed \033[0;0;0m')
    print('\033[1;45m installing sendemail \033[0;0;0m')
    os.system('sudo apt-get install sendemail -y')
    print('\033[1;45m sendemail is installed \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for the installation of xdg-utils
print('\033[1;45m Checking for the installation of xdg-utils \033[0;0;0m')
if os.system('xdg-open --version') == 0:
    print('\033[1;45m xdg-utils is installed \033[0;0;0m')
else:
    print('\033[1;45m xdg-utils is not installed \033[0;0;0m')
    print('\033[1;45m installing xdg-utils \033[0;0;0m')
    os.system('sudo apt-get install xdg-utils -y')
    print('\033[1;45m xdg-utils is installed \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for the installation of tabix
print('\033[1;45m Checking for the installation of tabix \033[0;0;0m')
if os.system('tabix --version') == 0:
    print('\033[1;45m tabix is installed \033[0;0;0m')
else:
    print('\033[1;45m tabix is not installed \033[0;0;0m')
    print('\033[1;45m installing tabix \033[0;0;0m')
    os.system('sudo apt-get install tabix -y')
    print('\033[1;45m tabix is installed \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#downloading SRA Toolkit, you can also try sudo apt install sra-toolkit (do not need to enter export path)
print('\033[1;45m Downloading SRA Toolkit \033[0;0;0m')
os.system('wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz')
os.system('tar -xvzf sratoolkit.current-ubuntu64.tar.gz')
os.system('rm sratoolkit.current-ubuntu64.tar.gz')
os.system('mv sratoolkit.2.10.9-ubuntu64 sratoolkit')
os.system('export PATH=$PATH:~/sratoolkit/bin')
print('\033[1;45m SRA Toolkit is downloaded \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install fastqc
print('\033[1;45m Ready to install FASTQC \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('sudo apt-get install -y fastqc')
print('\033[1;45m FASTQC is installed \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#checking for the installation of java
print('\033[1;45m Checking for the installation of java \033[0;0;0m')
if os.system('java -version') == 0:
    print('\033[1;45m java is installed \033[0;0;0m')
else:
    print('\033[1;45m java is not installed \033[0;0;0m')
    print('\033[1;45m installing java \033[0;0;0m')
    os.system('sudo apt-get install default-jre -y')
    os.system('sudo apt-get install default-jdk -y')
    print('\033[1;45m java is installed \033[0;0;0m')

input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install trimmomatic
print('\033[1;45m Ready to install Trimmomatic \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
os.system('wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip')
os.system('unzip Trimmomatic-0.39.zip')
os.system('rm Trimmomatic-0.39.zip')
os.system('mkdir -p local/bin')
os.system('sudo cp Trimmomatic-0.39/trimmomatic-0.39.jar $HOME/local/bin')
os.system('sudo chmod 755 $HOME/local/bin/trimmomatic-0.39.jar')
input('\033[1;45m Press enter to test Trimmomatic... \033[0;0;0m \n')
print('\033[1;45m Checking if Trimmomatic was properly installed \033[0;0;0m')
os.system('java -jar $HOME/local/bin/trimmomatic-0.39.jar')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install bwa
print('\033[1;45m Installing BWA \033[0;0;0m')
os.system('sudo apt-get install bwa -y')
print('\033[1;45m BWA is installed \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install Bowtie2
print('\033[1;45m Installing Bowtie2 \033[0;0;0m')
os.system('sudo apt-get install bowtie2 -y')
print('\033[1;45m Bowtie2 is installed \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install samtools
print('\033[1;45m Installing samtools \033[0;0;0m')
os.system('sudo apt-get install samtools -y')
print('\033[1;45m samtools is installed \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')

#install bcftools
print('\033[1;45m Installing bcftools \033[0;0;0m')
os.system('sudo apt-get install bcftools -y')
print('\033[1;45m bcftools is installed \033[0;0;0m')
input('\033[1;45m Press enter to continue... \033[0;0;0m \n')
print('\n')
print('\033[1;45m All tools have been installed. \033[0;0;0m')
print('033[1;45m Use the command "vdb-config -i" to configure your SRA Toolkit \033[0;0;0m')

