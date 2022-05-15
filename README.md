# SRA Analysis Tools

ABOUT
1) This is an attempt to simplify the installation and operation of Bioinformatics tools that can be used for DNA analysis.
2) As of now, all tool parameters are set for you.
3) This is a work in progress.

DOWNLOADING ALL NECESSARY TOOLS
Downloads

Linux and Windows Subsystem for Linux (WSL)

1) The file "install_linux_prog.py" must be copied and pasted into the terminal home directory.
2) To run the program, run the command "python3 install_linux_prog.py" on the terminal.
3) First set of commands updates the terminal and installs initial tools.
4) Python and Pip versions are then checked.
5) To install Python and/or Pip:
    a) For Python, visit website: https://www.python.org/
    b) For Pip, run the following terminal command: $ sudo apt install pip
6) The SRA Toolkit must be configured with the Quick Configuration Guide (instructions included in program).
7) The Java and JDK installation are required for Trimmomatic to work.
8) Our programs use BWA for indexing and Bowtie2 for mapping.

MacOS

1) The file "install_mac_prog.py" must be copied and pasted into the terminal home directory.
2) To run the program, run the command "python3 install_mac_prog.py" on the terminal.
3) First set of commands downloads and installs Homebrew, vim, and wget.
4) Python and Pip versions are then checked.
5) To install Python and/or Pip:
    a) For Python, visit website: https://www.python.org/
    b) For Pip, run the following terminal command: $ curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py
6) The SRA Toolkit must be configured with the Quick Configuration Guide (instructions included in program).
7) The Java and JDK installation are required for Trimmomatic to work.
8) Our programs use BWA for indexing and Bowtie2 for mapping. As of now, we have not been able to install BWA on the newer Macbooks that use an M1 chip.

REGULAR 

