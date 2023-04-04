# SNP-identifier-Python

ABOUT
1) This is an attempt to simplify the installation and operation of Bioinformatics tools that can be used indentify SNPs.
2) As of now, all tool parameters are set for you.
3) This is a work in progress.

DOWNLOADING ALL NECESSARY TOOLS
Downloads

Linux and Windows Subsystem for Linux (WSL)

1) The file "install_linux_prog_v??.py" must be copied and pasted into the terminal home directory.
2) To run the program, run the command "python3 install_linux_prog_v??.py" on the terminal.
3) First set of commands updates the terminal and installs initial tools.
4) Python and Pip versions are then checked and installed if necessary.
5) The SRA Toolkit must be configured with the Quick Configuration Guide (instructions here: https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration).
6) The Java and JDK installation are required for Trimmomatic to work.
8) Our programs use BWA for indexing and Bowtie2 for mapping.

MacOS

1) The file "install_mac_prog_v??.py" must be copied and pasted into the terminal home directory.
2) To run the program, run the command "python3 install_mac_prog_v??.py" on the terminal.
3) First set of commands downloads and installs Homebrew, vim, and wget.
4) Python and Pip versions are then checked.
5) The SRA Toolkit must be configured with the Quick Configuration Guide (instructions here: https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration).
6) The Java and JDK installation are required for Trimmomatic to work.
7) Our programs use BWA for indexing and Bowtie2 for mapping. As of now, we have not been able to install BWA on the newer Macbooks that use an M1 chip.

SRA Clinical Cohort Analysis 

1) After identifying a clinical cohort in the SRA database, download the accession list. NOTE: Our programs use the accession number "SRR" for detection and identification. Make sure "SRR" is not included in the accession list's name.
2) The accession list and untrimmed_analysis_tools_v?? must be in the same directory where you wish to run the analysis.

