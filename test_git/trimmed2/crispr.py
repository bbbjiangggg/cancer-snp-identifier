import os

# Define input parameters
bowtie2_path = input("Enter the path to Bowtie2: ")
ref_chrom = input("Enter the path to the reference chromosome: ")

# Loop through all directories starting with SRR or ERR
for dir in os.listdir():
    if dir.startswith("SRR") or dir.startswith("ERR"):
        fq_files = [file for file in os.listdir(dir) if file.endswith("trimmed.fq.gz")]
        if len(fq_files) > 0:
            print(f"Processing directory {dir}...")
            # Map reads using Bowtie2
            print("Mapping reads using Bowtie2...")
            os.system(f"{bowtie2_path} --very-fast-local -x {ref_chrom} {os.path.join(dir, fq_files[0])} -S {os.path.join(dir, 'now_mapped.sam')} && samtools view -S -b {os.path.join(dir, 'now_mapped.sam')} > {os.path.join(dir, 'now_mapped.bam')} && samtools sort {os.path.join(dir, 'now_mapped.bam')} > {os.path.join(dir, 'now_mapped.sorted.bam')}")
            # Summarize the base calls of aligned reads to the reference sequence (mpileup)
            print("Summarizing the base calls of aligned reads to the reference sequence (mpileup)...")
            os.system(f"bcftools mpileup -f {ref_chrom} {os.path.join(dir, 'now_mapped.sorted.bam')} | bcftools call -mv -Ob -o {os.path.join(dir, 'now_mapped.raw.bcf')} && bcftools view {os.path.join(dir, 'now_mapped.raw.bcf')} | vcfutils.pl varFilter - > {os.path.join(dir, 'now_mapped.var.-final.vcf')}")
            # Remove intermediate files
            print("Removing intermediate files...")
            os.remove(os.path.join(dir, 'now_mapped.sam'))
            os.remove(os.path.join(dir, 'now_mapped.bam'))
            os.remove(os.path.join(dir, 'now_mapped.sorted.bam'))
            os.remove(os.path.join(dir, 'now_mapped.raw.bcf'))
        else:
            print(f"No trimmed FASTQ files found in directory {dir}. Skipping...")

