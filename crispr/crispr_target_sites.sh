#!/bin/bash

# Define input parameters
output_folder="output"
num_threads=4
combined_file="combined_results.txt"

# Initialize combined file
echo "Directory\tGuide Sequence\tPotential Targets" > $combined_file

# Loop through all directories starting with SRR or ERR
for dir in SRR* ERR*; do
  echo "Processing directory $dir..."

  # Check if directory contains BAM files
  if [ -e "$dir/sorted.bam" ]; then
    # Identify potential targets using Bowtie2
    echo "Identifying potential targets using Bowtie2..."
    bowtie2 --quiet --no-unal --very-sensitive-local --no-hd --no-sq -x bowtie2_path -U $dir/now_trimmed.fq.gz -S $dir/now_mapped.sam && samtools view -S -b $dir/now_mapped.sam > $dir/now_mapped.bam && samtools sort $dir/now_mapped.bam > $dir/now_mapped.sorted.bam && samtools index $dir/now_mapped.sorted.bam
    python -c "import CRISPResso2; CRISPResso2.fastq2Crispresso('$dir', '$dir/$output_folder', '$dir/now_mapped.sorted.bam', '$dir/now_mapped.sorted.bam.bai', '$dir/now_mapped.sam', '$dir/now_trimmed.fq.gz', '', '', 'bowtie2', True, 50, 50, '', '', '', '', '', '', False, False, False, False, False, False, False, False, 50, 50, False, False, False, True, False, False, '', '', '', '', '', False, '', '', '', '', '', False, '', '', '', '', '', '')"

    # Get potential targets
    echo "Getting potential targets..."
    potential_targets=$(cat $dir/$output_folder/CRISPResso2_on_target_summary.txt | awk '{if ($4 >= 50) {print $1}}' | tr '\n' ',')

    # Save results to the combined file
    echo "Saving results to $combined_file..."
    echo -e "$(basename $dir)\t$guide_seq\t$potential_targets\n" >> $combined_file

    # Remove intermediate files
    echo "Removing intermediate files..."
    rm $dir/now_mapped.sam $dir/now_mapped.bam $dir/now_mapped.sorted.bam $dir/now_mapped.sorted.bam.bai
  else
    echo "No BAM files found in directory $dir. Skipping..."
  fi

done

echo "Analysis complete. Results are saved to $combined_file."

