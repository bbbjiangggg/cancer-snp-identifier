#!/bin/bash

# Define input parameters
guide_seq="CAGAGCCAGGCGCCCGAGGA"
output_folder="output"
num_threads=4
combined_file="combined_results.txt"

# Initialize combined file
echo "Guide Sequence\tTotal Reads\tModified Reads\tModified Frequency" > $combined_file

# Loop through all directories starting with SRR or ERR
for dir in SRR* ERR*; do
  echo "Processing directory $dir..."
  
  # Check if directory contains BAM files
  if [ -e "$dir/sorted.bam" ]; then
    # Mapping reads using Bowtie2
    echo "Mapping reads using Bowtie2..."
    bowtie2 --very-fast-local -x bowtie2_path $dir/sorted.bam -S $dir/now_mapped.sam && samtools view -S -b $dir/now_mapped.sam > $dir/now_mapped.bam && samtools sort $dir/now_mapped.bam > $dir/now_mapped.sorted.bam

    # Summarizing the base calls of aligned reads to the reference sequence (mpileup)
    echo "Summarizing the base calls of aligned reads to the reference sequence (mpileup)..."
    bcftools mpileup -f ref_chrom $dir/now_mapped.sorted.bam | bcftools call -mv -Ob -o $dir/now_mapped.raw.bcf && bcftools view $dir/now_mapped.raw.bcf | vcfutils.pl varFilter - > $dir/now_mapped.var.-final.vcf

    # Run CRISPResso2 analysis
    echo "Running CRISPResso2 analysis..."
    python -c "import CRISPResso2; CRISPResso2.run(alignment_file='$dir/sorted.bam', guide_seq='$guide_seq', guide_seq_pam_right='NGG', output_folder='$dir/$output_folder', n_threads=$num_threads, overwrite=True, quantify_pooled_edits=True, include_unmodified=True, no_plots=True, no_alignment_plots=True, count_mismatches_after_end=False, exclude_bp_from_left=0, exclude_bp_from_right=0, suppress_report=True, do_jit=True)"

    # Combine results to the combined file
    echo "Combining results to $combined_file..."
    echo -e "$(basename $dir)\t$(cat $dir/$output_folder/CRISPResso2_info.json | grep "total_count" | cut -d ":" -f 2 | tr -d ',' | tr -d ' ')\t$(cat $dir/$output_folder/CRISPResso2_info.json | grep "modified_count" | cut -d ":" -f 2 | tr -d ',' | tr -d ' ')\t$(cat $dir/$output_folder/CRISPResso2_info.json | grep "modified_frequency" | cut -d ":" -f 2 | tr -d ',' | tr -d ' ')\n" >> $combined_file

    # Remove intermediate files
    echo "Removing intermediate files..."
    rm $dir/now_mapped.sam $dir/now_mapped.bam $dir/now_mapped.sorted.bam $dir/now_mapped.raw.bcf
  else
    echo "No BAM files found in directory $dir. Skipping..."
  fi

done

echo "Analysis complete. Results are saved to $combined_file."

