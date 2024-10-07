from Bio import Entrez
import primer3
from termcolor import colored
import random
import re
import os
import csv

# Function to wrap text with ANSI codes correctly
def wrap_ansi(text, width):
    ansi_escape = re.compile(r'\x1b\[[0-9;]*[mK]')
    lines = []
    line = ''
    length = 0

    # Split text into ANSI and non-ANSI parts
    parts = ansi_escape.split(text)
    ansi_codes = ansi_escape.findall(text)

    # Reconstruct the text with markers
    tokens = []
    pos = 0
    for part in parts:
        tokens.append(('text', part))
        if pos < len(ansi_codes):
            tokens.append(('ansi', ansi_codes[pos]))
            pos += 1

    # Iterate over tokens to build lines
    for token_type, token_value in tokens:
        if token_type == 'text':
            for char in token_value:
                if char == '\n':
                    lines.append(line)
                    line = ''
                    length = 0
                else:
                    char_length = 1  # Each character counts as one
                    if length + char_length > width:
                        lines.append(line)
                        line = ''
                        length = 0
                    line += char
                    length += char_length
        else:
            # ANSI code, add without affecting length
            line += token_value

    if line:
        lines.append(line)

    return '\n'.join(lines)

# Display banner
print("\n===============================================")
print("Primer Design Using Primer3 - Sorted by Penalty")
print("===============================================\n")

# User input for accession or FASTA sequence
print("Enter the accession number or FASTA sequence (type 'EOF' on a new line to finish):")
lines = []
while True:
    line = input()
    if line.strip() == 'EOF':
        break
    lines.append(line)
accession_or_fasta = '\n'.join(lines).strip()

# User input for primer design parameters
print("\nEnter primer design parameter ranges (leave blank for default suggestions):")

# Enforce primer length between 18-22 bases
primer_length_min_input = input("Minimum primer length (default 18): ")
primer_length_min = int(primer_length_min_input) if primer_length_min_input else 18
if primer_length_min < 18 or primer_length_min > 22:
    print("Minimum primer length adjusted to 18 bases.")
    primer_length_min = 18

primer_length_max_input = input("Maximum primer length (default 22): ")
primer_length_max = int(primer_length_max_input) if primer_length_max_input else 22
if primer_length_max < 18 or primer_length_max > 22:
    print("Maximum primer length adjusted to 22 bases.")
    primer_length_max = 22

# Enforce GC content between 40-60%
gc_content_min_input = input("Minimum GC content (default 40%): ")
gc_content_min = float(gc_content_min_input) if gc_content_min_input else 40
if gc_content_min < 40 or gc_content_min > 60:
    print("Minimum GC content adjusted to 40%.")
    gc_content_min = 40

gc_content_max_input = input("Maximum GC content (default 60%): ")
gc_content_max = float(gc_content_max_input) if gc_content_max_input else 60
if gc_content_max < 40 or gc_content_max > 60:
    print("Maximum GC content adjusted to 60%.")
    gc_content_max = 60

# Enforce Tm between 55-75°C
tm_min_input = input("Minimum melting temperature (Tm) (default 55°C): ")
tm_min = float(tm_min_input) if tm_min_input else 55
if tm_min < 55 or tm_min > 75:
    print("Minimum Tm adjusted to 55°C.")
    tm_min = 55

tm_max_input = input("Maximum melting temperature (Tm) (default 75°C): ")
tm_max = float(tm_max_input) if tm_max_input else 75
if tm_max < 55 or tm_max > 75:
    print("Maximum Tm adjusted to 75°C.")
    tm_max = 75

product_size_min_input = input("Enter minimum product size (default 100): ")
product_size_min = int(product_size_min_input) if product_size_min_input else 100

product_size_max_input = input("Enter maximum product size (default 1000): ")
product_size_max = int(product_size_max_input) if product_size_max_input else 1000

# Process the input
header_line = ''
input_lines = accession_or_fasta.strip().split('\n')
if input_lines[0].startswith('>') or len(input_lines) > 1:
    # Assume it's a FASTA sequence
    header_line = input_lines[0] if input_lines[0].startswith('>') else ''
    sequence = ''.join(input_lines[1:]) if header_line else ''.join(input_lines)
    sequence = sequence.replace(' ', '').replace('\r', '').replace('\n', '')
    # Extract accession number and gene symbol from header if possible
    accession_number = ''
    gene_symbol = ''
    if header_line.startswith('>'):
        # Extract accession number
        parts = header_line[1:].split()
        if parts:
            accession_number = parts[0]
        # Extract gene symbol from parentheses
        match = re.search(r'\(([^)]+)\)', header_line)
        if match:
            gene_symbol = match.group(1)
        else:
            gene_symbol = ''
else:
    # Assume it's an accession number and fetch the sequence
    accession_number = accession_or_fasta
    print("\nFetching sequence from NCBI...")
    Entrez.email = "your.email@example.com"  # Replace with your email
    try:
        # Fetch the FASTA sequence
        handle = Entrez.efetch(db="nucleotide", id=accession_or_fasta, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        fasta_lines = fasta_data.strip().split('\n')
        header_line = fasta_lines[0]
        sequence = ''.join(fasta_lines[1:])
        sequence = sequence.replace(' ', '').replace('\r', '').replace('\n', '')
        # Fetch the gene symbol
        handle = Entrez.esummary(db="nucleotide", id=accession_or_fasta, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        # Try to extract gene symbol from the title
        title = records[0].get('Title', '')
        match = re.search(r'\(([^)]+)\)', title)
        if match:
            gene_symbol = match.group(1)
        else:
            gene_symbol = ''
    except Exception as e:
        print(f"Error fetching sequence: {e}")
        exit()

# Print sequence length for debugging
print(f"\nSequence length: {len(sequence)}")

# Check if sequence length is sufficient
if len(sequence) < product_size_min:
    print("Error: Sequence length is shorter than the minimum product size.")
    exit()

print("\nDesigning primers...")

# Set up primer3 input with rules enforced
# Increase PRIMER_NUM_RETURN to get more primers for selection
primer3_results = primer3.bindings.design_primers(
    {
        'SEQUENCE_ID': 'Target',
        'SEQUENCE_TEMPLATE': sequence,
    },
    {
        'PRIMER_OPT_SIZE': (primer_length_min + primer_length_max) // 2,
        'PRIMER_MIN_SIZE': primer_length_min,
        'PRIMER_MAX_SIZE': primer_length_max,
        'PRIMER_OPT_TM': (tm_min + tm_max) / 2,
        'PRIMER_MIN_TM': tm_min,
        'PRIMER_MAX_TM': tm_max,
        'PRIMER_MIN_GC': gc_content_min,
        'PRIMER_MAX_GC': gc_content_max,
        'PRIMER_PRODUCT_SIZE_RANGE': [[product_size_min, product_size_max]],
        'PRIMER_NUM_RETURN': 20,  # Increased to get more primers for selection
        'PRIMER_PAIR_MAX_COMPL_END': 3.00,
        'PRIMER_MAX_SELF_ANY': 8.00,
        'PRIMER_MAX_SELF_END': 3.00,
        'PRIMER_MAX_POLY_X': 3,
    }
)

# Check if primers were found
num_primers = primer3_results['PRIMER_PAIR_NUM_RETURNED']
if num_primers == 0:
    print("No primer pairs found with the given parameters.")
else:
    # Collect unique, non-overlapping primer pairs with penalty scores
    used_regions = set()
    accepted_primers = []
    for i in range(num_primers):
        left_start, left_len = primer3_results[f'PRIMER_LEFT_{i}']
        right_start, right_len = primer3_results[f'PRIMER_RIGHT_{i}']

        # Create sets of positions for the forward and reverse primers
        left_positions = set(range(left_start, left_start + left_len))
        right_positions = set(range(right_start, right_start + right_len))

        # Check if primers overlap with used regions
        if left_positions.isdisjoint(used_regions) and right_positions.isdisjoint(used_regions):
            # Retrieve penalty scores
            pair_penalty = primer3_results[f'PRIMER_PAIR_{i}_PENALTY']

            # Accept the primer pair
            accepted_primers.append({
                'index': i,
                'pair_penalty': pair_penalty
            })
            # Update used_regions
            used_regions.update(left_positions)
            used_regions.update(right_positions)
        else:
            # Skip this primer pair due to overlap
            continue

    if not accepted_primers:
        print("No unique, non-overlapping primer pairs found with the given parameters.")
    else:
        # Sort accepted primers by total pair penalty (ascending order)
        accepted_primers.sort(key=lambda x: x['pair_penalty'])

        print("\nUnique Primer Pairs Found (sorted by total pair penalty):")

        # Define available colors
        available_colors = ['red', 'green', 'yellow', 'blue', 'magenta', 'cyan']
        color_index = 0

        primer_regions = []
        primer_data_list = []  # For CSV output

        for idx, primer_info in enumerate(accepted_primers):
            i = primer_info['index']
            pair_penalty = primer_info['pair_penalty']
            color = available_colors[color_index % len(available_colors)]
            color_index += 1

            left_primer = primer3_results[f'PRIMER_LEFT_{i}_SEQUENCE']
            right_primer = primer3_results[f'PRIMER_RIGHT_{i}_SEQUENCE']
            left_start, left_len = primer3_results[f'PRIMER_LEFT_{i}']
            right_start, right_len = primer3_results[f'PRIMER_RIGHT_{i}']
            product_size = primer3_results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE']

            left_tm = primer3_results[f'PRIMER_LEFT_{i}_TM']
            right_tm = primer3_results[f'PRIMER_RIGHT_{i}_TM']
            left_gc = primer3_results[f'PRIMER_LEFT_{i}_GC_PERCENT']
            right_gc = primer3_results[f'PRIMER_RIGHT_{i}_GC_PERCENT']

            # Retrieve individual primer penalties
            left_penalty = primer3_results[f'PRIMER_LEFT_{i}_PENALTY']
            right_penalty = primer3_results[f'PRIMER_RIGHT_{i}_PENALTY']

            # Print primer pair information with color
            print(f"\n{colored(f'Primer Pair {idx + 1}:', color)}")
            print(f"Total Pair Penalty: {pair_penalty:.2f}")
            print(colored(f"Forward Primer: {left_primer}", color))
            print(f"  Start: {left_start}, Length: {left_len}")
            print(f"  Tm: {left_tm:.2f}°C, GC%: {left_gc:.2f}%")
            print(f"  Penalty: {left_penalty:.2f}")
            print(colored(f"Reverse Primer: {right_primer}", color))
            print(f"  Start: {right_start}, Length: {right_len}")
            print(f"  Tm: {right_tm:.2f}°C, GC%: {right_gc:.2f}%")
            print(f"  Penalty: {right_penalty:.2f}")
            print(f"Product Size: {product_size} bp")

            # Collect primer regions for visualization
            primer_regions.append({
                'start': left_start,
                'end': left_start + left_len,
                'color': color,
                'pair': idx + 1,
                'type': 'Forward'
            })
            primer_regions.append({
                'start': right_start,
                'end': right_start + right_len,
                'color': color,
                'pair': idx + 1,
                'type': 'Reverse'
            })

            # Collect data for CSV
            primer_data_list.append({
                'Primer Pair': idx + 1,
                'Total Pair Penalty': f"{pair_penalty:.2f}",
                'Forward Primer': left_primer,
                'Forward Start': left_start,
                'Forward Length': left_len,
                'Forward Tm': f"{left_tm:.2f}",
                'Forward GC%': f"{left_gc:.2f}",
                'Forward Penalty': f"{left_penalty:.2f}",
                'Reverse Primer': right_primer,
                'Reverse Start': right_start,
                'Reverse Length': right_len,
                'Reverse Tm': f"{right_tm:.2f}",
                'Reverse GC%': f"{right_gc:.2f}",
                'Reverse Penalty': f"{right_penalty:.2f}",
                'Product Size': product_size
            })

        # Prepare the sequence with colored primers
        seq_chars = list(sequence)

        # For each position, keep track of the colors applied
        position_colors = [None for _ in range(len(seq_chars))]

        # Apply colors to the sequence based on primer regions
        for region in primer_regions:
            for pos in range(region['start'], region['end']):
                if pos < len(seq_chars):
                    position_colors[pos] = region['color']

        # Reconstruct the sequence with colors
        colored_sequence = ''
        for i, base in enumerate(seq_chars):
            if position_colors[i]:
                base_colored = colored(base, position_colors[i])
                colored_sequence += base_colored
            else:
                colored_sequence += base

        # Include the header line
        if header_line:
            colored_sequence = header_line + '\n' + colored_sequence

        # Wrap the sequence for display using the custom function
        wrapped_sequence = wrap_ansi(colored_sequence, width=80)

        # Display the color legend
        print("\nColor Legend:")
        for idx, primer_info in enumerate(accepted_primers):
            color = available_colors[idx % len(available_colors)]
            print(colored(f"Primer Pair {idx + 1}", color))

        # Display the sequence once with all primers highlighted
        print(f"\nSequence with All Primers Highlighted:")
        print(wrapped_sequence)

        # Generate CSV file
        if accession_number and gene_symbol:
            filename = f"{accession_number}_{gene_symbol}_primers.csv"
        elif accession_number:
            filename = f"{accession_number}_primers.csv"
        else:
            filename = "primers.csv"

        fieldnames = [
            'Primer Pair', 'Total Pair Penalty',
            'Forward Primer', 'Forward Start', 'Forward Length',
            'Forward Tm', 'Forward GC%', 'Forward Penalty',
            'Reverse Primer', 'Reverse Start', 'Reverse Length',
            'Reverse Tm', 'Reverse GC%', 'Reverse Penalty',
            'Product Size'
        ]

        with open(filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            for primer_data in primer_data_list:
                writer.writerow(primer_data)

        print(f"\nPrimer information saved to {filename}")
