import requests
from tqdm import tqdm

# URLs for the human reference genome and annotation files in FASTA and GTF format from Ensembl
fasta_url = 'http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
gtf_url = 'http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz'

# Names of the output files
fasta_output_file = 'GRCh38_reference.fa.gz'
gtf_output_file = 'GRCh38_annotation.gtf.gz'

# Download the files and write them to disk with a progress bar
response = requests.get(fasta_url, stream=True)
with open(fasta_output_file, 'wb') as outfile:
    file_size = int(response.headers.get('content-length', 0))
    chunk_size = 1024
    progress_bar = tqdm(total=file_size, unit='B', unit_scale=True)
    for data in response.iter_content(chunk_size=chunk_size):
        progress_bar.update(len(data))
        outfile.write(data)
    progress_bar.close()

response = requests.get(gtf_url, stream=True)
with open(gtf_output_file, 'wb') as outfile:
    file_size = int(response.headers.get('content-length', 0))
    chunk_size = 1024
    progress_bar = tqdm(total=file_size, unit='B', unit_scale=True)
    for data in response.iter_content(chunk_size=chunk_size):
        progress_bar.update(len(data))
        outfile.write(data)
    progress_bar.close()

