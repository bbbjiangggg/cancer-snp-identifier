import requests

# URL for the human reference genome annotation file in GTF format from Ensembl
url = 'http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz'

# Name of the output file
output_file = 'human_reference_genome.gtf.gz'

# Download the file and write it to disk
response = requests.get(url)
with open(output_file, 'wb') as outfile:
    outfile.write(response.content)

