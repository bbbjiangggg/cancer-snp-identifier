from Bio import Entrez

def get_snp_accession_number(chromosome, position):
    Entrez.email = "your_email_address@your_domain.com"  # Enter your email address here
    handle = Entrez.esearch(db="snp", term=f"{chromosome}[CHR] AND {position}[CHRPOS]")
    record = Entrez.read(handle)
    handle.close()
    if int(record["Count"]) == 0:
        return None
    else:
        return record["IdList"][0]
