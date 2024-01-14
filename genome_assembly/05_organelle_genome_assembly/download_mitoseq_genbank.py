#!/bin/env python


from Bio import Entrez
from Bio import SeqIO
import argparse

def download_genbank(id_list, output_genbank):
    Entrez.email = "your_email@example.com" #replace with your email address
    handle = Entrez.efetch(db="nucleotide", id=id_list, rettype="gb", retmode="text")
    records = list(SeqIO.parse(handle, "gb"))
    handle.close()
    # write the records to a GenBank file
    SeqIO.write(records, output_genbank, "genbank")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Download genbank files from NCBI with input lists of IDs")
    parser.add_argument("-i", "--input", required=True, help="File containing the list of IDs, one ID per line", metavar="id_list.txt")
    parser.add_argument("-o", "--output", required=True, help="Path to the output genbank file", metavar="output.gb")

    args = parser.parse_args()
    with open(args.input) as f:
        id_list = f.read().splitlines()
    download_genbank(id_list, args.output)


