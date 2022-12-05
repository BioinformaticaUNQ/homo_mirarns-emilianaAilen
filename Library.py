from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import subprocess
import os
import argparse

parser = argparse.ArgumentParser(description='miARNs')

parser.add_argument('-p', '--path', type=str, required=True,
                    help='Path to file which contains fasta sequence or gen id')

parser.add_argument('-t', '--type', type=str,  required=True,
                    choices=['FASTA', 'MIRNA_FASTA', 'GENE_ID'], help='Sequence type or Gen ID')

parser.add_argument('-s', '--specie', type=str,
                    required=True, help='Target Specie')

parser.add_argument('-d', '--db', type=str,  nargs='?',
                    choices=['RUMIMIR', 'MIRNEST', 'MIRBASE'], help='miRNA db', default='MIRBASE')

parser.add_argument('-e', '--evalue', type=float,
                    nargs='?', help='Maximum e-value', default=0.05)

parser.add_argument('-i', '--percentage', type=float,
                    nargs='?', help='Minimun identity percentage', default=40.0)

args = parser.parse_args()

# Entrez config
Entrez.email = 'emiliana.ailen@hotmail.com'


# database file names
DATABASES = {
    "RUMIMIR": "/db/rumimir",
    "TarBase": "",
    "MIRNEST": "/db/mirnest",
    "MIRBASE": "/db/mirbase",
}


def get_sequence_by_(seq_id):
    handle = Entrez.efetch(db="nucleotide", id=seq_id,
                           rettype="fasta", retmode="text")
    record = handle.read()
    out_handle = open('sequenceFound.fasta', 'w')
    out_handle.write(record.rstrip('\n'))


def lookup_miRNAs(sequence_path, sequence_type, target_specie, selected_db, evalue, perc_identity):
    db = DATABASES[selected_db]
    match sequence_type:
        case "FASTA":
            get_counterparts(sequence_path, db)
        case "MIRNA_FASTA":
            get_miARNs(sequence_path, db, target_specie, evalue, perc_identity)
        case "GENE_ID":
            get_counterparts_from_gene_id(sequence_path, db)
        case _:
            raise Exception("Yo have to provide a correct sequence type")


def get_counterparts_from_gene_id(gene_id, db):
    get_sequence_by_(gene_id)
    get_counterparts("sequenceFound.fasta", db)


def get_counterparts(sequence_path, db, E_VALUE_THRESH=0.04):
    sequence = get_sequence_from_file(sequence_path)

    # Busco homologos
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)

    # Creo un fasta para los homologos
    fileOutput = open('./blast_out.fasta', "w")
    alignment = blast_record.alignments[0]
    hsp = alignment.hsps[0]
    first_organism = blast_record.descriptions[0]
    fileOutput.write(">" + str(alignment.title) + "\n")
    fileOutput.write(str(hsp.query) + "\n")

    fileOutput.close()

    # Llamo al otro metodo
    get_miARNs('blast_out.fasta', db)


def get_miARNs(sequence_path, db, specie, evalue, perc_identity):
    print(db)
    bashCommand = f'blastn -task blastn -query {sequence_path} -db {db} -evalue {evalue} -perc_identity {perc_identity} -outfmt "6 sseqid pident evalue" -out blast.txt'
    print(bashCommand)
    subprocess.run(bashCommand, shell=True)
    blast_result = open("blast.txt", "r")
    hits = blast_result.readlines()
    write_result_from_hits(hits, specie)
    blast_result.close()


def write_result_from_hits(blast_miARNs_hits, specie):
    result = open("result.txt", "w")
    result.write("Description - Identity percentage - E-value\n")
    for hit in blast_miARNs_hits:
        if hit.__contains__(specie):
            result.write(hit)
    result.close()


def get_sequence_from_file(file_path):
    if os.path.isfile(file_path):
        text_file = open(file_path, "r")
        data = text_file.read()
        text_file.close()
    return data


lookup_miRNAs(args.path, args.type, args.specie,
              args.db, args.evalue, args.percentage)
