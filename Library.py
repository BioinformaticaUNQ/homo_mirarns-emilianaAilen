from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import subprocess
import os
import argparse
import constants

parser = argparse.ArgumentParser(description='miARNs')

parser.add_argument('-p', '--path', type=str, required=True,
                    help='Path to file which contains fasta sequence or gen id')

parser.add_argument('-t', '--type', type=str,  required=True,
                    choices=['FASTA', 'MIRNA_FASTA', 'GENE_ID'], help='Sequence type or Gen ID')

parser.add_argument('-s', '--specie', type=str,
                    required=True, help='Target Specie')

parser.add_argument('-d', '--db', type=str,  nargs='?',
                    choices=['RUMIMIR', 'MIRNEST', 'MIRBASE', 'TARBASE'], help='miRNA db', default='MIRNEST')

parser.add_argument('-e', '--evalue', type=float,
                    nargs='?', help='Maximum e-value', default=0.05)

parser.add_argument('-i', '--percentage', type=float,
                    nargs='?', help='Minimun identity percentage', default=40.0)

parser.add_argument('-z', '--entrezdb', type=str,
                    nargs='?', help='Entrez db', default='nucleotide')

args = parser.parse_args()

# Entrez config
Entrez.email = 'emiliana.ailen@hotmail.com'


def get_sequence_by_(seq_id, entrez_db):
    handle = Entrez.efetch(db=entrez_db, id=seq_id,
                           rettype="fasta", retmode="text")
    record = handle.read()
    out_handle = open('sequenceFound.fasta', 'w')
    out_handle.write(record.rstrip('\n'))


def get_db_name(sequence_type, selected_db):
    if sequence_type == "FASTA" or sequence_type == "GENE_ID":
        db = constants.PLAIN_DATABASES[selected_db]
    else:
        db = constants.BLAST_DATABASES[selected_db]
    if not db:
        raise Exception(
            f"{selected_db} DB is not available for {sequence_type} searching method")
    return db


def lookup_miRNAs(sequence_path, sequence_type, target_specie, selected_db, evalue, perc_identity, entrez_db):
    db = get_db_name(sequence_type, selected_db)
    match sequence_type:
        case "FASTA":
            get_counterparts(sequence_path, db, target_specie,
                             evalue, perc_identity)
        case "MIRNA_FASTA":
            get_miARNs(sequence_path, db, target_specie, evalue, perc_identity)
        case "GENE_ID":
            get_counterparts_from_gene_id(
                sequence_path, db, target_specie, evalue, perc_identity, entrez_db)
        case _:
            raise Exception("You have to provide a correct sequence type")


def get_counterparts_from_gene_id(gene_id, db, target_specie, evalue, perc_identity, entrez_db):
    get_sequence_by_(gene_id, entrez_db)
    get_counterparts("sequenceFound.fasta", db,
                     target_specie, evalue, perc_identity)


def get_counterparts(sequence_path, db, target_specie, evalue, perc_identity):
    sequence = get_sequence_from_file(sequence_path)
    result_handle = NCBIWWW.qblast(
        "blastn", "nt", sequence, perc_ident=perc_identity)
    blast_record = NCBIXML.read(result_handle)
    alignment_found = get_best_alignment_ID(
        evalue, target_specie, blast_record.alignments)
    get_result_from_DB(db, alignment_found)


def get_result_from_DB(db, gene_id):
    input = open(db, "r")
    out_file = open("result.txt", "w")
    data = input.readlines()
    for line in data:
        splited = line.split(' ')
        mirna_position = 1
        if (db == constants.PLAIN_DATABASES["MIRNEST"]):
            mirna_position = 0
        mirna = splited[mirna_position]
        if line.__contains__(gene_id):
            out_file.write(mirna)


def get_best_alignment_ID(E_VALUE_THRESH, specie, alignments):
    for alignment in alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH and specie in alignment.title:
                return alignment.title.split("|")[1]


def get_miARNs(sequence_path, db, specie, evalue, perc_identity):
    bashCommand = f'blastn -task blastn -query {sequence_path} -db {db} -evalue {evalue} -perc_identity {perc_identity} -outfmt "6 stitle pident evalue" -out blast.txt'
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
              args.db, args.evalue, args.percentage, args.entrezdb)
