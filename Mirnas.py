from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez
import subprocess
import os
import constants


def get_sequence_by_(seq_id, entrez_db, entrezemail):
    # Entrez config
    Entrez.email = entrezemail
    try:
        handle = Entrez.efetch(db=entrez_db, id=seq_id,
                               rettype="fasta", retmode="text")
        record = handle.read()
        out_handle = open('sequenceFound.fasta', 'w')
        out_handle.write(record.rstrip('\n'))
    except IOError as e:
        if e.getcode() == 400:
            print(f"There isn't any result for {seq_id} id")
    except:
        print("Cannot retrieve results from Entrez. Please try again.")
    finally:
        handle.close()
        out_handle.close()


def get_db_name(sequence_type, selected_db):
    if sequence_type == "FASTA" or sequence_type == "GENE_ID":
        db = constants.PLAIN_DATABASES[selected_db]
    else:
        db = constants.BLAST_DATABASES[selected_db]
    if not db:
        raise Exception(
            f"{selected_db} DB is not available for {sequence_type} searching method")
    return db


def lookup_miRNAs(input, sequence_type, target_specie, selected_mirna_db, evalue, perc_identity, entrez_db, output_path, entrezemail, blastdb):
    mirna_db = get_db_name(sequence_type, selected_mirna_db)
    if (sequence_type == 'FASTA' or sequence_type == 'MIRNA_FASTA') and (not os.path.isfile(input)):
        raise ValueError(
            "The input file path is wrong. Please try again with a correct one")
    match sequence_type:
        case "FASTA":
            get_counterparts(input, mirna_db, target_specie,
                             evalue, perc_identity, output_path, blastdb)
        case "MIRNA_FASTA":
            get_miARNs(input, mirna_db, target_specie,
                       evalue, perc_identity, output_path)
        case "GENE_ID":
            get_counterparts_from_gene_id(
                input, mirna_db, target_specie, evalue, perc_identity, entrez_db, entrezemail, output_path, blastdb)
        case _:
            raise ValueError(
                'You have to provide a correct sequence type: FASTA, MIRNA_FASTA or GENE_ID')


def get_counterparts_from_gene_id(gene_id, mirna_db, target_specie, evalue, perc_identity, entrez_db, entrezemail, output_path, blastdb):
    get_sequence_by_(gene_id, entrez_db, entrezemail)
    get_counterparts("sequenceFound.fasta", mirna_db,
                     target_specie, evalue, perc_identity, output_path, blastdb)


def get_counterparts(sequence_path, mirna_db, target_specie, evalue, perc_identity, output_path, blastdb):
    sequence = get_sequence_from_file(sequence_path)
    try:
        result_handle = NCBIWWW.qblast(
            "blastn", blastdb, sequence, perc_ident=perc_identity)
    except:
        raise Exception(
            "The blast searching cannot be performed. Please verify your input parameters and sequence file.")
    blast_record = NCBIXML.read(result_handle)
    alignments = blast_record.alignments
    if len(alignments) == 0:
        raise ValueError("No alignments found")
    ids = get_best_alignment_IDs(
        evalue, target_specie, alignments)
    get_result_from_DB(mirna_db, ids, output_path)


def get_result_from_DB(db, gene_ids, output_path):
    input = open(db, "r")
    out_file = open(output_path, "w")
    data = input.readlines()
    MirnaFound = False
    for line in data:
        if any(word in line for word in gene_ids):
            MirnaFound = True
            splited = line.split()
            mirna_position = 1
            if (db == constants.PLAIN_DATABASES["MIRNEST"]):
                mirna_position = 0
            mirna = splited[mirna_position]
            out_file.write(mirna + "\n")
    input.close()
    if not MirnaFound:
        for id in gene_ids:
            out_file.write(id + "\n")
        out_file.close()
        raise ValueError("no gene id matches the ones in the given database, you can see the ids found in the output file")
    out_file.close()

def get_best_alignment_IDs(E_VALUE_THRESH, specie, alignments):
    gene_ids = [] 
    for alignment in alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH and specie in alignment.title:
                gen_id = alignment.title.split("|")[1]
                if not gene_ids.__contains__(gen_id):
                    gene_ids.append(gen_id)
    if len(gene_ids) == 0:
        raise ValueError("No matching gene Id found")
    return gene_ids


def get_miARNs(sequence_path, db, specie, evalue, perc_identity, output_path):
    bashCommand = f'blastn -task blastn -query {sequence_path} -db {db} -evalue {evalue} -perc_identity {perc_identity} -outfmt "6 stitle pident evalue" -out blast.txt'
    subprocess.run(bashCommand, shell=True)
    blast_result = open("blast.txt", "r")
    hits = blast_result.readlines()
    write_result_from_hits(hits, specie, output_path)
    blast_result.close()


def write_result_from_hits(blast_miARNs_hits, specie, output_path):
    result = open(output_path, "w")
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
        if not data:
            raise ValueError(
                f"There's no sequence on the following file: {file_path}")
    else:
        raise ValueError("The input file path is incorrect, please fix it.")
    return data
