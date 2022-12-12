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


def lookup_miRNAs(sequence_path, sequence_type, target_specie, selected_db, evalue, perc_identity, entrez_db, output_path, entrezemail):
    db = get_db_name(sequence_type, selected_db)
    match sequence_type:
        case "FASTA":
            get_counterparts(sequence_path, db, target_specie,
                             evalue, perc_identity, output_path)
        case "MIRNA_FASTA":
            get_miARNs(sequence_path, db, target_specie, evalue, perc_identity)
        case "GENE_ID":
            get_counterparts_from_gene_id(
                sequence_path, db, target_specie, evalue, perc_identity, entrez_db, entrezemail)
        case _:
            raise Exception("You have to provide a correct sequence type")


def get_counterparts_from_gene_id(gene_id, db, target_specie, evalue, perc_identity, entrez_db, entrezemail):
    get_sequence_by_(gene_id, entrez_db, entrezemail)
    get_counterparts("sequenceFound.fasta", db,
                     target_specie, evalue, perc_identity)


def get_counterparts(sequence_path, db, target_specie, evalue, perc_identity, output_path):
    sequence = get_sequence_from_file(sequence_path)
    result_handle = NCBIWWW.qblast(
        "blastn", "nt", sequence, perc_ident=perc_identity)
    blast_record = NCBIXML.read(result_handle)
    alignment_found = get_best_alignment_ID(
        evalue, target_specie, blast_record.alignments)
    get_result_from_DB(db, alignment_found, output_path)


def get_result_from_DB(db, gene_id, output_path):
    input = open(db, "r")
    out_file = open(output_path, "w")
    data = input.readlines()
    for line in data:
        splited = line.split(' ')
        mirna_position = 1
        if (db == constants.PLAIN_DATABASES["MIRNEST"]):
            mirna_position = 0
        mirna = splited[mirna_position]
        if line.__contains__(str(gene_id)):
            out_file.write(mirna)


def get_best_alignment_ID(E_VALUE_THRESH, specie, alignments):
    for alignment in alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH and specie in alignment.title:
                return alignment.title.split("|")[1]


def get_miARNs(sequence_path, db, specie, evalue, perc_identity, output_path):
    bashCommand = f'blastn -task blastn -query {sequence_path} -db {db} -evalue {evalue} -perc_identity {perc_identity} -outfmt "6 stitle pident evalue" -out {output_path}'
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
