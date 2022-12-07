from Bio.Blast import NCBIWWW, NCBIXML
from Bio import Entrez, SeqIO
from pyfaidx import Fasta
import subprocess
import os
import tarfile
import glob


# Entrez config
Entrez.email = 'emiliana.ailen@hotmail.com'


GZ_DATABASE_FILE_NAMES = {
    "TarBase": "tarbase_v8_data",
    "miRNEST": "",
}

# database file names
DATABASES = {
    "Rumir": "RumimiR.fasta",
    "TarBase": "tarbase_v8_data.fasta",
    "miRNEST": "",
    "mirBase": "miRBase.fa",
}


def get_txt_without_full_name(name_extract):
    i = "identifier"
    filePaths = glob.glob(os.path.join(
        "./", f'*{name_extract}*.txt'.format(i)))

    if filePaths:
        return open(filePaths[0], 'r')


# for tarbase and miRNEST databases
def get_fasta_from_gz(filename):
    # open tar.gz file
    file = tarfile.open(f'./{filename}.tar.gz')

    # extracting file to a txt
    file.extractall('./')

    file.close()

    # File input
    fileInput = get_txt_without_full_name(filename[:7])

    # File output
    fileOutput = open(f'./{filename}.fasta', "w")

    # Start count
    count = 1

    # Loop through each line in the input file to convert to fasta
    for strLine in fileInput:

        # Strip the endline character from each input line
        strLine = strLine.rstrip("\n")

        # Output the header
        fileOutput.write(">" + str(count) + "\n")
        fileOutput.write(strLine + "\n")

        count = count + 1

    # Close the input and output file
    fileInput.close()
    fileOutput.close()


def get_sequence_by_(seq_id):
    handle = Entrez.efetch(db="nucleotide", id=seq_id,
                           rettype="fasta", retmode="text")
    record = handle.read()
    out_handle = open('sequenceFound.fasta', 'w')
    out_handle.write(record.rstrip('\n'))


def lookup_miRNAs(search_prop, sequence_type, target_specie, selected_db):
    match sequence_type:
        case "FASTA":
            get_counterparts(search_prop, target_specie)
        case "MIRNA_FASTA":
            get_miARNs(search_prop, selected_db)
        case "GENE_ID":
            get_counterparts_from_gene_id(
                search_prop, selected_db, target_specie)
        case _:
            raise Exception("Yo have to provide a correct sequence type")


def get_counterparts_from_gene_id(gene_id, db, target_specie):
    get_sequence_by_(gene_id)
    get_counterparts("sequenceFound.fasta", db, target_specie)


def get_counterparts(sequence_path, db, specie, E_VALUE_THRESH=0.05):
    sequence = get_sequence_from_file(sequence_path)

    # Busco homologos
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)
    # obtengo gene id del que mas corresponda
    alignment_found = get_best_alignment_ID(
        E_VALUE_THRESH, specie, blast_record.alignments)
    # voy a la base de mirnas con ese gene id
    get_result_from_DB(db, alignment_found)


def get_result_from_DB(db, gene_id):
    if (db == "TarBase" or db == "miRNEST"):
        get_fasta_from_gz(GZ_DATABASE_FILE_NAMES[db])
    data = Fasta(DATABASES[db])
    keys = list(data.keys())
    out_file = open("result.txt", "w")
    for key in keys:
        # To do: define strategies to get mirna by db
        sequence = data[key][:].seq.split()
        mirna = sequence[2]
        seq_gene_Id = sequence[0]
        if gene_id == str(seq_gene_Id):
            out_file.write(mirna + "\n")


def get_best_alignment_ID(E_VALUE_THRESH, specie, alignments):
    for alignment in alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH and specie in alignment.title:
                print(alignment.title)
                return alignment.title.split("|")[1]


def get_miARNs(sequence_path, db):
    if (db == "TarBase" or db == "miRNEST"):
        get_fasta_from_gz(GZ_DATABASE_FILE_NAMES[db])
    bashCommand = f"blastn -task blastn -query {sequence_path} -db  C:/Users/emili/Documents/homo_mirarns-emilianaAilen/{DATABASES[db]} -out blastn.txt"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()


def get_sequence_from_file(file_path):
    if os.path.isfile(file_path):
        text_file = open(file_path, "r")
        data = text_file.read()
        text_file.close()
    return data


sequence_data = ">HO190899.1 Dn0d_1_1654 Dendrobium nobile un-vernalized cDNA library Dendrobium nobile cDNA 5', mRNA sequence  \n GGGGAGAGTAAAAGAATTGAGATTGAGACAAACAGCAGCAGTGGAAGGGGCATGCAGAGGAGCTCTATCAATTGGCAAGTGTGACAACGTCAAAGCTCCTATGTATGCCTCCTCCACTTCTGCTGCTAGCCTCCTATGCTACTGTTTGTTTTGATATTGTCAAGCGGCAACATCAAAGCTCTTAGTTTCTCTCGTTTATCTCTACTGTTAAGCATCTTAATCAGTGTTGGGTATGACTCTCGCTCCTATTTGATTAGATTTGAGAAACTTTTATATCGTTATACATATTATAGAAGTTTTATTTTAAATTTTTTTATAATCATGCGCGAATATAATAAAAGTTTTACCTACTTTTAATTTGTGATAATTTTAATTATCTTCTCTTGCTAGAAACTCTAATTTCATCACTGATCTTGACCTTTTGCTCATAAGCTCTCTGACTCGCTTTAAGGTTTTGTTTGAGACGGTTTTCTTATGCTTCTTAAAGCAGTTTTTTAGATTGAGAAGTTGATTTAAGAAACAGTAAGAATGACATCCCAAAGCCTAAGATTTCGTTTGGGATGGTTTTCTCACTGCTTTCTAGTAAAAAAATTCTGATTCAAGAAATTTTTATGCTGAAAAGCCGTTTTGTAAAGCAATAAGAAAGTTGTCTCAA"

example_ = ">47 \n NNNNNNNNNNNNNNNNNNNNNNNNGNNNNNCGGGTNNCCNNNNGCTGGNTGNNNTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGAAGGGGGGGGATTCCTACTGGCCAGGGTGGCTAATACCGGGTAACGTCGNNNGANCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCACATGTGCCCGGATGGGATTAGCTTGTTGGTGAGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCNCACTGGAACTGAGACACGGTCCCGACTCCTACGGGAGGCANCAGTGGGGAATATTGCTCTTGGGCGCAAGCCTGATGCAGCCATGCCGCGGGTATGAGGAAGGCCTTCGGTTTGTAAAGTACTTTCTCCGGGGAGGAAGGNGTNGTGGTGAATAACCGCTACANTTGANNCTNCCCGCNNAANAACCACCNGNTAACTCCNTGCNNNNNGCCGCGGTAATACGGANGGTGCAAGNGTTAATCGNANTTACTGNNTGTTGAGCGCACGNNGGCGGCCTGTCNNNTCTNATGTGAGATCCCCGGGCTCNCCCTGNNACCTGCATTCGNNNNNTGNNANGCTNGANTCTTGNNNNGNNGNGGNAGNAATTCCNNGTGTNNCGNNGAAATGCNNANAGATCTGNANANANNACNGGNGNCCAANGNNGNCCCCTGNTCTCNGACTGACGCNNGAGTGCTGAANNGTGNAGAGCGNACAGGATTANANNNCCNGNTAGNCCGNCNCCNCACACCNATGTCTACATGNGAGGNTNNNGNNNNNTGNGGCNNNNCNNTCCNNNAGCTNANGNNGTTNAANTANATCNNNCTNNNNCNNGCNNGGGGNCANNANGGGNNNAAANNTNNNNATNAATNTGACGGANNNNNCNNNNNNNNCNNNNNNANCATGNGGATTNANNNTNNNTNNNNNCNNNNNANAACCNNANNNNNNNNNNNNNNTNNNNNNNANNNTNNNNNNNTNNNNNNGNNNNCNNNNNNNNACTNNNNNNNCNNNNNNNNCANNGNNNNNNNNNNNANGNTNNNNNTGNNNNAANNNTGNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTNNNNNNNTNNNNNNNTNNNNTNNNNNNNNNNNNNNNNNNNNNNNTNNNNNNNNTCNNNNNN"

# lookup_miRNAs("37222967", "GENE_ID", "Paphiopedilum", "TarBase")
get_result_from_DB("TarBase", "Y54E10A.16")

# PRUEBAS
"""for alignment in blast_record.alignments:
print(alignment.title)
for hsp in alignment.hsps:
    if hsp.expect < E_VALUE_THRESH:
        print("****Alignment****")
        print("sequence:", alignment.title)
        print("length:", alignment.length)
        print("e value:", hsp.expect)
        print(hsp.query[0:75] + "...")
        print(hsp.match[0:75] + "...")
        print(hsp.sbjct[0:75] + "...")"""

"""# seq_id = '37222967'"""
"""get_result_from_DB("TarBase", "Y54E10A.16")"""
