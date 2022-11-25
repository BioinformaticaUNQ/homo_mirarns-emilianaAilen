from Bio.Blast import NCBIWWW, NCBIXML
import subprocess
import os
import tarfile
import glob

GZ_DATABASE_FILE_NAMES = {
    "TarBase": "tarbase_v8_data",
    "miRNEST": "",
}

#database file names
DATABASES = {
    "Rumir": "",
    "TarBase": "tarbase_v8_data.fasta",
    "miRNEST": "",
    "mirBase": "miRBase.fa",
}


def get_txt_without_full_name(name_extract):
    i = "identifier"
    filePaths = glob.glob(os.path.join("./", f'*{name_extract}*.txt'.format(i)))

    if filePaths:
        return open(filePaths[0], 'r')


#for tarbase and miRNEST databases
def get_fasta_from_gz(filename):
    # open tar.gz file
    file = tarfile.open(f'./{filename}.tar.gz')
  
    # extracting file to a txt
    file.extractall('./')
  
    file.close()

    #File input
    fileInput = get_txt_without_full_name(filename[:7])

    #File output
    fileOutput = open(f'./{filename}.fasta', "w")

    #Start count
    count = 1 

    #Loop through each line in the input file to convert to fasta
    for strLine in fileInput:

        #Strip the endline character from each input line
        strLine = strLine.rstrip("\n")

        #Output the header
        fileOutput.write(">" + str(count) + "\n")
        fileOutput.write(strLine + "\n")

        count = count + 1

    #Close the input and output file
    fileInput.close()
    fileOutput.close()


def lookup_miRNAs(sequence_path, sequence_type, target_specie, selected_db):
    if sequence_type == "FASTA":
        get_counterparts(sequence_path)
    elif sequence_type == "MIRNA_FASTA":
        get_miARNs(sequence_path, selected_db)
    else:
        "You're too young to party"


def get_counterparts(sequence_path, E_VALUE_THRESH=0.04):
    # Capturo la secuencia desde el archivo
    sequence = get_sequence_from_file(sequence_path)

    # Busco homologos
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)

    # Creo un fasta para los homologos
    new_fasta = "> " + \
        blast_record.alignments[0] + "\n" + \
        blast_record.alignments[0].hsps[0].query

    # Lo guardo en un archivo para llamar al otro metodo
    bashCommand = f"echo \"{new_fasta}\" > new_fasta.fa"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)

    # Llamo al otro metodo
    get_miARNs("new_fasta.fa")


def get_miARNs(sequence_path, db):
    if(db == "TarBase" or db == "miRNEST"):
        get_fasta_from_gz(GZ_DATABASE_FILE_NAMES[db])
    ##To do: define place to save databases or common route
    bashCommand = f"blastn -task blastn -query {sequence_path} -db /home/oem/new_folder/{DATABASES[db]} -out blastn.txt"
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

lookup_miRNAs("Salmonella.genome.fas", "MIRNA_FASTA", "asd", "TarBase")

# get_fasta_from_gz('tarbase_v8_data') tested

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
