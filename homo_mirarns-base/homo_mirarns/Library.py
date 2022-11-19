from Bio.Blast import NCBIWWW, NCBIXML


def lookup_miRNAs(sequence, sequence_type, target_specie):
    if sequence_type == "FATSA":
        get_counterparts(sequence)
    elif sequence_type == "MIRNA_FATSA":
        get_counterparts(sequence)
    else: 
        "You're too young to party"

def get_counterparts(sequence, E_VALUE_THRESH = 0.04):
    result_handle  = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                print("****Alignment****")
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("e value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")   
    return blast_record

sequence_data = ">HO190899.1 Dn0d_1_1654 Dendrobium nobile un-vernalized cDNA library Dendrobium nobile cDNA 5', mRNA sequence  \n GGGGAGAGTAAAAGAATTGAGATTGAGACAAACAGCAGCAGTGGAAGGGGCATGCAGAGGAGCTCTATCAATTGGCAAGTGTGACAACGTCAAAGCTCCTATGTATGCCTCCTCCACTTCTGCTGCTAGCCTCCTATGCTACTGTTTGTTTTGATATTGTCAAGCGGCAACATCAAAGCTCTTAGTTTCTCTCGTTTATCTCTACTGTTAAGCATCTTAATCAGTGTTGGGTATGACTCTCGCTCCTATTTGATTAGATTTGAGAAACTTTTATATCGTTATACATATTATAGAAGTTTTATTTTAAATTTTTTTATAATCATGCGCGAATATAATAAAAGTTTTACCTACTTTTAATTTGTGATAATTTTAATTATCTTCTCTTGCTAGAAACTCTAATTTCATCACTGATCTTGACCTTTTGCTCATAAGCTCTCTGACTCGCTTTAAGGTTTTGTTTGAGACGGTTTTCTTATGCTTCTTAAAGCAGTTTTTTAGATTGAGAAGTTGATTTAAGAAACAGTAAGAATGACATCCCAAAGCCTAAGATTTCGTTTGGGATGGTTTTCTCACTGCTTTCTAGTAAAAAAATTCTGATTCAAGAAATTTTTATGCTGAAAAGCCGTTTTGTAAAGCAATAAGAAAGTTGTCTCAA"

lookup_miRNAs(sequence_data, "MIRNA_FATSA", "asd")


""" 

result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
print('termino')
print(result_handle)

print(sequence_data)

lookup_miRNAs('>asdasd')

    if code.startswith('>'):
        return 'es una secuencia FASTA'
    else:
        return 'es un gene id'
 """

 