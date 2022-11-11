from Bio.Blast import NCBIWWW


def hello_world():
    print('hello world')


def lookup_miRNAs(code):
    if code.startswith('>'):
        return 'es una secuencia FASTA'
    else:
        return 'es un gene id'


sequence_data = "›SeqABCD [organism=Mus musculus] [strain=C57BL/6]"

result_handle = NCBIWWW.qblast("blastn", "nt", sequence_data)
print('termino')
print(result_handle)

print(sequence_data)

lookup_miRNAs('>asdasd')
