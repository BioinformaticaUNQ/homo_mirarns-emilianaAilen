import argparse
from Mirnas import lookup_miRNAs

parser = argparse.ArgumentParser(description='miARNs')

parser.add_argument('-p', '--path', type=str, required=True,
                    help='Path to file which contains fasta sequence or gen id')

parser.add_argument('-t', '--type', type=str,  required=True,
                    choices=['FASTA', 'MIRNA_FASTA', 'GENE_ID'], help='Sequence type or Gen ID')

parser.add_argument('-s', '--specie', type=str,
                    required=True, help='Target Specie')

parser.add_argument('-m', '--entrezemail', type=str,
                    nargs='?', required=True, help='Entrez account email')

parser.add_argument('-d', '--db', type=str,  nargs='?',
                    choices=['RUMIMIR', 'MIRNEST', 'MIRBASE', 'TARBASE'], help='miRNA db', default='MIRNEST')

parser.add_argument('-e', '--evalue', type=float,
                    nargs='?', help='Maximum e-value', default=0.05)

parser.add_argument('-i', '--percentage', type=float,
                    nargs='?', help='Minimun identity percentage', default=40.0)

parser.add_argument('-z', '--entrezdb', type=str,
                    nargs='?', help='Entrez db', default='nucleotide')

parser.add_argument('-o', '--output', type=str,
                    nargs='?', help='Output file path', default='result.txt')

args = parser.parse_args()

lookup_miRNAs(args.path, args.type, args.specie,
              args.db, args.evalue, args.percentage, args.entrezdb, args.output, args.entrezemail)
