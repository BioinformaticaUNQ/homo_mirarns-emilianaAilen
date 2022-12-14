# TP final Bioinformática

## How to install the cli

1. Update your local system's repository list by entering the following command:
   - `sudo apt update`
2. Install python
   - `sudo apt install python3`
3. Install pip.
   - `sudo apt install python3-pip`
4. Install bio python (version 1.79). The script will probably works with later versions, but we recommened installing the metioned one. Moreover, the whole package is not needed. The script only uses `Bio.Blast` and `Bio.Entrez`. So feel free to install only those two modules if you want.
   - `pip install biopython==1.79`
5. Install blast local -> https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
   - `sudo apt install ncbi-blast+`
6. Clone the github repo
   - `git clone https://github.com/BioinformaticaUNQ/homo_mirarns-emilianaAilen.git`
7. Extract the files from **db.zip** and **db2.zip** files into the repo root directory

## How to use the cli

1. Open a terminal within the repo root directory
2. Create a file with your fasta sequence / gene id number.
3. Type the following command: `python3 Library.py [PARAMETERS]`
4. After running the script a filed called "result.txt" will be generated in the same directory (this is the deafult outpout name, but it could be changed)

The parameters you can use are the followings:

| Parameter name   | Type | Required | Description |
|---|---|---|---|
|'-p', '--path'| String | True | Gene ID number or path to file which contains fasta sequence (depending on the serching type)|
|'-t', '--type'| String | True | Searching method. Must be one of these: 'FASTA', 'MIRNA_FASTA', 'GENE_ID'|
|'-s', '--specie'| String | True | Target specie |
|'-d', '--db'| String | False | miARNs db where the searching will be performed. Must be one  of these: 'RUMIMIR', 'MIRNEST', 'MIRBASE', 'TARBASE' |
|'-e', '--evalue| Float | False | Maximum e-value (default=0.05)|
|'-i', '--percentage' | String | False | Minimun identity percentage (default=40.0)|
|'-z', '--entrezdb'| String | False | Entrez db to perform gene id serachings (default='nucleotide')|
|'-m', '--entrezemail'| String | True | Email account for entrez configuration|
|'-o', '--output'| String | False | Output file path (default='result.txt')|
|'-b', '--blastdb'| String | False | Blast nucleotide target database. This parameter is used to perform a blast seraching when the script searching method is "FASTA" or "MIRNA_FASTA" (default='nt')|

## Use Cases

### MIRNA_FASTA case:

Command:

`python3 get_mirnas.py --path './chr1_2303_mature.fasta' --type MIRNA_FASTA --specie 'lin' --db RUMIMIR --evalue 0.005 --percentage 50.0 -m 'test@gmail.com'`

The input file should contain something like:

   <blockquote> > Rum-oar-00023#chr1_2303_mature@@hsa-miR-138-2-3p#Adipose_tissue <br>
   UGGUGUCGUGAUGCUCU</blockquote>

## GENE_ID case:
`python3 get_mirnas.py --path '241989491' --type GENE_ID --specie 'Oryza sativa' --db MIRNEST --evalue 0.005 --percentage 50.0 -m 'test@gmail.com'`

## FASTA case:
`python3 get_mirnas.py --path './secuence.genome.fas' --type FASTA --specie 'Oryza sativa' --db MIRNEST --evalue 0.005 --percentage 50.0 -m 'test@gmail.com'`

The input file should contain something like:
<blockquote> > ADN query <br>
GTGCATTAGACTACATGTTTAGTGATAACAGCCTTGTTACTGCTATTATCATTATTATTATTATAAGATTGGCAGATGGCACCCTTGGTTAGAACATTGGACAAACTGCCTTACAACACTTTCCCTACAATGTTCAAATCCCACTGGTGTCAACTTAACCTTTCATCCTTCTTGGCTCAATAAACATCTGATTGGCCCTTTACCACTAAAATATTTGGGGTTACATGGATTCACTACACAAATTTGTGACCTTAAGCCAATATAAGGAAACATTGTTTCTATTATGTTTATAATTAAATTCGCTGTTGCTGGTATTGTTATTTACTGCTGCAGTCTTCTATAAATTATATTATCTAATAATAATAATAATAATAATAATTAAGATATTAATAATGGTTTCACATTTTGGCACAAAGCCTACAAGTTCGTGGGAGGGGGTAAGTTGATAATATCAACCACAGTGCTCAACTAGTAATTGTTTTATCAACGCTGAAAGGATAAGTCAATCCCAGCAGCATCTGAACTCGGAATGTAAAGTTGGAAGAAATACTGCTAAGCATTTTGTCTGGAACGGTAATGATTCTGCCAGATTGCTGCCTTAATAATAATAATGATGATGATGATGATGATGATGATAATAATAATAAAAAATAATAATAATGATCATAATAACTAAGCTTATATATATATATACATAAATTTATTTATGTATAACTATTTATATTTATATTTAATTATTACAAACTACTTGGATGTAATTAAATTTCTTTTAAAATT</blockquote>

## Run tests

Command:
`python -m unittest discover .`
