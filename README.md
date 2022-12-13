# TP final Bioinformática

## How to install the cli

1. Update your local system's repository list by entering the following command:
   - `sudo apt update`
2. Install python
   - `sudo apt install python3`
3. Install pip.
   - `sudo apt install python3-pip`
4. Install bio python
   - `pip install biopython`
5. Install blast local -> https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
   - `sudo apt install ncbi-blast+`
6. Clone the github repo
   - `git clone https://github.com/BioinformaticaUNQ/homo_mirarns-emilianaAilen.git`
7. Extract the files from db and db files into the repo root directory

## How to use the cli

1. Open a terminal within the repo root directory
2. Type the following command: `python3 Library.py [PARAMETERS]`
3. After running the script a filed called "result.txt" will be generated in the same directory

The parameters you can use are the followings:

| Parameter name   | Type | Required | Description |
|---|---|---|---|
|'-p', '--path'| String | True | Path to file which contains fasta sequence or gen id|
|'-t', '--type'| String | True | Searching method. Must be one of these: 'FASTA', 'MIRNA_FASTA', 'GENE_ID'|
|'-s', '--specie'| String | True | Target specie |
|'-d', '--db'| String | False | miARNs db where the searching will be performed. Must be one  of these: 'RUMIMIR', 'MIRNEST', 'MIRBASE', 'TARBASE' |
|'-e', '--evalue| Float | False | Maximum e-value (default=0.05)|
|'-i', '--percentage' | String | False | Minimun identity percentage (default=40.0)|
|'-z', '--entrezdb'| String | False | Entrez db to perform gene id serachings (default='nucleotide')|
|'-m', '--entrezemail'| String | True | Email account for entrez configuration|
|'-o', '--output'| String | False | Output file path (default='result.txt')|

Examples:

MIRNA_FASTA case:
`python3 Library.py --path './Salmonella.genome.fas' --type MIRNA_FASTA --specie 'lin' --db RUMIMIR --evalue 0.005 --percentage 50.0 -m 'test@gmail.com'`

GENE_ID case:
`python3 get_mirnas.py --path '241989491' --type GENE_ID --specie 'Oryza sativa' --db MIRNEST --evalue 0.005 --percentage 50.0 -m 'test@gmail.com'`

FASTA case:
`python3 get_mirnas.py --path './secuence.genome.fas' --type FASTA --specie 'Oryza sativa' --db MIRNEST --evalue 0.005 --percentage 50.0 -m 'test@gmail.com'`
