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

The parameters you can use are the followings:

| Parameter name   | Type | Required | Description |
|---|---|---|---|
|'-p', '--path'| String | True | Path to file which contains fasta sequence or gen id|
|'-t', '--type'| String | True | Searching method. Must be one of these: 'FASTA', 'MIRNA_FASTA', 'GENE_ID'|
|'-s', '--specie'| String | True | 'Target Specie'|
|'-d', '--db'| String | False | miARNs db where the searching will be performed. Must be one  of these: 'RUMIMIR', 'MIRNEST', 'MIRBASE', 'TARBASE' |
|'-e', '--evalue| Float | False | Maximum e-value (default=0.05)|
|'-i', '--percentage' | String | False | Minimun identity percentage (default=40.0)|
|'-z', '--entrezdb'| String | False | Entrez db to perform gene id serachings (default='nucleotide')|

Example:

`python3 Library.py --path './Salmonella.genome.fas' --type MIRNA_FASTA --specie lin --db RUMIMIR --evalue 0.005 --percentage 50.0`
