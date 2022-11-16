# TP final Bioinformática

este es el readme de la librería.

1. install pip. 
    - sudo apt install python3-pip
2. install bio python
    - pip install biopython
3. install blast local -> https://ftp.ncbi.nlm.nih.gov/blast/executables/LATEST/
    - sudo apt install ncbi-blast+
4. install databases
    - perl /bin/update_blastdb --passive --decompress 16S_ribosomal_RNA
5. create a databases folder
    - mkdir $HOME/blastdb
    - export BLASTDB=$HOME/blastdb
