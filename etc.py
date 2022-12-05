import time
import tarfile


GZ_DATABASE_FILE_NAMES = {
    "TarBase": "tarbase_v8_data",
    "miRNEST": "",
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
    time.sleep(1)

    # File input
    # get_txt_without_full_name(filename[:7])
    fileInput = open(f"TarBase_v8_download.txt", "r")

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
