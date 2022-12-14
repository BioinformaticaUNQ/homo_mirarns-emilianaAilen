import unittest
import sys
import os
sys.path.append(r'..')
from Mirnas import lookup_miRNAs

"""
    These tests use miRNAs databases, and you need to uncompress db.zip and db2.zip in the root directory in order to be able to run them
"""


class TestGetMirnas(unittest.TestCase):

    def setUp(self):
        f = open("input_sequence.fasta", "w", encoding="utf-8")
        f.write("> Rum-chi-00001#chr1_811#Mammary_tissue\nAAUGGCUCAGGUCAGCUGCCAA")
        f.close()

    def tearDown(self):
        if os.path.isfile("input_sequence.fasta"):
            os.remove("input_sequence.fasta")
        if os.path.isfile("blast.txt"):
            os.remove("blast.txt")
        if os.path.isfile("blast.txt"):
            os.remove("result.txt")

    def test_get_mirarns_returns_the_blast_hits_results_which_met_the_conditions(self):
        lookup_miRNAs(input="input_sequence.fasta", sequence_type="MIRNA_FASTA", selected_mirna_db="RUMIMIR",
                      target_specie="chi", evalue=0.05, perc_identity=40, output_path="result.txt",
                      entrez_db="any", entrezemail="any@mail.com", blastdb="nt")

        f = open("./result.txt", "r")

        expected_result = "Description - Identity percentage - E-value\nRum-chi-00001#chr1_811#Mammary_tissue	100.000	5.17e-07\n"
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)

    def test_get_best_alignment_id_returns_none_if_all_the_aligments_have_an_evalue_above_the_maximum(self):
        lookup_miRNAs(input="input_sequence.fasta", sequence_type="MIRNA_FASTA", selected_mirna_db="RUMIMIR",
                      target_specie="chi", evalue=0.0000000001, perc_identity=40, output_path="result.txt",
                      entrez_db="any", entrezemail="any@mail.com", blastdb="nt")

        f = open("./result.txt", "r")

        expected_result = "Description - Identity percentage - E-value\n"
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)


if __name__ == '__main__':
    unittest.main()
