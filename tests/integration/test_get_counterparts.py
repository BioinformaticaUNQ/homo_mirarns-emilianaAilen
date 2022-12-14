from collections import namedtuple
import unittest
import sys
import os
from unittest.mock import patch
sys.path.append(r'..')
from Mirnas import lookup_miRNAs

"""
    These tests use miRNAs databases, and you need to uncompress db.zip and db2.zip in the root directory in order to be able to run them
"""


class TestGetCounterparts(unittest.TestCase):

    def setUp(self):
        f = open("input_sequence.fasta", "w", encoding="utf-8")
        f.write("> A_MOCKED_FASTA_SEQUENCE\n AAUGGCUCAGGUCAGCUGCCAA")
        f.close()

    def tearDown(self):
        if os.path.isfile("input_sequence.fasta"):
            os.remove("input_sequence.fasta")
        if os.path.isfile("blast.txt"):
            os.remove("blast.txt")
        if os.path.isfile("blast.txt"):
            os.remove("result.txt")

    @patch('Bio.Blast.NCBIWWW.qblast')
    @patch('Bio.Blast.NCBIXML.read')
    def test_get_counterparts_returns_the_mirnas_sequences_which_met_the_conditions(self, read, qblast):
        qblast.return_value = "a_blast record"
        read.return_value = namedtuple('Record', "alignments")(
            alignments=[namedtuple('Record', "hsps title")(hsps=[
                namedtuple('HSP', "expect")(expect=0.01),
                namedtuple('HSP',  "expect")(expect=0.02)], title="gi|188432925| Rangpur lime root")
            ])

        lookup_miRNAs(input="input_sequence.fasta", sequence_type="FASTA", selected_mirna_db="MIRNEST",
                      target_specie="Rangpur lime root", evalue=0.05, perc_identity=40, output_path="result.txt",
                      entrez_db="any", entrezemail="any@mail.com", blastdb="nt")

        f = open("./result.txt", "r")

        expected_result = "TCTTCCCTACTCCACCCATGC\nTCTTCCCTACTCCACCCAT\n"
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)

    @patch('Bio.Blast.NCBIWWW.qblast')
    @patch('Bio.Blast.NCBIXML.read')
    def test_get_counterparts_returns_none_if_all_the_aligments_have_an_evalue_above_the_maximum(self, read, qblast):
        qblast.return_value = "a_blast record"
        read.return_value = namedtuple('Record', "alignments")(
            alignments=[namedtuple('Record', "hsps title")(hsps=[
                namedtuple('HSP', "expect")(expect=0.10),
                namedtuple('HSP',  "expect")(expect=0.10)], title="gi|188432925| Rangpur lime root")
            ])

        lookup_miRNAs(input="input_sequence.fasta", sequence_type="FASTA", selected_mirna_db="MIRNEST",
                      target_specie="Rangpur lime root", evalue=0.05, perc_identity=40, output_path="result.txt",
                      entrez_db="any", entrezemail="any@mail.com", blastdb="nt")

        f = open("./result.txt", "r")

        expected_result = ""
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)


if __name__ == '__main__':
    unittest.main()
