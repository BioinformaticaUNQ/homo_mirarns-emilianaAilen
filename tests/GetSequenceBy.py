from Bio import Entrez
from collections import namedtuple
import unittest
from unittest.mock import MagicMock, patch
import io, os, sys
sys.path.append(r'..')
from Mirnas import get_sequence_by_


class GetBestAlignmentId(unittest.TestCase):

    def setUp(self):
        f = open("myfile.txt", "w", encoding="utf-8")
        f.write("123456")
        f.close()

    @patch('Bio.Entrez.efetch')
    def test_get_best_alignment_id_returns_none_if_the_record_name_does_not_contain_the_target_specie(self, efetch):
        # Efetch returns a TextIOWrapper object, so we need to mock it
        efetch_mocked_response = open("myfile.txt", "r", encoding="utf-8")
        efetch.return_value = efetch_mocked_response

        get_sequence_by_("a_gene_id", "nucleotide", "a_mail@mail.com")

        # Getting the function result
        f = open("sequenceFound.fasta", "r")
        obtained_gene_id = f.read()
        f.close()
        expected_gene_id= "123456"

        self.assertEqual(expected_gene_id, obtained_gene_id)

    def tearDown(self):
        os.remove('myfile.txt')
        os.remove('sequenceFound.fasta')

if __name__ == '__main__':
    unittest.main()
