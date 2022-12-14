from io import StringIO
import unittest
from unittest.mock import  patch
import os
import sys
sys.path.append(r'..')
from Mirnas import get_sequence_by_


class TestGetSequenceBy(unittest.TestCase):

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
        expected_gene_id = "123456"

        self.assertEqual(expected_gene_id, obtained_gene_id)

    @patch('Bio.Entrez.efetch', **{'return_value.raiseError.side_effect': IOError("asdasd", 400)})
    def test_when_get_best_alignment_id_fails_it_prints_an_appropriate_message(self, efetch):
        capturedOutput = StringIO()
        sys.stdout = capturedOutput

        get_sequence_by_("any_gene_id", "any_db", "any_mail@mail.com")

        expectedValue = "Cannot retrieve results from Entrez. Please try again.\n"
        obtained_value = capturedOutput.getvalue()

        self.assertEqual(expectedValue, obtained_value)

    @classmethod
    def tearDownClass(cls):
        sys.stdout = sys.__stdout__
        os.remove('myfile.txt')
        os.remove('sequenceFound.fasta')


if __name__ == '__main__':
    unittest.main()
