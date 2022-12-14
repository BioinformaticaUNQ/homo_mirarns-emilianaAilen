import unittest
import sys
import os
from unittest.mock import patch
sys.path.append(r'..')
from Mirnas import lookup_miRNAs


class TestLookupMirnas(unittest.TestCase):
    def setUp(self):
        # create an empty input file
        f = open("any_input", "w")
        f.close()

    def tearDown(self):
        if os.path.isfile("any_input"):
            os.remove("any_input")

    @patch('Mirnas.get_db_name')
    def test_lookup_mirnas_raises_an_exception_if_the_searching_method_is_wrong(self, get_db_name):
        non_existent_seraching_type = "none"
        get_db_name.return_value = "any_path_to_mirna_db"

        with self.assertRaises(ValueError) as context:
            lookup_miRNAs("any_input", non_existent_seraching_type, "any_target_specie", "any_mirna_db", "any_evalue",
                          "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb")

        expected_exception_message = 'You have to provide a correct sequence type: FASTA, MIRNA_FASTA or GENE_ID'
        obtained_exception_message = str(context.exception)

        self.assertEqual(expected_exception_message,
                         obtained_exception_message)

    @patch('Mirnas.get_db_name')
    def test_lookup_mirnas_raises_an_exception_if_the_input_file_does_not_exists_and_the_serching_method_is_fasta(self, get_db_name):
        get_db_name.return_value = "any_path_to_mirna_db"
        non_existent_input_file = "none"
        seraching_method = "FASTA"

        with self.assertRaises(ValueError) as context:
            lookup_miRNAs(non_existent_input_file, seraching_method, "any_target_specie", "any_mirna_db", "any_evalue",
                          "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb")

        expected_exception_message = 'The input file path is wrong. Please try again with a correct one'
        obtained_exception_message = str(context.exception)

        self.assertEqual(expected_exception_message,
                         obtained_exception_message)

    @patch('Mirnas.get_db_name')
    def test_lookup_mirnas_raises_an_exception_if_the_input_file_does_not_exists_and_the_serching_method_is_mirna_fasta(self, get_db_name):
        get_db_name.return_value = "any_path_to_mirna_db"
        non_existent_input_file = "none"
        seraching_method = "MIRNA_FASTA"

        with self.assertRaises(ValueError) as context:
            self.assertRaises(lookup_miRNAs(non_existent_input_file, seraching_method, "any_target_specie", "any_mirna_db", "any_evalue",
                                            "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb"))

        expected_exception_message = 'The input file path is wrong. Please try again with a correct one'
        obtained_exception_message = str(context.exception)

        self.assertEqual(expected_exception_message,
                         obtained_exception_message)

    @patch('Mirnas.get_counterparts')
    @patch('Mirnas.get_db_name')
    def test_get_counterparts_is_called_with_arguments_when_fasta_serching_type_is_selected(self, get_db_name, get_counterparts):
        get_db_name.return_value = "any_path_to_mirna_db"
        serching_type = 'FASTA'

        lookup_miRNAs("any_input", serching_type, "any_target_specie", "any_mirna_db", "any_evalue",
                      "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb")

        get_counterparts.assert_called_with('any_input', "any_path_to_mirna_db", 'any_target_specie',
                                            'any_evalue', 'any_perc_identity', 'any_output_path', 'any_blastdb')

    @patch('Mirnas.get_miARNs')
    @patch('Mirnas.get_db_name')
    def test_get_mirnas_is_called_with_arguments_when_mirna_fasta_serching_type_is_selected(self, get_db_name, get_miARNs):
        get_db_name.return_value = "any_path_to_mirna_db"
        serching_type = 'MIRNA_FASTA'

        lookup_miRNAs("any_input", serching_type, "any_target_specie", "any_mirna_db", "any_evalue",
                      "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb")

        get_miARNs.assert_called_with('any_input', "any_path_to_mirna_db",
                                      'any_target_specie', 'any_evalue', 'any_perc_identity', 'any_output_path')

    @patch('Mirnas.get_counterparts_from_gene_id')
    @patch('Mirnas.get_db_name')
    def test_get_counterparts_from_gene_id_is_called_with_arguments_when_gene_id_serching_type_is_selected(self, get_db_name, get_counterparts_from_gene_id):
        get_db_name.return_value = "any_path_to_mirna_db"
        serching_type = 'GENE_ID'

        lookup_miRNAs("any_input", serching_type, "any_target_specie", "any_mirna_db", "any_evalue",
                      "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb")

        get_counterparts_from_gene_id.assert_called_with('any_input', "any_path_to_mirna_db", 'any_target_specie',
                                                         'any_evalue', 'any_perc_identity', 'any_entrez_db', 'any_entrezemail', 'any_output_path', 'any_blastdb')


if __name__ == '__main__':
    unittest.main()
