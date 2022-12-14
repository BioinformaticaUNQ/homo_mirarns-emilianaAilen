import unittest
import sys
from unittest.mock import patch
sys.path.append(r'..')
from Mirnas import lookup_miRNAs


class TestLookupMirnas(unittest.TestCase):

    @patch('Mirnas.get_db_name')
    def test_get_best_alignment_id_returns_none_if_the_record_name_does_not_contain_the_target_specie(self, get_db_name):
        non_existent_seraching_type = "none"
        get_db_name.return_value = "any_path_to_mirna_db"

        with self.assertRaises(Exception):
            lookup_miRNAs("any_input", non_existent_seraching_type, "any_target_specie", "any_mirna_db", "any_evalue",
                          "any_perc_identity", "any_entrez_db", "any_output_path", "any_entrezemail", "any_blastdb")

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
