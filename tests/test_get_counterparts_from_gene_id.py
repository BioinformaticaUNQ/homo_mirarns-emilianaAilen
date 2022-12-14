import unittest
import sys
from unittest.mock import patch
sys.path.append(r'..')
from Mirnas import get_counterparts_from_gene_id


class TestLookupMirnas(unittest.TestCase):

    @patch('Mirnas.get_sequence_by_')
    @patch('Mirnas.get_counterparts')
    def test_get_sequence_by_and_get_counterpart_are_called_with_arguments(self, get_counterparts, get_sequence_by_):

        get_counterparts_from_gene_id("any_gene_id", "any_mirna_db", "any_target_specie", "any_evalue",
                                      "any_perc_identity", "any_entrez_db", "any_mail", "any_output_path", "any_blastdb")

        get_sequence_by_.assert_called_with(
            "any_gene_id", "any_entrez_db", "any_mail")
        get_counterparts.assert_called_with("sequenceFound.fasta", "any_mirna_db",
                                            "any_target_specie", "any_evalue", "any_perc_identity", "any_output_path", "any_blastdb")


if __name__ == '__main__':
    unittest.main()
