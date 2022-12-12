from Mirnas import get_db_name
import unittest
import sys
sys.path.append(r'..')


class GetDBName(unittest.TestCase):

    # Fasta method
    def test_get_db_name_fasta_rumimir(self):
        with self.assertRaises(Exception):
            self.assertRaises(get_db_name("FASTA", "RUMIMIR"))

    def test_get_db_name_fasta_mirbase(self):
        with self.assertRaises(Exception):
            self.assertRaises(get_db_name("FASTA", "MIRBASE"))

    def test_get_db_name_fasta_tarbase(self):
        expected_name = "tarbase"
        obtained_name = get_db_name("FASTA", "TARBASE")
        self.assertEqual(expected_name, obtained_name)

    def test_get_db_name_fasta_mirnest(self):
        expected_name = "mirnest_targets"
        obtained_name = get_db_name("FASTA", "MIRNEST")
        self.assertEqual(expected_name, obtained_name)

    # Gene id method

    def test_get_db_name_gene_id_rumimir(self):
        with self.assertRaises(Exception):
            self.assertRaises(get_db_name("GENE_ID", "RUMIMIR"))

    def test_get_db_name_gene_id_mirbase(self):
        with self.assertRaises(Exception):
            self.assertRaises(get_db_name("GENE_ID", "MIRBASE"))

    def test_get_db_name_gene_id_tarbase(self):
        expected_name = "tarbase"
        obtained_name = get_db_name("GENE_ID", "TARBASE")
        self.assertEqual(expected_name, obtained_name)

    def test_get_db_name_gene_id_mirnest(self):
        expected_name = "mirnest_targets"
        obtained_name = get_db_name("GENE_ID", "MIRNEST")
        self.assertEqual(expected_name, obtained_name)

    # Gene Mirna Fasta method

    def test_get_db_name_mirna_fasta_rumimir(self):
        expected_name = "rumimir"
        obtained_name = get_db_name("MIRNA_FASTA", "RUMIMIR")
        self.assertEqual(expected_name, obtained_name)

    def test_get_db_name_mirna_fasta_mirbase(self):
        expected_name = "mirbase"
        obtained_name = get_db_name("MIRNA_FASTA", "MIRBASE")
        self.assertEqual(expected_name, obtained_name)

    def test_get_db_name_mirna_fasta_tarbase(self):
        with self.assertRaises(Exception):
            self.assertRaises(get_db_name("MIRNA_FASTA", "TARBASE"))

    def test_get_db_name_mirna_fasta_mirnest(self):
        expected_name = "mirnest"
        obtained_name = get_db_name("MIRNA_FASTA", "MIRNEST")
        self.assertEqual(expected_name, obtained_name)


if __name__ == '__main__':
    unittest.main()
