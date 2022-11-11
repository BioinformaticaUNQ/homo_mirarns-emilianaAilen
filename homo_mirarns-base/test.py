import unittest

from homo_mirarns import Library

lookup_miRNAs = Library.lookup_miRNAs

class Testing(unittest.TestCase):
    def test_is_fasta(self):
        result = lookup_miRNAs('>asdasd')
        self.assertEqual(result, 'es una secuencia FASTA')

    def test_is_gene_id(self):
        result = lookup_miRNAs('id')
        self.assertEqual(result, 'es un gene id')

if __name__ == '__main__':
    unittest.main()
