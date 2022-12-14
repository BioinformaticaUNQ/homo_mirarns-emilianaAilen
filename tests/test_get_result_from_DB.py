import unittest
import os
import sys
sys.path.append(r'..')
from Mirnas import get_result_from_DB


class TestGetResultFromDb(unittest.TestCase):

    def tearDown(self):
        if os.path.isfile("./result.txt"):
            os.remove("result.txt")

    def test_get_result_from_db_does_not_writes_anything_in_the_output_file_when_the_gen_id_does_not_match(self):

        get_result_from_DB("tarbase", "non_existent_gene_id", "result.txt")

        f = open("./result.txt", "r")

        expected_result = ""
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)

    def test_get_result_from_db_writes_the_mirna_name_in_the_output_file_when_the_gen_id_matches_and_the_selected_db_is_tarbase(self):

        get_result_from_DB("tarbase", "0910001A06Rik", "result.txt")

        f = open("./result.txt", "r")

        expected_result = "mmu-miR-124-3p\n"
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)

    def test_get_result_from_db_writes_the_mirna_sequence_in_the_output_file_when_the_gen_id_matches_and_the_selected_db_is_mirnest(self):

        get_result_from_DB("mirnest_targets", "188432925", "result.txt")

        f = open("./result.txt", "r")

        expected_result = "TCTTCCCTACTCCACCCATGC\nTCTTCCCTACTCCACCCAT\n"
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)
