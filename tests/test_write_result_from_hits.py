import unittest
import os
import sys
sys.path.append(r'..')
from Mirnas import write_result_from_hits

class TestWriteResultFromHits(unittest.TestCase):
    def tearDown(self):
        if os.path.isfile("./result.txt"):
            os.remove("result.txt")

    def test_write_result_from_hits_only_takes_the_results_which_match_the_target_specie(self):
        mocked_hits = [
            "hit_1_specie_x 75% 0.05\n",
            "hit_2_specie_y 85% 0.07\n", 
            "hit_3_specie_x 79% 0.05\n"
        ]

        write_result_from_hits(mocked_hits, "specie_x", "result.txt")

        f = open("./result.txt", "r")

        expected_result = "Description - Identity percentage - E-value\nhit_1_specie_x 75% 0.05\nhit_3_specie_x 79% 0.05\n"
        obtained_result = f.read()

        f.close()

        self.assertEqual(expected_result, obtained_result)
    
    def test_write_result_from_hits_only_just_writes_the_header_if_the_hits_list_is_empty(self):
        mocked_hits = []

        write_result_from_hits(mocked_hits, "specie_x", "result.txt")

        f = open("./result.txt", "r")

        expected_result = "Description - Identity percentage - E-value\n"
        obtained_result = f.read()
        
        f.close()

        self.assertEqual(expected_result, obtained_result)

if __name__ == '__main__':
    unittest.main()
