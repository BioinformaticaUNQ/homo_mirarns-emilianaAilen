from collections import namedtuple
import unittest
import sys
sys.path.append(r'..')
from Mirnas import get_best_alignment_ID

""" def get_best_alignment_ID(E_VALUE_THRESH, specie, alignments):
    for alignment in alignments:
        for hsp in alignment.hsps:
            print(alignment.title)
            if hsp.expect < E_VALUE_THRESH and specie in alignment.title:
                return alignment.title.split("|")[1]

 """


class GetBestAlignmentId(unittest.TestCase):

    def test_get_best_alignment_id_returns_none_if_all_the_aligments_have_an_evalue_above_the_maximum(self):
        mocked_aligments = [
            namedtuple('Record', "hsps title")(hsps=[
                namedtuple('HSP', "expect")(expect=0.90),
                namedtuple('HSP',  "expect")(expect=0.95)], title="gi|0001| a_specie")
        ]

        maximum_evalue = 0.50
        specie = "a_specie"

        expected_result = None
        obtained_result = get_best_alignment_ID(
            maximum_evalue, specie, mocked_aligments)

        self.assertEquals(expected_result, obtained_result)
    
    def test_get_best_alignment_id_returns_none_if_the_record_name_does_not_contain_the_target_specie(self):
        mocked_aligments = [
            namedtuple('Record', "hsps title")(hsps=[
                namedtuple('HSP', "expect")(expect=0.90),
                namedtuple('HSP',  "expect")(expect=0.95)], title="gi|0001| another_specie")
        ]

        maximum_evalue = 0.97
        specie = "a_specie"

        expected_result = None
        obtained_result = get_best_alignment_ID(
            maximum_evalue, specie, mocked_aligments)

        self.assertEquals(expected_result, obtained_result)

    def test_get_best_alignment_id_returns_the_gene_id_when_all_conditions_are_met(self):
        mocked_aligments = [
            namedtuple('Record', "hsps title")(hsps=[
                namedtuple('HSP', "expect")(expect=0.90),
                namedtuple('HSP',  "expect")(expect=0.95)], title="gi|0001| a_specie")
        ]

        maximum_evalue = 1
        specie = "a_specie"

        expected_result = "0001"
        obtained_result = get_best_alignment_ID(
            maximum_evalue, specie, mocked_aligments)

        self.assertEquals(expected_result, obtained_result)


if __name__ == '__main__':
    unittest.main()