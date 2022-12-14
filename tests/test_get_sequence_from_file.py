import unittest
import os
import sys
sys.path.append(r'..')
from Mirnas import get_sequence_from_file


class TestGetSequenceFromFile(unittest.TestCase):
    def tearDown(self):
        if os.path.isfile("./test-file.txt"):
            os.remove("test-file.txt")

    def test_get_sequence_from_file_fails_if_the_file_does_not_exists(self):
        with self.assertRaises(ValueError) as context:
            get_sequence_from_file("./non-existent-path-file")

        expected_exception_message = 'The input file path is incorrect, please fix it.'
        obtained_exception_message = str(context.exception)

        self.assertEqual(expected_exception_message,
                         obtained_exception_message)

    def test_get_sequence_from_file_fails_if_the_file_does_not_contain_any_data(self):
        f = open("test-file.txt", "w", encoding="utf-8")
        f.close()

        with self.assertRaises(ValueError) as context:
            get_sequence_from_file("test-file.txt")

        expected_exception_message = "There's no sequence on the following file: test-file.txt"
        obtained_exception_message = str(context.exception)

        self.assertEqual(expected_exception_message,
                         obtained_exception_message)

    def test_get_sequence_from_file_returns_the_file_data(self):
        f = open("test-file.txt", "w", encoding="utf-8")
        f.write("> a_fasta_sequence \n TCTTCCCTACTCCACCCAT")
        f.close()

        expected_data = "> a_fasta_sequence \n TCTTCCCTACTCCACCCAT"
        obtained_data = get_sequence_from_file("test-file.txt")

        self.assertEqual(expected_data,
                         obtained_data)


if __name__ == '__main__':
    unittest.main()
