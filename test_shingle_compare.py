import unittest
import sys
import os
import tempfile
import shutil
from collections import Counter

# Add the directory containing the original scripts to Python path
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Import functions to test
from shingle import shingles, most_common_k_shingles
from compare import remove_punctuation, get_shingle_from_file, calc_jaccard


class TestShingleAndCompareModules(unittest.TestCase):
    """
    Comprehensive unit test suite for shingle and compare modules.

    This test class provides thorough testing for shingle generation,
    file processing, and similarity calculation functions.
    """


    def setUp(self):
        """
        Prepare test environment with actual text.txt content.
        """
        # Get the directory of the current script
        self.current_dir = os.path.dirname(os.path.abspath(__file__))
        self.test_file_path = os.path.join(self.current_dir, 'text.txt')

    def test_get_shingle_from_file_content(self):
        """
        Test reading and generating shingles from text.txt file,
        with and without punctuation removal.
        """
        # Test without punctuation removal
        shingles_list = get_shingle_from_file(self.test_file_path, 2, remove_punct=False)
        expected_shingles = [
            ('let', 'it'), ('it', 'be,'), ('be,', 'let'),
            ('let', 'it'), ('it', 'be,'), ('be,', 'let'),
            ('let', 'it'), ('it', 'be'), ('be', "can't"),
            ("can't", 'hold'), ('hold', 'it'), ('it', 'back'),
            ('back', 'further'), ('further', 'more')
        ]

        self.assertEqual(shingles_list, expected_shingles)

        # Test with punctuation removal
        shingles_list_no_punct = get_shingle_from_file(self.test_file_path, 2, remove_punct=True)
        expected_shingles_no_punct = [
            ('let', 'it'), ('it', 'be'), ('be', 'let'),
            ('let', 'it'), ('it', 'be'), ('be', 'let'),
            ('let', 'it'), ('it', 'be'), ('be', "cant"),
            ("cant", 'hold'), ('hold', 'it'), ('it', 'back'),
            ('back', 'further'), ('further', 'more')
        ]
        self.assertEqual(shingles_list_no_punct, expected_shingles_no_punct)

    def test_shingles_basic(self):
        """Test basic k-shingle generation from a list of tokens."""
        tokens = ['a', 'b', 'c', 'd', 'e']
        result = shingles(tokens, 2)
        expected = [('a', 'b'), ('b', 'c'), ('c', 'd'), ('d', 'e')]
        self.assertEqual(result, expected)

    def test_shingles_edge_cases(self):
        """
        Test edge cases for shingle generation including
        empty list, oversized shingles, and zero-sized shingles.
        """
        self.assertEqual(shingles([], 2), [])
        self.assertEqual(shingles(['a', 'b'], 3), [])
        self.assertEqual(shingles(['a', 'b', 'c'], 0), [])

    def test_most_common_k_shingles(self):
        """
        Verify finding most common k-shingles with multiple occurrences.
        """
        tokens = ['let', 'it', 'be', 'let', 'it', 'be', 'let', 'it']
        result = most_common_k_shingles(tokens, 2, 3)

        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], (('let', 'it'), 3))
        self.assertEqual(result[1], (('it', 'be'), 2))

    def test_remove_punctuation(self):
        """Test removal of punctuation from text."""
        text = "Hello, world! How are you?"
        cleaned_text = remove_punctuation(text)
        self.assertEqual(cleaned_text, "Hello world How are you")

    def test_calc_jaccard(self):
        """
        Test Jaccard similarity calculation for various scenarios:
        No overlap, partial overlap, and complete overlap.
        """
        test_cases = [
            {
                'query': [('a', 'b'), ('b', 'c')],
                'target': [('x', 'y'), ('y', 'z')],
                'expected': '0/4'
            },
            {
                'query': [('a', 'b'), ('b', 'c'), ('a', 'b')],
                'target': [('b', 'c'), ('c', 'd'), ('a', 'b')],
                'expected': '2/4'
            },
            {
                'query': [('a', 'b'), ('b', 'c')],
                'target': [('a', 'b'), ('b', 'c')],
                'expected': '2/2'
            }
        ]

        for case in test_cases:
            result = calc_jaccard(case['query'], case['target'])
            self.assertEqual(result, case['expected'])


if __name__ == '__main__':
    unittest.main(argv=[''], verbosity=2, exit=False)