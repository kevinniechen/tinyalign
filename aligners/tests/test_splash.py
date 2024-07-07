import unittest
import sys
import os

# Add the parent directory to the Python path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from splash import run_splash, count_targets

class TestSplash(unittest.TestCase):
    def setUp(self):
        self.toy_data = {
            "sample1": [
                "AAAAABBBBBCCCCCDDDDDEEEEE",
                "AAAAABBBBBDDDDDEEEEEFFFF",
                "AAAAABBBBBCCCCCDDDDDEEEEE"
            ],
            "sample2": [
                "AAAAABBBBBDDDDDEEEEEFFFF",
                "AAAAABBBBBCCCCCDDDDDEEEEE",
                "AAAAABBBBBDDDDDEEEEEFFFF"
            ]
        }
        self.anchor_length = 5
        self.target_length = 5
        self.offset = 5

    def test_count_targets(self):
        count_tables = count_targets(self.toy_data, self.anchor_length, self.target_length, self.offset)
        
        # Test the number of unique anchors
        self.assertEqual(len(count_tables), 15, "Should have 15 unique anchors")

        # Test specific count tables
        self.assertEqual(count_tables["AAAAA"], 
                         {"CCCCC": {"sample1": 2, "sample2": 1},
                          "DDDDD": {"sample1": 1, "sample2": 2}},
                         "Incorrect count table for AAAAA")

        self.assertEqual(count_tables["BBBBB"], 
                         {"DDDDD": {"sample1": 2, "sample2": 1},
                          "EEEEE": {"sample1": 1, "sample2": 2}},
                         "Incorrect count table for BBBBB")

        self.assertEqual(count_tables["CCCCC"], 
                         {"EEEEE": {"sample1": 2, "sample2": 1}},
                         "Incorrect count table for CCCCC")

    def test_run_splash(self):
        # Capture the output of run_splash
        import io
        import sys
        captured_output = io.StringIO()
        sys.stdout = captured_output
        
        run_splash(self.toy_data, self.anchor_length, self.target_length, self.offset)
        
        sys.stdout = sys.__stdout__  # Reset redirect.
        output = captured_output.getvalue()

        # Check if the output contains expected information
        self.assertIn("Total number of unique anchors: 15", output)
        self.assertIn("Count table for anchor: AAAAA", output)
        self.assertIn("Count table for anchor: BBBBB", output)
        self.assertIn("Count table for anchor: CCCCC", output)

if __name__ == '__main__':
    unittest.main()
