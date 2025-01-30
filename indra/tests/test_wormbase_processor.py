from unittest import TestCase
from unittest.mock import patch
from io import StringIO
import sys
import os

# Get the path dynamically (e.g., if relative to the current script)
script_dir = os.path.dirname(os.path.abspath(__file__))  # Current script directory
project_root = os.path.dirname(script_dir)  # Navigate one level up (Aging_Project)
path_to_sources = os.path.join(project_root, 'sources')  # Adjust as per your structure
sys.path.append(path_to_sources)

from wormbase.processor import WormBaseProcessor

class TestWormBaseProcessor(TestCase):
    def test_download_wormbase_data_real(self):
        """Test the `_download_wormbase_data` function with the actual WormBase URL."""
        # Define the real WormBase data URL
        url = 'https://fms.alliancegenome.org/download/INTERACTION-GEN_WB.tsv.gz'

        print("\nStarting test for _download_wormbase_data with real WormBase data...")

        try:
            # Invoke the actual function
            csv_reader = WormBaseProcessor._download_wormbase_data(url)

            # Read the rows into a list
            rows = list(csv_reader)

            # Print some of the rows for visual inspection
            print("\nSample rows from WormBase data:")
            for i, row in enumerate(rows[:5]):  # Print the first 5 rows
                print(f"Row {i + 1}: {row}")

            # Assert there are rows (basic validation)
            print(f"Number of rows found in the WormBase file: {len(rows)}")
            self.assertGreater(len(rows), 0, "No data rows found in the WormBase file.")

            # Assert the structure of the rows
            expected_columns = 42
            self.assertEqual(
                len(rows[0]),
                expected_columns,
                f"Expected {expected_columns} columns; found {len(rows[0])}."
            )

            print("\nTest passed: WormBase data downloaded and validated successfully.")

        except Exception as e:
            # Fail the test if an exception occurs
            self.fail(f"Test failed due to an exception: {e}")

    def test_print_agent_data_in_init(self):
        """Test the __init__ method to inspect name_agent_a, wormbase_id_agent_a, and entrez_id_agent_a."""

        # Define the WormBase data URL
        url = 'https://fms.alliancegenome.org/download/INTERACTION-GEN_WB.tsv.gz'

        try:
            # Download the WormBase data
            csv_reader = WormBaseProcessor._download_wormbase_data(url)
            rows = list(csv_reader)

            # Ensure there is at least one data row
            self.assertGreater(len(rows), 0, "No data rows found in the WormBase file.")

            # Use the first row in the downloaded data
            first_row = rows[0]

            # Convert the first row into a test file (optional)
            test_file_path = os.path.join(os.path.dirname(__file__), "test_wormbase_data.tsv")
            with open(test_file_path, "w") as f:
                f.write("\t".join(first_row) + "\n")

            # Capture the output of the print statements
            with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
                # Instantiate WormBaseProcessor with the test file
                processor = WormBaseProcessor(wormbase_file=test_file_path)

                # Fetch and store the captured output
                output = mock_stdout.getvalue()

            # Print the captured output for inspection
            print("\nCaptured Output:\n", output)

            # Add assertions to check for specific print content (optional)
            if len(first_row) >= 5:  # Example assumption on available columns for aliases
                self.assertIn(first_row[0], output, f"Expected {first_row[0]} in output for name_agent_a.")
            if len(first_row) >= 42:  # Example assumption WormBase schema
                self.assertIn(first_row[41], output, f"Expected {first_row[41]} in output for identification method.")

        except Exception as e:
            self.fail(f"Test failed due to an exception: {e}")
        finally:
            # Clean up the test file
            if os.path.exists(test_file_path):
                os.remove(test_file_path)


if __name__ == "__main__":
    import unittest

    unittest.main()
