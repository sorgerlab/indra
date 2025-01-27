from unittest import TestCase
# from unittest.mock import patch
# from io import BytesIO
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