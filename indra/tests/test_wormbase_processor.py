from unittest import TestCase
import sys
import os

# Get path to wormbase module
script_dir = os.path.dirname(os.path.abspath(__file__))  # Current script directory
project_root = os.path.dirname(script_dir)  # Navigate one level up (indra)
path_to_sources = os.path.join(project_root, 'sources')
sys.path.append(path_to_sources)

from wormbase.processor import WormBaseProcessor, wormbase_file_url

class TestWormBaseProcessor(TestCase):
    def setUp(self):
        self.test_file_path = os.path.join(project_root, "tests/wormbase_tests_data/INTERACTION-GEN_WB_test.tsv")

    def test_download_wormbase_data(self):
        """Test the `_download_wormbase_data` function with the actual
        WormBase URL."""

        #  URL to WormBase data
        url = wormbase_file_url

        print("\nStarting test for _download_wormbase_data...")

        try:
            # Invoke the actual function
            rows = WormBaseProcessor._download_wormbase_data(url)

            # Assert there are rows (basic validation)
            print(f"Number of rows found in the test file: {len(rows)}")
            self.assertGreater(len(rows), 0, "No data rows found "
                                             "in the test file.")

            # Print some of the rows for visual inspection
            print("\nSample rows from test data:")
            for i, row in enumerate(rows[:5]):  # Print the first 5 rows
                print(f"Row {i + 1}: {row}")

            # Assert the structure of the rows
            expected_columns = 42
            self.assertEqual(
                len(rows[0]),
                expected_columns,
                f"Expected {expected_columns} columns; found {len(rows[0])}."
            )

            print("\nTest passed: Data downloaded and validated "
                  "successfully.")

        except Exception as e:
            # Fail the test if an exception occurs
            self.fail(f"Test failed due to an exception: {e}")

    def test_read_wormbase_data(self):
        """Test the `_read_wormbase_data` function with a locally
        downloaded WormBase tsv file named INTERACTION-GEN_WB_test.tsv.
        """
        print("\nStarting test for _read_wormbase_data...")

        # test_file_path = os.path.join(path_to_wb_test_data_dir, "INTERACTION-GEN_WB_test.tsv")

        try:
            rows = WormBaseProcessor(wormbase_file=self.test_file_path)._read_wormbase_data()

            # Assert there are rows (basic validation)
            print(f"Number of rows found in the test file: {len(rows)}")
            self.assertGreater(len(rows), 0, "No data rows found in "
                                             "the test file.")

            # Print some of the rows for visual inspection
            print("\nSample rows from test data:")
            for i, row in enumerate(rows[:5]):  # Print the first 5 rows
                print(f"Row {i + 1}: {row}")

            # Assert the structure of the rows
            expected_columns = 42
            self.assertEqual(
                len(rows[0]),
                expected_columns,
                f"Expected {expected_columns} columns; found {len(rows[0])}."
            )

            print("\nTest passed: Data read and validated successfully.")

        except Exception as e:
            # Fail the test if an exception occurs
            self.fail(f"Test failed due to an exception: {e}")

    def test_agent_extraction(self):
        """Test extraction of agent data from a row in WormBaseProcessor."""

        print("\nStarting test for Agent extraction "
              "(WormBaseProcessor.__init__)...")

        try:
            processor = WormBaseProcessor(wormbase_file=self.test_file_path)
            rows = processor._read_wormbase_data()

            sample = rows[0]

            # aliases_interactor_a
            name_info_agent_a = processor._alias_conversion(sample[4])
            name_agent_a = name_info_agent_a.get("public_name")

            # ids_interactor_a
            db_id_info_agent_a = processor._id_conversion(sample[0])
            worm_base_id_agent_a = db_id_info_agent_a.get("wormbase")
            entrez_id_agent_a = db_id_info_agent_a.get("entrez gene/locuslink")

            # Assertions for value validation
            self.assertIsInstance(name_info_agent_a, dict)
            # self.assertIn("public_name", name_info_agent_a)
            self.assertIsNotNone(name_info_agent_a)
            print(f"name_info_agent_a: {name_info_agent_a}")

            self.assertIsNotNone(name_agent_a)
            print(f"name_agent_a: {name_agent_a}")

            self.assertIsInstance(db_id_info_agent_a, dict)
            self.assertIn("wormbase", db_id_info_agent_a)
            print(f"db_id_info_agent_a: {db_id_info_agent_a}")

            self.assertIsNotNone(worm_base_id_agent_a)
            print(f"wormbase_id_agent_a: {worm_base_id_agent_a}")

            # Entrez ID is None in the first row of the test data
            self.assertIsNone(entrez_id_agent_a)
            print(f"entrez_id_agent_a: {entrez_id_agent_a}")

            print("\nTest passed: Agent extraction successful.")


        except Exception as e:
            # Fail the test if an exception occurs
            self.fail(f"Test failed due to an exception: {e}")


if __name__ == "__main__":
    import unittest

    unittest.main()
