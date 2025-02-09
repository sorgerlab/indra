from collections import Counter
from unittest import TestCase
import os
import tqdm
from indra.sources.wormbase.processor import WormBaseProcessor


this_dir = os.path.dirname(__file__)
test_file_gen = os.path.join(this_dir, 'wormbase_tests_data/INTERACTION-GEN_WB_test.tsv')
test_file_mol = os.path.join(this_dir, 'wormbase_tests_data/INTERACTION-MOL_WB_3_test.tsv')


class TestWormBaseProcessor(TestCase):

    def test_statements_creation(self):
        """Test attributes and properties of the `statements` list in WormBaseProcessor."""

        try:
            # Create a WormBaseProcessor instance with test files
            processor = WormBaseProcessor(test_file_gen, test_file_mol)
            # processor = WormBaseProcessor()

            # Ensure `statements` has been populated
            self.assertGreater(len(processor.statements), 0, "The statements list should not be empty.")

            # print(len(processor.statements))
            # interaction_type_counts = Counter(type(statement).__name__ for statement in processor.statements)
            # print(interaction_type_counts)

        except Exception as e:
            self.fail(f"Test failed due to an exception: {e}")


if __name__ == "__main__":
    import unittest

    unittest.main()
