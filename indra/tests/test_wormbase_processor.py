from unittest import TestCase
import os
from indra.sources.wormbase.processor import WormBaseProcessor, _WormBaseRow, columns #, wormbase_file_url
from indra.statements import Evidence, Activation, Association, Inhibition


this_dir = os.path.dirname(__file__)
test_file = os.path.join(this_dir, 'wormbase_tests_data/INTERACTION-GEN_WB_test.tsv')

class TestWormBaseProcessor(TestCase):
    # def test_download_wormbase_data(self):
    #     """Test the `_download_wormbase_data` function with the actual
    #     WormBase URL."""
    #
    #     #  URL to WormBase data
    #     url = wormbase_file_url
    #
    #     print("\nStarting test for _download_wormbase_data...")
    #
    #     try:
    #         # Invoke the actual function
    #         rows = WormBaseProcessor._download_wormbase_data(url)
    #
    #         # Assert there are rows (basic validation)
    #         print(f"Number of rows found in the test file: {len(rows)}")
    #         self.assertGreater(len(rows), 0, "No data rows found "
    #                                          "in the test file.")
    #
    #         # Print some of the rows for visual inspection
    #         print("\nSample rows from test data:")
    #         for i, row in enumerate(rows[:5]):  # Print the first 5 rows
    #             print(f"Row {i + 1}: {row}")
    #
    #         # Assert the structure of the rows
    #         expected_columns = 42
    #         self.assertEqual(
    #             len(rows[0]),
    #             expected_columns,
    #             f"Expected {expected_columns} columns; found {len(rows[0])}."
    #         )
    #
    #         print("\nTest passed: Data downloaded and validated "
    #               "successfully.")
    #
    #     except Exception as e:
    #         # Fail the test if an exception occurs
    #         self.fail(f"Test failed due to an exception: {e}")

    def test_read_wormbase_data(self):
        """Test the `_read_wormbase_data` function with a locally
        downloaded WormBase tsv file named INTERACTION-GEN_WB_test.tsv.
        """
        try:
            rows = WormBaseProcessor(wormbase_file=test_file)._read_wormbase_data()

            # Assert there are rows (basic validation)
            self.assertEqual(len(rows), 5)

            # Assert the structure of the rows
            expected_columns = 42
            self.assertEqual(
                len(rows[0]),
                expected_columns)

        except Exception as e:
            # Fail the test if an exception occurs
            self.fail(f"Test failed due to an exception: {e}")

    def test_agent_extraction(self):
        """Test extraction of agent data from a row in WormBaseProcessor."""

        try:
            processor = WormBaseProcessor(wormbase_file=test_file)
            rows = processor._read_wormbase_data()

            sample = rows[0]
            filt_sample = [None if item == '-' else item
                        for item in sample][:len(columns)]
            sample = _WormBaseRow(*filt_sample)

            # aliases_interactor_a
            name_info_agent_a = processor._alias_conversion(sample.aliases_interactor_a)
            name_agent_a = name_info_agent_a.get("public_name")

            # aliases_interactor_b
            name_info_agent_b = processor._alias_conversion(sample.aliases_interactor_b)
            name_agent_b = name_info_agent_b.get("public_name")

            # ids_interactor_a
            db_id_info_agent_a = processor._id_conversion(sample.ids_interactor_a)
            wormbase_id_agent_a = db_id_info_agent_a.get("wormbase")
            entrez_id_agent_a = db_id_info_agent_a.get("entrez gene/locuslink")

            # ids_interactor_b
            db_id_info_agent_b = processor._id_conversion(sample.ids_interactor_b)
            wormbase_id_agent_b = db_id_info_agent_b.get("wormbase")
            entrez_id_agent_b = db_id_info_agent_b.get("entrez gene/locuslink")

            # Ground agent
            agent_a = processor._make_agent(name_agent_a, wormbase_id_agent_a, entrez_id_agent_a)
            agent_b = processor._make_agent(name_agent_b, wormbase_id_agent_b, entrez_id_agent_b)

            # Skip any agents with neither HGNC grounding or string name
            if agent_a is None or agent_b is None:
                return
            # Get evidence
            pmid = processor._id_conversion(sample.publication_identifiers).get('pubmed')
            doi = processor._id_conversion(sample.publication_identifiers).get('doi')
            text_refs = {}
            if pmid:
                text_refs['PMID'] = pmid
            elif doi:
                text_refs['DOI'] = doi

            source_id = processor._id_conversion(sample.interaction_identifiers).get('wormbase')
            interaction_annotations = processor._id_conversion(sample.interaction_annotations).get('wormbase')
            ev = Evidence(source_api='wormbase',
                          source_id=source_id,
                          pmid=text_refs.get('PMID'),
                          text_refs=text_refs,
                          # annotations=interaction_annotations
                          annotations=dict(sample._asdict())
                          )
            # Make statement
            interaction_type = processor._interaction_type_conversion(sample.interaction_types).get('psi-mi')
            if 'enhancement' in interaction_type:
                s = Activation([agent_a, agent_b], evidence=ev)
            elif 'suppression' in interaction_type:
                s = Inhibition([agent_a, agent_b], evidence=ev)
            else:
                s = Association([agent_a, agent_b], evidence=ev)
            processor.statements.append(s)

            # Assertions
            self.assertDictEqual(
                name_info_agent_a,
                {"public_name": "dpy-21", "sequence_name": "Y59A8B.1"})
            self.assertEqual("dpy-21", name_agent_a)

            self.assertDictEqual(
                db_id_info_agent_a,
                {"wormbase": "WBGene00001080"})

            # Entrez ID is None in the first row of the test data
            self.assertIsNone(entrez_id_agent_a)

            print(agent_a)
            print(agent_a.name)
            print(agent_a.db_refs)

            # print(processor.statements)
            print(processor.statements[0])
            print(processor.statements[0].evidence)



        except Exception as e:
            # Fail the test if an exception occurs
            self.fail(f"Test failed due to an exception: {e}")


if __name__ == "__main__":
    import unittest

    unittest.main()
