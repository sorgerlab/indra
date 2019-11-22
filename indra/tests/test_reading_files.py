from os import path
from indra.tools.reading.read_files import read_files, get_reader_classes

from nose.plugins.attrib import attr


@attr('slow', 'nonpublic', 'notravis')
def test_read_files():
    "Test that the system can read files."
    # Create the test files.
    example_files = []

    # Get txt content
    abstract_txt = ("This is a paper that contains the phrase: MEK "
                    "phosphorylates ERK.")
    with open('test_abstract.txt', 'w') as f:
        f.write(abstract_txt)
    example_files.append('test_abstract.txt')

    # Get nxml content
    pmc_test_fpath = path.join(path.dirname(path.abspath(__file__)),
                               'pmc_cont_example.nxml')
    if path.exists(pmc_test_fpath):
        example_files.append(pmc_test_fpath)

    assert len(example_files), "No content available to test."

    # Now read them.
    reader_classes = get_reader_classes()
    readers = []
    for rc in reader_classes:
        readers.append(rc())
    outputs = read_files(example_files, readers)
    N_out = len(outputs)
    N_exp = len(readers)*len(example_files)
    assert N_out == N_exp, "Expected %d outputs, got %d." % (N_exp, N_out)
