from nose.plugins.attrib import attr
from indra.sources import hypothesis
from indra.sources.hypothesis.processor import HypothesisProcessor


@attr('nonpublic')
def test_process_indra_annnotations():
    hp = hypothesis.process_annotations()
    assert hp.statements
    print(hp.statements)
