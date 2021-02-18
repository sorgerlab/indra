from nose.plugins.attrib import attr
from gilda import ground
from indra.sources import hypothesis
from indra.sources import trips
from indra.statements import *
from indra.sources.hypothesis.processor import HypothesisProcessor, \
    parse_context_entry, parse_grounding_entry, get_text_refs
from indra.sources.hypothesis.annotator import statement_to_annotations, \
    evidence_to_annotation, get_annotation_text


@attr('nonpublic', 'slow', 'notravis')
def test_process_indra_annnotations():
    hp = hypothesis.process_annotations(reader=trips.process_text)
    assert hp.statements
    for stmt in hp.statements:
        print(stmt)
        print(stmt.evidence[0])


def test_grounding_annotation():
    hp = HypothesisProcessor(annotations=[grounding_annot_example])
    hp.extract_groundings()
    assert hp.groundings['HCQ'] == {'CHEBI': 'CHEBI:5801'}
    assert hp.groundings['Plaquenil'] == {'CHEBI': 'CHEBI:5801'}


@attr('slow')
def test_statement_annotation():
    hp = HypothesisProcessor(annotations=[statement_annot_example],
                             reader=trips.process_text)
    hp.extract_statements()
    assert len(hp.statements) == 1
    stmt = hp.statements[0]
    assert stmt.subj.name == 'AMPK'
    assert stmt.obj.name == 'STAT3'
    context = stmt.evidence[0].context
    assert context.location.name == 'nucleus', context
    assert context.location.db_refs == {'GO': 'GO:0005634', 'TEXT': 'nucleus'}
    assert context.organ.name == 'Liver', context
    assert context.organ.db_refs == {'MESH': 'D008099', 'TEXT': 'liver'}


def test_get_text_refs_pmid():
    url = 'https://www.ncbi.nlm.nih.gov/pubmed/32196952'
    refs = get_text_refs(url)
    assert refs.get('PMID') == '32196952', refs
    assert refs.get('URL') == url, refs


def test_get_text_refs_pmcid():
    url = 'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7071777/'
    refs = get_text_refs(url)
    assert refs.get('PMCID') == 'PMC7071777', refs
    assert refs.get('URL') == url, refs


def test_get_text_refs_biorxiv():
    url = 'https://www.biorxiv.org/content/10.1101/2020.04.16.044016v1'
    refs = get_text_refs(url)
    assert refs.get('URL') == url, refs
    assert refs.get('DOI') == '10.1101/2020.04.16.044016', refs
    url = 'https://www.biorxiv.org/content/10.1101/2020.04.16.044016v1.full'
    refs = get_text_refs(url)
    assert refs.get('URL') == url, refs
    assert refs.get('DOI') == '10.1101/2020.04.16.044016', refs


def test_parse_grounding_entry():
    entry = '[a and b] -> CHEBI:CHEBI:1234|PUBCHEM:5678'
    grounding = parse_grounding_entry(entry)
    assert grounding == {'a and b': {'CHEBI': 'CHEBI:1234',
                                     'PUBCHEM': '5678'}}, grounding


def test_parse_invalid_grounding_entry():
    entries = ['xxx', '[xxx]->a', 'xxx -> a', 'xxx -> a:1&b:4']
    for entry in entries:
        assert parse_grounding_entry(entry) is None


def test_parse_context_entry():
    context_dict = parse_context_entry('Cell type: antigen presenting cells',
                                       ground, 'antigen presenting cells')
    assert len(context_dict) == 1
    assert 'cell_type' in context_dict
    ref_context = context_dict['cell_type']
    assert ref_context.name == 'Antigen-Presenting Cells', ref_context
    assert ref_context.db_refs.get('MESH') == 'D000938'
    assert ref_context.db_refs.get('TEXT') == 'antigen presenting cells'


def test_parse_invalid_context_entry():
    entries = ['xxx: yyy', 'Disease:something', 'xxx']
    for entry in entries:
        assert parse_context_entry(entry, ground) is None


def test_parse_ungrounded_context_entry():
    entry = 'Cell type: CD4+ T-cells'
    context_dict = parse_context_entry(entry, ground)
    assert len(context_dict['cell_type'].db_refs) == 1, \
        context_dict['cell_type'].db_refs
    assert context_dict['cell_type'].db_refs['TEXT'] == \
        'CD4+ T-cells', context_dict['cell_type'].db_refs


grounding_annot_example = {
 'uri': 'https://en.wikipedia.org/wiki/Hydroxychloroquine',
 'text': '[Plaquenil] -> CHEBI:CHEBI:5801\n\n[HCQ] -> CHEBI:CHEBI:5801',
 'tags': ['gilda'],
 'target': [{'source': 'https://en.wikipedia.org/wiki/Hydroxychloroquine'}],
 'document': {'title': ['Hydroxychloroquine - Wikipedia']},
}


statement_annot_example = {
 'id': '4nBYAmqwEeq1ujf13__Y-w',
 'uri': 'https://www.ncbi.nlm.nih.gov/pubmed/32190173',
 'text': 'AMPK activates STAT3\nOrgan: liver\nLocation: nucleus',
 'tags': [],
}


def test_get_annotation_text():
    # Test statement with multiple grounded agents
    stmt = Inhibition(
        Agent('vemurafenib', db_refs={'CHEBI': 'CHEBI:63637'}),
        Agent('BRAF', db_refs={'HGNC': '1097'})
    )
    annot_text = get_annotation_text(stmt, annotate_agents=True)
    assert annot_text == \
        '[vemurafenib](https://identifiers.org/CHEBI:63637) inhibits ' \
        '[BRAF](https://identifiers.org/hgnc:1097).', annot_text
    annot_text = get_annotation_text(stmt, annotate_agents=False)
    assert annot_text == 'Vemurafenib inhibits BRAF.', annot_text

    # Test statement with ungrounded and None agents
    stmt = Phosphorylation(None, Agent('X'))
    annot_text = get_annotation_text(stmt, annotate_agents=True)
    assert annot_text == 'X is phosphorylated.', annot_text
    annot_text = get_annotation_text(stmt, annotate_agents=False)
    assert annot_text == 'X is phosphorylated.', annot_text


def test_evidence_to_annot():
    # No evidence text
    ev = Evidence(source_api='reach')
    assert evidence_to_annotation(ev) is None

    # No text refs
    ev = Evidence(source_api='reach', text='Some text')
    assert evidence_to_annotation(ev) is None

    # Various text refs
    ev = Evidence(source_api='reach', text='Some text',
                  pmid='12345')
    annot = evidence_to_annotation(ev)
    assert annot == {'url': 'https://pubmed.ncbi.nlm.nih.gov/12345/',
                     'target_text': 'Some text',
                     'tags': ['reach']}, annot

    ev = Evidence(source_api='reach', text='Some text',
                  pmid=None, text_refs={'PMCID': '12345'})
    annot = evidence_to_annotation(ev)
    assert annot['url'] == 'https://www.ncbi.nlm.nih.gov/pmc/articles/12345/'


    ev = Evidence(source_api='reach', text='Some text',
                  pmid=None, text_refs={'URL': 'https://wikipedia.org'})
    annot = evidence_to_annotation(ev)
    assert annot['url'] == 'https://wikipedia.org'


def test_statement_to_annotations():
    evs = [
        # This will get filtered out
        Evidence(source_api='reach'),
        # This will get added as an annotation
        Evidence(source_api='sparser', text='some text 1',
                 pmid='12345'),
    ]
    stmt = Dephosphorylation(None, Agent('X'), evidence=evs)
    annots = statement_to_annotations(stmt)
    assert len(annots) == 1
    assert annots[0]['target_text'] == 'some text 1'
