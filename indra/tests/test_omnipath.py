import requests
from indra.sources.omnipath import OmniPathModificationProcessor,\
    OmniPathLiganReceptorProcessor
from indra.sources.omnipath.api import op_url
from indra.statements import Agent, Phosphorylation, Complex
from indra.preassembler.grounding_mapper import GroundingMapper
from nose.plugins.attrib import attr

BRAF_UPID = 'P15056'
JAK2_UPID = 'O60674'
BRAF_AG = Agent(None, db_refs={'UP': BRAF_UPID})
GroundingMapper.standardize_agent_name(BRAF_AG)
JAK2_AG = Agent(None, db_refs={'UP': JAK2_UPID})
GroundingMapper.standardize_agent_name(JAK2_AG)


def test_omnipath_web_api():
    query_url = '%s/queries' % op_url
    res = requests.get(query_url)
    assert res.status_code == 200


def test_mods_from_web():
    params = {'format': 'json', 'substrates': JAK2_UPID,
              'fields': ['sources', 'references']}
    ptm_url = '%s/ptms' % op_url
    res = requests.get(ptm_url, params=params)
    assert res.status_code == 200
    assert res.text
    ptm_json = res.json()
    assert ptm_json[0]['substrate'] == JAK2_UPID, ptm_json[0]['substrate']
    stmts = OmniPathModificationProcessor(ptm_json).statements
    assert JAK2_AG.name in [a.name for a in stmts[0].agent_list()],\
        stmts[0].agent_list()
    assert 'omnipath' == stmts[0].evidence[0].source_api,\
        stmts[0].evidence[0].source_api


def test_pypath_import():
    # Import package
    try:
        import pypath
    except ImportError:
        pypath = None

    assert pypath, 'PyPath is not avaialble'

    # Import of main
    try:
        from pypath import main as pypath_main
    except ImportError:
        pypath_main = None
    assert pypath_main, 'Could not import pypath.main'

    # Data formats
    try:
        from pypath import data_formats
    except ImportError:
        data_formats = None
    assert data_formats, 'Could not import pypath.data_formats'


@attr('no-travis')
def test_lr_pypath_network():
    try:
        from pypath import main as pypath_main, data_formats
    except ImportError:
        pypath_main = None
        data_formats = None
    assert pypath_main and data_formats, 'Failed to import pypath'
    pa = pypath_main.PyPath()
    pa.init_network({
        'hpmr': data_formats.ligand_receptor['hpmr']
    })
    stmts = OmniPathLiganReceptorProcessor(pa).statements
    assert len(stmts) > 0, 'len(stmts) = %d' % len(stmts)
    stmt = stmts[0]
    assert isinstance(stmt, Complex)
    ev = stmt.evidence[0]
    assert 'omnipath' == ev.source_api, ev.source_api
    assert ev.pmid or ev.text_refs, 'pmid=%s, ev.text_refs=%s' % \
                                    (ev.pmid, ev.text_refs)
    assert 'source_sub_id' in ev.annotations, print(ev.annotations.keys())
    assert 'hpmr' == ev.annotations['source_sub_id'].lower(), \
        ev.annotations['source_sub_id']
