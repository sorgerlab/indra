import re
import logging
from indra.databases import uniprot_client
from indra.statements import Agent, Complex, Evidence
from indra.preassembler.grounding_mapper import standardize_agent_name


logger = logging.getLogger(__name__)


class VirhostnetProcessor:
    """A processor that takes a pandas DataFrame and extracts INDRA Statements.

    Parameters
    ----------
    df : pandas.DataFrame
        A pandas DataFrame representing VirHostNet interactions.

    Attributes
    ----------
    df : pandas.DataFrame
        A pandas DataFrame representing VirHostNet interactions.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements extracted from the DataFrame.
    """
    def __init__(self, df, up_web_fallback=False):
        self.df = df
        self.up_web_fallback = up_web_fallback
        self.statements = []

    def extract_statements(self):
        for _, row in self.df.iterrows():
            stmt = process_row(row, up_web_fallback=self.up_web_fallback)
            if stmt:
                self.statements.append(stmt)


def process_row(row, up_web_fallback=False):
    """Process one row of the DataFrame into an INDRA Statement."""
    host_agent = get_agent_from_grounding(row['host_grounding'],
                                          up_web_fallback=up_web_fallback)
    vir_agent = get_agent_from_grounding(row['vir_grounding'],
                                         up_web_fallback=up_web_fallback)

    # There's a column that is always a - character
    assert row['dash'] == '-', row['dash']

    exp_method_id, exp_method_name = parse_psi_mi(row['exp_method'])
    int_type__id, int_type_name = parse_psi_mi(row['int_type'])

    assert row['host_tax'].startswith('taxid:'), row['host_tax']
    _, host_tax = row['host_tax'].split(':')
    assert row['vir_tax'].startswith('taxid:'), row['vir_tax']
    _, vir_tax = row['vir_tax'].split(':')
    assert row['score'].startswith('virhostnet-miscore:'), row['score']
    _, score = row['score'].split(':')
    score = float(score)

    source_ids = parse_source_ids(row['source_id'])

    annotations = {
        'exp_method': {'id': exp_method_id, 'name': exp_method_name},
        'int_type': {'id': int_type__id, 'name': int_type_name},
        'host_tax': host_tax,
        'vir_tax': vir_tax,
        'score': score,
        **source_ids,
    }

    text_refs = parse_text_refs(row['publication'])

    ev = Evidence(source_api='virhostnet', annotations=annotations,
                  text_refs=text_refs, pmid=text_refs.get('PMID'),
                  source_id=source_ids.get('virhostnet-rid'))

    stmt = Complex([host_agent, vir_agent], evidence=[ev])
    return stmt


def get_agent_from_grounding(grounding, up_web_fallback=False):
    """Return an INDRA Agent based on a grounding annotation."""
    db_ns, db_id = grounding.split(':')
    # Assume UniProt or RefSeq IDs
    assert db_ns in {'uniprotkb', 'refseq', 'ddbj/embl/genbank'}, db_ns
    if db_ns == 'uniprotkb':
        if '-' in db_id:
            up_id, feat_id = db_id.split('-')
            # Assume it's a feature ID
            assert feat_id.startswith('PRO'), feat_id
            db_refs = {'UP': up_id, 'UPPRO': feat_id}
        else:
            db_refs = {'UP': db_id}
    elif db_ns == 'refseq':
        db_refs = {'REFSEQ_PROT': db_id}
    else:
        db_refs = {'GENBANK': db_id}
    agent = Agent(db_id, db_refs=db_refs)
    standardized = standardize_agent_name(agent)
    if up_web_fallback:
        # Handle special case of unreviewed UP entries
        if not standardized and 'UP' in db_refs:
            name = uniprot_client.get_gene_name(db_refs['UP'],
                                                web_fallback=True)
            if name:
                agent.name = name
    return agent


def parse_psi_mi(psi_mi_str):
    """Parse a PSI-MI annotation into an ID and name pair."""
    # Example: psi-mi:"MI:0018"(two hybrid)
    match = re.match(r'psi-mi:"(.+)"\((.+)\)', psi_mi_str)
    mi_id, name = match.groups()
    return mi_id, name


def parse_text_refs(text_ref_str):
    """Parse a text reference annotation into a text_refs dict."""
    tr_ns, tr_id = text_ref_str.split(':')
    assert tr_ns == 'pubmed', text_ref_str
    if re.match(r'^\d+$', tr_id):
        return {'PMID': tr_id}
    else:
        match = re.match(r'^https\(//doi.org/(.+)\)$', tr_id)
        if not match:
            logger.warning('Failed to parse text ref: %s' % text_ref_str)
            return {}
        doi = match.groups()[0]
        return {'DOI': doi}


def parse_source_ids(source_id_str):
    """Parse VirHostNet source id annotations into a dict."""
    ids = source_id_str.split('|')
    assert len(ids) == 2
    ids_dict = {id.split(':')[0]: id.split(':')[1] for id in ids}
    return ids_dict
