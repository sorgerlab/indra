import sys
import json
import logging
from bottle import route, run, request, default_app, response
from indra.sources import trips, reach, bel, biopax
from indra.databases import hgnc_client
from indra.statements import *
from indra.assemblers import PysbAssembler, CxAssembler, GraphAssembler,\
    CyJSAssembler, SifAssembler, EnglishAssembler
import indra.tools.assemble_corpus as ac
from indra.databases import cbio_client
from indra.sources.indra_db_rest import get_statements

logger = logging.getLogger('rest_api')
logger.setLevel(logging.DEBUG)


#   ALLOW CORS   #
def allow_cors(func):
    """This is a decorator which enable CORS for the specified endpoint."""
    def wrapper(*args, **kwargs):
        response.headers['Access-Control-Allow-Origin'] = '*'
        response.headers['Access-Control-Allow-Methods'] = \
            'PUT, GET, POST, DELETE, OPTIONS'
        response.headers['Access-Control-Allow-Headers'] = \
            'Origin, Accept, Content-Type, X-Requested-With, X-CSRF-Token'
        return func(*args, **kwargs)
    return wrapper

#     INPUT PROCESSING     #


#   TRIPS   #
@route('/trips/process_text', method=['POST', 'OPTIONS'])
@allow_cors
def trips_process_text():
    """Process text with TRIPS and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    text = body.get('text')
    tp = trips.process_text(text)
    if tp and tp.statements:
        stmts = stmts_to_json(tp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


@route('/trips/process_xml', method=['POST', 'OPTIONS'])
@allow_cors
def trips_process_xml():
    """Process TRIPS EKB XML and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    xml_str = body.get('xml_str')
    tp = trips.process_xml(xml_str)
    if tp and tp.statements:
        stmts = stmts_to_json(tp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res
################


#   REACH   #
@route('/reach/process_text', method=['POST', 'OPTIONS'])
@allow_cors
def reach_process_text():
    """Process text with REACH and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    text = body.get('text')
    offline = True if body.get('offline') else False
    rp = reach.process_text(text, offline=offline)
    if rp and rp.statements:
        stmts = stmts_to_json(rp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


@route('/reach/process_json', method=['POST', 'OPTIONS'])
@allow_cors
def reach_process_json():
    """Process REACH json and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    json_str = body.get('json')
    rp = reach.process_json_str(json_str)
    if rp and rp.statements:
        stmts = stmts_to_json(rp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


@route('/reach/process_pmc', method=['POST', 'OPTIONS'])
@allow_cors
def reach_process_pmc():
    """Process PubMedCentral article and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    pmcid = body.get('pmcid')
    rp = reach.process_pmc(pmcid)
    if rp and rp.statements:
        stmts = stmts_to_json(rp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res
##################


#   BEL   #
@route('/bel/process_pybel_neighborhood', method=['POST', 'OPTIONS'])
@allow_cors
def bel_process_pybel_neighborhood():
    """Process BEL Large Corpus neighborhood and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    genes = body.get('genes')
    bp = bel.process_pybel_neighborhood(genes)
    if bp and bp.statements:
        stmts = stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


@route('/bel/process_belrdf', method=['POST', 'OPTIONS'])
@allow_cors
def bel_process_belrdf():
    """Process BEL RDF and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    belrdf = body.get('belrdf')
    bp = bel.process_belrdf(belrdf)
    if bp and bp.statements:
        stmts = stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


#   BioPAX   #
@route('/biopax/process_pc_pathsbetween', method=['POST', 'OPTIONS'])
@allow_cors
def biopax_process_pc_pathsbetween():
    """Process PathwayCommons paths between genes, return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    genes = body.get('genes')
    bp = biopax.process_pc_pathsbetween(genes)
    if bp and bp.statements:
        stmts = stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


@route('/biopax/process_pc_pathsfromto', method=['POST', 'OPTIONS'])
@allow_cors
def biopax_process_pc_pathsfromto():
    """Process PathwayCommons paths from-to genes, return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    source = body.get('source')
    target = body.get('target')
    bp = biopax.process_pc_pathsfromto(source, target)
    if bp and bp.statements:
        stmts = stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res


@route('/biopax/process_pc_neighborhood', method=['POST', 'OPTIONS'])
@allow_cors
def biopax_process_pc_neighborhood():
    """Process PathwayCommons neighborhood, return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    genes = body.get('genes')
    bp = biopax.process_pc_neighborhood(genes)
    if bp and bp.statements:
        stmts = stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res

#   OUTPUT ASSEMBLY   #


#   PYSB   #
@route('/assemblers/pysb', method=['POST', 'OPTIONS'])
@allow_cors
def assemble_pysb():
    """Assemble INDRA Statements and return PySB model string."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    export_format = body.get('export_format')
    stmts = stmts_from_json(stmts_json)
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model()
    if not export_format:
        model_str = pa.print_model()
    else:
        try:
            model_str = pa.export_model(format=export_format)
        except Exception as e:
            logger.exception(e)
            model_str = ''
    res = {'model': model_str}
    return res


#   CX   #
@route('/assemblers/cx', method=['POST', 'OPTIONS'])
@allow_cors
def assemble_cx():
    """Assemble INDRA Statements and return CX network json."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    ca = CxAssembler(stmts)
    model_str = ca.make_model()
    res = {'model': model_str}
    return res


#   SHARE CX   #
@route('/share_model', method=['POST', 'OPTIONS'])
@allow_cors
def share_model():
    """Upload the model to NDEX"""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    cyjs_model_str = body.get('cyjs_model')
    stmts = stmts_from_json(stmts_json)
    ca = CxAssembler(stmts)
    ca.cx['networkAttributes'].append({'n': 'cyjs_model',
                                       'v': cyjs_model_str,
                                       'd': 'string'})
    ca.make_model()
    network_id = ca.upload_model(private=False)
    return {'network_id': network_id}


#  GRAPH   #
@route('/assemblers/graph', method=['POST', 'OPTIONS'])
@allow_cors
def assemble_graph():
    """Assemble INDRA Statements and return Graphviz graph dot string."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    ga = GraphAssembler(stmts)
    model_str = ga.make_model()
    res = {'model': model_str}
    return res


#   CyJS   #
@route('/assemblers/cyjs', method=['POST', 'OPTIONS'])
@allow_cors
def assemble_cyjs():
    """Assemble INDRA Statements and return Cytoscape JS network."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    cja = CyJSAssembler()
    cja.add_statements(stmts)
    cja.make_model(grouping=True)
    model_str = cja.print_cyjs_graph()
    return model_str


#   English   #
@route('/assemblers/english', method=['POST', 'OPTIONS'])
@allow_cors
def assemble_english():
    """Assemble each statement into """
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    sentences = {}
    for st in stmts:
        enga = EnglishAssembler()
        enga.add_statements([st])
        model_str = enga.make_model()
        sentences[st.uuid] = model_str
    res = {'sentences': sentences}
    return res


@route('/assemblers/sif/loopy', method=['POST', 'OPTIONS'])
@allow_cors
def assemble_loopy():
    """Assemble INDRA Statements into a Loopy model using SIF Assembler."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    sa = SifAssembler(stmts)
    sa.make_model(use_name_as_key=True)
    model_str = sa.print_loopy(as_url=True)
    res = {'loopy_url': model_str}
    return res


#   CCLE mRNA amounts   #
@route('/databases/cbio/get_ccle_mrna', method=['POST', 'OPTIONS'])
@allow_cors
def get_ccle_mrna_levels():
    """Get CCLE mRNA amounts using cBioClient"""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    gene_list = body.get('gene_list')
    cell_lines = body.get('cell_lines')
    mrna_amounts = cbio_client.get_ccle_mrna(gene_list, cell_lines)
    res = {'mrna_amounts': mrna_amounts}
    return res


#   CCLE CNA   #
@route('/databases/cbio/get_ccle_cna', method=['POST', 'OPTIONS'])
@allow_cors
def get_ccle_cna():
    """Get CCLE CNA
    -2 = homozygous deletion
    -1 = hemizygous deletion
     0 = neutral / no change
     1 = gain
     2 = high level amplification
    """
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    gene_list = body.get('gene_list')
    cell_lines = body.get('cell_lines')
    cna = cbio_client.get_ccle_cna(gene_list, cell_lines)
    res = {'cna': cna}
    return res


@route('/databases/cbio/get_ccle_mutations', method=['POST', 'OPTIONS'])
@allow_cors
def get_ccle_mutations():
    """Get CCLE mutations
    returns the amino acid changes for a given list of genes and cell lines
    """
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    gene_list = body.get('gene_list')
    cell_lines = body.get('cell_lines')
    mutations = cbio_client.get_ccle_mutations(gene_list, cell_lines)
    res = {'mutations': mutations}
    return res


@route('/preassembly/map_grounding', method=['POST', 'OPTIONS'])
@allow_cors
def map_grounding():
    """Map grounding on a list of INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.map_grounding(stmts)
    if stmts_out:
        stmts_json = stmts_to_json(stmts_out)
        res = {'statements': stmts_json}
        return res
    else:
        res = {'statements': []}
    return res


@route('/preassembly/map_sequence', method=['POST', 'OPTIONS'])
@allow_cors
def map_sequence():
    """Map sequence on a list of INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.map_sequence(stmts)
    if stmts_out:
        stmts_json = stmts_to_json(stmts_out)
        res = {'statements': stmts_json}
        return res
    else:
        res = {'statements': []}
    return res


@route('/preassembly/run_preassembly', method=['POST', 'OPTIONS'])
@allow_cors
def run_preassembly():
    """Run preassembly on a list of INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.run_preassembly(stmts)
    if stmts_out:
        stmts_json = stmts_to_json(stmts_out)
        res = {'statements': stmts_json}
        return res
    else:
        res = {'statements': []}
    return res


@route('/preassembly/filter_by_type', method=['POST', 'OPTIONS'])
@allow_cors
def filter_by_type():
    """Filter to a given INDRA Statement type."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmt_type_str = body.get('type')
    stmt_type_str = stmt_type_str.capitalize()
    stmt_type = getattr(sys.modules[__name__], stmt_type_str)
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.filter_by_type(stmts, stmt_type)
    if stmts_out:
        stmts_json = stmts_to_json(stmts_out)
        res = {'statements': stmts_json}
        return res
    else:
        res = {'statements': []}
    return res


@route('/preassembly/filter_grounded_only', method=['POST', 'OPTIONS'])
@allow_cors
def filter_grounded_only():
    """Filter to grounded Statements only."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.filter_grounded_only(stmts)
    if stmts_out:
        stmts_json = stmts_to_json(stmts_out)
        res = {'statements': stmts_json}
        return res
    else:
        res = {'statements': []}
    return res


@route('/indra_db_rest/get_evidence', method=['POST', 'OPTIONS'])
@allow_cors
def get_evidence_for_stmts():
    if request.method == 'OPTIONS':
        return {}
    req = request.body.read().decode('utf-8')
    body = json.loads(req)
    stmt_json = body.get('statement')
    stmt = Statement._from_json(stmt_json)

    def _get_agent_ref(agent):
        """Get the preferred ref for an agent for db web api."""
        if agent is None:
            return None
        ag_hgnc_id = hgnc_client.get_hgnc_id(agent.name)
        if ag_hgnc_id is not None:
            return ag_hgnc_id + "@HGNC"
        db_refs = agent.db_refs
        for namespace in ['HGNC', 'FPLX', 'CHEBI', 'TEXT']:
            if namespace in db_refs.keys():
                # TODO: Remove. This is a temporary workaround.
                if namespace == 'FPLX':
                    return '%s@%s' % (db_refs[namespace], 'BE')
                return '%s@%s' % (db_refs[namespace], namespace)
        return '%s@%s' % (agent.name, 'TEXT')

    def _get_matching_stmts(stmt_ref):
        # Filter by statement type.
        stmt_type = stmt_ref.__class__.__name__
        agent_name_list = [_get_agent_ref(ag) for ag in stmt_ref.agent_list()]
        non_binary_statements = (Complex, SelfModification, ActiveForm)
        # TODO: We should look at more than just the agent name.
        # Doing so efficiently may require changes to the web api.
        if isinstance(stmt_ref, non_binary_statements):
            agent_list = [ag_name for ag_name in agent_name_list
                          if ag_name is not None]
            kwargs = {}
        else:
            agent_list = []
            kwargs = {k: v for k, v in zip(['subject', 'object'],
                                           agent_name_list)}
            if not any(kwargs.values()):
                return []
            print(agent_list)
        stmts = get_statements(agents=agent_list, stmt_type=stmt_type,
                               **kwargs)
        return stmts

    stmts_out = _get_matching_stmts(stmt)
    agent_name_list = [ag.name for ag in stmt.agent_list()]
    stmts_out = stmts = ac.filter_concept_names(stmts_out, agent_name_list, 'all')
    if stmts_out:
        stmts_json = stmts_to_json(stmts_out)
        res = {'statements': stmts_json}
        return res
    else:
        res = {'statements': []}
    return res



app = default_app()

if __name__ == '__main__':
    run(app, host='0.0.0.0', port='8080')
