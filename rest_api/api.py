import os
import sys
import json
import base64
import logging
from bottle import route, run, request, default_app, response, static_file
from indra import get_config
from indra.sources import trips, reach, bel, biopax
from indra.sources import eidos, hume, cwms, sofia
from indra.databases import hgnc_client
from indra.statements import *
from indra.assemblers.pysb import PysbAssembler
import indra.assemblers.pysb.assembler as pysb_assembler
from indra.assemblers.cx import CxAssembler
from indra.assemblers.graph import GraphAssembler
from indra.assemblers.cyjs import CyJSAssembler
from indra.assemblers.sif import SifAssembler
from indra.assemblers.english import EnglishAssembler
import indra.tools.assemble_corpus as ac
from indra.databases import cbio_client
from indra.sources.indra_db_rest import get_statements
from indra.sources.ndex_cx.api import process_ndex_network
from indra.sources.reach.api import reach_nxml_url, reach_text_url
from indra.belief.wm_scorer import get_eidos_scorer
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap
from indra.pipeline import AssemblyPipeline


logger = logging.getLogger('rest_api')
logger.setLevel(logging.DEBUG)


def _return_stmts(stmts):
    if stmts:
        stmts_json = stmts_to_json(stmts)
        res = {'statements': stmts_json}
    else:
        res = {'statements': []}
    return res


def _stmts_from_proc(proc):
    if proc and proc.statements:
        stmts = stmts_to_json(proc.statements)
        res = {'statements': stmts}
    else:
        res = {'statements': []}
    return res


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

@route('/', method=['GET', 'OPTIONS'])
@allow_cors
def root():
    """API root."""
    if request.method == 'OPTIONS':
        return {}
    return {'This is the INDRA REST API. See documentation at '
            'http://www.indra.bio/rest_api/docs.'}

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
    return _stmts_from_proc(tp)


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
    return _stmts_from_proc(tp)
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
    given_url = body.get('url')
    config_url = get_config('REACH_TEXT_URL', failure_ok=True)
    # Order: URL given as an explicit argument in the request. Then any URL
    # set in the configuration. Then, unless offline is set, use the default
    # REACH web service URL.
    if 'url' in body:  # This is to take None if explicitly given
        url = given_url
    elif config_url:
        url = config_url
    elif not offline:
        url = reach_text_url
    else:
        url = None
    # If a URL is set, prioritize it over the offline setting
    if url:
        offline = False
    rp = reach.process_text(text, offline=offline, url=url)
    return _stmts_from_proc(rp)


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
    return _stmts_from_proc(rp)


@route('/reach/process_pmc', method=['POST', 'OPTIONS'])
@allow_cors
def reach_process_pmc():
    """Process PubMedCentral article and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    pmcid = body.get('pmcid')
    offline = True if body.get('offline') else False
    given_url = body.get('url')
    config_url = get_config('REACH_NXML_URL', failure_ok=True)
    # Order: URL given as an explicit argument in the request. Then any URL
    # set in the configuration. Then, unless offline is set, use the default
    # REACH web service URL.
    if 'url' in body:  # This is to take None if explicitly given
        url = given_url
    elif config_url:
        url = config_url
    elif not offline:
        url = reach_nxml_url
    else:
        url = None
    # If a URL is set, prioritize it over the offline setting
    if url:
        offline = False
    rp = reach.process_pmc(pmcid, offline=offline, url=url)
    return _stmts_from_proc(rp)

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
    return _stmts_from_proc(bp)


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
    return _stmts_from_proc(bp)


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
    return _stmts_from_proc(bp)


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
    return _stmts_from_proc(bp)


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
    return _stmts_from_proc(bp)


# TODO: document this in Swagger once webservice is available
@route('/eidos/process_text', method=['POST', 'OPTIONS'])
@allow_cors
def eidos_process_text():
    """Process text with EIDOS and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    req = request.body.read().decode('utf-8')
    body = json.loads(req)
    text = body.get('text')
    webservice = body.get('webservice')
    grounding_ns = body.get('grounding_ns')
    if not webservice:
        response.status = 400
        response.content_type = 'application/json'
        return json.dumps({'error': 'No web service address provided.'})
    ep = eidos.process_text(text, webservice=webservice,
                            grounding_ns=grounding_ns)
    return _stmts_from_proc(ep)


@route('/eidos/process_jsonld', method=['POST', 'OPTIONS'])
@allow_cors
def eidos_process_jsonld():
    """Process an EIDOS JSON-LD and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    eidos_json = body.get('jsonld')
    grounding_ns = body.get('grounding_ns')
    ep = eidos.process_json_str(eidos_json, grounding_ns=grounding_ns)
    return _stmts_from_proc(ep)


@route('/cwms/process_text', method=['POST', 'OPTIONS'])
@allow_cors
def cwms_process_text():
    """Process text with CWMS and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    text = body.get('text')
    cp = cwms.process_text(text)
    return _stmts_from_proc(cp)


@route('/hume/process_jsonld', method=['POST', 'OPTIONS'])
@allow_cors
def hume_process_jsonld():
    """Process Hume JSON-LD and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    jsonld_str = body.get('jsonld')
    jsonld = json.loads(jsonld_str)
    hp = hume.process_jsonld(jsonld)
    return _stmts_from_proc(hp)


# TODO: the transfer of the binary file uploaded here will need to be
# implemented
#@route('/sofia/process_table', method=['POST', 'OPTIONS'])
#@allow_cors
#def sofia_process_table():
#    """Process Sofia table and return INDRA Statements."""
#    if request.method == 'OPTIONS':
#        return {}
#    response = request.body.read().decode('utf-8')
#    body = json.loads(response)
#    table = body.get('table')
#    sp = sofia.process_table(table)
#    return _stmts_from_proc(sp)
#

@route('/sofia/process_text', method=['POST', 'OPTIONS'])
@allow_cors
def sofia_process_text():
    """Process text with Sofia and return INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    text = body.get('text')
    auth = body.get('auth')
    sp = sofia.process_text(text, auth=auth)
    return _stmts_from_proc(sp)




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
    try:
        for m in pa.model.monomers:
            pysb_assembler.set_extended_initial_condition(pa.model, m, 0)
    except Exception as e:
        logger.exception(e)

    if not export_format:
        model_str = pa.print_model()
    elif export_format in ('kappa_im', 'kappa_cm'):
        fname = 'model_%s.png' % export_format
        root = os.path.dirname(os.path.abspath(fname))
        graph = pa.export_model(format=export_format, file_name=fname)
        with open(fname, 'rb') as fh:
            data = 'data:image/png;base64,%s' % \
                base64.b64encode(fh.read()).decode() 
            return {'image': data}
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
@route('/share_model_ndex', method=['POST', 'OPTIONS'])
@allow_cors
def share_model_ndex():
    """Upload the model to NDEX"""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_str = body.get('stmts')
    stmts_json = json.loads(stmts_str)
    stmts = stmts_from_json(stmts_json["statements"])
    ca = CxAssembler(stmts)
    for n, v in body.items():
        ca.cx['networkAttributes'].append({'n': n, 'v': v, 'd': 'string'})
    ca.make_model()
    network_id = ca.upload_model(private=False)
    return {'network_id': network_id}


@route('/fetch_model_ndex', method=['POST', 'OPTIONS'])
@allow_cors
def fetch_model_ndex():
    """Download model and associated pieces from NDEX"""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    network_id = body.get('network_id')
    cx = process_ndex_network(network_id)
    network_attr = [x for x in cx.cx if x.get('networkAttributes')]
    network_attr = network_attr[0]['networkAttributes']
    keep_keys = ['txt_input', 'parser',
                 'model_elements', 'preset_pos', 'stmts',
                 'sentences', 'evidence', 'cell_line', 'mrna', 'mutations']
    stored_data = {}
    for d in network_attr:
        if d['n'] in keep_keys:
            stored_data[d['n']] = d['v']
    return stored_data


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
    return _return_stmts(stmts_out)


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
    return _return_stmts(stmts_out)


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
    scorer = body.get('scorer')
    return_toplevel = body.get('return_toplevel')
    if scorer == 'wm':
        belief_scorer = get_eidos_scorer()
    else:
        belief_scorer = None
    stmts_out = ac.run_preassembly(stmts, belief_scorer=belief_scorer,
                                   return_toplevel=return_toplevel)
    return _return_stmts(stmts_out)


@route('/preassembly/map_ontologies', method=['POST', 'OPTIONS'])
@allow_cors
def map_ontologies():
    """Run ontology mapping on a list of INDRA Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    stmts = stmts_from_json(stmts_json)
    om = OntologyMapper(stmts, wm_ontomap, scored=True, symmetric=False)
    om.map_statements()
    return _return_stmts(stmts)


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
    return _return_stmts(stmts_out)


@route('/preassembly/filter_grounded_only', method=['POST', 'OPTIONS'])
@allow_cors
def filter_grounded_only():
    """Filter to grounded Statements only."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    score_threshold = body.get('score_threshold')
    if score_threshold is not None:
        score_threshold = float(score_threshold)
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.filter_grounded_only(stmts, score_threshold=score_threshold)
    return _return_stmts(stmts_out)


@route('/preassembly/filter_belief', method=['POST', 'OPTIONS'])
@allow_cors
def filter_belief():
    """Filter to beliefs above a given threshold."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    belief_cutoff = body.get('belief_cutoff')
    if belief_cutoff is not None:
        belief_cutoff = float(belief_cutoff)
    stmts = stmts_from_json(stmts_json)
    stmts_out = ac.filter_belief(stmts, belief_cutoff)
    return _return_stmts(stmts_out)


@route('/preassembly/pipeline', method=['POST', 'OPTIONS'])
@allow_cors
def run_pipeline():
    """Run an assembly pipeline for a list of Statements."""
    if request.method == 'OPTIONS':
        return {}
    response = request.body.read().decode('utf-8')
    body = json.loads(response)
    stmts_json = body.get('statements')
    pipeline = body.get('pipeline')
    stmts = stmts_from_json(stmts_json)
    ap = AssemblyPipeline(pipeline)
    stmts_out = ap.run(stmts)
    return _return_stmts(stmts_out)


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
                               simple_response=True, **kwargs)
        return stmts

    stmts_out = _get_matching_stmts(stmt)
    agent_name_list = [ag.name for ag in stmt.agent_list()]
    stmts_out = stmts = ac.filter_concept_names(stmts_out, agent_name_list, 'all')
    return _return_stmts(stmts_out)


app = default_app()

if __name__ == '__main__':
    run(app, host='0.0.0.0', port='8080')
