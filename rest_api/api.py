import json
import logging
from bottle import route, run, request, default_app
from indra import trips, reach, bel, biopax
from indra.statements import *
from indra.assemblers import PysbAssembler, CxAssembler, GraphAssembler,\
                             CyJSAssembler

logger = logging.getLogger('rest_api')
logger.setLevel(logging.DEBUG)

###### INPUT PROCESSING #######

### TRIPS ###
@route('/trips/process_text', method='POST')
def trips_process_text():
    body = json.load(request.body)
    text = body.get('text')
    tp = trips.process_text(text)
    if tp and tp.statements:
        stmts = _stmts_to_json(tp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res
################

### REACH ###
@route('/reach/process_text', method='POST')
def reach_process_text():
    body = json.load(request.body)
    text = body.get('text')
    rp = reach.process_text(text)
    if rp and rp.statements:
        stmts = _stmts_to_json(rp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res

@route('/reach/process_pmc', method='POST')
def reach_process_pmc():
    body = json.load(request.body)
    pmcid = body.get('pmcid')
    rp = reach.process_pmc(pmcid)
    if rp and rp.statements:
        stmts = _stmts_to_json(rp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res
##################

### BEL ###
@route('/bel/process_ndex_neighborhood', method='POST')
def bel_process_ndex_neighborhood():
    body = json.load(request.body)
    genes = body.get('genes')
    logger.debug(genes)
    bp = bel.process_ndex_neighborhood(genes)
    if bp and bp.statements:
        stmts = _stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res

### BioPAX ###
@route('/biopax/process_pc_pathsbetween', method='POST')
def biopax_process_pc_pathsbetween():
    body = json.load(request.body)
    genes = body.get('genes')
    logger.debug(genes)
    bp = bel.process_ndex_neighborhood(genes)
    if bp and bp.statements:
        stmts = _stmts_to_json(bp.statements)
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res

### OUTPUT ASSEMBLY ####################

### PYSB ###

@route('/assemblers/pysb', method='POST')
def assemble_pysb():
    body = json.load(request.body)
    stmts_str = body.get('statements')
    stmts = _stmts_from_json(stmts_str)
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model()
    model_str = pa.print_model()
    res = {'model': model_str}
    return res

### CX ###

@route('/assemblers/cx', method='POST')
def assemble_cx():
    body = json.load(request.body)
    stmts_str = body.get('statements')
    stmts = _stmts_from_json(stmts_str)
    ca = CxAssembler(stmts)
    model_str = ca.make_model()
    res = {'model': model_str}
    return res

### GRAPH ###

@route('/assemblers/graph', method='POST')
def assemble_grapg():
    body = json.load(request.body)
    stmts_str = body.get('statements')
    stmts = _stmts_from_json(stmts_str)
    ga = GraphAssembler(stmts)
    model_str = ga.make_model()
    res = {'model': model_str}
    return res

### CyJS ###

@route('/assemblers/cyjs', method='POST')
def assemble_cyjs():
    body = json.load(request.body)
    stmts_str = body.get('statements')
    stmts = _stmts_from_json(stmts_str)
    cja = CyJSAssembler()
    cja.add_statements(stmts)
    cja.make_model(grouping=True,
                   drop_virtual_edges=False,
                   add_edge_weights=True)
    line = body.get('line')
    if line is not None:
        cja.set_context(cell_type = line,
                        bin_expression = True,
                        n_bins = 9)
    model_str = cja.print_cyjs()
    return model_str

def _stmts_to_json(stmts):
    stj = json.dumps([json.loads(st.to_json()) for st in stmts])
    return stj

def _stmts_from_json(stmts_str):
    stmts_json = json.loads(stmts_str)
    stmts = [Statement.from_json(json.dumps(stj)) for stj in stmts_json]
    return stmts

if __name__ == '__main__':
    app = default_app()
    run(app)

