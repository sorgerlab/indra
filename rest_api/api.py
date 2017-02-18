import json
from bottle import route, run, request, post, default_app
from indra import trips, reach, bel, biopax
from indra.statements import *
from indra.assemblers import PysbAssembler

###### INPUT PROCESSING #######

### TRIPS ###
@route('/trips/process_text', method='POST')
def trips_process_text():
    body = json.load(request.body)
    text = body.get('text')
    tp = trips.process_text(text)
    if tp and tp.statements:
        stmts = json.dumps([json.loads(st.to_json()) for st
                            in tp.statements])
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
        stmts = json.dumps([json.loads(st.to_json()) for st
                            in rp.statements])
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
        stmts = json.dumps([json.loads(st.to_json()) for st
                            in rp.statements])
        res = {'statements': stmts}
        return res
    else:
        res = {'statements': []}
    return res
##################

### OUTPUT ASSEMBLY ####################

### PYSB ###

@route('/assemblers/pysb', method='POST')
def assemble_pysb():
    body = json.load(request.body)
    stmts_str = body.get('statements')
    stmts_json = json.loads(stmts_str)
    stmts = [Statement.from_json(json.dumps(stj)) for stj in stmts_json]
    pa = PysbAssembler()
    pa.add_statements(stmts)
    pa.make_model()
    model_str = pa.print_model()
    res = {'model': model_str}
    return res

if __name__ == '__main__':
    app = default_app()
    run(app)
