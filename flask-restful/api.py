import json
from bottle import route, run, request, post, default_app
from indra import trips, reach, bel, biopax
from indra.statements import *


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


if __name__ == '__main__':
    app = default_app()
    run(app)
