import sys
import json
import logging
from functools import wraps

from flask import Flask, request, abort, Response
from flask_compress import Compress
from flask_cors import CORS

from indra.db.client import get_statement_jsons_from_agents, \
    get_statement_jsons_from_hashes, get_statement_jsons_from_papers, has_auth
from indra.statements import make_statement_camel
from indra.databases import hgnc_client
from indra.util import batch_iter

logger = logging.getLogger("db-api")

app = Flask(__name__)
Compress(app)
CORS(app)


MAX_STATEMENTS = int(1e3)


class DbAPIError(Exception):
    pass


def __process_agent(agent_param):
    """Get the agent id and namespace from an input param."""
    if not agent_param.endswith('TEXT'):
        param_parts = agent_param.split('@')
        if len(param_parts) == 2:
            ag, ns = param_parts
        elif len(param_parts) == 1:
            ag = agent_param
            ns = 'HGNC-SYMBOL'
        else:
            raise DbAPIError('Unrecognized agent spec: \"%s\"' % agent_param)
    else:
        ag = agent_param[:-5]
        ns = 'TEXT'

    if ns == 'HGNC-SYMBOL':
        original_ag = ag
        ag = hgnc_client.get_hgnc_id(original_ag)
        if ag is None and 'None' not in agent_param:
            raise DbAPIError('Invalid agent name: \"%s\"' % original_ag)
        ns = 'HGNC'

    return ag, ns


def get_source(ev_json):
    notes = ev_json.get('annotations')
    if notes is None:
        return
    src = notes.get('content_source')
    if src is None:
        return
    return src.lower()


REDACT_MESSAGE = '[MISSING/INVALID API KEY: limited to 200 char for Elsevier]'


def _query_wrapper(f):
    @wraps(f)
    def decorator():
        logger.info("Got query for %s!" % f.__name__)

        query_dict = request.args.copy()
        offs = query_dict.pop('offset', None)
        ev_limit = query_dict.pop('ev_limit', 10)
        best_first = query_dict.pop('best_first', True)
        do_stream_str = query_dict.pop('stream', 'false')
        do_stream = True if do_stream_str == 'true' else False

        api_key = query_dict.pop('api-key', None)

        result = f(query_dict, offs, MAX_STATEMENTS, ev_limit, best_first)

        # Redact elsevier content for those without permission.
        if api_key is None or not has_auth(api_key):
            for stmt_json in result['statements'].values():
                for ev_json in stmt_json['evidence']:
                    if get_source(ev_json) == 'elsevier':
                        text = ev_json['text']
                        if len(text) > 200:
                            ev_json['text'] = text[:200] + REDACT_MESSAGE
        result['offset'] = offs
        result['evidence_limit'] = ev_limit
        result['statement_limit'] = MAX_STATEMENTS

        if do_stream:
            # Returning a generator should stream the data.
            resp_json_bts = json.dumps(result)
            gen = batch_iter(resp_json_bts, 10000)
            resp = Response(gen, mimetype='application/json')
        else:
            resp = Response(json.dumps(result), mimetype='application/json')
        logger.info("Exiting with %d statements with %d evidence of size %f MB."
                    % (len(result['statements']), result['total_evidence'],
                       sys.getsizeof(resp.data)/1e6))
        return resp
    return decorator


@app.route('/')
def welcome():
    logger.info("Got request for welcome info.")
    return Response("Welcome the the INDRA database webservice!\n"
                    "\n"
                    "Use modes:\n"
                    "----------\n"
                    "/            - (you are here) Welcome page.\n"
                    "/statements  - Get detailed instructions for querying "
                    "statements.\n"
                    "/statements/?<query_string> - Get a list of statement "
                    "jsons."
                    "\n")


@app.route('/statements', methods=['GET'])
def get_statements_query_format():
    return Response('To get a list of statements, include a query after '
                    '/statements/ with the following keys:\n\n'
                    'type : the type of interaction (e.g. Phosphorylation)\n'
                    'namespace : select the namespace in which agents are '
                    'identified.\n'
                    '[subject, object, agent] : the agents, indicated by '
                    'their role. Note that at least one agent is needed in '
                    'a query. If agent is use, that agent will be matched to '
                    'any argument in the statement.\n\n'
                    'For example: /statements/?subject=MAP2K1&object=MAPK1'
                    '&type=Phosphorylation'
                    'Most statements have a subject and an object, but unary '
                    'and n-ary statements should have agents specified by '
                    '\"other\".')


@app.route('/statements/', methods=['GET'])
@_query_wrapper
def get_statements(query_dict, offs, max_stmts, ev_limit, best_first):
    """Get some statements constrained by query."""
    logger.info("Getting query details.")
    try:
        # Get the agents without specified locations (subject or object).
        free_agents = [__process_agent(ag)
                       for ag in query_dict.poplist('agent')]
        ofaks = {k for k in query_dict.keys() if k.startswith('agent')}
        free_agents += [__process_agent(query_dict.pop(k)) for k in ofaks]

        # Get the agents with specified roles.
        roled_agents = {role: __process_agent(query_dict.pop(role))
                        for role in ['subject', 'object']
                        if query_dict.get(role) is not None}
    except DbAPIError as e:
        logger.exception(e)
        abort(Response('Failed to make agents from names: %s\n' % str(e), 400))
        return

    # Get the raw name of the statement type (we allow for variation in case).
    act_raw = query_dict.pop('type', None)

    # Fix the case, if we got a statement type.
    act = None if act_raw is None else make_statement_camel(act_raw)

    # If there was something else in the query, there shouldn't be, so someone's
    # probably confused.
    if query_dict:
        abort(Response("Unrecognized query options; %s." % list(query_dict.keys()),
                       400))
        return

    # Make sure we got SOME agents. We will not simply return all
    # phosphorylations, or all activations.
    if not any(roled_agents.values()) and not free_agents:
        logger.error("No agents.")
        abort(Response(("No agents. Must have 'subject', 'object', or "
                        "'other'!\n"), 400))

    # Check to make sure none of the agents are None.
    assert None not in roled_agents.values() and None not in free_agents, \
        "None agents found. No agents should be None."

    # Now find the statements.
    logger.info("Getting statements...")
    agent_iter = [(role, ag_dbid, ns)
                  for role, (ag_dbid, ns) in roled_agents.items()]
    agent_iter += [(None, ag_dbid, ns) for ag_dbid, ns in free_agents]

    result = \
        get_statement_jsons_from_agents(agent_iter, stmt_type=act, offset=offs,
                                        max_stmts=max_stmts, ev_limit=ev_limit,
                                        best_first=best_first)
    return result


@app.route('/statements/from_hashes', methods=['POST'])
@_query_wrapper
def get_statements_by_hash(query_dict, offs, max_stmts, ev_limit, best_first):
    hashes = request.json.get('hashes')
    if not hashes:
        logger.error("No hashes provided!")
        abort(Response("No hashes given!", 400))
    if len(hashes) > max_stmts:
        logger.error("Too many hashes given!")
        abort(Response("Too many hashes given, %d allowed." % max_stmts,
                       400))

    result = get_statement_jsons_from_hashes(hashes, max_stmts=max_stmts,
                                             offset=offs, ev_limit=ev_limit,
                                             best_first=best_first)
    return result


@app.route('/papers/', methods=['GET'])
@_query_wrapper
def get_paper_statements(query_dict, offs, max_stmts, ev_limit, best_first):
    """Get and preassemble statements from a paper given by pmid."""
    # Get the paper id.
    id_val = query_dict.get('id')
    if id_val is None:
        logger.error("No id provided!")
        abort(Response("No id in request!", 400))

    # Get the id type, if given.
    id_type = query_dict.get('type', 'pmid')

    # Now get the statements.
    logger.info('Getting statements for %s=%s...' % (id_type, id_val))
    if id_type in ['trid', 'tcid']:
        id_val = int(id_val)
    result = get_statement_jsons_from_papers([(id_type, id_val)],
                                             max_stmts=max_stmts,
                                             offset=offs, ev_limit=ev_limit,
                                             best_first=best_first)
    return result


if __name__ == '__main__':
    app.run()
