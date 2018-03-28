import re
import sys
import logging
import json

from flask import Flask, request, abort, jsonify, Response
from flask_compress import Compress

from indra.db.util import get_statements_by_gene_role_type
from indra.statements import make_statement_camel
from indra.databases import hgnc_client

logger = logging.getLogger("db-api")

app = Flask(__name__)
Compress(app)


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
        ag = agent_param[:-5]
        ns = 'TEXT'

    if ns == 'HGNC-SYMBOL':
        ag = hgnc_client.get_hgnc_id(ag)
        ns = 'HGNC'

    return ag, ns


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
                    'namespace : select the namespace in which agents are identified.\n'
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
def get_statments():
    """Get some statements constrained by query."""
    logger.info("Got query for statements!")
    query_dict = request.args.copy()

    logger.info("Getting query details.")
    free_agents = [__process_agent(ag) for ag in query_dict.poplist('agent')]
    agents = {role: __process_agent(query_dict.pop(role, None)) for role
              in ['subject', 'object'] if query_dict.get(role) is not None}
    act_uncamelled = query_dict.pop('type', None)
    if act_uncamelled is not None:
        act = make_statement_camel(act_uncamelled)
    else:
        act = None
    if query_dict:
        abort(Response("Unrecognized query options; %s." %
                       list(query_dict.keys()),
                       400))

    if not any(agents.values()) and not free_agents:
        logger.error("No agents.")
        abort(Response(("No agents. Must have 'subject', 'object', or "
                        "'other'!\n"), 400))

    stmts = []
    logger.info("Getting statements...")
    for role, (agent, ns) in agents.items():
        logger.debug("Checking agent %s in namespace %s." % (agent, ns))
        if not stmts:
            # Get an initial list
            stmts = get_statements_by_gene_role_type(agent_id=agent,
                                                     agent_ns=ns,
                                                     role=role.upper(),
                                                     stmt_type=act,
                                                     do_stmt_count=False)
        elif role.lower() == 'subject':
            stmts = [s for s in stmts if len(s.agent_list()) > 0
                     and s.agent_list()[0].db_refs.get(ns) is not None
                     and s.agent_list()[0].db_refs.get(ns) == agent]
        elif role.lower() == 'object':
            stmts = [s for s in stmts if len(s.agent_list()) > 1
                     and s.agent_list()[1].db_refs.get(ns) is not None
                     and s.agent_list()[1].db_refs.get(ns) == agent]
    for agent, ns in free_agents:
        logger.debug("Checking agent %s in namespace %s." % (agent, ns))
        if not stmts:
            # Get an initial list
            stmts = get_statements_by_gene_role_type(agent_id=agent,
                                                     agent_ns=ns,
                                                     stmt_type=act,
                                                     do_stmt_count=False)
        else:
            stmts = [
                s for s in stmts if len(s.agent_list()) > 0 and agent in [
                    ag for other_agent in s.agent_list() if other_agent is not None
                    for ag in other_agent.db_refs.values()
                    ] + [
                    ag for ag in s.agent_list() if ag is None
                    ]
                ]

    resp = jsonify([stmt.to_json() for stmt in stmts])
    logger.info("Exiting with %d statements of nominal size %f MB."
                % (len(stmts), sys.getsizeof(resp.data)/1e6))
    return resp

if __name__ == '__main__':
    app.run()
