from flask import Flask, request, abort, jsonify, Response
from indra.db.util import get_statements_by_gene_role_type
import logging
import re
from indra.statements import make_statement_camel

logger = logging.getLogger("db-api")

app = Flask(__name__)

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
    free_agents = query_dict.poplist('agent')
    agents = {role: query_dict.pop(role, None) for role
              in ['subject', 'object']}
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
    for role, agent in agents.items():
        if agent is not None:
            if not stmts:
                # Get an initial list
                stmts = get_statements_by_gene_role_type(agent_id=agent,
                                                         role=role.upper(),
                                                         stmt_type=act,
                                                         do_stmt_count=False)
            elif role.lower() == 'subject':
                stmts = [s for s in stmts if len(s.agent_list()) > 0
                         and s.agent_list()[0].name == agent]
            elif role.lower() == 'object':
                stmts = [s for s in stmts if len(s.agent_list()) > 1
                         and s.agent_list()[1].name == agent]
    for agent in free_agents:
        if not stmts and agent is not None:
            # Get an initial list
            stmts = get_statements_by_gene_role_type(agent_id=agent,
                                                     stmt_type=act,
                                                     do_stmt_count=False)
        else:
            stmts = [s for s in stmts if len(s.agent_list()) > 0
                     and agent in [ag.name for ag in s.agent_list()
                                   if ag is not None]]

    logger.info("Exiting with %d statements." % len(stmts))
    return jsonify([stmt.to_json() for stmt in stmts])

if __name__ == '__main__':
    app.run()