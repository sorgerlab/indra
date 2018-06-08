import sys
import logging

from flask import Flask, request, abort, jsonify, Response, make_response
from flask_compress import Compress

from indra.db.client import get_statements_by_gene_role_type, \
    get_statements_by_paper, get_evidence
from indra.statements import make_statement_camel
from indra.databases import hgnc_client

logger = logging.getLogger("db-api")

app = Flask(__name__)
Compress(app)


MAX_STATEMENTS = int(1e4)


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


def _filter_statements(stmts_in, ns, ag_id, role=None):
    """Return statements filtered to ones where agent is at given position."""
    # Make sure the role is good.
    assert role in ['SUBJECT', 'OBJECT', None], \
        "The value of role must be either 'SUBJECT', 'OBJECT', or None."

    # Map the role to an index:
    agent_pos = 0 if role == 'SUBJECT' else 1 if role == 'OBJECT' else None

    # Filter the statements.
    stmts_out = []
    for stmt in stmts_in:
        # Make sure the statement has enough agents to get the one at the
        # position of interest e.g. has only 1 agent but the agent_pos is not 0
        agents = stmt.agent_list()
        if agent_pos is not None:
            if len(agents) <= agent_pos:
                continue
            # Get the agent at the position of interest and make sure it's an
            # actual Agent
            agent = agents[agent_pos]
            if agent is not None:
                # Check if the db_refs for the namespace of interest matches the
                # value
                if agent.db_refs.get(ns) == ag_id:
                    stmts_out.append(stmt)
        else:
            # Search through all the agents looking for an agent with a matching
            # db ref.
            for agent in agents:
                if agent is not None and agent.db_refs.get(ns) == ag_id:
                    stmts_out.append(stmt)
                    break
    return stmts_out


def _get_relevant_statements(stmts, ag_id, ns, stmt_type, role=None):
    """Get statements that are relevant to the criteria included.

    If stmts is an empty list or None (bool evaluates to False), then get a
    matching set of statements from the database. Otherwise, filter down the
    existing list of statements.
    """
    logger.debug("Checking agent %s in namespace %s." % (ag_id, ns))
    # TODO: This is a temporary measure, remove ASAP.
    if role:
        role = role.upper()

    if not stmts:
        # Get an initial list
        stmts = get_statements_by_gene_role_type(agent_id=ag_id, agent_ns=ns,
                                                 role=role, stmt_type=stmt_type,
                                                 do_stmt_count=False,
                                                 with_evidence=False,
                                                 with_support=False)
    else:
        stmts = _filter_statements(stmts, ns, ag_id, role)

    return stmts


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
def get_statements():
    """Get some statements constrained by query."""
    logger.info("Got query for statements!")
    query_dict = request.args.copy()
    limit_behavior = query_dict.pop('on_limit', 'sample')
    if limit_behavior not in ['sample', 'truncate', 'error']:
        abort(Response('Unrecognized value for on_limit: %s' % limit_behavior,
                       400))

    logger.info("Getting query details.")
    try:
        # Get the agents without specified locations (subject or object).
        free_agents = [__process_agent(ag)
                       for ag in query_dict.poplist('agent')]

        # Get the agents with specified roles.
        roled_agents = {role: __process_agent(query_dict.pop(role))
                        for role in ['subject', 'object']
                        if query_dict.get(role) is not None}
    except DbAPIError as e:
        logger.exception(e)
        abort(Response('Failed to make agents from names: %s\n' % str(e), 400))
        return

    # Get the raw name of the statement type (we allow for variation in case).
    act_uncamelled = query_dict.pop('type', None)

    # Fix the case, if we got a statement type.
    if act_uncamelled is not None:
        act = make_statement_camel(act_uncamelled)
    else:
        act = None

    # If there was something else in the query, there shouldn't be, so someone's
    # probably confused.
    if query_dict:
        abort(Response("Unrecognized query options; %s." %
                       list(query_dict.keys()),
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
    stmts = []
    logger.info("Getting statements...")

    # First look for statements matching the role'd agents.
    for role, (ag_id, ns) in roled_agents.items():
        stmts = _get_relevant_statements(stmts, ag_id, ns, act, role)
        if not len(stmts):
            break
    else:
        # Second (if we're still looking), look for statements from free agents.
        for ag_id, ns in free_agents:
            stmts = _get_relevant_statements(stmts, ag_id, ns, act)
            if not len(stmts):
                break

    # Check to make sure we didn't get too many statements
    hit_limit = (len(stmts) > MAX_STATEMENTS)
    if hit_limit:
        if limit_behavior == 'error':
            # Divide the statements up by type, and get counts of each type.
            stmt_counts = {}
            for stmt in stmts:
                stmt_type = type(stmt).__name__
                if stmt_type not in stmt_counts:
                    stmt_counts[stmt_type] = {'count': 0}
                stmt_counts[stmt_type]['count'] += 1

            return jsonify({'limited': True, 'statements': stmt_counts}), 413
        elif limit_behavior == 'truncate':
            stmts = stmts[:MAX_STATEMENTS]
        elif limit_behavior == 'sample':
            from random import sample
            stmts = sample(stmts, MAX_STATEMENTS)

    # Retrieve the evidence for the statements.
    logger.info("Getting evidence for the %d resulting statements."
                % len(stmts))
    if stmts:
        get_evidence(stmts, fix_refs=False)

    # Create the json response, and send off.
    resp = jsonify({'limited': hit_limit,
                    'statements':[stmt.to_json() for stmt in stmts]})
    logger.info("Exiting with %d statements of nominal size %f MB."
                % (len(stmts), sys.getsizeof(resp.data)/1e6))
    return resp


@app.route('/papers/', methods=['GET'])
def get_paper_statements():
    """Get and preassemble statements from a paper given by pmid."""
    logger.info("Got query for statements from a paper!")
    query_dict = request.args.copy()

    # Get the paper id.
    id_val = query_dict.get('id')
    if id_val is None:
        logger.error("No id provided!")
        abort(Response("No id in request!", 400))

    # Get the id type, if given.
    id_type = query_dict.get('type', 'pmid')

    # Now get the statements.
    logger.info('Getting statements for %s=%s...' % (id_type, id_val))
    stmts = get_statements_by_paper(id_val, id_type, do_stmt_count=False)
    if stmts is None:
        msg = "Invalid or unavailable id %s=%s!" % (id_type, id_val)
        logger.error(msg)
        abort(Response(msg, 404))
    logger.info("Found %d statements." % len(stmts))

    resp = jsonify({'statements': [stmt.to_json() for stmt in stmts]})
    logger.info("Exiting with %d statements." % len(stmts))
    return resp


if __name__ == '__main__':
    app.run()
