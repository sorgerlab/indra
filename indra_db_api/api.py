from flask import Flask, request, abort, jsonify, Response
from indra.db.util import get_statements_by_gene_role_type
import logging
import re

logger = logging.getLogger("db-api")

app = Flask(__name__)

@app.route('/')
def welcome():
    logger.info("Got request for welcome info.")
    return ("Welcome the the INDRA database webservice!\n"
            "\n"
            "Use modes:\n"
            "----------\n"
            "/            - (you are here) Welcome page.\n"
            "/statements  - Get statement content from INDRA's database.\n"
            "\n")

@app.route('/statements', methods=['GET'])
def get_statements_query_format():
    return ('To get a list of statements, use ?[object,subject,action]=<val>, '
            'for example ?object=MAP2K1?subject=MAPK1?action=Phosphorylation')

@app.route('/statements/<query_str>', methods=['GET'])
def get_statments(query_str):
    """Get some statements constrained by query."""
    logger.info("Got query for statements!")
    arg_patt = re.compile('\?(\w+)=(\w+)')
    query_dict = dict(arg_patt.findall(query_str))

    logger.info("Getting query details.")
    obj = query_dict.pop('object', None)
    act = query_dict.pop('action', None)
    sub = query_dict.pop('subject', None)
    if query_dict:
        abort(Response("Unrecognized query options; %s." %
                       list(query_dict.keys()),
                       400))

    if sub is None and obj is None:
        logger.error("No subject or object.")
        abort(Response("No subject or object!\n", 400))

    stmts = []
    logger.info("Getting statements...")
    if sub is not None:
        logger.info("Getting statements by subject...")
        stmts = get_statements_by_gene_role_type(agent_id=sub,
                                                 role='SUBJECT',
                                                 stmt_type=act,
                                                 do_stmt_count=False)
        if obj is not None:
            logger.info("Filtering by object...")
            stmts = [s for s in stmts if len(s.agent_list()) > 1
                     and s.agent_list()[1].name == obj]
    elif obj is not None:
        logger.info("Getting statements by object...")
        stmts = get_statements_by_gene_role_type(agent_id=obj,
                                                 role='OBJECT',
                                                 stmt_type=act,
                                                 do_stmt_count=False)

    logger.info("Exiting with %d statements." % len(stmts))
    return jsonify([stmt.to_json() for stmt in stmts])

if __name__ == '__main__':
    app.run()