from flask import Flask, request, abort, jsonify, Response
from indra.db.util import get_statements_by_gene_role_type
import logging

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

@app.route('/statements', methods=['POST'])
def get_statments():
    """Get some statements constrained by query."""
    logger.info("Got request for statements!")
    headers = request.headers
    if 'Content-Type' not in headers.keys() or headers['Content-Type'] != 'application/json':
        abort(Response("Content-Type not set to application/json.\n", 400))
    json_req = None
    try:
        json_req = request.get_json()
        logger.info("Got json data: %s" % json_req)
    except:
        logger.error("Could not load json request data.")
        abort(Response("Failed to get json data.\n"), 400)
    if not json_req:
        abort(Response("No data provided!\n", 400))
    logger.info("Getting query details.")
    obj = json_req.get('object')
    act = json_req.get('action')
    sub = json_req.get('subject')

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