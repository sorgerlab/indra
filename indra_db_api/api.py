from flask import Flask, request, abort, jsonify
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
    json_req = None
    try:
        json_req = request.get_json()
        logger.info("Got json data: %s" % json_req)
    except:
        logger.error("Could not load json request data.")
        return "Failed to get json data."
    if not json_req:
        return "No data provided!\n"
    logger.info("Getting query details.")
    obj = json_req.get('object')
    act = json_req.get('action')
    sub = json_req.get('subject')

    if sub is None and obj is None:
        logger.error("No subject or object.")
        return "No subject or object!\n"

    stmts = []
    logger.info("Getting statements...")
    if sub is not None:
        logger.debug("Getting statements by subject...")
        stmts = get_statements_by_gene_role_type(agent_id=sub,
                                                 role='SUBJECT',
                                                 stmt_type=act,
                                                 do_stmt_count=False)
        if obj is not None:
            logger.debug("Filtering by object...")
            stmts = [s for s in stmts if len(s.agent_list()) > 1
                     and s.agent_list()[1].name == obj]
    elif obj is not None:
        logger.debug("Getting statements by object...")
        stmts = get_statements_by_gene_role_type(agent_id=obj,
                                                 role='OBJECT',
                                                 stmt_type=act,
                                                 do_stmt_count=False)

    logger.debug("Exiting with %d statements." % len(stmts))
    return jsonify([stmt.to_json() for stmt in stmts])

if __name__ == '__main__':
    app.run()