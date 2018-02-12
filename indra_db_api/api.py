from flask import Flask, request, abort, jsonify
from indra.db.util import get_statements_by_gene_role_type

app = Flask(__name__)

@app.route('/statements', methods=['GET'])
def get_statments():
    """Get some statements constrained by query."""
    json_req = request.get_json()
    if not json_req:
        abort(400)
    obj = json_req.get('object')
    act = json_req.get('action')
    sub = json_req.get('subject')

    if sub is None and obj is None:
        abort(400)

    stmts = []
    if sub is not None:
        stmts = get_statements_by_gene_role_type(agent_id=sub,
                                                 role='SUBJECT',
                                                 stmt_type=act)
        if obj is not None:
            stmts = [s for s in stmts if s.agent_list()[1].name == obj]
    elif obj is not None:
        stmts = get_statements_by_gene_role_type(agent_id=obj,
                                                 role='OBJECT',
                                                 stmt_type=act)

    return jsonify([stmt.to_json() for stmt in stmts])