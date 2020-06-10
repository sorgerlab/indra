from flask import Flask, request, jsonify
from indra.ontology.bio import bio_ontology

app = Flask(__name__)


bio_ontology.initialize()
ontologies = {'bio': bio_ontology}


@app.route('/child_rel', methods=['GET'])
def child_rel():
    ont = request.json.get('ontology')
    ontology = ontologies.get(ont)
    kwargs = ('ns', 'id', 'rel_types')
    return jsonify(list(ontology.child_rel(
        **{k: v for k, v in request.json.items() if k in kwargs})))


@app.route('/parent_rel', methods=['GET'])
def parent_rel():
    ont = request.json.get('ontology')
    ontology = ontologies.get(ont)
    kwargs = ('ns', 'id', 'rel_types')
    return jsonify(list(ontology.parent_rel(
        **{k: v for k, v in request.json.items() if k in kwargs})))


@app.route('/get_node_property', methods=['GET'])
def get_node_property():
    ont = request.json.get('ontology')
    ontology = ontologies.get(ont)
    kwargs = ('ns', 'id', 'property')
    return jsonify(ontology.get_node_property(
        **{k: v for k, v in request.json.items() if k in kwargs}))


if __name__ == '__main__':
    app.run(host='0.0.0.0', port=6667)