"""This service implements all base functions of the ontology graph as a REST
service. The three key functions that most ontology methods rely on are
child_rel, parent_rel, and get_node_property. There are a few other bookkeeping
functions that also need to be implemented here since they access ontology
attributes directly."""
import argparse
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


@app.route('/get_id_from_name', methods=['GET'])
def get_id_from_name():
    ont = request.json.get('ontology')
    ontology = ontologies.get(ont)
    kwargs = ('ns', 'name')
    return jsonify(ontology.get_id_from_name(
        **{k: v for k, v in request.json.items() if k in kwargs}))


@app.route('/get_xrefs', methods=['GET'])
def get_xrefs():
    ont = request.json.get('ontology')
    ontology = ontologies.get(ont)
    kwargs = ('ns', 'id')
    return jsonify(ontology.get_mappings(
        **{k: v for k, v in request.json.items() if k in kwargs}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Run the INDRA Ontology service.')
    parser.add_argument('--port', help='The port to run the server on.',
                        default=8082)
    args = parser.parse_args()
    app.run(host='0.0.0.0', port=args.port)
