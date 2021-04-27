"""This is a Python-based web server that can be run to
read with Eidos. To run the server, do

    python -m indra.sources.eidos.server

and then submit POST requests to the `localhost:5000/process_text` endpoint
with JSON content as `{'text': 'text to read'}`. The response will be the
Eidos JSON-LD output. Another endpoint for regrounding entity texts
is also available on the `reground` endpoint.
"""
import sys
import json
from flask import Flask, request
from indra.sources.eidos.reader import EidosReader

app = Flask(__name__)


@app.route('/process_text', methods=['POST'])
def process_text():
    text = request.json.get('text')
    if not text:
        return {}
    res = er.process_text(text)
    return json.dumps(res)


@app.route('/reground', methods=['POST'])
def reground():
    text = request.json.get('text')
    ont_yml = request.json.get('ont_yml')
    if not ont_yml:
        from indra_world.ontology import world_ontology
        ont_yml = world_ontology.dump_yml_str()
    topk = request.json.get('topk', 10)
    is_canonicalized = request.json.get('is_canonicalized', False)
    if not text:
        return []
    if isinstance(text, str):
        text = [text]
    res = er.reground_texts(text, ont_yml, topk=topk,
                            is_canonicalized=is_canonicalized)
    return json.dumps(res)


if __name__ == '__main__':
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 6666
    er = EidosReader()
    er.process_text('hello')  # This is done to initialize the system
    app.run(host='0.0.0.0', port=port)
