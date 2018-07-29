import json
from flask import Flask, request
from indra.sources.eidos.eidos_reader import EidosReader


app = Flask(__name__)


@app.route('/process_text', methods=['POST'])
def process_text():
    text = request.json.get('text')
    if not text:
        return {}
    res = er.process_text(text, 'json_ld')
    return json.dumps(res)


if __name__ == '__main__':
    er = EidosReader()
    er.process_text('hello', 'json_ld')
    app.run(host='0.0.0.0')
