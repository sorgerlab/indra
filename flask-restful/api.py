from flask import Flask, request
from flask_restful import Resource, Api, reqparse
from indra import reach
from indra.statements import *
import json

app = Flask(__name__)
api = Api(app)
parser = reqparse.RequestParser()
parser.add_argument('txt')

class InputText(Resource):
    def post(self):
        args = parser.parse_args()
        txt = args['txt']
	rp = reach.process_text(txt, offline=False)
	st = rp.statements
	json_statements = {}
        json_statements['statements'] = []
	for s in st:
            s_json = s.to_json()
            json_statements['statements'].append(s_json)
        return json_statements, 201

api.add_resource(InputText, '/')

if __name__ == '__main__':
    app.run(debug=True)
