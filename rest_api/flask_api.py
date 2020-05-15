import argparse
import json
from functools import partial
from flask import Flask, request
from flask_restful import Resource, Api, reqparse
from indra.pipeline import AssemblyPipeline, pipeline_functions
from indra.tools.assemble_corpus import *
from indra.statements import stmts_from_json


app = Flask(__name__)
api = Api(app)
parser = reqparse.RequestParser()
parser.add_argument('statements')


def _return_stmts(stmts):
    if stmts:
        stmts_json = stmts_to_json(stmts)
        res = {'statements': stmts_json}
    else:
        res = {'statements': []}
    return res


class PreassembleStatements(Resource):
    def post(self, func_name):
        args = request.form.to_dict()
        # args = parser.parse_args()
        stmts = stmts_from_json(json.loads(args['statements']))
        stmts_out = pipeline_functions[func_name](stmts)
        return _return_stmts(stmts_out)


api.add_resource(PreassembleStatements, '/preassembly/<func_name>')


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('Run the INDRA REST API')
    argparser.add_argument('--host', default='0.0.0.0')
    argparser.add_argument('--port', default=8081, type=int)
    args = argparser.parse_args()
    app.run(host=args.host, port=args.port)
