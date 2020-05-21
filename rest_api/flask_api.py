import argparse
import json
from functools import partial
from flask import Flask, request
from flask_restful import Resource, Api, reqparse
from indra.pipeline import AssemblyPipeline, pipeline_functions
from indra.tools.assemble_corpus import *
from indra.statements import stmts_from_json


boolean_args = [
    'do_rename', 'use_adeft', 'do_methionine_offset', 'do_orthology_mapping',
    'do_isoform_mapping', 'use_cache', 'return_toplevel', 'flatten_evidence',
    'normalize_equivalences', 'normalize_opposites', 'invert', 'remove_bound',
    'specific_only', 'allow_families', 'match_suffix', 'update_belief']
list_args = [
    'gene_list', 'name_list', 'values', 'source_apis', 'uuids', 'curations',
    'correct_tags', 'ignores']
dict_args = [
    'grounding_map', 'misgrounding_map', 'hierarchies', 'whitelist',
    'mutations']
float_args = ['score_threshold', 'belief_cutoff']
int_args = ['poolsize', 'size_cutoff']
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


def add_arg_to_parser(arg, parser):
    if arg in list_args:
        dtype = list
    elif arg in boolean_args:
        dtype = inputs.boolean
    elif arg in int_args:
        dtype = int
    elif arg in float_args:
        dtype = float
    elif arg in dict_args:
        dtype = dict
    else:
        dtype = str
    parser.add_argument(arg, location='json')


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
