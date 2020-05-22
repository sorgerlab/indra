import argparse
import inspect
from collections import namedtuple
from flask import Flask, request
from flask_restx import Api, Resource, reqparse, Namespace, inputs
from indra.pipeline import AssemblyPipeline, pipeline_functions
from indra.tools.assemble_corpus import *
from indra.statements import stmts_from_json, get_statement_by_name
from indra.belief.wm_scorer import get_eidos_scorer
from indra.preassembler.custom_preassembly import *


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
api = Api(app, title='INDRA REST API')
preassembly = Namespace('preassembly', path='/preassembly/')
api.add_namespace(preassembly)
parsers = {}


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
    parser.add_argument(arg, type=dtype, location='json')


class PreassembleStatements(Resource):
    func_name = None

    def process_args(self, args_json):
        for arg in args_json:
            if arg == 'stmt_type':
                args_json[arg] = get_statement_by_name(args_json[arg])
            elif arg in ['matches_fun', 'refinement_fun']:
                args_json[arg] = pipeline_functions[args_json[arg]]
            elif arg == 'curations':
                Curation = namedtuple(
                    'Curation', ['pa_hash', 'source_hash', 'tag'])
                args_json[arg] = [
                    Curation(cur['pa_hash'], cur['source_hash'], cur['tag'])
                    for cur in args_json[arg]]
            elif arg == 'belief_scorer':
                if args_json[arg] == 'wm' or \
                        args_json[arg] == 'get_eidos_scorer':
                    args_json[arg] = get_eidos_scorer()
            elif arg == 'whitelist' or arg == 'mutations':
                args_json[arg] = {
                    gene: [tuple(mod) for mod in mods]
                    for gene, mods in args_json[arg].items()}
        return args_json

    def post(self):
        args = self.process_args(request.json)
        stmts = stmts_from_json(args.pop('statements'))
        print(args)
        stmts_out = pipeline_functions[self.func_name](stmts, **args)
        return _return_stmts(stmts_out)


for func_name, func in pipeline_functions.items():
    if func.__module__ == 'indra.tools.assemble_corpus':
        doc = ''
        if func.__doc__:
            doc = func.__doc__.split('\n')[0]
        parser = reqparse.RequestParser()
        parser.add_argument('statements', type=list, location='json')
        for arg in inspect.getargspec(func).args:
            if arg != 'stmts_in':
                add_arg_to_parser(arg, parser)
        parsers[func_name] = parser

        @preassembly.expect(parsers[func_name])
        @preassembly.route(f'/{func_name}', doc={'description': doc})
        class NewFunction(PreassembleStatements):
            func_name = func_name


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('Run the INDRA REST API')
    argparser.add_argument('--host', default='0.0.0.0')
    argparser.add_argument('--port', default=8081, type=int)
    argparserargs = argparser.parse_args()
    app.run(host=argparserargs.host, port=argparserargs.port)
