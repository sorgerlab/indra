import argparse
import inspect
from collections import namedtuple
from flask import Flask, request
from flask_restx import Api, Resource, Namespace, inputs, fields

from indra import get_config
from indra.sources import trips, reach, bel, biopax, eidos, hume, cwms, sofia
from indra.databases import hgnc_client
from indra.statements import stmts_from_json, get_statement_by_name
from indra.assemblers.pysb import PysbAssembler
import indra.assemblers.pysb.assembler as pysb_assembler
from indra.assemblers.cx import CxAssembler
from indra.assemblers.graph import GraphAssembler
from indra.assemblers.cyjs import CyJSAssembler
from indra.assemblers.sif import SifAssembler
from indra.assemblers.english import EnglishAssembler
from indra.tools.assemble_corpus import *
from indra.databases import cbio_client
from indra.sources.indra_db_rest import get_statements
from indra.sources.ndex_cx.api import process_ndex_network
from indra.sources.reach.api import reach_nxml_url, reach_text_url
from indra.belief.wm_scorer import get_eidos_scorer
from indra.preassembler.ontology_mapper import OntologyMapper, wm_ontomap
from indra.pipeline import AssemblyPipeline, pipeline_functions
from indra.preassembler.custom_preassembly import *


# Create Flask app, api, namespaces, and models
app = Flask(__name__)
api = Api(app, title='INDRA REST API', description='REST API for INDRA webservice')

preassembly = Namespace('Preassembly', path='/preassembly/')
sofia = Namespace('Sofia', path='/sofia/')
eidos = Namespace('Eidos', path='/eidos/')
hume = Namespace('Hume', path='/hume/')
bel = Namespace('BEL', path='/bel/')
trips = Namespace('TRIPS', path='/trips/')
reach = Namespace('REACH', path='/reach/')
cwms = Namespace('CWMS', path='/cwms/')
biopax = Namespace('BioPAX', path='/biopax/')
assemblers = Namespace('Assemblers', path='/assemblers/')
api.add_namespace(preassembly)
api.add_namespace(sofia)
api.add_namespace(eidos)
api.add_namespace(hume)
api.add_namespace(bel)
api.add_namespace(trips)
api.add_namespace(reach)
api.add_namespace(cwms)
api.add_namespace(biopax)
api.add_namespace(assemblers)

stmt_model = api.model('Statement', {})
stmts_model = api.model('Statements', {
    'statements': fields.List(fields.Nested(stmt_model))
})
cur_model = api.model('Curation', {})
cur_dict = api.model('dict', {'cur': fields.Nested(cur_model)})

# Store the arguments by type
int_args = ['poolsize', 'size_cutoff']
float_args = ['score_threshold', 'belief_cutoff']
boolean_args = [
    'do_rename', 'use_adeft', 'do_methionine_offset', 'do_orthology_mapping',
    'do_isoform_mapping', 'use_cache', 'return_toplevel', 'flatten_evidence',
    'normalize_equivalences', 'normalize_opposites', 'invert', 'remove_bound',
    'specific_only', 'allow_families', 'match_suffix', 'update_belief']
list_args = [
    'gene_list', 'name_list', 'values', 'source_apis', 'uuids', 'curations',
    'correct_tags', 'ignores', 'deletions']
dict_args = [
    'grounding_map', 'misgrounding_map', 'whitelist', 'mutations']


def _return_stmts(stmts):
    if stmts:
        stmts_json = stmts_to_json(stmts)
        res = {'statements': stmts_json}
    else:
        res = {'statements': []}
    return res


def _stmts_from_proc(proc):
    if proc and proc.statements:
        stmts = stmts_to_json(proc.statements)
        res = {'statements': stmts}
    else:
        res = {'statements': []}
    return res


# Create Resources in Preassembly Namespace
class PreassembleStatements(Resource):
    """Parent Resource for Preassembly resources."""
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
        stmts_out = pipeline_functions[self.func_name](stmts, **args)
        return _return_stmts(stmts_out)


def make_preassembly_model(func):
    """Create new Flask model with function arguments."""
    args = inspect.signature(func).parameters
    # We can reuse Staetments model if only stmts_in or stmts and **kwargs are
    # arguments of the function
    if ((len(args) == 1 and ('stmts_in' in args or 'stmts' in args)) or
            (len(args) == 2 and 'kwargs' in args and
                ('stmts_in' in args or 'stmts' in args))):
        return stmts_model
    # Inherit a model if there are other arguments
    model_fields = {}
    for arg in args:
        if arg != 'stmts_in' and arg != 'stmts' and arg != 'kwargs':
            example = None
            if args[arg].default is not inspect.Parameter.empty:
                example = args[arg].default
            if arg in boolean_args:
                model_fields[arg] = fields.Boolean(example=example)
            elif arg in int_args:
                model_fields[arg] = fields.Integer(example=example)
            elif arg in float_args:
                model_fields[arg] = fields.Float(example=example)
            elif arg in list_args:
                if arg == 'curations':
                    model_fields[arg] = fields.List(fields.Nested(cur_model))
                else:
                    model_fields[arg] = fields.List(
                        fields.String, example=example)
            elif arg in dict_args:
                item_model = api.model(arg, {})
                model_fields[arg] = fields.Nested(item_model, example=example)
            else:
                model_fields[arg] = fields.String(example=example)
    new_model = api.inherit(
        ('%s_input' % func.__name__), stmts_model, model_fields)
    return new_model


# Create resources for each of assembly_corpus functions
for func_name, func in pipeline_functions.items():
    if func.__module__ == 'indra.tools.assemble_corpus':
        doc = ''
        # Get the function description from docstring
        if func.__doc__:
            doc = func.__doc__.split('\n')[0]

        new_model = make_preassembly_model(func)

        @preassembly.expect(new_model)
        @preassembly.route(('/%s' % func_name), doc={'description': doc})
        class NewFunction(PreassembleStatements):
            func_name = func_name


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('Run the INDRA REST API')
    argparser.add_argument('--host', default='0.0.0.0')
    argparser.add_argument('--port', default=8081, type=int)
    argparserargs = argparser.parse_args()
    app.run(host=argparserargs.host, port=argparserargs.port)
