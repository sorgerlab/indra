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

preassembly_ns = Namespace('Preassembly', path='/preassembly/')
sofia_ns = Namespace('Sofia', path='/sofia/')
eidos_ns = Namespace('Eidos', path='/eidos/')
hume_ns = Namespace('Hume', path='/hume/')
bel_ns = Namespace('BEL', path='/bel/')
trips_ns = Namespace('TRIPS', path='/trips/')
reach_ns = Namespace('REACH', path='/reach/')
cwms_ns = Namespace('CWMS', path='/cwms/')
biopax_ns = Namespace('BioPAX', path='/biopax/')
assemblers_ns = Namespace('Assemblers', path='/assemblers/')
api.add_namespace(preassembly_ns)
api.add_namespace(sofia_ns)
api.add_namespace(eidos_ns)
api.add_namespace(hume_ns)
api.add_namespace(bel_ns)
api.add_namespace(trips_ns)
api.add_namespace(reach_ns)
api.add_namespace(cwms_ns)
api.add_namespace(biopax_ns)
api.add_namespace(assemblers_ns)

# Models that can be inherited and reused in different namespaces
stmt_model = api.model('Statement', {})
stmts_model = api.model('Statements', {
    'statements': fields.List(fields.Nested(stmt_model))
})
dict_model = api.model('dict', {})
bio_text_model = api.model('BioText', {
    'text': fields.String(example='GRB2 binds SHC.')})
wm_text_model = api.model('WMText', {
    'text': fields.String(example='Rainfall causes floods.')})

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
            default = None
            if args[arg].default is not inspect.Parameter.empty:
                default = args[arg].default
            if arg in boolean_args:
                model_fields[arg] = fields.Boolean(default=default)
            elif arg in int_args:
                model_fields[arg] = fields.Integer(default=default)
            elif arg in float_args:
                model_fields[arg] = fields.Float(default=default)
            elif arg in list_args:
                if arg == 'curations':
                    model_fields[arg] = fields.List(fields.Nested(dict_model))
                else:
                    model_fields[arg] = fields.List(
                        fields.String, default=default)
            elif arg in dict_args:
                model_fields[arg] = fields.Nested(dict_model, default=default)
            else:
                model_fields[arg] = fields.String(default=default)
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

        @preassembly_ns.expect(new_model)
        @preassembly_ns.route(('/%s' % func_name), doc={'description': doc})
        class NewFunction(PreassembleStatements):
            func_name = func_name

# Create models and resources for REACH namespace
reach_text_model = api.inherit('ReachText', bio_text_model, {
    'offline': fields.Boolean(default=False),
    'url': fields.String
})
reach_json_model = api.model('ReachJSON', {'json': fields.Nested(dict_model)})
reach_pmc_model = api.model('ReachPMC', {'pmcid': fields.String})


@reach_ns.expect(reach_text_model)
@reach_ns.route('/process_text')
class ReachProcessText(Resource):
    def post(self):
        args = request.json
        text = args.get('text')
        offline = True if args.get('offline') else False
        given_url = args.get('url')
        config_url = get_config('REACH_TEXT_URL', failure_ok=True)
        # Order: URL given as an explicit argument in the request. Then any URL
        # set in the configuration. Then, unless offline is set, use the
        # default REACH web service URL.
        if 'url' in args:  # This is to take None if explicitly given
            url = given_url
        elif config_url:
            url = config_url
        elif not offline:
            url = reach_text_url
        else:
            url = None
        # If a URL is set, prioritize it over the offline setting
        if url:
            offline = False
        rp = reach.process_text(text, offline=offline, url=url)
        return _stmts_from_proc(rp)


@reach_ns.expect(reach_json_model)
@reach_ns.route('/process_json')
class ReachProcessJson(Resource):
    def post(self):
        args = request.json
        json_str = args.get('json')
        rp = reach.process_json_str(json_str)
        return _stmts_from_proc(rp)


@reach_ns.expect(reach_pmc_model)
@reach_ns.route('/process_pmc')
class ReachProcessPmc(Resource):
    def post(self):
        args = request.json
        pmcid = args.get('pmcid')
        offline = True if args.get('offline') else False
        given_url = args.get('url')
        config_url = get_config('REACH_NXML_URL', failure_ok=True)
        # Order: URL given as an explicit argument in the request. Then any URL
        # set in the configuration. Then, unless offline is set, use the
        # default REACH web service URL.
        if 'url' in args:  # This is to take None if explicitly given
            url = given_url
        elif config_url:
            url = config_url
        elif not offline:
            url = reach_nxml_url
        else:
            url = None
        # If a URL is set, prioritize it over the offline setting
        if url:
            offline = False
        rp = reach.process_pmc(pmcid, offline=offline, url=url)
        return _stmts_from_proc(rp)


# Create models and resources for TRIPS namespace
xml_model = api.model('XML', {'xml_str': fields.String})


@trips_ns.expect(bio_text_model)
@trips_ns.route('/process_text')
class TripsProcessText(Resource):
    def post(self):
        args = request.json
        text = args.get('text')
        tp = trips.process_text(text)
        return _stmts_from_proc(tp)


@trips_ns.expect(xml_model)
@trips_ns.route('/process_xml')
class TripsProcessText(Resource):
    def post(self):
        args = request.json
        xml_str = args.get('xml_str')
        tp = trips.process_xml(xml_str)
        return _stmts_from_proc(tp)


# Create models and resources for Sofia namespace
text_auth_model = api.inherit('TextAuth', wm_text_model, {
    'auth': fields.List(fields.String, example=['USER', 'PASS'])})


@sofia_ns.expect(text_auth_model)
@sofia_ns.route('/process_text')
class SofiaProcessText(Resource):
    def post(self):
        args = request.json
        text = args.get('text')
        auth = args.get('auth')
        sp = sofia.process_text(text, auth=auth)
        return _stmts_from_proc(sp)


# Create models and resources for Eidos namespace
eidos_text_model = api.inherit('EidosText', wm_text_model, {
    'webservice': fields.String,
    'grounding_ns': fields.String
})
jsonld_model = api.model('jsonld', {
    'jsonld': fields.Nested(dict_model)})


@eidos_ns.expect(eidos_text_model)
@eidos_ns.route('/process_text')
class EidosProcessText(Resource):
    def post(self):
        args = request.json
        text = args.get('text')
        webservice = args.get('webservice')
        grounding_ns = args.get('grounding_ns')
        # if not webservice:
        #     response.status = 400
        #     response.content_type = 'application/json'
        #     return json.dumps({'error': 'No web service address provided.'})
        ep = eidos.process_text(text, webservice=webservice,
                                grounding_ns=grounding_ns)
        return _stmts_from_proc(ep)


# Create models and resources for Hume namespace

if __name__ == '__main__':
    argparser = argparse.ArgumentParser('Run the INDRA REST API')
    argparser.add_argument('--host', default='0.0.0.0')
    argparser.add_argument('--port', default=8081, type=int)
    argparserargs = argparser.parse_args()
    app.run(host=argparserargs.host, port=argparserargs.port)
