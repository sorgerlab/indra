import argparse
import inspect
import logging
import json
import base64

from docstring_parser import parse
from flask import Flask, request
from flask_restx import Api, Resource, fields, abort
from flask_cors import CORS

from indra import get_config
from indra.sources import trips, reach, bel, biopax, eidos
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
from indra.ontology.bio import bio_ontology
from indra.pipeline import AssemblyPipeline, pipeline_functions
from indra.preassembler.custom_preassembly import *


logger = logging.getLogger('rest_api')
logger.setLevel(logging.DEBUG)


# Create Flask app, api, namespaces, and models
app = Flask(__name__)
api = Api(
    app, title='INDRA REST API', description='REST API for INDRA webservice')
CORS(app)

preassembly_ns = api.namespace(
    'Preassembly', 'Preassemble INDRA Statements', path='/preassembly/')
sources_ns = api.namespace(
    'Sources', 'Get INDRA Statements from various sources', path='/')
assemblers_ns = api.namespace(
    'Assemblers', 'Assemble INDRA Statements into models', path='/assemblers/')
ndex_ns = api.namespace('NDEx', 'Use NDEx service', path='/')
indra_db_rest_ns = api.namespace(
    'INDRA DB REST', 'Use INDRA DB REST API', path='/indra_db_rest/')
databases_ns = api.namespace(
    'Databases', 'Access external databases', path='/databases/')

# Models that can be inherited and reused in different namespaces
dict_model = api.model('dict', {})

stmts_model = api.model('Statements', {
    'statements': fields.List(fields.Nested(dict_model), example=[{
        "id": "acc6d47c-f622-41a4-8ae9-d7b0f3d24a2f",
        "type": "Complex",
        "members": [
            {"db_refs": {"TEXT": "MEK", "FPLX": "MEK"}, "name": "MEK"},
            {"db_refs": {"TEXT": "ERK", "FPLX": "ERK"}, "name": "ERK"}
        ],
        "sbo": "https://identifiers.org/SBO:0000526",
        "evidence": [{"text": "MEK binds ERK", "source_api": "trips"}]
        }])})
bio_text_model = api.model('BioText', {
    'text': fields.String(example='GRB2 binds SHC.')})
jsonld_model = api.model('jsonld', {
    'jsonld': fields.String(example='{}')})
genes_model = api.model('Genes', {
    'genes': fields.List(fields.String, example=['BRAF', 'MAP2K1'])})

# Store the arguments by type
int_args = ['members_allowed', 'protocol', 'poolsize', 'size_cutoff']
float_args = ['score_threshold', 'belief_cutoff']
boolean_args = [
    'do_rename', 'use_adeft', 'do_methionine_offset', 'do_orthology_mapping',
    'do_isoform_mapping', 'use_cache', 'return_toplevel', 'flatten_evidence',
    'normalize_equivalences', 'normalize_opposites', 'invert', 'remove_bound',
    'specific_only', 'allow_families', 'match_suffix', 'update_belief',
    'in_place', 'print_report_before', 'print_report_after',
    'prior_hash_annots']
list_args = [
    'gene_list', 'name_list', 'values', 'source_apis', 'uuids', 'curations',
    'correct_tags', 'ignores', 'deletions']
dict_args = [
    'grounding_map', 'misgrounding_map', 'whitelist', 'mutations']

str_args = {'stmt_type': "Modification", 'policy': 'all'}


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

# Manually add preassembly resources not based on assembly corpus functions
pipeline_model = api.inherit('Pipeline', stmts_model, {
    'pipeline': fields.List(fields.Nested(dict_model), example=[
        {'function': 'filter_grounded_only'},
        {'function': 'run_preassembly', 'kwargs': {'return_toplevel': False}}
    ])
})

# There's an extra blank line between parameters here and in all the following
# docstrings for better visualization in Swagger
@preassembly_ns.expect(pipeline_model)
@preassembly_ns.route('/pipeline')
class RunPipeline(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Run an assembly pipeline for a list of Statements.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to run the pipeline.

        pipeline : list[dict]
            A list of dictionaries representing steps in the pipeline. Each
            step should have a 'function' key and, if appropriate, 'args' and
            'kwargs' keys. For more documentation and examples, see
            https://indra.readthedocs.io/en/latest/modules/pipeline.html

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            The list of INDRA Statements resulting from running the pipeline
            on the list of input Statements.
        """
        args = request.json
        stmts = stmts_from_json(args.get('statements'))
        pipeline_steps = args.get('pipeline')
        ap = AssemblyPipeline(pipeline_steps)
        stmts_out = ap.run(stmts)
        return _return_stmts(stmts_out)


# Dynamically generate resources for assembly corpus functions
class PreassembleStatements(Resource):
    """Parent Resource for Preassembly resources."""
    func_name = None

    def process_args(self, args_json):
        for arg in args_json:
            if arg == 'stmt_type':
                args_json[arg] = get_statement_by_name(args_json[arg])
            elif arg in ['matches_fun', 'refinement_fun']:
                args_json[arg] = pipeline_functions[args_json[arg]]
            elif arg == 'belief_scorer':
                # Here we could handle various string values of args_json[arg]
                # but there currently aren't any specific options
                args_json[arg] = None
            elif arg == 'ontology':
                # Here we could handle various string values of args_json[arg]
                # but there currently aren't any specific options
                args_json[arg] = bio_ontology
            elif arg == 'whitelist' or arg == 'mutations':
                args_json[arg] = {
                    gene: [tuple(mod) for mod in mods]
                    for gene, mods in args_json[arg].items()}
        return args_json

    @api.doc(False)
    def options(self):
        return {}

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
            # Need to use default for boolean and example for other types
            if arg in boolean_args:
                model_fields[arg] = fields.Boolean(default=default)
            elif arg in int_args:
                model_fields[arg] = fields.Integer(example=default)
            elif arg in float_args:
                model_fields[arg] = fields.Float(example=0.7)
            elif arg in list_args:
                if arg == 'curations':
                    model_fields[arg] = fields.List(
                        fields.Nested(dict_model),
                        example=[{'pa_hash': '1234', 'source_hash': '2345',
                                  'tag': 'wrong_relation'}])
                else:
                    model_fields[arg] = fields.List(
                        fields.String, example=default)
            elif arg in dict_args:
                model_fields[arg] = fields.Nested(dict_model)
            elif arg in str_args.keys():
                model_fields[arg] = fields.String(example=str_args[arg])
            else:
                model_fields[arg] = fields.String(example=default)
    new_model = api.inherit(
        ('%s_input' % func.__name__), stmts_model, model_fields)
    return new_model


def update_docstring(func):
    doc = func.__doc__
    docstring = parse(doc)
    new_doc = docstring.short_description + '\n\n'
    if docstring.long_description:
        new_doc += (docstring.long_description + '\n\n')
    new_doc += ('Parameters\n----------\n')
    for param in docstring.params:
        if param.arg_name in ['save', 'save_unique']:
            continue
        elif param.arg_name in ['stmts', 'stmts_in']:
            param.arg_name = 'statements'
            param.type_name = 'list[indra.statements.Statement.to_json()]'
        elif param.arg_name == 'belief_scorer':
            param.type_name = 'Optional[str] or None'
            param.description = (
                'Type of BeliefScorer to use in calculating Statement '
                'probabilities. If None is provided (default), then the '
                'default scorer is used (good for biology use case). '
                'For WorldModelers use case belief scorer should be set '
                'to "wm".')
        elif param.arg_name == 'ontology':
            param.type_name = 'Optional[str] or None'
            param.description = (
                'Type of ontology to use for preassembly ("bio" or "wm"). '
                'If None is provided (default), then the bio ontology is used.'
                'For WorldModelers use case ontology should be set to "wm".')
        elif param.arg_name in ['matches_fun', 'refinement_fun']:
            param.type_name = 'str'
        elif param.arg_name == 'curations':
            param.type_name = 'list[dict]'
            param.description = (
                'A list of dictionaries representing curations. Each '
                'dictionary must have "pa_hash" (preassembled statement hash)'
                ', "source_hash", (evidence hash) and "tag" (e.g. "correct", '
                '"wrong_relation", etc.) keys.')
        new_doc += (param.arg_name + ' : ' + param.type_name + '\n' +
                    param.description + '\n\n')
    new_doc += 'Returns\n----------\n'
    new_doc += 'statements : list[indra.statements.Statement.to_json()]\n'
    new_doc += 'A list of processed INDRA Statements'
    return docstring.short_description, new_doc


# Create resources for each of assembly_corpus functions
for func_name, func in pipeline_functions.items():
    if func.__module__ == 'indra.tools.assemble_corpus':
        doc = ''
        short_doc = ''
        # Get the function description from docstring
        if func.__doc__:
            short_doc, doc = update_docstring(func)
        new_model = make_preassembly_model(func)

        @preassembly_ns.expect(new_model)
        @preassembly_ns.route(('/%s' % func_name),
                              doc={'summary': short_doc})
        class NewFunction(PreassembleStatements):
            func_name = func_name

            def post(self):
                return super().post()

            post.__doc__ = doc


# Create resources for Sources namespace

# REACH
reach_text_model = api.inherit('ReachText', bio_text_model, {
    'offline': fields.Boolean(default=False),
    'url': fields.String(example=reach_text_url)
})
reach_json_model = api.model('ReachJSON', {'json': fields.String(example='{}')})
reach_pmc_model = api.model('ReachPMC', {
    'pmc_id': fields.String(example='PMC8511698'),
    'offline': fields.Boolean(default=False),
    'url': fields.String(example=reach_nxml_url)
})


@sources_ns.expect(reach_text_model)
@sources_ns.route('/reach/process_text')
class ReachProcessText(Resource):
    # TODO: REACH web service is down. Need to use a local Reach
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process text with REACH and return INDRA Statements.

        Parameters
        ----------
        text : str
            The text to be processed.

        offline : Optional[bool]
            If set to True, the REACH system is run offline via a JAR file.
            Otherwise (by default) the web service is called. Default: False

        url : Optional[str]
            URL for a REACH web service instance, which is used for reading if
            provided. If not provided but offline is set to False (its default
            value), REACH_TEXT_URL set in configuration will be used. If not
            provided in configuration, the Arizona REACH web service is called
            (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
            Default: None

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
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


@sources_ns.expect(reach_json_model)
@sources_ns.route('/reach/process_json')
class ReachProcessJson(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process REACH json and return INDRA Statements.

        Parameters
        ----------
        json : str
            The json string to be processed.

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        json_str = args.get('json')
        rp = reach.process_json_str(json_str)
        return _stmts_from_proc(rp)


@sources_ns.expect(reach_pmc_model)
@sources_ns.route('/reach/process_pmc')
class ReachProcessPmc(Resource):
    #TODO: REACH web service is down. Need to use a local Reach
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process PubMedCentral article and return INDRA Statements.

        Parameters
        ----------
        pmc_id : str
            The ID of a PubmedCentral article. The string may start with PMC
            but passing just the ID also works.
            Examples: 8511698, PMC8511698
            https://www.ncbi.nlm.nih.gov/pmc/

        offline : Optional[bool]
            If set to True, the REACH system is run offline via a JAR file.
            Otherwise (by default) the web service is called. Default: False

        url : Optional[str]
            URL for a REACH web service instance, which is used for reading if
            provided. If not provided but offline is set to False (its default
            value), REACH_NXML_URL set in configuration will be used. If not
            provided in configuration, the Arizona REACH web service is called
            (http://agathon.sista.arizona.edu:8080/odinweb/api/help).
            Default: None

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        pmcid = args.get('pmc_id')
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


# TRIPS
xml_model = api.model('XML', {'xml_str': fields.String})


@sources_ns.expect(bio_text_model)
@sources_ns.route('/trips/process_text')
class TripsProcessText(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process text with TRIPS and return INDRA Statements.

        Parameters
        ----------
        text : str
            The text to be processed.

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        text = args.get('text')
        tp = trips.process_text(text)
        return _stmts_from_proc(tp)


@sources_ns.expect(xml_model)
@sources_ns.route('/trips/process_xml')
class TripsProcessText(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process TRIPS EKB XML and return INDRA Statements.

        Parameters
        ----------
        xml_string : str
            A TRIPS extraction knowledge base (EKB) string to be processed.
            http://trips.ihmc.us/parser/api.html

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        xml_str = args.get('xml_str')
        tp = trips.process_xml(xml_str)
        return _stmts_from_proc(tp)


# Eidos
eidos_text_model = api.inherit('EidosText', bio_text_model, {
    'webservice': fields.String
})


# Hide docs until webservice is available
@sources_ns.expect(eidos_text_model)
@sources_ns.route('/eidos/process_text', doc=False)
class EidosProcessText(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process text with EIDOS and return biology INDRA Statements.

        Parameters
        ----------
        text : str
            The text to be processed.

        webservice : Optional[str]
            An Eidos reader web service URL to send the request to.
            If None, the reading is assumed to be done with the Eidos JAR
            rather than via a web service. Default: None

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        text = args.get('text')
        webservice = args.get('webservice')
        if not webservice:
            abort(400, 'No web service address provided.')
        ep = eidos.process_text_bio(text, webservice=webservice)
        return _stmts_from_proc(ep)


@sources_ns.expect(jsonld_model)
@sources_ns.route('/eidos/process_jsonld')
class EidosProcessJsonld(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process an EIDOS JSON-LD and return biology INDRA Statements.

        Parameters
        ----------
        jsonld : str
            The JSON-LD string to be processed.

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        eidos_json = args.get('jsonld')
        jj = json.loads(eidos_json)
        ep = eidos.process_json_bio(jj)
        return _stmts_from_proc(ep)


# BEL
bel_rdf_model = api.model('BelRdf', {'belrdf': fields.String})


@sources_ns.expect(genes_model)
@sources_ns.route('/bel/process_pybel_neighborhood')
class BelProcessNeighborhood(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process BEL Large Corpus neighborhood and return INDRA Statements.

        Parameters
        ----------
        genes : list[str]
            A list of entity names (e.g., gene names) which will be used as the
            basis of filtering the result. If any of the Agents of an extracted
            INDRA Statement has a name appearing in this list, the Statement is
            retained in the result.

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        genes = args.get('genes')
        bp = bel.process_pybel_neighborhood(genes)
        return _stmts_from_proc(bp)


@sources_ns.expect(bel_rdf_model)
@sources_ns.route('/bel/process_belrdf')
class BelProcessBelRdf(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process BEL RDF and return INDRA Statements.

        Parameters
        ----------
        belrdf : str
            A BEL/RDF string to be processed. This will usually come from
            reading a .rdf file.

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        belrdf = args.get('belrdf')
        bp = bel.process_belrdf(belrdf)
        return _stmts_from_proc(bp)


# BioPax
source_target_model = api.model('SourceTarget', {
    'source': fields.List(fields.String, example=['BRAF', 'RAF1', 'ARAF']),
    'target': fields.List(fields.String, example=['MAP2K1', 'MAP2K2'])
})


@sources_ns.expect(genes_model)
@sources_ns.route('/biopax/process_pc_pathsbetween')
class BiopaxPathsBetween(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """
        Process PathwayCommons paths between genes, return INDRA Statements.

        Parameters
        ----------
        genes : list
            A list of HGNC gene symbols to search for paths between.
            Examples: ['BRAF', 'MAP2K1']

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        genes = args.get('genes')
        bp = biopax.process_pc_pathsbetween(genes)
        return _stmts_from_proc(bp)


@sources_ns.expect(source_target_model)
@sources_ns.route('/biopax/process_pc_pathsfromto')
class BiopaxPathsFromTo(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """
        Process PathwayCommons paths from-to genes, return INDRA Statements.

        Parameters
        ----------
        source : list
            A list of HGNC gene symbols that are the sources of paths being
            searched for.
            Examples: ['BRAF', 'RAF1', 'ARAF']

        target : list
            A list of HGNC gene symbols that are the targets of paths being
            searched for.
            Examples: ['MAP2K1', 'MAP2K2']

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        source = args.get('source')
        target = args.get('target')
        bp = biopax.process_pc_pathsfromto(source, target)
        return _stmts_from_proc(bp)


@sources_ns.expect(genes_model)
@sources_ns.route('/biopax/process_pc_neighborhood')
class BiopaxNeighborhood(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Process PathwayCommons neighborhood, return INDRA Statements.

        Parameters
        ----------
        genes : list
            A list of HGNC gene symbols to search the neighborhood of.
            Examples: ['BRAF'], ['BRAF', 'MAP2K1']

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of extracted INDRA Statements.
        """
        args = request.json
        genes = args.get('genes')
        bp = biopax.process_pc_neighborhood(genes)
        return _stmts_from_proc(bp)


# Create resources for Assemblers namespace
pysb_stmts_model = api.inherit('PysbStatements', stmts_model, {
    'export_format': fields.String(example='kappa')
})


@assemblers_ns.expect(pysb_stmts_model)
@assemblers_ns.route('/pysb')
class AssemblePysb(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Assemble INDRA Statements and return PySB model string.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        export_format : str
            The format to export into, for instance "kappa", "bngl",
            "sbml", "matlab", "mathematica", "potterswheel". See
            http://pysb.readthedocs.io/en/latest/modules/export/index.html
            for a list of supported formats. In addition to the formats
            supported by PySB itself, this method also provides "sbgn"
            output.

        Returns
        -------
        image or model
            Assembled exported model. If export_format is kappa_im or kappa_cm,
            image is returned. Otherwise model string is returned.
        """
        args = request.json
        stmts_json = args.get('statements')
        export_format = args.get('export_format')
        stmts = stmts_from_json(stmts_json)
        pa = PysbAssembler()
        pa.add_statements(stmts)
        pa.make_model()
        try:
            for m in pa.model.monomers:
                pysb_assembler.set_extended_initial_condition(pa.model, m, 0)
        except Exception as e:
            logger.exception(e)

        if not export_format:
            model_str = pa.print_model()
        elif export_format in ('kappa_im', 'kappa_cm'):
            fname = 'model_%s.png' % export_format
            root = os.path.dirname(os.path.abspath(fname))
            graph = pa.export_model(format=export_format, file_name=fname)
            with open(fname, 'rb') as fh:
                data = 'data:image/png;base64,%s' % \
                    base64.b64encode(fh.read()).decode()
                return {'image': data}
        else:
            try:
                model_str = pa.export_model(format=export_format)
            except Exception as e:
                logger.exception(e)
                model_str = ''
        res = {'model': model_str}
        return res


@assemblers_ns.expect(stmts_model)
@assemblers_ns.route('/cx')
class AssembleCx(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Assemble INDRA Statements and return CX network json.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        Returns
        -------
        model
            Assembled model string.
        """
        args = request.json
        stmts_json = args.get('statements')
        stmts = stmts_from_json(stmts_json)
        ca = CxAssembler(stmts)
        model_str = ca.make_model()
        res = {'model': model_str}
        return res


@assemblers_ns.expect(stmts_model)
@assemblers_ns.route('/graph')
class AssembleGraph(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Assemble INDRA Statements and return Graphviz graph dot string.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        Returns
        -------
        model
            Assembled model string.
        """
        args = request.json
        stmts_json = args.get('statements')
        stmts = stmts_from_json(stmts_json)
        ga = GraphAssembler(stmts)
        model_str = ga.make_model()
        res = {'model': model_str}
        return res


@assemblers_ns.expect(stmts_model)
@assemblers_ns.route('/cyjs')
class AssembleCyjs(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Assemble INDRA Statements and return Cytoscape JS network.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        Returns
        -------
        json_model : dict
            Json dictionary containing graph information.
        """
        args = request.json
        stmts_json = args.get('statements')
        stmts = stmts_from_json(stmts_json)
        cja = CyJSAssembler(stmts)
        cja.make_model(grouping=True)
        model_str = cja.print_cyjs_graph()
        return json.loads(model_str)


@assemblers_ns.expect(stmts_model)
@assemblers_ns.route('/english')
class AssembleEnglish(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Assemble each statement into English sentence.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        Returns
        -------
        sentences : dict
            Dictionary mapping Statement UUIDs with English sentences.
        """
        args = request.json
        stmts_json = args.get('statements')
        stmts = stmts_from_json(stmts_json)
        sentences = {}
        for st in stmts:
            enga = EnglishAssembler()
            enga.add_statements([st])
            model_str = enga.make_model()
            sentences[st.uuid] = model_str
        res = {'sentences': sentences}
        return res


@assemblers_ns.expect(stmts_model)
@assemblers_ns.route('/sif/loopy')
class AssembleLoopy(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Assemble INDRA Statements into a Loopy model using SIF Assembler.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        Returns
        -------
        loopy_url : str
            Assembled Loopy model string.
        """
        args = request.json
        stmts_json = args.get('statements')
        stmts = stmts_from_json(stmts_json)
        sa = SifAssembler(stmts)
        sa.make_model(use_name_as_key=True)
        model_str = sa.print_loopy(as_url=True)
        res = {'loopy_url': model_str}
        return res


# Create resources for NDEx namespace
network_model = api.model('Network', {'network_id': fields.String})


@ndex_ns.expect(stmts_model)
@ndex_ns.route('/share_model_ndex')
class ShareModelNdex(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Upload the model to NDEX.

        Parameters
        ----------
        statements : list[indra.statements.Statement.to_json()]
            A list of INDRA Statements to assemble.

        Returns
        -------
        network_id : str
            ID of uploaded NDEx network.
        """
        args = request.json
        stmts_json = args.get('statements')
        stmts = stmts_from_json(stmts_json)
        ca = CxAssembler(stmts)
        for n, v in args.items():
            ca.cx['networkAttributes'].append({'n': n, 'v': v, 'd': 'string'})
        ca.make_model()
        network_id = ca.upload_model(private=False)
        return {'network_id': network_id}


@ndex_ns.expect(network_model)
@ndex_ns.route('/fetch_model_ndex')
class FetchModelNdex(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Download model and associated pieces from NDEX.

        Parameters
        ----------
        network_id : str
            ID of NDEx network to fetch.

        Returns
        -------
        stored_data : dict
            Dictionary representing the network.
        """
        args = request.json
        network_id = args.get('network_id')
        cx = process_ndex_network(network_id)
        network_attr = [x for x in cx.cx if x.get('networkAttributes')]
        network_attr = network_attr[0]['networkAttributes']
        keep_keys = ['txt_input', 'parser',
                     'model_elements', 'preset_pos', 'stmts',
                     'sentences', 'evidence', 'cell_line', 'mrna', 'mutations']
        stored_data = {}
        for d in network_attr:
            if d['n'] in keep_keys:
                stored_data[d['n']] = d['v']
        return stored_data


# Create resources for INDRA DB REST namespace
stmt_model = api.model('Statement', {'statement': fields.Nested(dict_model)})


@indra_db_rest_ns.expect(stmt_model)
@indra_db_rest_ns.route('/get_evidence')
class GetEvidence(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Get all evidence for a given INDRA statement.

        Parameters
        ----------
        statements : indra.statements.Statement.to_json()
            An INDRA Statement to get evidence for.

        Returns
        -------
        statements : list[indra.statements.Statement.to_json()]
            A list of retrieved INDRA Statements with evidence.
        """
        args = request.json
        stmt_json = args.get('statement')
        stmt = Statement._from_json(stmt_json)

        def _get_agent_ref(agent):
            """Get the preferred ref for an agent for db web api."""
            if agent is None:
                return None
            ag_hgnc_id = hgnc_client.get_hgnc_id(agent.name)
            if ag_hgnc_id is not None:
                return ag_hgnc_id + "@HGNC"
            db_refs = agent.db_refs
            for namespace in ['HGNC', 'FPLX', 'CHEBI', 'TEXT']:
                if namespace in db_refs.keys():
                    return '%s@%s' % (db_refs[namespace], namespace)
            return '%s@%s' % (agent.name, 'TEXT')

        def _get_matching_stmts(stmt_ref):
            # Filter by statement type.
            stmt_type = stmt_ref.__class__.__name__
            agent_name_list = [
                _get_agent_ref(ag) for ag in stmt_ref.agent_list()]
            non_binary_statements = (Complex, SelfModification, ActiveForm)
            # TODO: We should look at more than just the agent name.
            # Doing so efficiently may require changes to the web api.
            if isinstance(stmt_ref, non_binary_statements):
                agent_list = [ag_name for ag_name in agent_name_list
                              if ag_name is not None]
                kwargs = {}
            else:
                agent_list = []
                kwargs = {k: v for k, v in zip(['subject', 'object'],
                                               agent_name_list)}
                if not any(kwargs.values()):
                    return []
                print(agent_list)
            ip = get_statements(agents=agent_list, stmt_type=stmt_type,
                                **kwargs)
            return ip.statements

        stmts_out = _get_matching_stmts(stmt)
        agent_name_list = [ag.name for ag in stmt.agent_list()]
        stmts_out = stmts = filter_concept_names(
            stmts_out, agent_name_list, 'all')
        return _return_stmts(stmts_out)


# Create resources for Databases namespace
cbio_model = api.model('Cbio', {
    'gene_list': fields.List(fields.String, example=["FOSL1", "GRB2"]),
    'cell_lines': fields.List(fields.String, example=['COLO679_SKIN'])
})


@databases_ns.expect(cbio_model)
@databases_ns.route('/cbio/get_ccle_mrna')
class CbioMrna(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Get CCLE mRNA amounts using cBioClient

        Parameters
        ----------
        gene_list : list[str]
            A list of HGNC gene symbols to get mRNA amounts for.

        cell_lines : list[str]
            A list of CCLE cell line names to get mRNA amounts for.

        Returns
        -------
        mrna_amounts : dict[dict[float]]
            A dict keyed to cell lines containing a dict keyed to genes
            containing float
        """
        args = request.json
        gene_list = args.get('gene_list')
        cell_lines = args.get('cell_lines')
        mrna_amounts = cbio_client.get_ccle_mrna(gene_list, cell_lines)
        res = {'mrna_amounts': mrna_amounts}
        return res


@databases_ns.expect(cbio_model)
@databases_ns.route('/cbio/get_ccle_cna')
class CbioCna(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Get CCLE CNA

        -2 = homozygous deletion
        -1 = hemizygous deletion
        0 = neutral / no change
        1 = gain
        2 = high level amplification

        Parameters
        ----------
        gene_list : list[str]
            A list of HGNC gene symbols to get mutations in.

        cell_lines : list[str]
            A list of CCLE cell line names to get mutations for.

        Returns
        -------
        cna : dict[dict[int]]
            A dict keyed to cases containing a dict keyed to genes
            containing int
        """
        args = request.json
        gene_list = args.get('gene_list')
        cell_lines = args.get('cell_lines')
        cna = cbio_client.get_ccle_cna(gene_list, cell_lines)
        res = {'cna': cna}
        return res


@databases_ns.expect(cbio_model)
@databases_ns.route('/cbio/get_ccle_mutations')
class CbioMutations(Resource):
    @api.doc(False)
    def options(self):
        return {}

    def post(self):
        """Get CCLE mutations

        Parameters
        ----------
        gene_list : list[str]
            A list of HGNC gene symbols to get mutations in

        cell_lines : list[str]
            A list of CCLE cell line names to get mutations for.

        Returns
        -------
        mutations : dict
            The result from cBioPortal as a dict in the format
            {cell_line : {gene : [mutation1, mutation2, ...] }}
        """
        args = request.json
        gene_list = args.get('gene_list')
        cell_lines = args.get('cell_lines')
        mutations = cbio_client.get_ccle_mutations(gene_list, cell_lines)
        res = {'mutations': mutations}
        return res


if __name__ == '__main__':
    argparser = argparse.ArgumentParser('Run the INDRA REST API')
    argparser.add_argument('--host', default='0.0.0.0')
    argparser.add_argument('--port', default=8080, type=int)
    argparserargs = argparser.parse_args()
    app.run(host=argparserargs.host, port=argparserargs.port)
