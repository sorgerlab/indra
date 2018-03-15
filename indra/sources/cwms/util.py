import os
import rdflib
import logging
from rdflib.tools import rdf2dot
from indra.sources.trips import trips_client as tc

try:
    import pydot
    HAVE_PYDOT = True
except ImportError:
    print("Could not import module pydot. Will not be able to create pdf "
          "cag graphs.")
    HAVE_PYDOT = False

logger = logging.getLogger('cwms_util')


def make_cag_graphs(cag_rdf_dicts, skip_if_no_dot=True):
    name_fmt = 'cag_%s_sentence_%d'
    cag_graphs = {}
    for cag_label, rdf_dict in cag_rdf_dicts.items():
        cag_graphs[cag_label] = {}
        for i_sent, (sentence, rdf_str) in enumerate(rdf_dict.items()):
            g_rdf, g_dot = make_graphs(rdf_str, name_fmt % (cag_label, i_sent))
            cag_graphs[cag_label][sentence] = (g_rdf, g_dot)
            if g_dot is None:
                if skip_if_no_dot:
                    del cag_graphs[cag_label]
                    logger.warning("No dot graph produced. Aborting %s."
                                   % cag_label)
                    break
                else:
                    continue
            g_dot.set_label(sentence)
            g_dot.write_pdf('sources/cwms/example_pdfs/%s.pdf'
                            % (name_fmt % (cag_label, i_sent)))
    return cag_graphs


def make_graphs(rdf_str, fname_base):
    try:
        rdf_fname = 'sources/cwms/example_rdfs/' + fname_base + '.rdf'
        with open(rdf_fname, 'w') as f:
            f.write(rdf_str)
        g_rdf = rdflib.Graph()
        g_rdf.parse(rdf_fname, publicID=fname_base)
        dot_fname = 'sources/cwms/example_dots/' + fname_base + '.dot'
        with open(dot_fname, 'w') as f:
            rdf2dot.rdf2dot(g_rdf, f)
        if HAVE_PYDOT:
            g_dot = pydot.graph_from_dot_file(dot_fname)[0]
            logger.warning("pydot is not available, so pydot graph was not "
                           "produced. To get pydot, run `pip install pydot`.")
        else:
            g_dot = None
    except Exception as e:
        logger.error("Got exeption for rdf string: %s" % rdf_str)
        logger.exception(e)
        return (None, None)
    return (g_rdf, g_dot)


def get_example_extractions(fname):
    "Get extractions from one of the examples in `cag_examples`."
    with open(fname, 'r') as f:
        sentences = f.read().splitlines()
    rdf_xml_dict = {}
    for sentence in sentences:
        logger.info("Reading \"%s\"..." % sentence)
        html = tc.send_query(sentence, 'cwms')
        try:
            rdf_xml_dict[sentence] = tc.get_xml(html, 'rdf:RDF',
                                                fail_if_empty=True)
        except AssertionError as e:
            logger.error("Got error for %s." % sentence)
            logger.exception(e)
    return rdf_xml_dict


def make_example_graphs():
    "Make graphs from all the examples in cag_examples."
    cag_example_rdfs = {}
    for i, fname in enumerate(os.listdir('cag_examples')):
        cag_example_rdfs[i+1] = get_example_extractions(fname)
    return make_cag_graphs(cag_example_rdfs)
