import re

from docutils import nodes, utils
from docutils.parsers.rst import roles

pubmed_uri_pattern = "http://www.ncbi.nlm.nih.gov/pubmed/%i"
doi_uri_pattern = "http://dx.doi.org/%s"

def pmid_reference_role(role, rawtext, text, lineno, inliner,
                        options={}, content=[]):
    try:
        pmid = int(text)
        if pmid <= 0:
            raise ValueError
    except ValueError:
        msg = inliner.reporter.error(
            'pmid number must be a number greater than or equal to 1; '
            '"%s" is invalid.' % text, line=lineno)
        prb = inliner.problematic(rawtext, rawtext, msg)
        return [prb], [msg]
    ref = pubmed_uri_pattern % pmid
    nodelist = []
    nodelist.append(nodes.inline(text='PMID:'))
    nodelist.append(nodes.reference(rawtext, utils.unescape(text), refuri=ref,
                                    **options))
    return nodelist, []

def doi_reference_role(role, rawtext, text, lineno, inliner,
                       options={}, content=[]):
    ref = doi_uri_pattern % text
    nodelist = []
    nodelist.append(nodes.inline(text='doi:'))
    nodelist.append(nodes.reference(rawtext, utils.unescape(text), refuri=ref,
                                    **options))
    return nodelist, []

def setup(app):
    app.add_role('pmid', pmid_reference_role)
    app.add_role('doi', doi_reference_role)
