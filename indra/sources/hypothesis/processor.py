import re
import logging
from indra.statements import BioContext, RefContext
from indra.preassembler.grounding_mapper.standardize import \
    standardize_db_refs, name_from_grounding

logger = logging.getLogger(__name__)


class HypothesisProcessor:
    """Processes hypothes.is annotations into INDRA Statements or groundings.

    Parameters
    ----------
    annotations : list[dict]
        A list of annotations fetched from hypothes.is in JSON-deserialized
        form represented as a list of dicts.
    reader : Optiona[function]
        A handle for a function which takes a single str argument
        (text to process) and returns a processor object with a statements
        attribute containing INDRA Statements. By default, the REACH reader's
        process_text function is used with default parameters. Note that
        if the function requires extra parameters other than the input text,
        functools.partial can be used to set those.
    grounder : Optional[function]
        A handle for a function which takes a positional str argument (entity
        text to ground) and an optional context key word argument and returns
        a list of objects matching the structure of gilda.grounder.ScoredMatch.
        By default, Gilda's ground function is used for grounding.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements extracted from the given annotations.
    groundings : dict
        A dict of entity text keys with an associated dict of grounding
        references.
    """
    def __init__(self, annotations, reader=None, grounder=None):
        self.annotations = annotations
        self.statements = []
        self.groundings = {}
        if reader is None:
            from indra.sources import reach
            self.reader = reach.process_text
        if grounder is None:
            from gilda import ground
            self.grounder = ground

    def extract_statements(self):
        """Sets statements attribute to list of extracted INDRA Statements."""
        for annotation in self.annotations:
            tags = annotation.get('tags')
            # Allow no tags or indra as a tag
            if not tags or 'indra' in tags:
                stmts = self.stmts_from_annotation(annotation)
                if stmts:
                    self.statements += stmts

    def extract_groundings(self):
        """Sets groundings attribute to list of extracted groundings."""
        for annotation in self.annotations:
            tags = annotation.get('tags')
            if tags and 'gilda' in tags:
                groundings = self.groundings_from_annotation(annotation)
                if groundings:
                    for txt, refs in groundings.items():
                        if txt in self.groundings and \
                                (self.groundings[txt] != refs):
                            logger.info(
                                'There is already a curation for %s: %s, '
                                'overwriting with %s' % (txt,
                                                         str(groundings[txt]),
                                                         str(refs)))

                        self.groundings[txt] = refs

    @staticmethod
    def groundings_from_annotation(annotation):
        """Return a dict of groundings from a single annotation."""
        text = annotation.get('text')
        if not text:
            return {}
        parts = [t for t in text.split('\n') if t]
        groundings = {}
        for entry in parts:
            grounding = parse_grounding_entry(entry)
            if grounding:
                groundings.update(grounding)
        return groundings

    def stmts_from_annotation(self, annotation):
        """Return a list of Statements extracted from a single annotation."""
        text = annotation.get('text')
        if not text:
            return []
        parts = [t for t in text.split('\n') if t]
        text = parts[0]
        rp = self.reader(text)
        if not rp or not rp.statements:
            logger.warning('Could not extract any statements from %s'
                           % text)
            return []

        contexts = {}
        # We assume that all other parts are related to context
        for part in parts[1:]:
            context_dict = get_context_entry(part, self.grounder, text)
            if context_dict:
                contexts.update(context_dict)
        bio_context = BioContext(**contexts) if contexts else None
        text_refs = get_text_refs(annotation['uri'])
        # In case we got multiple statements out, we apply the same
        # annotations to each
        for stmt in rp.statements:
            # There is expected to be exactly one evidence in all cases
            # but this is still a good way to work with it
            for ev in stmt.evidence:
                ev.source_api = 'hypothes.is'
                ev.text = text
                ev.text_refs = text_refs
                if 'PMID' in text_refs:
                    ev.pmid = text_refs['PMID']
                ev.annotations['hypothes.is'] = annotation
                ev.context = bio_context
        return rp.statements


def get_context_entry(entry, grounder, sentence):
    """Returnn a dict of context type and object processed from an entry."""
    match = re.match(r'(.*): (.*)', entry)
    if not match:
        return None
    context_type, context_txt = match.groups()
    if context_type not in allowed_contexts:
        logger.warning('Unknown context type %s' % context_type)
        return None

    terms = grounder(context_txt, context=sentence)
    if not terms:
        logger.warning('Could not ground %s context: %s'
                       % (context_type, context_txt))
    db_refs = {}
    if terms:
        db_refs = standardize_db_refs({terms[0].term.db:
                                       terms[0].term.id})
    db_refs['TEXT'] = context_txt
    standard_name = None
    if terms:
        standard_name = name_from_grounding(terms[0].term.db,
                                            terms[0].term.id)
    name = standard_name if standard_name else context_txt
    context = RefContext(name=name, db_refs=db_refs)
    return {allowed_contexts[context_type]: context}


def parse_grounding_entry(entry):
    """Return a dict representing single grounding curation entry string."""
    entry = entry.strip()
    # We now try to match the standard pattern for grounding curation
    match = re.match(r'^\[(.*)\] -> ([^ ]+)$', entry)
    # We log any instances of curations that don't match the pattern
    if not match:
        logger.warning('"%s" by %s does not match the grounding curation '
                       'pattern.' % entry)
        return None
    txt, dbid_str = match.groups()
    # We now get a dict of curated mappings to return
    try:
        dbid_entries = [entry.split(':', maxsplit=1)
                        for entry in dbid_str.split('|')]
        dbids = {k: v for k, v in dbid_entries}
    except Exception as e:
        logger.warning('Could not interpret DB IDs: %s for %s' %
                       (dbid_str, txt))
        return None
    return {txt: dbids}


def get_text_refs(url):
    """Return the parsed out text reference dict from an URL."""
    text_refs = {'URL': url}
    match = re.match(r'https://www.ncbi.nlm.nih.gov/pubmed/(\d+)', url)
    if match:
        text_refs['PMID'] = match.groups()[0]
    match = re.match(r'https://www.ncbi.nlm.nih.gov/pmc/articles/PMC(\d+)/',
                     url)
    if match:
        text_refs['PMCID'] = match.groups()[0]
    return text_refs


allowed_contexts = {
    'Location': 'location',
    'Cell line': 'cell_line',
    'Cell type': 'cell_type',
    'Organ': 'organ',
    'Disease': 'disease',
    'Species': 'species'
}
