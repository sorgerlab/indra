import logging
import requests

from indra.statements import Phosphorylation, Evidence, Agent

from .phosphoelm_mapping import phosphoelm_mapping

gilda_url = 'http://grounding.indra.bio/ground'
logger = logging.getLogger(__name__)


def _gilda_grounder(entity_str):
    # Try to find a namespace for the enzyme entity string, return the
    # string that provided the match
    res = requests.post(gilda_url, json={'text': entity_str})
    if res.status_code == 200 and res.json():
        db_ns = res.json()[0]['term']['db']
        db_id = res.json()[0]['term']['id']
        return entity_str, db_ns, db_id
    else:
        if res.status_code != 200:
            logger.warning('Gilda service responded with status code %d' %
                           res.status_code)
        return entity_str, None, None


class PhosphoElmProcessor(object):
    """Processes data dumps from the phospho.ELM database.

    See http://phospho.elm.eu.org/dataset.html

    Parameters
    ----------
    phosphoelm_data : list(dict)
        JSON compatible list of entries from a phospho.ELM data dump

    Attributes
    ----------
    statements : list(indra.statement.Phosphorylation)
        A list of the phosphorylation statements produced by the entries
        in phosphoelm_data
    """
    def __init__(self, phosphoelm_data):
        self.statements = []
        self._phosphoelm_data = phosphoelm_data

    def process_phosphorylations(self, keep_empty=False):
        """Create Phosphorylation statements from phosphoelm_data

        Parameters
        ----------
        keep_empty : bool
            Default: False. If true, also create statements when upstream
            kinases in entry['kinases'] are not known.
        """
        for entry in self._phosphoelm_data:
            if not keep_empty and not entry['kinases'] or\
                    not entry['species'].lower() == 'homo sapiens':
                # Skip entries without any kinases or if species is other
                # than human.
                continue
            # Entries:
            # 'acc': '<UP ID>', <-- substrate
            # 'sequence': '<protein sequence>',
            # 'position': '<sequence position>',
            # 'code': '<phosphorylated residue>',
            # 'pmids': '<pmid>',
            # 'kinases': '<responsible kinase>', <-- enzyme
            # 'source': 'HTP|LTP',
            # 'species': '<species name in latin>',
            # 'entry_date': 'yyyy-mm-dd HH:MM:SS.mmmmmm'
            substrate = Agent(None, db_refs={'UP': entry['acc']})
            used_name, enz = self._get_enzyme(entry['kinases']) if\
                entry.get('kinases') else '', None

            # Skip if we hit one of the special cases
            if enz is None:
                continue

            # Build evidence, add statement
            evidence = Evidence(
                source_api='phospho.ELM',
                pmid=entry['pmids'],
                annotations={
                    'data_source': 'High-ThroughPut' if
                    entry['source'].lower == 'htp' else (
                        'Low-ThroughPut' if entry['source'].lower() == 'ltp'
                        else None),
                    'phosphoelm_substrate_name': entry['acc'],
                    'phosphoelm_kinase_name': entry.get('kinases', None),
                    'used_kinse_name': used_name,
                    'entry_date': entry['entry_date'],
                    'sequence': entry['sequence']
                }
            )
            self.statements.append(Phosphorylation(
                enz=enz,
                sub=substrate,
                residue=entry['code'],
                position=entry['position'],
                evidence=evidence)
            )

    @staticmethod
    def _get_enzyme(upstream_kinase):
        """Handle the upstream kinases

        Parameters
        ----------
        upstream_kinase : str
            The string occuring in the entry 'upstream_kinases'

        Returns
        -------
        used_str, ag : tuple(str, indra.statements.Agent)
            A tuple containing a string and an agent. The string is the
            string from within 'upstream_kinases' that was used to create
            the agent.
        """
        # Check if kinase is in the hard coded mapping
        if upstream_kinase in phosphoelm_mapping:
            ns, _id = phosphoelm_mapping[upstream_kinase]
            if ns is None and _id is None:
                return upstream_kinase, None
            return upstream_kinase, Agent(None, db_refs={ns: _id})

        strip_words = ['_group', '_drome', '_Caeel']
        # Pre process: strip 'strip words' and any trailing space
        for word in strip_words:
            upstream_kinase = upstream_kinase.replace(word, '').rstrip()
        used_str, ns, _id = _gilda_grounder(upstream_kinase)

        # Split on '/'
        if ns is None and _id is None and '/' in used_str:
            used_str = used_str.split('/')[0]
            used_str, ns, _id = _gilda_grounder(used_str)

        # Replace '_' with '-'
        if ns is None and _id is None and '_' in used_str:
            used_str = used_str.replace('_', '-')
            used_str, ns, _id = _gilda_grounder(used_str)

        if ns is None and _id is None:
            ns = 'TEXT'
            _id = used_str

        ag = Agent(None, db_refs={ns: _id})
        return used_str, ag
