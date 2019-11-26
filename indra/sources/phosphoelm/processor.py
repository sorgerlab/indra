import logging
import requests

from indra.statements import Phosphorylation, Evidence, Agent

gilda_url = 'http://grounding.indra.bio/ground'
logger = logging.getLogger(__name__)


def _gilda_grounder(entity_str):
    # Try to find a namespace for the enzyme entity string, return the string
    # that provided the match
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
    def __init__(self, phosphoelm_data=None):
        self.statements = []
        self._phosphoelm_data = phosphoelm_data

    def process_phosphorylations(self, keep_empty=False):
        """Create Phosphorylation statements from the json entries

        Parameters
        ----------
        keep_empty : bool
            If true, also create statements when upstream kinases in
            entry['kinases'] are not known.

        Returns
        -------
        statements : list[indra.statement.Phosphorylation]
            A list of the phosphorylation statements produced by the entries
            in the json
        """
        if self._phosphoelm_data is None:
            return
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
        strip_words = ['_group', '_drome', '_Caeel']
        # Pre process: strip 'strip words' and any trailing space
        for word in strip_words:
            upstream_kinase = upstream_kinase.replace(word, '').rstrip()
        used_str, ns, id = _gilda_grounder(upstream_kinase)

        # Split on '_'
        if ns is None and id is None and '_' in used_str:
            used_str, suffix = used_str.split('_')
            used_str, ns, id = _gilda_grounder(used_str)

        # Split on '/'
        if ns is None and id is None and '/' in used_str:
            used_str = used_str.split('/')[0]
            used_str, ns, id = _gilda_grounder(used_str)

        if ns is None and id is None:
            ns = 'TEXT'
            id = used_str

        ag = Agent(None, db_refs={ns: id})
        return used_str, ag
