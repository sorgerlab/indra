import requests
import logging

from indra.statements import Phosphorylation, Evidence, Agent

gilda_url = 'http://grounding.indra.bio/ground'
logger = logging.getLogger(__file__)


def _gilda_grounder(entity_str):
    # If match found, return the string that provided the match
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


class PhosphoELMPRocessor(object):
    def __init__(self, file_dump_json=None, keep_empty=False):
        self.statements = []
        self.statements.extend(self._from_file_dump_json(file_dump_json,
                                                         keep_empty))

    def _from_file_dump_json(self, fd_json, keep_empty=False):
        """Structuring the json entry to Phosphorylation statements

        fd_json : list(json)
            JSON comatible list of entries
        keep_empty : bool
            If true, also create statements when upstream kinases
            (in entry['kinases']) are not known.

        Returns
        -------
        """
        if fd_json is None:
            return []
        statements = []
        for entry in fd_json:
            if not keep_empty and not entry['kinases'] or\
                    not entry['species'].lower() == 'homo sapiens':
                # Skip entries without any kinases or if species is other
                # than human.
                continue
            # Entries:
            # 'acc': '<UP ID>', <-- substrate
            # 'sequence': '<protein sequence>',
            # 'position': '<position>',
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
                source_api='phosphoELM',
                pmid=entry['pmids'],
                annotations={
                    'data_source': 'High-ThroughPut' if
                    entry['source'].lower == 'htp' else 'Low-ThroughPut',
                    'phosphoelm_substrate_name': entry['acc'],
                    'phosphoelm_kinase_name': entry.get('kinases', ''),
                    'used_kinse_name': used_name,
                    'entry_date': entry['entry_date'],
                    'sequence': entry['sequence']
                }
            )
            statements.append(Phosphorylation(
                enz=enz,
                sub=substrate,
                residue=entry['code'],
                position=entry['position'],
                evidence=evidence)
            )
        return statements

    @staticmethod
    def _get_enzyme(upstream_kinase):
        """Handle the upstream kinases

        Parameters
        ----------
        upstream_kinase : str
            The string occuring in the entry 'upstream_kinases'

        Returns
        -------
        kinases : indra.statements.Agent
            The agents contained in 'upstream_kinases'
        """
        strip_words = ['_group', 'kinase', '_drome', '_Caeel']
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
