import os
import json
import logging
import re
from copy import deepcopy

from indra.statements import *
from indra.preassembler.grounding_mapper import load_grounding_map,\
        GroundingMapper

logger = logging.getLogger('isi')


class IsiProcessor(object):
    """Processes the output of the ISI reader.

    Attributes
    ----------
    output_dir: str
        The output directory of the ISI reader, with one file per input file.
    verbs: set[str]
        A list of verbs that have appeared in the processed ISI output
    pmids: dict
        Dictionary mapping preprocessed text file ids to pmids
    extra_annotations: dict
        Dictionary mapping preprocessed text file ids to a associated
        annotations to be included with each statement from the document
    statements: list[indra.statements.Statement]
        Extracted statements
    """
    def __init__(self, output_dir, preprocessor):
        self.pmids = preprocessor.pmids
        self.extra_annotations_each_doc = preprocessor.extra_annotations

        # Load grounding information
        path_this = os.path.dirname(os.path.abspath(__file__))
        gm_fname = os.path.join(path_this, '../../resources/',
                                'extracted_reach_grounding_map.csv')
        try:
            gm = load_grounding_map(gm_fname)
        except BaseException:
            raise Exception('Could not load the grounding map from ' +
                            gm_fname)
        mapper = GroundingMapper(gm)

        # Extract statements
        self.statements = []
        self.verbs = set()
        if output_dir is not None:
            self.process_directory(output_dir)

        # Ground statements
        self.statements = mapper.map_agents(self.statements)

    def process_directory(self, dir_name):
        """Recursively extracts statements from all ISI output files in the
        given directory and subdirectories.

        Parameters
        ----------
        dir_name: str
            The directory to traverse
        """
        for entry in os.listdir(dir_name):
            full_entry_path = os.path.join(dir_name, entry)

            if os.path.isdir(full_entry_path):
                logger.warning('ISI processor: did not expect any ' +
                               'subdirectories in the output directory.')
                self.process_directory(full_entry_path)
            elif entry.endswith('.json'):
                # Extract the corresponding file id
                m = re.match('([0-9]+)\.json', entry)
                if m is None:
                    logger.warning('ISI processor:', entry, ' does not ' +
                                   ' match expected format for output files.')
                    pmid = None
                    extra_annotations = {}
                else:
                    doc_id = int(m.group(1))
                    pmid = self.pmids.get(doc_id)
                    extra_annotations = self.extra_annotations_each_doc.get(
                            doc_id)
                self.process_file(full_entry_path, pmid, extra_annotations)
            else:
                logger.warning('ISI processor: did not expect any non-json ' +
                               'files in the output directory')

    def process_file(self, filename, pmid, extra_annotations):
        """Extracts statements from the given ISI output file.

        Parameters
        ----------
        filename: str
            The ISI output file from which to extract statements
        pmid: int
            The PMID of the document being preprocessed, or None if not
            specified
        extra_annotations: dict
            Extra annotations to be added to each statement from this document
            (can be the empty dictionary)
        """
        print('Extracting from', filename)
        with open(filename, 'r') as f:
            j = json.load(f)
            for k in j:
                text = j[k]['text']
                interactions = j[k]['interactions']
                for interaction in interactions:
                    self.process_interaction(k, interaction, text, pmid,
                                             extra_annotations)

    def process_interaction(self, source_id, interaction, text, pmid,
                            extra_annotations):
        """Process an interaction JSON tuple from the ISI output, and adds up
        to one statement to the list of extracted statements.

        Parameters
        ----------
        source_id: str
            the JSON key corresponding to the sentence in the ISI output
            interaction: the JSON list with subject/verb/object information
            about the event in the ISI output
        text: str
            the text of the sentence
        pmid: int
            the PMID of the article from which the information was extracted
        extra_annotations: dict
            Additional annotations to add to the statement's evidence,
            potentially containing metadata about the source. Annotations
            with the key "interaction" will be overridden by the JSON
            interaction tuple from the ISI output
        """
        verb = interaction[0].lower()
        subj = interaction[-2]
        obj = interaction[-1]

        # Make ungrounded agent objects for the subject and object
        # Grounding will happen after all statements are extracted in __init__
        subj = self._make_agent(subj)
        obj = self._make_agent(obj)

        # Make an evidence object
        annotations = deepcopy(extra_annotations)
        annotations['interaction'] = interaction
        ev = Evidence(source_api='isi',
                      source_id=source_id,
                      pmid=pmid,
                      text=text,
                      annotations=annotations)

        # For binding time interactions, it is said that a catayst might be
        # specified. We don't use this for now, but extract in case we want
        # to in the future
        cataylst_specified = False
        if len(interaction) == 4:
            catalyst = interaction[1]
            if catalyst is not None:
                cataylst_specified = True
        self.verbs.add(verb)

        statement = None
        if verb == 'transcription':
            pass
        elif verb == 'binding':
            statement = Complex([subj, obj], evidence=ev)
        elif verb == 'complexed':
            statement = Complex([subj, obj], evidence=ev)
        elif verb == 'binds':
            statement = Complex([subj, obj], evidence=ev)
        elif verb == 'formation':
            pass
        elif verb == 'dimer':
            pass
        elif verb == 'dissociation':
            pass
        elif verb == 'form':
            pass
        elif verb == 'bound':
            statement = Complex([subj, obj], evidence=ev)
        elif verb == 'complex':
            statement = Complex([subj, obj], evidence=ev)
        elif verb == 'association':
            pass
        elif verb == 'forms':
            pass
        elif verb == 'interacts':
            pass
        elif verb == 'associate':
            pass
        elif verb == 'phosphorylation':
            statement = Phosphorylation(subj, obj, evidence=ev)
        elif verb == 'phosphorylates':
            statement = Phosphorylation(subj, obj, evidence=ev)
        elif verb == 'complexes':
            pass
        elif verb == 'inhibition':
            pass
        elif verb == 'bind':
            statement = Complex([subj, obj], evidence=ev)
        elif verb == 'expressing':
            pass
        elif verb == 'associated':
            pass
        elif verb == 'hydrolysis':
            pass
        elif verb == 'phosphorylated':
            statement = Phosphorylation(subj, obj, evidence=ev)
        elif verb == 'interact':
            pass
        elif verb == 'interaction':
            pass
        else:
            logger.error('Do not know how to process verb:', verb)

        if statement is not None:
            self.statements.append(statement)

    def _make_agent(self, agent_str):
        """Makes an ungrounded Agent object from a string specifying an
        entity.

        Parameters
        ----------
        agent_str: str
            A string specifying the agent

        Returns
        -------
        agent: indra.statements.Agent
            An ungrounded Agent object referring to the specified text
        """
        return Agent(agent_str, db_refs={'TEXT': agent_str})


if __name__ == '__main__':
    f = '/Users/daniel/workspace/isi/output/test1.json'
    d = '/Users/daniel/workspace/isi/output'
    ip = IsiProcessor(d)
    print('Verbs:', ip.verbs)
    print(ip.statements)
