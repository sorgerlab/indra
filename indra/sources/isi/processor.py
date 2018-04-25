import os
import json
import logging

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
    statements: list[indra.statements.Statement]
        Extracted statements
    """
    def __init__(self, output_dir):
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
                self.process_directory(full_entry_path)
            elif entry.endswith('.json'):
                self.process_file(full_entry_path)

    def process_file(self, filename):
        """Extracts statements from the given ISI output file.

        Parameters
        ----------
        filename: str
            The ISI output file from which to extract statements
        """
        print('Extracting from', filename)
        with open(filename, 'r') as f:
            j = json.load(f)
            for k in j:
                text = j[k]['text']
                interactions = j[k]['interactions']
                for interaction in interactions:
                    self.process_interaction(k, interaction, text)

    def process_interaction(self, source_id, interaction, text):
        verb = interaction[0].lower()
        subj = interaction[-2]
        obj = interaction[-1]

        # Make ungrounded agent objects for the subject and object
        # Grounding will happen after all statements are extracted in __init__
        subj = self._make_agent(subj)
        obj = self._make_agent(obj)

        # For binding time interactions, it is said that a catayst might be
        # specified. We don't use this for now, but extract in case we want
        # to in the future
        cataylst_specified = False
        if len(interaction) == 4:
            catalyst = interaction[1]
            if catalyst is not None:
                cataylst_specified=True
                print(catalyst)
        self.verbs.add(verb)
        #print(subj, verb, obj)

        statement = None
        if verb == 'transcription':
            pass
        elif verb == 'binding':
            statement = Complex([subj, obj])
        elif verb == 'complexed':
            statement = Complex([subj, obj])
        elif verb == 'binds':
            statement = Complex([subj, obj])
        elif verb == 'formation':
            pass
        elif verb == 'dimer':
            pass
        elif verb == 'dissociation':
            pass
        elif verb == 'form':
            pass
        elif verb == 'bound':
            statement = Complex([subj, obj])
        elif verb == 'complex':
            statement = Complex([subj, obj])
        elif verb == 'association':
            pass
        elif verb == 'forms':
            pass
        elif verb == 'interacts':
            pass
        elif verb == 'associate':
            pass
        elif verb == 'phosphorylation':
            statement = Phosphorylation(subj, obj)
        elif verb == 'phosphorylates':
            statement = Phosphorylation(subj, obj)
        elif verb == 'complexes':
            pass
        elif verb == 'inhibition':
            pass
        elif verb == 'bind':
            statement = Complex([subj, obj])
        elif verb == 'expressing':
            pass
        elif verb == 'associated':
            pass
        elif verb == 'hydrolysis':
            pass
        elif verb == 'phosphorylated':
            statement = Phosphorylation(subj, obj)
        elif verb == 'interact':
            pass
        elif verb == 'interaction':
            pass
        else:
            logger.error('Do not know how to process verb:', verb)

        if statement is not None:
            self.statements.append(statement)
    
    def _make_agent(self, agent_str):
        return Agent(agent_str, db_refs={'TEXT': agent_str})

if __name__ == '__main__':
    f = '/Users/daniel/workspace/isi/output/test1.json'
    d = '/Users/daniel/workspace/isi/output'
    ip = IsiProcessor(d)
    print('Verbs:', ip.verbs)
    print(ip.statements)
