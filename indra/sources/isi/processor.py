import os
import json
import logging
import re
from copy import deepcopy

import indra.statements
from indra.statements import *
from indra.preassembler.grounding_mapper import load_grounding_map, \
        GroundingMapper

logger = logging.getLogger('isi')

# Load the mapping between ISI verb and INDRA statement type
def _build_verb_statement_mapping():
    """Build the mapping between ISI verb strings and INDRA statement classes.

    Looks up the INDRA statement class name, if any, in a resource file,
    and resolves this class name to a class.

    Returns
    -------
    verb_to_statement_type : dict
        Dictionary mapping verb name to an INDRA statment class
    """
    path_this = os.path.dirname(os.path.abspath(__file__))
    map_path = os.path.join(path_this, 'isi_verb_to_indra_statement_type.tsv')
    with open(map_path, 'r') as f:
        first_line = True
        verb_to_statement_type = {}
        for line in f:
            if not first_line:
                line = line[:-1]
                tokens = line.split('\t')

                if len(tokens) == 2 and len(tokens[1]) > 0:
                    verb = tokens[0]
                    s_type = tokens[1]
                    try:
                        statement_class = getattr(indra.statements, s_type)
                        verb_to_statement_type[verb] = statement_class
                    except Exception:
                        pass
            else:
                first_line = False
    return verb_to_statement_type


verb_to_statement_type = _build_verb_statement_mapping()


class IsiProcessor(object):
    """Processes the output of the ISI reader.

    Attributes
    ----------
    output_dir : str
        The output directory of the ISI reader, with one file per input file.
    verbs : set[str]
        A list of verbs that have appeared in the processed ISI output
    pmids : dict
        Dictionary mapping preprocessed text file ids to pmids
    extra_annotations : dict
        Dictionary mapping preprocessed text file ids to a associated
        annotations to be included with each statement from the document
    statements : list[indra.statements.Statement]
        Extracted statements
    """
    def __init__(self, output_dir, pmids=None, extra_annotations=None):
        self.pmids = pmids if pmids is not None else {}
        self.extra_annotations_each_doc = extra_annotations if \
            extra_annotations is not None else {}

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
        dir_name : str
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
        filename : str
            The ISI output file from which to extract statements
        pmid : int
            The PMID of the document being preprocessed, or None if not
            specified
        extra_annotations : dict
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
        source_id : str
            the JSON key corresponding to the sentence in the ISI output
            interaction: the JSON list with subject/verb/object information
            about the event in the ISI output
        text : str
            the text of the sentence
        pmid : int
            the PMID of the article from which the information was extracted
        extra_annotations : dict
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
        if 'interaction' in extra_annotations:
            logger.warning("'interaction' key of extra_annotations ignored" +
                           " since this is reserved for storing the raw ISI " +
                           "input.")
        annotations['interaction'] = interaction
        ev = Evidence(source_api='isi',
                      source_id=source_id,
                      pmid=pmid,
                      text=text.rstrip(),
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
        if verb in verb_to_statement_type:
            statement_class = verb_to_statement_type[verb]

            if statement_class == Complex:
                statement = Complex([subj, obj], evidence=ev)
            else:
                statement = statement_class(subj, obj, evidence=ev)

        if statement is not None:
            # For Complex statements, the ISI reader produces two events:
            # binds(A, B) and binds(B, A)
            # We want only one Complex statement for each sentence, so check
            # to see if we already have a Complex for this source_id with the
            # same members
            already_have = False
            if type(statement) == Complex:
                for old_s in self.statements:
                    old_id = statement.evidence[0].source_id
                    new_id = old_s.evidence[0].source_id
                    if type(old_s) == Complex and old_id == new_id:
                        old_statement_members = \
                                [m.db_refs['TEXT'] for m in old_s.members]
                        old_statement_members = sorted(old_statement_members)

                        new_statement_members = [m.db_refs['TEXT']
                                                 for m in statement.members]
                        new_statement_members = sorted(new_statement_members)

                        if old_statement_members == new_statement_members:
                            already_have = True
                            break

            if not already_have:
                self.statements.append(statement)

    def _make_agent(self, agent_str):
        """Makes an ungrounded Agent object from a string specifying an
        entity.

        Parameters
        ----------
        agent_str : str
            A string specifying the agent

        Returns
        -------
        agent : indra.statements.Agent
            An ungrounded Agent object referring to the specified text
        """
        return Agent(agent_str, db_refs={'TEXT': agent_str})
