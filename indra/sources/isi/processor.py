import os
import json
import logging
from copy import deepcopy
import indra.statements as ist
from indra.preassembler.grounding_mapper.gilda import ground_statements

logger = logging.getLogger(__name__)


class IsiProcessor(object):
    """Processes the output of the ISI reader.

    Parameters
    ----------
    reader_output : json
        The output JSON of the ISI reader as a json object.
    pmid : Optional[str]
        The PMID to assign to the extracted Statements
    extra_annotations : Optional[dict]
        Annotations to be included with each extracted Statement
    add_grounding : Optional[bool]
        If True, Gilda is used as a service to ground the Agents in
        the extracted Statements.

    Attributes
    ----------
    verbs : set[str]
        A list of verbs that have appeared in the processed ISI output
    statements : list[indra.statements.Statement]
        Extracted statements
    """
    def __init__(self, reader_output, pmid=None, extra_annotations=None,
                 add_grounding=False):
        self.reader_output = reader_output
        self.pmid = pmid
        self.extra_annotations = extra_annotations if \
            extra_annotations is not None else {}
        self.verbs = set()

        self.statements = []
        self.add_grounding = add_grounding

    def get_statements(self):
        """Process reader output to produce INDRA Statements."""
        for k, v in self.reader_output.items():
            for interaction in v['interactions']:
                self._process_interaction(k, interaction, v['text'], self.pmid,
                                          self.extra_annotations)
        if self.add_grounding:
            ground_statements(self.statements)

    def _process_interaction(self, source_id, interaction, text, pmid,
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
        pmid : str
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
        annotations['source_id'] = source_id
        annotations['interaction'] = interaction
        ev = ist.Evidence(source_api='isi',
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

            if statement_class == ist.Complex:
                statement = ist.Complex([subj, obj], evidence=ev)
            else:
                statement = statement_class(subj, obj, evidence=ev)

        if statement is not None:
            # For Complex statements, the ISI reader produces two events:
            # binds(A, B) and binds(B, A)
            # We want only one Complex statement for each sentence, so check
            # to see if we already have a Complex for this source_id with the
            # same members
            already_have = False
            if type(statement) == ist.Complex:
                for old_s in self.statements:
                    old_id = statement.evidence[0].source_id
                    new_id = old_s.evidence[0].source_id
                    if type(old_s) == ist.Complex and old_id == new_id:
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

    @staticmethod
    def _make_agent(agent_str):
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
        return ist.Agent(agent_str, db_refs={'TEXT': agent_str})


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
                        statement_class = getattr(ist, s_type)
                        verb_to_statement_type[verb] = statement_class
                    except Exception:
                        pass
            else:
                first_line = False
    return verb_to_statement_type


verb_to_statement_type = _build_verb_statement_mapping()
