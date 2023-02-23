from typing import Any, Callable, Mapping, Optional

from indra.statements import Agent
from indra.statements import Activation, Inhibition
from indra.ontology.standardize import standardize_agent_name
from .processor import EidosProcessor

GrounderResult = Mapping[str, str]
Grounder = Callable[[str, Optional[str]], GrounderResult]


class EidosBioProcessor(EidosProcessor):
    """Class to extract biology-oriented INDRA statements from Eidos output
    in a way that agents are grounded to biomedical ontologies."""

    def __init__(self, json_dict, grounder: Optional[Grounder] = None):
        super().__init__(json_dict)
        if grounder:
            self.grounder = grounder
        else:
            self.grounder = default_grounder_wrapper

    def get_regulate_activity(self, stmt):
        context = stmt.evidence[0].text
        subj = self.get_agent_bio(stmt.subj.concept, context=context)
        obj = self.get_agent_bio(stmt.obj.concept, context=context)
        if not subj or not obj:
            return None
        pol = stmt.overall_polarity()
        stmt_type = Activation if pol == 1 or not pol else Inhibition
        bio_stmt = stmt_type(subj, obj, evidence=stmt.evidence)
        return bio_stmt

    def extract_statements(self):
        self.extract_causal_relations()
        bio_stmts = []
        for stmt in self.statements:
            bio_stmt = self.get_regulate_activity(stmt)
            if bio_stmt:
                bio_stmts.append(bio_stmt)
        self.statements = bio_stmts

    def get_agent_bio(self, concept, context=None):
        return get_agent_bio(concept, context=context, grounder=self.grounder)


def get_agent_bio(concept, context=None, grounder: Optional[Grounder] = None):
    if not grounder:
        grounder = default_grounder_wrapper
    # Note that currently concept.name is the canonicalized entity text
    # whereas db_refs['TEXT'] is the unaltered original entity text
    raw_txt = concept.db_refs['TEXT']
    norm_txt = concept.name
    # We ground first the raw entity text and if that cannot be grounded,
    # the normalized entity text. The agent name is chosen based on the
    # first text that was successfully grounded, or if no grounding was
    # obtained, is chosen as the normalized text
    for txt in (raw_txt, norm_txt):
        gr = grounder(txt, context=context)
        if gr:
            name = txt
            break
    else:
        gr = {}
        name = norm_txt
    # We take whatever grounding and name are available and then
    # standardize the agent.
    agent = Agent(name, db_refs={'TEXT_NORM': norm_txt,
                                 'TEXT': raw_txt, **gr})
    standardize_agent_name(agent, standardize_refs=True)
    return agent


def default_grounder_wrapper(text: str, context: Optional[str]) -> GrounderResult:
    # Import here to avoid this when working in INDRA World context
    from indra.preassembler.grounding_mapper.gilda import get_grounding
    grounding, _ = get_grounding(text, context=context, mode='local')
    return grounding
