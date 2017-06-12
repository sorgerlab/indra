from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import copy
import logging
import lxml.etree
import lxml.builder
from indra.statements import *

logger = logging.getLogger('sbgn_assembler')

sbgn_ns = 'http://sbgn.org/libsbgn/pd/0.1'
emaker = lxml.builder.ElementMaker(nsmap={None: sbgn_ns})

class SBGNAssembler(object):
    def __init__(self, statements=None):
        if not statements:
            self.statements = []
        else:
            self.statements = statements
        self.agent_set = None
        self.id_counter = 0
        self.sbgn = None
        self._map = None
        self.agent_ids = {}

    def make_model(self):
        self.sbgn = emaker.sbgn()
        self._map = emaker.map()
        self.sbgn.append(self._map)
        for stmt in self.statements:
            if isinstance(stmt, Modification):
                self._assemble_modification(stmt)
            elif isinstance(stmt, RegulateActivity):
                self._assemble_regulateactivity(stmt)
            elif isinstance(stmt, RegulateAmount):
                self._assemble_regulateamount(stmt)
            elif isinstance(stmt, Complex):
                self._assemble_complex(stmt)
            else:
                logger.warning("Unhandled Statement type %s" % type(s))
                continue
        return lxml.etree.tostring(self.sbgn, pretty_print=True)

    def _assemble_modification(self, stmt):
        if not stmt.enz:
            return
        # Make glyph for enz
        enz_glyph = self._agent_glyph(stmt.enz)
        mc_changed = stmt._get_mod_condition()
        # Make glyphs for sub
        sub_changed = copy.deepcopy(stmt.sub)
        sub_changed.mods.append(mc_changed)
        sub_in, sub_out = \
            (stmt.sub, sub_changed) if isinstance(stmt, AddModification) else \
            (sub_changed, stmt.sub)
        sub_in_glyph = self._agent_glyph(sub_in)
        sub_out_glyph = self._agent_glyph(sub_out)
        # Make the process glyph
        process_glyph = self._process_glyph('process')
        # Add the arcs
        self._arc('consumption', sub_in_glyph, process_glyph)
        self._arc('production', process_glyph, sub_out_glyph)
        self._arc('catalysis', enz_glyph, process_glyph)

    def _assemble_regulateactivity(self, stmt):
        # Make glyph for subj
        subj_glyph = self._agent_glyph(stmt.subj)
        # Make glyphs for obj
        obj_act = copy.deepcopy(stmt.obj)
        obj_inact = copy.deepcopy(stmt.obj)
        obj_act.activity = ActivityCondition(stmt.obj_activity, True)
        obj_inact.activity = ActivityCondition(stmt.obj_activity, False)
        obj_in, obj_out = (obj_inact, obj_act) if stmt.is_activation else \
                          (obj_act, obj_inact)
        obj_in_glyph = self._agent_glyph(obj_in)
        obj_out_glyph = self._agent_glyph(obj_out)
        # Make the process glyph
        process_glyph = self._process_glyph('process')
        # Add the arcs
        self._arc('consumption', obj_in_glyph, process_glyph)
        self._arc('production', process_glyph, obj_out_glyph)
        self._arc('catalysis', subj_glyph, process_glyph)

    def _assemble_regulateamount(self, stmt):
        # Make glyphs for obj
        obj_glyph = self._agent_glyph(stmt.obj)
        obj_none_glyph = self._none_glyph()
        obj_in_glyph, obj_out_glyph = \
            (obj_none_glyph, obj_glyph) if \
            isinstance(stmt, IncreaseAmount) else \
            (obj_glyph, obj_none_glyph)
        # Make the process glyph
        process_glyph = self._process_glyph('process')
        # Add the arcs
        self._arc('consumption', obj_in_glyph, process_glyph)
        self._arc('production', process_glyph, obj_out_glyph)
        # Make glyph for subj and add arc if needed
        if stmt.subj:
            subj_glyph = self._agent_glyph(stmt.subj)
            self._arc('catalysis', subj_glyph, process_glyph)

    def _assemble_complex(self, stmt):
        # Make glyph for individual members
        member_glyphs = [self._agent_glyph(m) for m in stmt.members]
        # Make glyph for complex
        # First we need to unroll all members and their bound conditions
        # into a single list with a single prime agent
        all_members = []
        for i, member in enumerate(stmt.members):
            member_tmp = copy.deepcopy(member)
            bound = [bc.agent for bc in member_tmp.bound_conditions
                     if bc.is_bound]
            member_tmp.bound_conditions = []
            if i==0:
                prime_agent = member_tmp
            else:
                all_members.append(member_tmp)
            all_members += bound
        # Now we set all the other members as bound conditions on the prime
        # agent
        prime_agent.bound_conditions = \
            [BoundCondition(m, True) for m in all_members]
        complex_glyph = self._agent_glyph(prime_agent)
        process_glyph = self._process_glyph('association')
        for member_glyph in member_glyphs:
            self._arc('consumption', member_glyph, process_glyph)
        self._arc('production', process_glyph, complex_glyph)


    def _arc(self, class_name, source, target):
        arc_id = self._make_id()
        arc = emaker.arc(class_(class_name), source=source, target=target,
                         id=arc_id)
        self._map.append(arc)

    def _process_glyph(self, class_name):
        process_id = self._make_id()
        process_glyph = emaker.glyph(emaker.bbox(x='0', y='0', w='20', h='20'),
                                     class_(class_name), id=process_id)
        self._map.append(process_glyph)
        return process_id

    def _none_glyph(self):
        glyph_id = self._make_id()
        none_glyph = emaker.glyph(emaker.bbox(x='0', y='0', w='20', h='20'),
                                  class_('source and sink'), id=glyph_id)
        self._map.append(none_glyph)
        return glyph_id

    def _agent_glyph(self, agent, append=True):
        # Make the main glyph for the agent
        # TODO: handle bound conditions
        agent_id = self._make_agent_id(agent)
        agent_type = _get_agent_type(agent)
        glyph = emaker.glyph(emaker.label(text=agent.name),
                             emaker.bbox(x='0', y='0', w='140', h='60'),
                             class_(agent_type), id=agent_id)

        # Make a glyph for the agent type
        # TODO: handle other agent types
        type_glyph = emaker.glyph(emaker.label(text='mt:prot'),
                                  class_('unit of information'),
                                  emaker.bbox(x='0', y='0', w='53', h='18'),
                                  id=self._make_id())
        glyph.append(type_glyph)

        # Make glyphs for agent state
        # TODO: handle location, mutation
        for m in agent.mods:
            if m.residue is not None:
                mod = m.residue
            else:
                mod = abbrevs[m.mod_type]
            mod_pos = m.position if m.position is not None else ''
            variable = '%s%s' % (mod, mod_pos)
            value = states[m.mod_type][1 if m.is_modified else 0]
            state = emaker.state(variable=variable, value=value)
            state_glyph = \
                emaker.glyph(state, emaker.bbox(x='1', y='1', w='70', h='30'),
                             class_('state variable'), id=self._make_id())
            glyph.append(state_glyph)
        if agent.activity:
            value = 'a' if agent.activity.is_active else 'i'
            state = {'variable': abbrevs[agent.activity.activity_type],
                     'value': value}
            state_glyph = \
                emaker.glyph(state, emaker.bbox(x='1', y='1', w='70', h='30'),
                             class_('state variable'), id=self._make_id())

        # Handle bound conditions as complexes
        if agent.bound_conditions:
            members = [glyph]
            for bc in agent.bound_conditions:
                if bc.is_bound:
                    member_glyph = self._agent_glyph(bc.agent, append=False)
                    members.append(member_glyph)
            # Exclude the case where only negative bound conditions
            # are given and so members has only 1 element
            if len(members) > 1:
                complex_glyph = \
                    emaker.glyph(emaker.bbox(x='1', y='1', w='70', h='30'),
                                 class_('complex'), id=self._make_id())
                for member in members:
                    complex_glyph.append(member)
                glyph = complex_glyph
        if append:
            self._map.append(glyph)
            return agent_id
        return glyph

    def _make_id(self):
        element_id = 'id_%d' % self.id_counter
        self.id_counter += 1
        return element_id

    def _make_agent_id(self, agent):
        key = agent.matches_key()
        mapped_id = self.agent_ids.get(key)
        if mapped_id:
            return mapped_id
        new_id = self._make_id()
        self.agent_ids[key] = new_id
        return new_id

    def statement_exists(self, stmt):
        for s in self.statements:
            if stmt.matches(s):
                return True
        return False

    def add_statements(self, stmts):
        for stmt in stmts:
            if not self.statement_exists(stmt):
                self.statements.append(stmt)

def _get_agent_type(agent):
    if agent.db_refs.get('UP') or agent.db_refs.get('HGNC') or \
        agent.db_refs.get('BE') or agent.db_refs.get('PF'):
        return 'macromolecule'
    elif agent.db_refs.get('CHEBI') or agent.db_refs.get('PUBCHEM'):
        return 'simple chemical'
    elif agent.db_refs.get('GO'):
        return 'phenotype'
    return 'unspecified entity'

def class_(name):
    return {'class': name}

abbrevs = {
    'phosphorylation': 'phospho',
    'ubiquitination': 'ub',
    'farnesylation': 'farnesyl',
    'hydroxylation': 'hydroxyl',
    'acetylation': 'acetyl',
    'sumoylation': 'sumo',
    'glycosylation': 'glycosyl',
    'methylation': 'methyl',
    'modification': 'mod',
    'activity': 'act',
    'kinase': 'kin'
}

states = {
    'phosphorylation': ['u', 'p'],
    'ubiquitination': ['n', 'y'],
    'farnesylation': ['n', 'y'],
    'hydroxylation': ['n', 'y'],
    'acetylation': ['n', 'y'],
    'sumoylation': ['n', 'y'],
    'glycosylation': ['n', 'y'],
    'methylation': ['n', 'y'],
    'modification': ['n', 'y'],
}

