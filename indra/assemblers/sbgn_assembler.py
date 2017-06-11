from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import itertools
import copy
import logging
import collections
import lxml.builder
import lxml.etree
from indra.trips import trips_api
from indra.statements import *

logger = logging.getLogger('sbgn_assembler')

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

    def _agent_glyph(self, agent):
        # Make the main glyph for the agent
        # TODO: handle other agent types
        # TODO: use a dict to map matches keys to simple ids
        # TODO: handle bound conditions
        agent_id = agent.matches_key()
        glyph = emaker.glyph(
            emaker.label(text=agent.name),
            emaker.bbox(x='0', y='0', w='140', h='60'),
            class_('macromolecule'), id=agent_id,
            )

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
        self._map.append(glyph)
        return agent_id

        def glyph_for_complex(agent):
            glyph = E.glyph(
                E.bbox(x='0', y='0', w='120', h='60'),
                class_('complex'), id=agent_ids[agent.matches_key()],
                )
            for component in complex_components(agent):
               glyph.append(glyph_for_monomer(component, in_complex=True))
            return glyph

        E = lxml.builder.ElementMaker(nsmap={None: 'http://sbgn.org/libsbgn/pd/0.1'})
        root = E.sbgn()
        map = E.map()
        root.append(map)
        base = agents_for_statements(self.statements)
        transformed = transformed_agents(self.statements)
        agents = distinct_agents(base + transformed)
        agent_ids = {a.matches_key(): make_id() for a in agents}
        for a in agents:
            if not a.bound_conditions:
                glyph = glyph_for_monomer(a)
            else:
                glyph = glyph_for_complex(a)
            map.append(glyph)
        for s in self.statements:
            if isinstance(s, Modification) or \
                isinstance(s, RegulateActivity) or \
                isinstance(s, RegulateAmount) or \
                isinstance(s, ActiveForm):
                class_name = 'process'
            elif isinstance(s, Complex):
                class_name = 'association'
            else:
                logger.warning("WARNING: skipping %s" % type(s))
                continue
            consumed = statement_consumed(s)
            st_prod = statement_product(s)
            if st_prod is None:
                logger.warning("WARNING: skipping %s" % type(s))
                continue
            produced = [st_prod]
            map.append(process_glyph)
            for c in consumed:
                map.append(
                    E.arc(class_('consumption'),
                          source=agent_ids[c.matches_key()],
                          target=pg_id,
                          id=make_id(),
                          )
                    )
            for p in produced:
                map.append(
                    E.arc(class_('production'),
                          source=pg_id,
                          target=agent_ids[p.matches_key()],
                          id=make_id(),
                          )
                    )
            if isinstance(s, Modification):
                map.append(
                    E.arc(class_('catalysis'),
                          source=agent_ids[s.enz.matches_key()],
                          target=pg_id,
                          id=make_id(),
                          )
                    )
            if isinstance(s, RegulateActivity):
                map.append(
                    E.arc(class_('catalysis'),
                          source=agent_ids[s.subj.matches_key()],
                          target=pg_id,
                          id=make_id(),
                          )
                    )

    def _make_id(self):
        element_id = 'id_%d' % self.id_counter
        self.id_counter += 1
        return element_id

    def statement_exists(self, stmt):
        for s in self.statements:
            if stmt.matches(s):
                return True
        return False

    def add_statements(self, stmts):
        for stmt in stmts:
            if any([a is None for a in stmt.agent_list()]):
                continue
            stmt = copy.deepcopy(stmt)
            if not self.statement_exists(stmt):
                self.statements.append(stmt)



def agents_for_statements(statements):
    return [a for stmt in statements for a in stmt.agent_list()]

def transformed_agents(statements):
    # Following filter not needed once all statement types are implemented.
    agents = []
    for s in statements:
        cs = statement_consumed(s)
        if cs is not None:
            agents += cs
    agents += [statement_product(s) for s in statements]
    return [a for a in agents if a is not None]

def remove_agent_mod(agent, mc):
    agent.mods = []
    for mod in agent.mods:
        if not mod.matches(mc):
            agent.mods.append(mc)
    return agent

def statement_product(stmt):
    if isinstance(stmt, Modification):
        modtype = modclass_to_modtype[stmt.__class__]
        if isinstance(stmt, RemoveModification):
            modtype = modtype_to_inverse[modtype]
        mc = ModCondition(modtype, stmt.residue, stmt.position)
        product = copy.deepcopy(stmt.sub)
        product.mods.append(mc)
    elif isinstance(stmt, Complex):
        product = copy.deepcopy(stmt.members[0])
        for member in stmt.members[1:]:
            bc = BoundCondition(member, True)
            product.bound_conditions.append(bc)
    elif isinstance(stmt, RegulateActivity):
        product = copy.deepcopy(stmt.obj)
        if stmt.is_activation:
            mc = ModCondition(stmt.obj_activity)
            product.mods.append(mc)
    elif isinstance(stmt, IncreaseAmount):
        product = copy.deepcopy(stmt.obj)
    elif isinstance(stmt, ActiveForm):
        product = copy.deepcopy(stmt.agent)
        mc = ModCondition(stmt.activity)
        product.mods.append(mc)
    else:
        logger.warning("WARNING: skipping %s" % type(stmt))
        product = None
    return product

def statement_consumed(stmt):
    if isinstance(stmt, AddModification):
        consumed = [copy.deepcopy(stmt.sub)]
    elif isinstance(stmt, RemoveModification):
        modtype = modclass_to_modtype[stmt.__class__]
        if isinstance(stmt, RemoveModification):
            modtype = modtype_to_inverse[modtype]
        mc = ModCondition(modtype, stmt.residue, stmt.position)
        consumed1 = copy.deepcopy(stmt.sub)
        consumed1.mods.append(mc)
        consumed = [consumed1]
    elif isinstance(stmt, Complex):
        consumed = stmt.members
    elif isinstance(stmt, RegulateActivity):
        consumed1 = copy.deepcopy(stmt.obj)
        if not stmt.is_activation:
            mc = ModCondition(stmt.obj_activity)
            consumed1.mods.append(mc)
        consumed = [consumed1]
    elif isinstance(stmt, DecreaseAmount):
        consumed = copy.deepcopy(stmt.obj)
    elif isinstance(stmt, ActiveForm):
        consumed = [copy.deepcopy(stmt.agent)]
    else:
        logger.warning("WARNING: skipping %s" % type(stmt))
        consumed = None
    return consumed

def distinct_agents(agents):
    agents = sorted(agents, key=Agent.matches_key)
    gb = itertools.groupby(agents, Agent.matches_key)
    distinct = [next(g[1]) for g in gb]
    return distinct

def complex_components(agent):
    agent_copy = copy.copy(agent)
    agent_copy.bound_conditions = []
    agents = [agent_copy]
    for bc in agent.bound_conditions:
        agents += complex_components(bc.agent)
    return agents

def class_(name):
    return {'class': name}
