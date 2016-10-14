from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import itertools
import copy
import logging
import collections
import lxml.builder
import lxml.etree
from indra.trips import trips_api
from indra import statements as ist

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
    'activity': ['i', 'a'],
    'kinase': ['i', 'a']
}

class SBGNAssembler(object):

    def __init__(self, policies=None):
        self.statements = []
        self.agent_set = None

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
            uppercase_agents(stmt)
            if not self.statement_exists(stmt):
                self.statements.append(stmt)

    def make_model(self):

        def make_id(_counter=[0]):
            id_ = 'id_%d' % _counter[0]
            _counter[0] += 1
            return id_

        def class_(name):
            return {'class': name}

        def glyph_for_monomer(agent, in_complex=False):
            if in_complex:
                agent_id = make_id()
            else:
                agent_id = agent_ids[agent.matches_key()]
            glyph = E.glyph(
                E.label(text=agent.name),
                E.bbox(x='0', y='0', w='140', h='60'),
                class_('macromolecule'), id=agent_id,
                )
            glyph.append(
                E.glyph(
                    E.label(text='mt:prot'),
                    class_('unit of information'),
                    E.bbox(x='0', y='0', w='53', h='18'),
                    id=make_id())
                )
            for st in sbgn_states_for_agent(agent):
                glyph.append(
                    E.glyph(
                        E.state(**st._asdict()),
                        E.bbox(x='1', y='1', w='70', h='30'),
                        class_('state variable'), id=make_id(),
                        )
                    )
            return glyph

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
            if isinstance(s, ist.Modification) or \
                isinstance(s, ist.Activation) or \
                isinstance(s, ist.ActiveForm):
                class_name = 'process'
            elif isinstance(s, ist.Complex):
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
            pg_id = make_id()
            process_glyph = E.glyph(E.bbox(x='0', y='0', w='20', h='20'),
                                    class_(class_name), id=pg_id)
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
            if isinstance(s, ist.Modification):
                map.append(
                    E.arc(class_('catalysis'),
                          source=agent_ids[s.enz.matches_key()],
                          target=pg_id,
                          id=make_id(),
                          )
                    )
            if isinstance(s, ist.Activation):
                map.append(
                    E.arc(class_('catalysis'),
                          source=agent_ids[s.subj.matches_key()],
                          target=pg_id,
                          id=make_id(),
                          )
                    )
        return lxml.etree.tostring(root, pretty_print=True)

SBGNState = collections.namedtuple('SBGNState', 'variable value')

def sbgn_states_for_agent(agent):
    agent_states = []
    for m in agent.mods:
        if m.residue is not None:
            mod = m.residue
        else:
            mod = abbrevs[m.mod_type]
        mod_pos = m.position if m.position is not None else ''
        variable = '%s%s' % (mod, mod_pos)
        value = states[m.mod_type][1].upper()
        agent_states.append(SBGNState(variable, value))
    return agent_states

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
    if isinstance(stmt, ist.Phosphorylation):
        product = copy.deepcopy(stmt.sub)
        mc = ist.ModCondition('phosphorylation', stmt.residue, stmt.position)
        product.mods.append(mc)
    elif isinstance(stmt, ist.Dephosphorylation):
        product = copy.deepcopy(stmt.sub)
        mc = ist.ModCondition('phosphorylation', stmt.residue, stmt.position)
        product = remove_agent_mod(product, mc)
    elif isinstance(stmt, ist.Ubiquitination):
        product = copy.deepcopy(stmt.sub)
        mc = ist.ModCondition('ubiquitination', stmt.residue, stmt.position)
        product.mods.append(mc)
    elif isinstance(stmt, ist.Deubiquitination):
        product = copy.deepcopy(stmt.sub)
        mc = ist.ModCondition('ubiquitination', stmt.residue, stmt.position)
        product = remove_agent_mod(product, mc)
    elif isinstance(stmt, ist.Acetylation):
        product = copy.deepcopy(stmt.sub)
        mc = ist.ModCondition('acetylation', stmt.residue, stmt.position)
        product.mods.append(mc)
    elif isinstance(stmt, ist.Deacetylation):
        product = copy.deepcopy(stmt.sub)
        mc = ist.ModCondition('acetylation', stmt.residue, stmt.position)
        product = remove_agent_mod(product, mc)
    elif isinstance(stmt, ist.Complex):
        product = copy.deepcopy(stmt.members[0])
        for member in stmt.members[1:]:
            bc = ist.BoundCondition(member, True)
            product.bound_conditions.append(bc)
    elif isinstance(stmt, ist.Activation):
        product = copy.deepcopy(stmt.obj)
        if stmt.is_activation:
            mc = ist.ModCondition(stmt.obj_activity)
            product.mods.append(mc)
    elif isinstance(stmt, ist.ActiveForm):
        product = copy.deepcopy(stmt.agent)
        mc = ist.ModCondition(stmt.activity)
        product.mods.append(mc)
    else:
        logger.warning("WARNING: skipping %s" % type(stmt))
        product = None
    return product

def statement_consumed(stmt):
    if isinstance(stmt, ist.Phosphorylation):
        consumed = [copy.deepcopy(stmt.sub)]
    elif isinstance(stmt, ist.Ubiquitination):
        consumed = [copy.deepcopy(stmt.sub)]
    elif isinstance(stmt, ist.Acetylation):
        consumed = [copy.deepcopy(stmt.sub)]
    elif isinstance(stmt, ist.Dephosphorylation):
        stmt_mc = ist.ModCondition('phosphorylation',
                                   stmt.residue, stmt.position)
        consumed1 = copy.deepcopy(stmt.sub)
        consumed1.mods.append(stmt_mc)
        consumed = [consumed1]
    elif isinstance(stmt, ist.Complex):
        consumed = stmt.members
    elif isinstance(stmt, ist.Activation):
        consumed1 = copy.deepcopy(stmt.obj)
        if not stmt.is_activation:
            mc = ist.ModCondition(stmt.obj_activity)
            consumed1.mods.append(mc)
        consumed = [consumed1]
    elif isinstance(stmt, ist.ActiveForm):
        consumed = [copy.deepcopy(stmt.agent)]
    else:
        logger.warning("WARNING: skipping %s" % type(stmt))
        consumed = None
    return consumed

def distinct_agents(agents):
    agents = sorted(agents, key=ist.Agent.matches_key)
    gb = itertools.groupby(agents, ist.Agent.matches_key)
    distinct = [next(g[1]) for g in gb]
    return distinct

def complex_components(agent):
    agent_copy = copy.copy(agent)
    agent_copy.bound_conditions = []
    agents = [agent_copy]
    for bc in agent.bound_conditions:
        agents += complex_components(bc.agent)
    return agents

def uppercase_agents(statement):

    def uppercase_single(agent):
        agent.name = agent.name.upper()
        for bc in agent.bound_conditions:
            uppercase_single(bc.agent)

    for agent in statement.agent_list():
        uppercase_single(agent)

def text_to_sbgn(text=None, trips_xml=None):
    if ((text is None and trips_xml is None) or
        (text is not None and trips_xml is not None)):
        raise ValueError("Must provide ONE of 'text' or 'trips_xml'")
    elif text is not None:
        tp = trips_api.process_text(text)
    elif trips_xml is not None:
        tp = trips_api.process_xml(trips_xml)
    else:
        raise RuntimeError("Unexpected or impossible combination of arguments")
    sa = SBGNAssembler()
    sa.add_statements(tp.statements)
    sbgn_output = sa.make_model()
    return sbgn_output
