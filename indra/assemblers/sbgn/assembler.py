from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import copy
import logging
import lxml.etree
import lxml.builder
from indra.statements import *
from indra.assemblers.pysb_assembler import PysbPreassembler

logger = logging.getLogger('sbgn_assembler')

sbgn_ns = 'http://sbgn.org/libsbgn/pd/0.1'
emaker = lxml.builder.ElementMaker(nsmap={None: sbgn_ns})


class SBGNAssembler(object):
    """This class assembles an SBGN model from a set of INDRA Statements.

    The Systems Biology Graphical Notation (SBGN) is a widely used
    graphical notation standard for systems biology models.
    This assembler creates SBGN models following the Process Desctiption (PD)
    standard, documented at:
    https://github.com/sbgn/process-descriptions/blob/master/UserManual/sbgn_PD-level1-user-public.pdf.
    For more information on SBGN, see: http://sbgn.github.io/sbgn/

    Parameters
    ----------
    stmts : Optional[list[indra.statements.Statement]]
        A list of INDRA Statements to be assembled.

    Attributes
    ----------
    statements : list[indra.statements.Statement]
        A list of INDRA Statements to be assembled.
    sbgn : lxml.etree.ElementTree
        The structure of the SBGN model that is assembled, represented as an
        XML ElementTree.
    """

    process_style = {'x': '0', 'y': '0', 'w': '10', 'h': '10'}
    source_sink_style = {'x': '0', 'y': '0', 'w': '10', 'h': '10'}
    monomer_style = {'x': '0', 'y': '0', 'w': '60', 'h': '30'}
    complex_style = {'x': '1', 'y': '1', 'w': '60', 'h': '65'}
    entity_type_style = {'x': '0', 'y': '0', 'w': '30', 'h': '12'}
    entity_state_style = {'x': '1', 'y': '1', 'w': '28', 'h': '12'}

    def __init__(self, statements=None):
        if not statements:
            self.statements = []
        else:
            self.statements = statements
        self.sbgn = emaker.sbgn()
        self._map = emaker.map()
        self.sbgn.append(self._map)
        self._id_counter = 0
        self._agent_ids = {}

    def add_statements(self, stmts):
        """Add INDRA Statements to the assembler's list of statements.

        Parameters
        ----------
        stmts : list[indra.statements.Statement]
            A list of :py:class:`indra.statements.Statement`
            to be added to the statement list of the assembler.
        """
        for stmt in stmts:
            if not self.statement_exists(stmt):
                self.statements.append(stmt)

    def make_model(self):
        """Assemble the SBGN model from the collected INDRA Statements.

        This method assembles an SBGN model from the set of INDRA Statements.
        The assembled model is set as the assembler's sbgn attribute (it is
        represented as an XML ElementTree internally). The model is returned
        as a serialized XML string.

        Returns
        -------
        sbgn_str : str
            The XML serialized SBGN model.
        """
        ppa = PysbPreassembler(self.statements)
        ppa.replace_activities()
        self.statements = ppa.statements
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
            elif isinstance(stmt, ActiveForm):
                #self._assemble_activeform(stmt)
                pass
            else:
                logger.warning("Unhandled Statement type %s" % type(stmt))
                continue
        sbgn_str = self.print_model()
        return sbgn_str

    def print_model(self, pretty=True, encoding='utf8'):
        """Return the assembled SBGN model as an XML string.

        Parameters
        ----------
        pretty : Optional[bool]
            If True, the SBGN string is formatted with indentation (for human
            viewing) otherwise no indentation is used. Default: True

        Returns
        -------
        sbgn_str : bytes (str in Python 2)
            An XML string representation of the SBGN model.
        """
        return lxml.etree.tostring(self.sbgn, pretty_print=pretty,
                                   encoding=encoding)

    def save_model(self, file_name='model.sbgn'):
        """Save the assembled SBGN model in a file.

        Parameters
        ----------
        file_name : Optional[str]
            The name of the file to save the SBGN network to.
            Default: model.sbgn
        """
        model = self.print_model()
        with open(file_name, 'wb') as fh:
            fh.write(model)

    def _assemble_modification(self, stmt):
        if not stmt.enz:
            return
        # Make glyph for enz
        enz_glyph = self._agent_glyph(stmt.enz)
        mc_changed = stmt._get_mod_condition()
        mc_unchanged = stmt._get_mod_condition()
        mc_unchanged.is_modified = not mc_unchanged.is_modified
        # Make glyphs for sub
        sub_changed = copy.deepcopy(stmt.sub)
        sub_changed.mods.append(mc_changed)
        sub_unchanged = copy.deepcopy(stmt.sub)
        sub_unchanged.mods.append(mc_unchanged)
        sub_in, sub_out = \
            (sub_unchanged, sub_changed) if isinstance(stmt, AddModification) else \
            (sub_changed, sub_unchanged)
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
        # Make the process glyph
        process_glyph = self._process_glyph('process')
        # Add the arcs
        if isinstance(stmt, DecreaseAmount):
            self._arc('consumption', obj_glyph, process_glyph)
        else:
            self._arc('production', process_glyph, obj_glyph)
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
            if i == 0:
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

    def _assemble_activeform(self, stmt):
        agent_glyph = self._agent_glyph(stmt.agent)
        agent_active = copy.deepcopy(stmt.agent)
        agent_active.activity = ActivityCondition(stmt.activity,
                                                  stmt.is_active)
        agent_active_glyph = self._agent_glyph(agent_active)
        process_glyph = self._process_glyph('process')
        self._arc('consumption', agent_glyph, process_glyph)
        self._arc('production', process_glyph, agent_active_glyph)

    def _arc(self, class_name, source, target):
        arc_id = self._make_id()
        arc = emaker.arc(class_(class_name), source=source, target=target,
                         id=arc_id)
        self._map.append(arc)

    def _process_glyph(self, class_name):
        process_id = self._make_id()
        process_glyph = emaker.glyph(emaker.bbox(**self.process_style),
                                     class_(class_name), id=process_id)
        self._map.append(process_glyph)
        return process_id

    def _none_glyph(self):
        glyph_id = self._make_id()
        none_glyph = emaker.glyph(emaker.bbox(**self.source_sink_style),
                                  class_('source and sink'), id=glyph_id)
        self._map.append(none_glyph)
        return glyph_id

    def _agent_glyph(self, agent, append=True):
        # Make the main glyph for the agent
        # TODO: handle bound conditions
        agent_id = self._make_agent_id(agent)
        agent_type = _get_agent_type(agent)
        glyph = emaker.glyph(emaker.label(text=agent.name),
                             emaker.bbox(**self.monomer_style),
                             class_(agent_type), id=agent_id)

        # Temporarily remove
        # Make a glyph for the agent type
        # TODO: handle other agent types
        #type_glyph = emaker.glyph(emaker.label(text='mt:prot'),
        #                          class_('unit of information'),
        #                          emaker.bbox(**self.entity_type_style),
        #                          id=self._make_id())
        #glyph.append(type_glyph)

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
                emaker.glyph(state, emaker.bbox(**self.entity_state_style),
                             class_('state variable'), id=self._make_id())
            glyph.append(state_glyph)
        if agent.activity:
            value = 'a' if agent.activity.is_active else 'i'
            state = emaker.state(variable=abbrevs[agent.activity.activity_type],
                                 value=value)
            state_glyph = \
                emaker.glyph(state, emaker.bbox(**self.entity_state_style),
                             class_('state variable'), id=self._make_id())
            glyph.append(state_glyph)

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
                    emaker.glyph(emaker.bbox(**self.complex_style),
                                 class_('complex'), id=self._make_id())
                for member in members:
                    complex_glyph.append(member)
                glyph = complex_glyph
        if append:
            self._map.append(glyph)
            return agent_id
        return glyph

    def _glyph_for_complex_pattern(self, pattern):
        """Add glyph and member glyphs for a PySB ComplexPattern."""
        # Make the main glyph for the agent
        monomer_glyphs = []
        for monomer_pattern in pattern.monomer_patterns:
            glyph = self._glyph_for_monomer_pattern(monomer_pattern)
            monomer_glyphs.append(glyph)

        if len(monomer_glyphs) > 1:
            pattern.matches_key = lambda: str(pattern)
            agent_id = self._make_agent_id(pattern)
            complex_glyph = \
                emaker.glyph(emaker.bbox(**self.complex_style),
                             class_('complex'), id=agent_id)
            for glyph in monomer_glyphs:
                glyph.attrib['id'] = agent_id + glyph.attrib['id']
                complex_glyph.append(glyph)
            return complex_glyph
        return monomer_glyphs[0]

    def _glyph_for_monomer_pattern(self, pattern):
        """Add glyph for a PySB MonomerPattern."""
        pattern.matches_key = lambda: str(pattern)
        agent_id = self._make_agent_id(pattern)
        # Handle sources and sinks
        if pattern.monomer.name in ('__source', '__sink'):
            return None
        # Handle molecules
        glyph = emaker.glyph(emaker.label(text=pattern.monomer.name),
                             emaker.bbox(**self.monomer_style),
                             class_('macromolecule'), id=agent_id)
        # Temporarily remove this
        # Add a glyph for type
        #type_glyph = emaker.glyph(emaker.label(text='mt:prot'),
        #                          class_('unit of information'),
        #                          emaker.bbox(**self.entity_type_style),
        #                          id=self._make_id())
        #glyph.append(type_glyph)
        for site, value in pattern.site_conditions.items():
            if value is None or isinstance(value, int):
                continue
            # Make some common abbreviations
            if site == 'phospho':
                site = 'p'
            elif site == 'activity':
                site = 'act'
                if value == 'active':
                    value = 'a'
                elif value == 'inactive':
                    value = 'i'
            state = emaker.state(variable=site, value=value)
            state_glyph = \
                emaker.glyph(state, emaker.bbox(**self.entity_state_style),
                             class_('state variable'), id=self._make_id())
            glyph.append(state_glyph)
        return glyph

    def _make_id(self):
        element_id = 'id_%d' % self._id_counter
        self._id_counter += 1
        return element_id

    def _make_agent_id(self, agent):
        key = agent.matches_key()
        mapped_id = self._agent_ids.get(key)
        if mapped_id:
            return mapped_id
        new_id = self._make_id()
        self._agent_ids[key] = new_id
        return new_id

    def statement_exists(self, stmt):
        for s in self.statements:
            if stmt.matches(s):
                return True
        return False


def _get_agent_type(agent):
    if agent.db_refs.get('UP') or agent.db_refs.get('HGNC') or \
        agent.db_refs.get('FPLX') or agent.db_refs.get('PF'):
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
