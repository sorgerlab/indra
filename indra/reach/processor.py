from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import re
import logging
import objectpath
from indra.statements import *
from indra.util import read_unicode_csv
from indra.databases import hgnc_client
import indra.databases.uniprot_client as up_client

logger = logging.getLogger('reach')

class ReachProcessor(object):
    """The ReachProcessor extracts INDRA Statements from REACH parser output.

    Parameters
    ----------
    json_dict : dict
        A JSON dictionary containing the REACH extractions.
    pmid : Optional[str]
        The PubMed ID associated with the extractions. This can be passed
        in case the PMID cannot be determined from the extractions alone.`

    Attributes
    ----------
    tree : objectpath.Tree
        The objectpath Tree object representing the extractions.
    statements : list[indra.statements.Statement]
        A list of INDRA Statements that were extracted by the processor.
    citation : str
        The PubMed ID associated with the extractions.
    all_events : dict[str, str]
        The frame IDs of all events by type in the REACH extraction.
    """
    def __init__(self, json_dict, pmid=None):
        self.tree = objectpath.Tree(json_dict)
        self.statements = []
        self.citation = pmid
        if pmid is None:
            if self.tree is not None:
                self.citation =\
                    self.tree.execute("$.events.object_meta.doc_id")
        self.get_all_events()

    def print_event_statistics(self):
        """Print the number of events in the REACH output by type."""
        logger.info('All events by type')
        logger.info('-------------------')
        for k, v in self.all_events.items():
            logger.info('%s, %s' % (k, len(v)))
        logger.info('-------------------')

    def get_all_events(self):
        """Gather all event IDs in the REACH output by type.

        These IDs are stored in the self.all_events dict.
        """
        self.all_events = {}
        events = self.tree.execute("$.events.frames")
        if events is None:
            return
        for e in events:
            event_type = e.get('type')
            frame_id = e.get('frame_id')
            try:
                self.all_events[event_type].append(frame_id)
            except KeyError:
                self.all_events[event_type] = [frame_id]

    def print_regulations(self):
        qstr = "$.events.frames[(@.type is 'regulation')]"
        res = self.tree.execute(qstr)
        if res is None:
            return
        for r in res:
            print(r['subtype'])
            for a in r['arguments']:
                print(a['type'], '/', a['argument-type'], ':', a['text'])

    def get_modifications(self):
        """Extract Modification INDRA Statements."""
        qstr = "$.events.frames[(@.type is 'protein-modification')]"
        res = self.tree.execute(qstr)
        if res is None:
            return
        for r in res:
            modification_type = r.get('subtype')
            epistemics = self._get_epistemics(r)
            if epistemics.get('negative'):
                continue
            context = self._get_context(r)
            frame_id = r['frame_id']
            args = r['arguments']
            site = None
            theme = None

            for a in args:
                if self._get_arg_type(a) == 'theme':
                    theme = a['arg']
                elif self._get_arg_type(a) == 'site':
                    site = a['text']
            theme_agent = self._get_agent_from_entity(theme)
            if site is not None:
                residue, pos = self._parse_site_text(site)
            else:
                residue = None
                pos = None
            qstr = "$.events.frames[(@.type is 'regulation') and " + \
                   "(@.arguments[0].arg is '%s')]" % frame_id
            reg_res = self.tree.execute(qstr)
            reg_res = list(reg_res)
            for reg in reg_res:
                controller_agent = None
                for a in reg['arguments']:
                    if self._get_arg_type(a) == 'controller':
                        controller = a.get('arg')
                        if controller is not None:
                            controller_agent = \
                                self._get_agent_from_entity(controller)
                            break

                sentence = reg['verbose-text']
                ev = Evidence(source_api='reach', text=sentence,
                              annotations=context, pmid=self.citation,
                              epistemics=epistemics)
                args = [controller_agent, theme_agent, residue, pos, ev]

                # Here ModStmt is a sub-class of Modification
                ModStmt = modtype_to_modclass.get(modification_type)
                if ModStmt is None:
                    logger.warning('Unhandled modification type: %s' %
                                   modification_type)
                else:
                    # Handle this special case here because only
                    # enzyme argument is needed
                    if modification_type == 'autophosphorylation':
                        args = [theme_agent, residue, pos, ev]
                    self.statements.append(ModStmt(*args))

    def get_regulate_amounts(self):
        """Extract RegulateAmount INDRA Statements."""
        qstr = "$.events.frames[(@.type is 'transcription')]"
        res = self.tree.execute(qstr)
        all_res = []
        if res is not None:
            all_res += list(res)
        qstr = "$.events.frames[(@.type is 'amount')]"
        res = self.tree.execute(qstr)
        if res is not None:
            all_res += list(res)

        for r in all_res:
            subtype = r.get('subtype')
            epistemics = self._get_epistemics(r)
            if epistemics.get('negative'):
                continue
            context = self._get_context(r)
            frame_id = r['frame_id']
            args = r['arguments']
            theme = None
            for a in args:
                if self._get_arg_type(a) == 'theme':
                    theme = a['arg']
                    break
            if theme is None:
                continue
            theme_agent = self._get_agent_from_entity(theme)
            qstr = "$.events.frames[(@.type is 'regulation') and " + \
                   "(@.arguments[0].arg is '%s')]" % frame_id
            reg_res = self.tree.execute(qstr)
            for reg in reg_res:
                controller_agent = None
                for a in reg['arguments']:
                    if self._get_arg_type(a) == 'controller':
                        controller = a.get('arg')
                        if controller is not None:
                            controller_agent = \
                                    self._get_agent_from_entity(controller)
                            break
                sentence = reg['verbose-text']

                ev = Evidence(source_api='reach', text=sentence,
                              annotations=context, pmid=self.citation,
                              epistemics=epistemics)
                args = [controller_agent, theme_agent, ev]
                subtype = reg.get('subtype')
                if subtype == 'positive-regulation':
                    st = IncreaseAmount(*args)
                else:
                    st = DecreaseAmount(*args)
                self.statements.append(st)


    def get_complexes(self):
        """Extract INDRA Complex Statements."""
        qstr = "$.events.frames[@.type is 'complex-assembly']"
        res = self.tree.execute(qstr)
        if res is None:
            return
        for r in res:
            epistemics = self._get_epistemics(r)
            if epistemics.get('negative'):
                continue
            context = self._get_context(r)
            args = r['arguments']
            sentence = r['verbose-text']
            members = []
            for a in args:
                agent = self._get_agent_from_entity(a['arg'])
                members.append(agent)
            ev = Evidence(source_api='reach', text=sentence,
                          annotations=context, pmid=self.citation,
                          epistemics=epistemics)
            self.statements.append(Complex(members, ev))

    def get_activation(self):
        """Extract INDRA Activation Statements."""
        qstr = "$.events.frames[@.type is 'activation']"
        res = self.tree.execute(qstr)
        if res is None:
            return
        for r in res:
            epistemics = self._get_epistemics(r)
            if epistemics.get('negative'):
                continue
            sentence = r['verbose-text']
            context = self._get_context(r)
            ev = Evidence(source_api='reach', text=sentence,
                          pmid=self.citation, annotations=context,
                          epistemics=epistemics)
            args = r['arguments']
            for a in args:
                if self._get_arg_type(a) == 'controller':
                    controller = a.get('arg')
                    # When the controller is not a simple entity
                    if controller is None:
                        if a['argument-type'] == 'complex':
                            controllers = list(a.get('args').values())
                            controller_agent =\
                                self._get_agent_from_entity(controllers[0])
                            bound_agents = [self._get_agent_from_entity(c) 
                                            for c in controllers[1:]]
                            bound_conditions = [BoundCondition(ba, True) for
                                                ba in bound_agents]
                            controller_agent.bound_conditions = \
                                    bound_conditions
                    else:
                        controller_agent =\
                            self._get_agent_from_entity(controller)
                if self._get_arg_type(a) == 'controlled':
                    controlled = a['arg']
            controlled_agent = self._get_agent_from_entity(controlled)
            if r['subtype'] == 'positive-activation':
                st = Activation(controller_agent, controlled_agent,
                                evidence=ev)
            else:
                st = Inhibition(controller_agent, controlled_agent,
                                evidence=ev)
            self.statements.append(st)

    def get_translocation(self):
        """Extract INDRA Translocation Statements."""
        qstr = "$.events.frames[@.type is 'translocation']"
        res = self.tree.execute(qstr)
        if res is None:
            return
        for r in res:
            epistemics = self._get_epistemics(r)
            if epistemics.get('negative'):
                continue
            sentence = r['verbose-text']
            context = self._get_context(r)
            ev = Evidence(source_api='reach', text=sentence,
                          pmid=self.citation, annotations=context,
                          epistemics=epistemics)
            args = r['arguments']
            from_location = None
            to_location = None
            for a in args:
                if self._get_arg_type(a) == 'theme':
                    agent = self._get_agent_from_entity(a['arg'])
                    if agent is None:
                        continue
                elif self._get_arg_type(a) == 'source':
                    from_location = self._get_location_by_id(a['arg'])
                elif self._get_arg_type(a) == 'destination':
                    to_location = self._get_location_by_id(a['arg'])
            st = Translocation(agent, from_location, to_location,
                               evidence=ev)
            self.statements.append(st)

    def _get_location_by_id(self, loc_id):
        qstr = "$.entities.frames[(@.frame_id is \'%s\')]" % loc_id
        res = self.tree.execute(qstr)
        if res is None:
            return None
        try:
            entity_term = next(res)
        except StopIteration:
            logger.debug(' %s is not an entity' % entity_id)
            return None
        name = entity_term.get('text')
        go_id = None
        for xr in entity_term['xrefs']:
            ns = xr['namespace']
            if ns == 'go':
                go_id = xr['id']
        # Try to get valid location based on GO id
        if go_id is not None:
            try:
                loc = get_valid_location(go_id)
                return loc
            except InvalidLocationError:
                pass
        # See if the raw name is a valid cellular component
        try:
            loc = get_valid_location(name.lower())
            return loc
        except InvalidLocationError:
            pass
        return None

    def _get_agent_from_entity(self, entity_id):
        qstr = "$.entities.frames[(@.frame_id is \'%s\')]" % entity_id
        res = self.tree.execute(qstr)
        if res is None:
            return None
        try:
            entity_term = next(res)
        except StopIteration:
            logger.debug(' %s is not an entity' % entity_id)
            return None
        # This is the default name, which can be overwritten 
        # below for specific database entries
        agent_name = self._get_valid_name(entity_term['text'])
        db_refs = {}
        for xr in entity_term['xrefs']:
            ns = xr['namespace']
            if ns == 'uniprot':
                up_id = xr['id']
                db_refs['UP'] = up_id
                # Look up official names in UniProt
                gene_name = up_client.get_gene_name(up_id)
                if gene_name is not None:
                    agent_name = self._get_valid_name(gene_name)
                    # If the gene name corresponds to an HGNC ID, add it to the
                    # db_refs
                    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                    if hgnc_id:
                        db_refs['HGNC'] = hgnc_id
            elif ns == 'hgnc':
                hgnc_id = xr['id']
                db_refs['HGNC'] = hgnc_id
                # Look up the standard gene symbol and set as name
                hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
                if hgnc_name:
                    agent_name = hgnc_name
                # Look up the corresponding uniprot id
                up_id = hgnc_client.get_uniprot_id(hgnc_id)
                if up_id:
                    db_refs['UP'] = up_id
            elif ns == 'pfam':
                be_id = bioentities_map.get(('PF', xr['id']))
                if be_id:
                    db_refs['BE'] = be_id
                    agent_name = be_id
                db_refs['PF'] = xr['id']
            elif ns == 'interpro':
                be_id = bioentities_map.get(('IP', xr['id']))
                if be_id:
                    db_refs['BE'] = be_id
                    agent_name = be_id
                db_refs['PF'] = xr['id']
            elif ns == 'chebi':
                db_refs['CHEBI'] = xr['id']
            elif ns == 'pubchem':
                db_refs['PUBCHEM'] = 'PUBCHEM:%s' % xr['id']
            elif ns == 'go':
                db_refs['GO'] = xr['id']
            elif ns == 'mesh':
                db_refs['MESH'] = xr['id']
            elif ns == 'hmdb':
                db_refs['HMDB'] = xr['id']
            elif ns == 'simple_chemical':
                if xr['id'].startswith('HMDB'):
                    db_refs['HMDB'] = xr['id']
            elif ns == 'be':
                db_refs['BE'] = xr['id']
                agent_name = db_refs['BE']
            # These name spaces are ignored
            elif ns in ['uaz']:
                pass
            else:
                logger.warning('Unhandled xref namespace: %s' % ns)
        db_refs['TEXT'] = entity_term['text']

        mod_terms = entity_term.get('modifications')
        mods = []
        muts = []
        if mod_terms is not None:
            for m in mod_terms:
                if m['type'].lower() == 'mutation':
                    # Evidence is usualy something like "V600E"
                    # We could parse this to get the amino acid
                    # change that happened.
                    mutation_str = m.get('evidence')
                    # TODO: sometimes mutation_str is "mutant", "Mutant",
                    # "mutants" - this indicates that there is a mutation
                    # but not the specific type. We should encode this
                    # somehow as a "blank" mutation condition
                    mut = self._parse_mutation(mutation_str)
                    if mut is not None:
                        muts.append(mut)
                else:
                    mc = self._get_mod_condition(m)
                    if mc is not None:
                        mods.append(mc)

        agent = Agent(agent_name, db_refs=db_refs, mods=mods, mutations=muts)
        return agent

    def _get_mod_condition(self, mod_term):
        site = mod_term.get('site')
        if site is not None:
            mod_res, mod_pos = self._parse_site_text(site)
        else:
            mod_res = None
            mod_pos = None
        mod_type_str = mod_term['type'].lower()
        mod_state = agent_mod_map.get(mod_type_str)
        if mod_state is not None:
            mc = ModCondition(mod_state[0], residue=mod_res, position=mod_pos,
                              is_modified=mod_state[1])
            return mc
        logger.warning('Unhandled entity modification type: %s' % mod_type_str)
        return None

    def _get_context(self, frame_term):
        context = {}
        context['found_by'] = frame_term['found_by']
        try:
            context_id = frame_term['context']
        except KeyError:
            return context
        # For backwards compatibility with older versions
        # of REACH
        if isinstance(context_id, dict):
            context_term = context_id
            species = context_term.get('Species')
            cell_type = context_term.get('CellType')
            cell_line = None
            location = None
            tissue = None
            organ = None
        else:
            qstr = "$.entities.frames[(@.frame_id is \'%s\')]" % context_id[0]
            res = self.tree.execute(qstr)
            if res is None:
                return context
            context_frame = next(res)
            facets = context_frame['facets']
            cell_line = facets.get('cell-line')
            cell_type = facets.get('cell-type')
            species = facets.get('organism')
            location = facets.get('location')
            tissue = facets.get('tissue_type')
            organ = facets.get('organ')
        context['species'] = species
        context['cell_type'] = cell_type
        context['cell_line'] = cell_line
        context['location'] = location
        context['tissue'] = tissue
        context['organ'] = organ
        return context

    def _get_epistemics(self, event):
        epistemics = {}
        # Check whether information is negative
        neg = event.get('is_negated')
        if neg is True:
            epistemics['negative'] = True
        # Check if it is a hypothesis
        hyp = event.get('is_hypothesis')
        if hyp is True:
            epistemics['hypothesis'] = True
        # Check if it is direct
        if 'is_direct' in event:
            direct = event['is_direct']
            epistemics['direct'] = direct
        # Get the section of the paper it comes from
        section = self._get_section(event)
        epistemics['section_type'] = section
        return epistemics

    _section_list = ['title', 'abstract', 'introduction', 'background',
                     'results', 'methods', 'discussion', 'conclusion',
                     'supplementary', 'figure']

    def _get_section(self, event):
        """Get the section of the paper that the event is from."""
        sentence_id = event.get('sentence')
        section = None
        if sentence_id:
            qstr = "$.sentences.frames[(@.frame_id is \'%s\')]" % sentence_id
            res = self.tree.execute(qstr)
            if res:
                sentence_frame = list(res)[0]
                passage_id = sentence_frame.get('passage')
                if passage_id:
                    qstr = "$.sentences.frames[(@.frame_id is \'%s\')]" % \
                            passage_id
                    res = self.tree.execute(qstr)
                    if res:
                        passage_frame = list(res)[0]
                        section = passage_frame.get('section-id')
        # If the section is in the standard list, return as is
        if section in self._section_list:
            return section
        # Next, handle a few special cases that come up in practice
        elif section.startswith('fig'):
            return 'figure'
        elif section.startswith('supm'):
            return 'supplementary'
        elif section == 'article-title':
            return 'title'
        elif section in ['subjects|methods', 'methods|subjects']:
            return 'methods'
        elif section == 'conclusions':
            return 'conclusion'
        elif section == 'intro':
            return 'introduction'
        else:
            return None

    @staticmethod
    def _get_arg_type(arg):
        """Return the type of the argument with backwards compatibility."""
        if arg.get('argument_label') is not None:
            return arg.get('argument_label')
        else:
            return arg.get('type')

    @staticmethod
    def _get_valid_name(txt):
        """Produce valid agent name from string."""
        name = ''.join(ch if ch.isalnum() else '_' for ch in txt)
        if name and name[0].isdigit():
            name = 'p' + name
        return name

    @staticmethod
    def _parse_mutation(s):
        m = re.match(r'([A-Z])([0-9]+)([A-Z])', s)
        if m is not None:
            parts = [str(g) for g in m.groups()]
            residue_from = get_valid_residue(parts[0])
            residue_to = get_valid_residue(parts[2])
            position = parts[1]
            mut = MutCondition(position, residue_from, residue_to)
            return mut
        return None

    @staticmethod
    def _parse_site_text(s):
        for p in (_site_pattern1, _site_pattern2, _site_pattern3):
            m = re.match(p, s.upper())
            if m is not None:
                residue = get_valid_residue(m.groups()[0])
                site = m.groups()[1]
                return residue, site
        m = re.match(_site_pattern4, s.upper())
        if m is not None:
            site = m.groups()[0]
            residue = m.groups()[1]
            return residue, site
        for p in (_site_pattern5, _site_pattern6, _site_pattern7):
            m = re.match(p, s.upper())
            if m is not None:
                residue = get_valid_residue(m.groups()[0])
                site = None
                return residue, site
        m = re.match(_site_pattern8, s.upper())
        if m is not None:
            site = m.groups()[0]
            residue = None
            return residue, site
        logger.warning('Could not parse site text %s' % s)
        return None, None

_site_pattern1 = '([' + ''.join(list(amino_acids.keys())) + '])[-]?([0-9]+)$'
_site_pattern2 = '(' + '|'.join([v['short_name'].upper() for
                                 v in amino_acids.values()]) + \
                        ')[- ]?([0-9]+)$'
_site_pattern3 = '(' + '|'.join([v['indra_name'].upper() for
                                 v in amino_acids.values()]) + \
                        ')[^0-9]*([0-9]+)$'
_site_pattern4 = '([0-9]+)([' + ''.join(list(amino_acids.keys())) + '])$'
_site_pattern5 = '^([' + ''.join(list(amino_acids.keys())) + '])$'
_site_pattern6 = '^(' + '|'.join([v['short_name'].upper() for
                                 v in amino_acids.values()]) + ')$'
_site_pattern7 = '.*(' + '|'.join([v['indra_name'].upper() for
                                 v in amino_acids.values()]) + ').*'
_site_pattern8 = '([0-9]+)$'

# Subtypes that exist but we don't handle: methylation, hydrolysis
agent_mod_map = {
    'phosphorylation': ('phosphorylation', True),
    'phosphorylated': ('phosphorylation', True),
    'dephosphorylation': ('phosphorylation', False),
    'acetylation': ('acetylation', True),
    'deacetylation': ('acetylation', False),
    'ubiquitination': ('ubiquitination', True),
    'deubiquitination': ('ubiquitination', False),
    'hydroxylation': ('hydroxylation', True),
    'dehydroxylation': ('hydroxylation', False),
    'sumoylation': ('sumoylation', True),
    'desumoylation': ('sumoylation', False),
    'glycosylation': ('glycosylation', True),
    'deglycosylation': ('glycosylation', False),
    'farnesylation': ('farnesylation', True),
    'defarnesylation': ('farnesylation', False),
    'ribosylation': ('ribosylation', True),
    'deribosylation': ('ribosylation', False),
    'methylation': ('methylation', True),
    'demethylation': ('methylation', False),
    'unknown': ('modification', True),
}

def _read_bioentities_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../resources/bioentities_map.tsv')
    bioentities_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        source_ns = row[0]
        source_id = row[1]
        be_id = row[2]
        bioentities_map[(source_ns, source_id)] = be_id
    return bioentities_map

bioentities_map = _read_bioentities_map()
