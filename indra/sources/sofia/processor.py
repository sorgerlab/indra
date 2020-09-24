import itertools
from indra.statements import Influence, Concept, Event, Evidence, \
    WorldContext, TimeContext, RefContext, QualitativeDelta


pos_rels = ['provide', 'led', 'lead', 'driv', 'support', 'enabl', 'develop',
            'increas', 'ris']
neg_rels = ['restrict', 'worsen', 'declin', 'limit', 'constrain',
            'decreas', 'hinder', 'deplet', 'reduce', 'hamper']
neu_rels = ['affect', 'impact', 'due', 'caus', 'because']
first_gen_ont_nodes = ['concept', 'process', 'property', 'entity', 'time']
bad_grnd = {'event1'}


class SofiaProcessor(object):
    def __init__(self, score_cutoff=None):
        self._entities = {}
        self._events = {}
        self._score_cutoff = score_cutoff

    @staticmethod
    def process_event(event_dict):
        mappings = [
            ('Event_Type', 'Event_Type'),
            ('Relation', 'Relation'),
            ('Location', 'Location'),
            ('Time', 'Time'),
            ('Source', 'Source_File'),
            ('Score', 'Score'),
            ('Text', 'Sentence'),
            ('Agent_index', 'Agent Index'),
            ('Patient_index', 'Patient Index'),
            ('Agent', 'Agent'),
            ('Patient', 'Patient'),
            ('Event Index', 'Event Index'),
            ('Span', 'Span')
        ]
        return {k: event_dict.get(v) for k, v in mappings}

    @staticmethod
    def process_entity(ent_dict):
        key_mappings = [
            ('Source', 'Source_File'),
            ('Query', 'Query'),
            ('Score', 'Score'),
            ('Entity Index', 'Entity Index'),
            ('Entity', 'Entity'),
            ('Entity_Type', 'Entity_Type'),
            ('Indicator', 'Indicator'),
            ('Qualifier', 'Qualifier'),
            ('Text', 'Sentence'),
            ('Span', 'Span')
        ]
        return {k: ent_dict.get(v) for k, v in key_mappings}

    def _build_influences(self, rel_dict):
        stmt_list = []
        cause_entries = rel_dict.get('Cause Index')
        effect_entries = rel_dict.get('Effect Index')

        # FIXME: Handle cases in which there is a missing cause/effect
        if not cause_entries or not effect_entries:
            return []
        causes = [c.strip() for c in cause_entries.split(', ')]
        effects = [e.strip() for e in effect_entries.split(', ')]
        rel = rel_dict.get('Relation')
        if _in_rels(rel, pos_rels):
            pol = 1
        elif _in_rels(rel, neg_rels):
            pol = -1
        elif _in_rels(rel, neu_rels):
            pol = None
        # If we don't recognize this relation, we don't get any
        # statements
        else:
            return []

        text = rel_dict.get('Sentence')
        annot_keys = ['Relation']
        annots = {k: rel_dict.get(k) for k in annot_keys}
        ref = rel_dict.get('Source_File')

        for cause_idx, effect_idx in itertools.product(causes, effects):
            cause = self._events.get(cause_idx)
            effect = self._events.get(effect_idx)
            if not cause or not effect:
                continue
            subj = self.get_event(cause)
            obj = self.get_event(effect)

            ev = Evidence(source_api='sofia', pmid=ref,
                          annotations=annots, text=text)
            stmt = Influence(subj, obj, evidence=[ev])
            # Use the polarity of the events, if object does not have a
            # polarity, use overall polarity
            if stmt.obj.delta.polarity is None:
                stmt.obj.delta.set_polarity(pol)

            stmt_list.append(stmt)
        return stmt_list

    def get_event(self, event_entry):
        # Get get compositional grounding
        comp_name, comp_grnd = self.get_compositional_grounding(event_entry)
        if comp_name is not None and \
                comp_grnd[0] is not None and \
                comp_grnd[0][0] is not None:
            concept = Concept(comp_name, db_refs={
                'TEXT': comp_name,
                'WM': [comp_grnd]
            })
        # If not try to get old style Sofia grounding
        else:
            name = event_entry['Relation']
            concept = Concept(name, db_refs={'TEXT': name})
            if event_entry['Event_Type']:
                concept.db_refs['SOFIA'] = event_entry['Event_Type']

        context = WorldContext()
        time = event_entry.get('Time')
        if time:
            context.time = TimeContext(text=time.strip())
        loc = event_entry.get('Location')
        if loc:
            context.geo_location = RefContext(name=loc)

        text = event_entry.get('Text')
        ref = event_entry.get('Source')
        agent = event_entry.get('Agent')
        patient = event_entry.get('Patient')
        anns = {}
        if agent:
            anns['agent'] = agent
        if patient:
            anns['patient'] = patient
        ev = Evidence(source_api='sofia', pmid=ref, text=text,
                      annotations=anns, source_id=event_entry['Event Index'])
        pol = event_entry.get('Polarity')
        event = Event(concept, context=context, evidence=[ev],
                      delta=QualitativeDelta(polarity=pol, adjectives=None))

        return event

    def get_relation_events(self, rel_dict):
        # Save indexes of events that are part of causal relations
        relation_events = []
        cause_entries = rel_dict.get('Cause Index')
        effect_entries = rel_dict.get('Effect Index')
        causes = [c.strip() for c in cause_entries.split(',')]
        effects = [e.strip() for e in effect_entries.split(',')]
        for ci in causes:
            relation_events.append(ci)
        for ei in effects:
            relation_events.append(ei)
        return relation_events

    def get_meaningful_events(self, raw_event_dict):
        # Only keep meaningful events and extract polarity information from
        # events showing change
        processed_event_dict = {}
        for event_index, event_info in raw_event_dict.items():
            ai = event_info['Agent_index']
            agent_index = [] if not ai else ai.split(', ')
            pi = event_info['Patient_index']
            patient_index = [] if not pi else pi.split(', ')

            if _in_rels(event_info['Relation'], pos_rels):
                pol = 1
            elif _in_rels(event_info['Relation'], neg_rels):
                pol = -1
            else:
                pol = None

            agent_embedded_events = agent_index and \
                all(a in raw_event_dict for a in agent_index)
            patient_embedded_events = patient_index and \
                all(a in raw_event_dict for a in patient_index)
            # If the agent is itself an event, we use that as the reference
            if agent_embedded_events:
                for event_ix in agent_index:
                    processed_event_dict[event_ix] = raw_event_dict[event_ix]
                    processed_event_dict[event_ix]['Polarity'] = pol
            # If the patient is itself an event, we use that as the reference
            elif patient_embedded_events:
                for event_ix in patient_index:
                    processed_event_dict[event_ix] = raw_event_dict[event_ix]
                    processed_event_dict[event_ix]['Polarity'] = pol
            # Otherwise we take the event itself
            else:
                processed_event_dict[event_index] = raw_event_dict[event_index]
        return processed_event_dict

    def get_compositional_grounding(self, event_entry):
        # According to the Sofia team:
        # Process is the "Event_Type" of event
        # Theme is the "Entity_Type" of patient of event, a wm-concept
        # Property is the "Qualifier" of patient of event, a wm-property
        #
        # Decision tree:
        #                     Theme
        #                   Y/     \N
        #                  Proc    continue
        #                Y/    \N
        #            Prop        Prop
        #          Y/    \N      Y/    \N
        #   (,-,P,P) (,-,P,-)  (,P,-,-) (,-,-,-) <- Grounding results
        #
        # ToDo:
        #  - How do we pick among multiple agents or patients? Are they
        #    assumed to be referring to the same grounded entity?
        #    + Assuming YES for now
        #  - If we have process and process property available, but the
        #    process gets removed because it's not grounded properly,
        #    should the property be remove or put to the theme property?
        name = None
        theme, theme_prop, theme_proc, theme_proc_prop = (None, )*4

        # First, try to get the theme
        event_patients = event_entry['Patient_index'].strip().split(', ')
        if event_patients != ['']:
            name, theme = self._get_entity_grounding(event_patients)

        # See if we have a theme process
        proc = self._clean_grnd_filter(event_entry['Event_Type'],
                                       float(event_entry['Score']) or 0.7,
                                       'process')

        # Next, see if we have a theme property
        if event_patients != ['']:
            prop = self._get_theme_prop(event_patients)
        else:
            prop = None

        # Set correct combination of groundings:
        # If no theme, see if we have a process that can be promoted to
        # theme, and add theme property if it exists
        if not theme or theme[0] is None:
            if proc and proc[0] is not None:
                name = event_entry['Relation']  # Match name to theme (process)
                theme = proc
                if prop and prop[0] is not None:
                    theme_prop = prop
            # We don't have a grounding, return nothing
            else:
                return None, (None, ) * 4
        # If we have a theme and process
        elif proc and proc[0] is not None:
            theme_proc = proc
            # If we have a process, add property as process property
            if prop and prop[0] is not None:
                theme_proc_prop = prop
        # If we have theme and property, but no process
        elif prop and prop[0] is not None:
            theme_prop = prop

        # Return 4-tuple of:
        # Theme, Theme Property, Theme Process, Theme Process Property
        assert not all(co == (None, 0.0) for co in
                       [theme, theme_prop, theme_proc, theme_proc_prop])
        return name, (theme, theme_prop, theme_proc, theme_proc_prop)

    def _get_theme_prop(self, entity_inds):
        qualifiers = [
            (self._entities[ai]['Qualifier'],
             float(self._entities[ai]['Score']) or 0.7)  # ignore zero scores
            for ai in entity_inds if self._entities[ai]['Qualifier']
        ]
        if qualifiers:
            # Sort by highest score first
            qualifiers.sort(key=lambda t: t[1], reverse=True)
            for qlfr in qualifiers:
                qlfr = self._clean_grnd_filter(
                    *qlfr, grnd_type='property',
                    score_cutoff=self._score_cutoff or None)
                if qlfr[0] is not None:
                    return qlfr
        return None, 0.0

    def _get_entity_grounding(self, entity_inds):
        # Get name and grounding
        grnd_ent_list = [
            (self._entities[ai]['Entity'] or None,  # Set None if ''
             (self._entities[ai]['Entity_Type'],
              float(self._entities[ai]['Score']) or 0.7))  # ignore zero scores
            for ai in entity_inds if self._entities[ai]['Entity_Type']
        ]
        if grnd_ent_list:
            # Sort by highest score first
            grnd_ent_list.sort(key=lambda t: t[1][1], reverse=True)
            for grnd_ent in grnd_ent_list:
                name, ent = grnd_ent
                ent = self._clean_grnd_filter(
                    *ent, grnd_type='concept',
                    score_cutoff=self._score_cutoff or None
                )
                if ent[0] is not None:
                    return name, ent
        return None, (None, 0.0)

    @staticmethod
    def _clean_grnd_filter(grnd, score, grnd_type, score_cutoff=0.0):
        assert isinstance(score, float)
        if not grnd or '/' not in grnd:
            return None, 0.0
        # Filter low scores if provided
        if score_cutoff and score < score_cutoff:
            return None, 0.0
        # Remove initial slash
        if grnd.startswith('/'):
            grnd = grnd[1:]
        # Add initial wm
        if grnd and not grnd.startswith('wm'):
            grnd = f'wm_compositional/{grnd_type}/{grnd}'

        # Remove special misgrounding
        if any(mg in grnd for mg in bad_grnd):
            grnd = '/'.join([g for g in grnd.split('/') if g not in bad_grnd])

        grnd.replace('//', '/')
        return grnd, score


class SofiaJsonProcessor(SofiaProcessor):
    def __init__(self, jd):
        super().__init__()
        self._entities = self.process_entities(jd)
        self._events = self.process_events(jd)
        self.statements = []
        self.relation_subj_obj_ids = []

    def process_entities(self, jd):
        ent_dict = {}
        entity_list = jd['entities']
        for entity in entity_list:
            ent_index = entity['Entity Index']
            ent_dict[ent_index] = self.process_entity(entity)

        return ent_dict

    def process_events(self, jd):
        event_dict = {}
        # First get all events from reader output
        events = jd['events']
        for event in events:
            event_index = event.get('Event Index')
            event_dict[event_index] = self.process_event(event)
        processed_event_dict = self.get_meaningful_events(event_dict)
        return processed_event_dict

    def extract_relations(self, jd):
        stmts = []
        for rel_dict in jd['causal']:
            # Save indexes of causes and effects
            current_relation_events = self.get_relation_events(rel_dict)
            for ix in current_relation_events:
                self.relation_subj_obj_ids.append(ix)
            # Make Influence Statements
            stmt_list = self._build_influences(rel_dict)
            if not stmt_list:
                continue
            stmts = stmts + stmt_list
        for stmt in stmts:
            self.statements.append(stmt)

    def extract_events(self, json_list):
        # First confirm we have extracted event information and relation events
        if not self._events:
            self.process_events(json_list)
        if not self.relation_subj_obj_ids:
            self.extract_relations(json_list)
        # Only make Event Statements from standalone events
        for event_index in self._events:
            if event_index not in self.relation_subj_obj_ids:
                event = self.get_event(self._events[event_index])
                self.statements.append(event)


class SofiaExcelProcessor(SofiaProcessor):
    def __init__(self, relation_rows, event_rows, entity_rows):
        super().__init__()
        self._events = self.process_events(event_rows)
        self.statements = []
        self.relation_subj_obj_ids = []

    def process_events(self, event_rows):
        header = [cell.value for cell in next(event_rows)]
        event_dict = {}
        for row in event_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            event_index = row_dict.get('Event Index')
            event_dict[event_index] = self.process_event(row_dict)
        processed_event_dict = self.get_meaningful_events(event_dict)
        return processed_event_dict

    def extract_relations(self, relation_rows):
        header = [cell.value for cell in next(relation_rows)]
        stmts = []
        for row in relation_rows:
            row_values = [r.value for r in row]
            row_dict = {h: v for h, v in zip(header, row_values)}
            # Save indexes of causes and effects
            current_relation_events = self.get_relation_events(row_dict)
            for ix in current_relation_events:
                self.relation_subj_obj_ids.append(ix)
            # Make Influence Statements
            stmt_list = self._build_influences(row_dict)
            if not stmt_list:
                continue
            stmts = stmts + stmt_list
        for stmt in stmts:
            self.statements.append(stmt)

    def extract_events(self, event_rows, relation_rows):
        # First confirm we have extracted event information and relation events
        if not self._events:
            self.process_events(event_rows)
        if not self.relation_subj_obj_ids:
            self.extract_relations(relation_rows)
        # Only make Event Statements from standalone events
        for event_index in self._events:
            if event_index not in self.relation_subj_obj_ids:
                event = self.get_event(self._events[event_index])
                self.statements.append(event)


def _in_rels(value, rels):
    for rel in rels:
        if value is not None:
            if value.lower().startswith(rel):
                return True
            if rel in value.lower():
                return True
    return False
