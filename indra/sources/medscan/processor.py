import re
import os
import codecs
import lxml.etree
import collections
from indra.statements import *
from indra.databases.chebi_client import get_chebi_id_from_cas
from indra.databases.hgnc_client import get_hgnc_from_entrez, get_uniprot_id, \
        get_hgnc_name
from indra.util import read_unicode_csv
from indra.sources.reach.processor import ReachProcessor, Site

logger = logging.getLogger('medscan')


MedscanEntity = collections.namedtuple('MedscanEntity', ['name', 'urn', 'type',
                                                         'properties'])


MedscanProperty = collections.namedtuple('MedscanProperty',
                                         ['type', 'name', 'urn'])

def _read_famplex_map():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         '../../resources/famplex_map.tsv')
    famplex_map = {}
    csv_rows = read_unicode_csv(fname, delimiter='\t')
    for row in csv_rows:
        source_ns = row[0]
        source_id = row[1]
        be_id = row[2]
        famplex_map[(source_ns, source_id)] = be_id
    return famplex_map


famplex_map = _read_famplex_map()


def is_statement_in_list(statement, statement_list):
    """Return True of given statement is equivalent to on in a list

    Determines whether the statement is equivalent to any statement in the
    given list of statements, with equivalency determined by Statement's
    equals method.

    Parameters
    ----------
    statement : indra.statements.Statement
        The statement to compare with
    statement_list : list[indra.statements.Statement]
        The statement list whose entries we compare with statement

    Returns
    -------
    in_list : bool
        True if statement is equivalent to any statements in the list
    """
    for s in statement_list:
        if s.equals(statement):
            return True
    return False


class ProteinSiteInfo(object):
    """Represent a site on a protein, extracted from a StateEffect event.

    Parameters
    ----------
    site_text : str
        The site as a string (ex. S22)
    object_text : str
        The protein being modified, as the string that appeared in the original
        sentence
    """
    def __init__(self, site_text, object_text):
        self.site_text = site_text
        self.object_text = object_text

    def get_sites(self):
        """Parse the site-text string and return a list of sites.

        Returns
        -------
        sites : list[Site]
            A list of position-residue pairs corresponding to the site-text
        """
        st = self.site_text
        suffixes = [' residue', ' residues', ',', '/']
        for suffix in suffixes:
            if st.endswith(suffix):
                st = st[:-len(suffix)]
        assert(not st.endswith(','))

        # Strip parentheses
        st = st.replace('(', '')
        st = st.replace(')', '')
        st = st.replace(' or ', ' and ')  # Treat end and or the same

        sites = []
        parts = st.split(' and ')
        for part in parts:
            if part.endswith(','):
                part = part[:-1]
            if len(part.strip()) > 0:
                sites.extend(ReachProcessor._parse_site_text(part.strip()))
        return sites


class MedscanProcessor(object):
    """Processes Medscan data into INDRA statements.

    The special StateEffect event conveys information about the binding
    site of a protein modification. Sometimes this is paired with additional
    event information in a seperate SVO. When we encounter a StateEffect, we
    don't process into an INDRA statement right away, but instead store
    the site information and use it if we encounter a ProtModification
    event within the same sentence.

    Attributes
    ----------
    statements : list<str>
        A list of extracted INDRA statements
    sentence_statements : list<str>
        A list of statements for the sentence we are currently processing.
        Deduplicated and added to the main statement list when we finish
        processing a sentence.
    num_entities : int
        The total number of subject or object entities the processor attempted
        to resolve
    num_entities_not_found : int
        The number of subject or object IDs which could not be resolved by
        looking in the list of entities or tagged phrases.
    last_site_info_in_sentence : SiteInfo
        Stored protein site info from the last StateEffect event within the
        sentence, allowing us to combine information from StateEffect and
        ProtModification events within a single sentence in a single INDRA
        statement. This is reset at the end of each sentence
    """
    def __init__(self):
        self.statements = []
        self.sentence_statements = []
        self.num_entities_not_found = 0
        self.num_entities = 0
        self.relations = []
        self.last_site_info_in_sentence = None

    def process_csxml_from_file_handle(self, f, num_documents):
        """Processes a filehandle to MedScan csxml input into INDRA
        statements.

        The CSXML format consists of a top-level `<batch>` root element
        containing a series of `<doc>` (document) elements, in turn containing
        `<sec>` (section) elements, and in turn containing `<sent>` (sentence)
        elements.

        Within the `<sent>` element, a series of additional elements appear in
        the following order:

        * `<toks>`, which contains a tokenized form of the sentence in its text
          attribute
        * `<textmods>`, which describes any preprocessing/normalization done to
          the underlying text
        * `<match>` elements, each of which contains one of more `<entity>`
          elements, describing entities in the text with their identifiers.
          The local IDs of each entities are given in the `msid` attribute of
          this element; these IDs are then referenced in any subsequent SVO
          elements.
        * `<svo>` elements, representing subject-verb-object triples. SVO
          elements with a `type` attribute of `CONTROL` represent normalized
          regulation relationships; they often represent the normalized
          extraction of the immediately preceding (but unnormalized SVO
          element). However, in some cases there can be a "CONTROL" SVO
          element without its parent immediately preceding it.

        Parameters
        ----------
        f : file object
            A filehandle to a source of MedScan csxml data
        num_documents : int
            The number of documents to process, or None to process all
            documents in the input stream
        """
        pmid = None
        sec = None
        tagged_sent = None
        svo_list = []
        doc_counter = 0
        entities = {}
        match_text = None
        in_prop = False
        last_relation = None
        property_entities = []
        property_name = None

        self.log_entities = collections.defaultdict(int)

        # Go through the document again and extract statements
        for event, elem in lxml.etree.iterparse(f, events=('start', 'end'),
                                                encoding='utf-8',
                                                recover=True):
            # If opening up a new doc, set the PMID
            if event == 'start' and elem.tag == 'doc':
                pmid = elem.attrib.get('uri')
            # If getting a section, set the section type
            elif event == 'start' and elem.tag == 'sec':
                sec = elem.attrib.get('type')
            # Set the sentence context
            elif event == 'start' and elem.tag == 'sent':
                tagged_sent = elem.attrib.get('msrc')

                # Reset last_relation between sentences, since we will only be
                # interested in the relation immediately preceding a CONTROL
                # statement but within the same sentence.
                last_relation = None

                entities = {}
            elif event == 'end' and elem.tag == 'sent':
                # End of sentence; deduplicate and copy statements from this
                # sentence to the main statements list

                # Make a list of those statments in sentence_statements sans
                # duplicates
                statements_to_add = []
                for s in self.sentence_statements:
                    if not is_statement_in_list(s, statements_to_add):
                        statements_to_add.append(s)

                # Add deduplicated statements to the main statements list
                self.statements.extend(statements_to_add)

                # Reset sentence statements list to prepare for processing the
                # next sentence
                self.sentence_statements = []

                # Reset site info
                self.last_site_info_in_sentence = None
            elif event == 'start' and elem.tag == 'match':
                match_text = elem.attrib.get('chars')
            elif event == 'start' and elem.tag == 'entity':
                if not in_prop:
                    ent_id = elem.attrib['msid']
                    ent_urn = elem.attrib.get('urn')
                    ent_type = elem.attrib['type']
                    entities[ent_id] = MedscanEntity(match_text, ent_urn,
                                                     ent_type, {})

                    tuple_key = (ent_type, elem.attrib.get('name'), ent_urn)
                    if ent_type == 'Complex' or ent_type == 'FunctionalClass':
                        self.log_entities[tuple_key] = \
                                self.log_entities[tuple_key] + 1
                else:
                    ent_type = elem.attrib['type']
                    ent_urn = elem.attrib['urn']
                    ent_name = elem.attrib['name']
                    property_entities.append(MedscanEntity(ent_name, ent_urn,
                                                           ent_type, None))
            elif event == 'start' and elem.tag == 'svo':
                subj = elem.attrib.get('subj')
                verb = elem.attrib.get('verb')
                obj = elem.attrib.get('obj')
                svo_type = elem.attrib.get('type')
                svo = {'uri': pmid,
                       'sec': sec,
                       'text': tagged_sent,
                       'entities': entities}
                svo.update(elem.attrib)

                # Aggregate information about the relation
                relation = MedscanRelation(
                                       uri=pmid,
                                       sec=sec,
                                       tagged_sentence=tagged_sent,
                                       entities=entities,
                                       subj=subj,
                                       verb=verb,
                                       obj=obj,
                                       svo_type=svo_type,
                                      )
                self.relations.append(relation)
                if svo_type == 'CONTROL':
                    self.process_relation(relation, last_relation)
                else:
                    # Sometimes a CONTROL SVO can be after an unnormalized SVO
                    # that is a more specific but less uniform version of the
                    # same extracted statement.
                    last_relation = relation
            elif event == 'start' and elem.tag == 'prop':
                in_prop = True
                property_name = elem.attrib.get('name')
                property_entities = []
            elif event == 'end' and elem.tag == 'prop':
                in_prop = False
                entities[ent_id].properties[property_name] = property_entities
            elif event == 'end' and elem.tag == 'doc':
                doc_counter += 1
                # Give a status update
                if doc_counter % 100 == 0:
                    logger.info("Processed %d documents" % doc_counter)
                if num_documents is not None and doc_counter >= num_documents:
                    break


    def process_relation(self, relation, last_relation):
        """Process a relation into an INDRA statement.

        Parameters
        ----------
        relation : MedscanRelation
            The relation to process (a CONTROL svo with normalized verb)
        last_relation : MedscanRelation
            The relation immediately proceding the relation to process within
            the same sentence, or None if there are no preceding relations
            within the same sentence. This proceeding relation, if available,
            will refer to the same interaction but with an unnormalized
            (potentially more specific) verb, and is used when processing
            protein modification events.
        """
        subj = self.agent_from_entity(relation, relation.subj)
        obj = self.agent_from_entity(relation, relation.obj)
        if subj is None or obj is None:
            # Don't extract a statement if the subject or object cannot
            # be resolved
            return

        # Make evidence object
        untagged_sentence = _untag_sentence(relation.tagged_sentence)
        source_id = relation.uri
        m = re.match('info:pmid/([0-9]+)', source_id)
        if m is not None:
            # Extract the pmid from the URI if the URI refers to a pmid
            pmid = m.group(1)
        if last_relation:
            last_verb = last_relation.verb
        else:
            last_verb = None
        annotations = {'verb': relation.verb, 'last_verb': last_verb}
        epistemics = dict()
        epistemics['direct'] = False  # Overridden later if needed
        ev = [Evidence(source_api='medscan', source_id=source_id, pmid=pmid,
                       text=untagged_sentence, annotations=annotations,
                       epistemics=epistemics)]

        # These normalized verbs are mapped to IncreaseAmount statements
        increase_amount_verbs = ['ExpressionControl-positive',
                                 'MolSynthesis-positive',
                                 'CellExpression',
                                 'QuantitativeChange-positive',
                                 'PromoterBinding']

        # These normalized verbs are mapped to DecreaseAmount statements
        decrease_amount_verbs = ['ExpressionControl-negative',
                                 'MolSynthesis-negative',
                                 'miRNAEffect-negative',
                                 'QuantitativeChange-negative']

        # These normalized verbs are mapped to Activation statements (indirect)
        activation_verbs = ['UnknownRegulation-positive',
                            'Regulation-positive']
        # These normalized verbs are mapped to Activation statements (direct)
        d_activation_verbs = ['DirectRegulation-positive',
                              'DirectRegulation-positive--direct interaction']
        # All activation verbs
        all_activation_verbs = list(activation_verbs)
        all_activation_verbs.extend(d_activation_verbs)

        # These normalized verbs are mapped to Inhibition statements (indirect)
        inhibition_verbs = ['UnknownRegulation-negative',
                            'Regulation-negative']
        # These normalized verbs are mapped to Inhibition statements (direct)
        d_inhibition_verbs = ['DirectRegulation-negative',
                              'DirectRegulation-negative--direct interaction']
        # All inhibition verbs
        all_inhibition_verbs = list(inhibition_verbs)
        all_inhibition_verbs.extend(d_inhibition_verbs)

        if relation.verb in increase_amount_verbs:
            # If the normalized verb corresponds to an IncreaseAmount statement
            # then make one
            self.sentence_statements.append(
                                   IncreaseAmount(subj, obj, evidence=ev)
                                  )
        elif relation.verb in decrease_amount_verbs:
            # If the normalized verb corresponds to a DecreaseAmount statement
            # then make one
            self.sentence_statements.append(
                                   DecreaseAmount(subj, obj, evidence=ev)
                                  )
        elif relation.verb in all_activation_verbs:
            # If the normalized verb corresponds to an Activation statement,
            # then make one
            if relation.verb in d_activation_verbs:
                ev[0].epistemics['direction'] = True
            self.sentence_statements.append(
                    Activation(subj, obj, evidence=ev)
                    )
        elif relation.verb in all_inhibition_verbs:
            # If the normalized verb corresponds to an Inhibition statement,
            # then make one
            if relation.verb in d_inhibition_verbs:
                ev[0].epistemics['direct'] = True
            self.sentence_statements.append(
                    Inhibition(subj, obj, evidence=ev)
                    )

        elif relation.verb == 'ProtModification':
            # The normalized verb 'ProtModification' is too vague to make
            # an INDRA statement. We look at the unnormalized verb in the
            # previous svo element, if available, to decide what type of
            # INDRA statement to construct.

            if last_relation is None:
                # We cannot make a statement unless we have more fine-grained
                # information on the relation type from a preceding
                # unnormalized SVO
                return

            # Map the unnormalized verb to an INDRA statement type
            statement_type = None
            if last_relation.verb == 'TK{phosphorylate}':
                statement_type = Phosphorylation
            elif last_relation.verb == 'TK{dephosphorylate}':
                statement_type = Dephosphorylation
            elif last_relation.verb == 'TK{ubiquitinate}':
                statement_type = Ubiquitination
            elif last_relation.verb == 'TK{acetylate}':
                statement_type = Acetylation
            elif last_relation.verb == 'TK{methylate}':
                statement_type = Methylation
            elif last_relation.verb == 'TK{deacetylate}':
                statement_type = Deacetylation
            elif last_relation.verb == 'TK{demethylate}':
                statement_type = Demethylation
            elif last_relation.verb == 'TK{hyperphosphorylate}':
                statement_type = Phosphorylation
            elif last_relation.verb == 'TK{hydroxylate}':
                statement_type = Hydroxylation
            elif last_relation.verb == 'TK{sumoylate}':
                statement_type = Sumoylation
            elif last_relation.verb == 'TK{palmitoylate}':
                statement_type = Palmitoylation
            elif last_relation.verb == 'TK{glycosylate}':
                statement_type = Glycosylation
            elif last_relation.verb == 'TK{ribosylate}':
                statement_type = Ribosylation
            elif last_relation.verb == 'TK{deglycosylate}':
                statement_type = Deglycosylation
            elif last_relation.verb == 'TK{myristylate}':
                statement_type = Myristoylation
            elif last_relation.verb == 'TK{farnesylate}':
                statement_type = Farnesylation
            elif last_relation.verb == 'TK{desumoylate}':
                statement_type = Desumoylation
            elif last_relation.verb == 'TK{geranylgeranylate}':
                statement_type = Geranylgeranylation
            elif last_relation.verb == 'TK{deacylate}':
                statement_type = Deacetylation
            else:
                # This unnormalized verb is not handled, do not extract an
                # INDRA statement
                return

            obj_text = obj.db_refs['TEXT']
            last_info = self.last_site_info_in_sentence
            if last_info is not None and obj_text == last_info.object_text:
                for site in self.last_site_info_in_sentence.get_sites():
                    r = site.residue
                    p = site.position

                    s = statement_type(subj, obj, residue=r, position=p,
                                       evidence=ev)
                    self.sentence_statements.append(s)
            else:
                self.sentence_statements.append(statement_type(subj, obj,
                                                evidence=ev))

        elif relation.verb == 'Binding':
            # The Binding normalized verb corresponds to the INDRA Complex
            # statement.
            self.sentence_statements.append(
                                   Complex([subj, obj], evidence=ev)
                                  )
        elif relation.verb == 'ProtModification-negative':
            pass  # TODO? These occur so infrequently so maybe not worth it
        elif relation.verb == 'Regulation-unknown':
            pass  # TODO? These occur so infrequently so maybe not worth it
        elif relation.verb == 'StateEffect-positive':
            pass
            # self.sentence_statements.append(
            #                       ActiveForm(subj, obj, evidence=ev)
            #                      )
            # TODO: disabling for now, since not sure whether we should set
            # the is_active flag
        elif relation.verb == 'StateEffect':
            self.last_site_info_in_sentence = \
                    ProteinSiteInfo(site_text=subj.name,
                                    object_text=obj.db_refs['TEXT'])

    def agent_from_entity(self, relation, entity_id):
        """Create a (potentially grounded) INDRA Agent object from a given a
        Medscan entity describing the subject or object.

        Uses helper functions to convert a Medscan URN to an INDRA db_refs
        grounding dictionary.

        If the entity has properties indicating that it is a protein with
        a mutation or modification, then constructs the needed ModCondition
        or MutCondition.

        Parameters
        ----------
        relation : MedscanRelation
            The current relation being processed
        entity_id : str
            The ID of the entity to process

        Returns
        -------
        agent : indra.statements.Agent
            A potentially grounded INDRA agent representing this entity
        """
        # Extract sentence tags mapping ids to the text. We refer to this
        # mapping only if the entity doesn't appear in the grounded entity
        # list
        tags = _extract_sentence_tags(relation.tagged_sentence)

        if entity_id is None:
            return None
        self.num_entities = self.num_entities + 1

        entity_id = _extract_id(entity_id)

        if entity_id not in relation.entities and \
                entity_id not in tags:
            # Could not find the entity in either the list of grounded
            # entities of the items tagged in the sentence. Happens for
            # a very small percentage of the dataset.
            self.num_entities_not_found = self.num_entities_not_found + 1
            return None

        if entity_id not in relation.entities:
            # The entity is not in the grounded entity list
            # Instead, make an ungrounded entity, with TEXT corresponding to
            # the words with the given entity id tagged in the sentence.
            entity_text = tags[entity_id]
            db_refs = {'TEXT': entity_text}
            return Agent(normalize_medscan_name(db_refs['TEXT']),
                         db_refs=db_refs)
        else:
            entity = relation.entities[entity_id]

            prop = entity.properties
            if len(prop.keys()) == 2 and 'Protein' in prop \
               and 'Mutation' in prop:
                # Handle the special case where the entity is a protein
                # with a mutation or modification, with those details
                # described in the entity properties
                protein = prop['Protein']
                assert(len(protein) == 1)
                protein = protein[0]

                mutation = prop['Mutation']
                assert(len(mutation) == 1)
                mutation = mutation[0]

                db_refs, db_name = _urn_to_db_refs(protein.urn)

                if db_refs is None:
                    return None
                db_refs['TEXT'] = protein.name

                if db_name is None:
                    agent_name = db_refs['TEXT']
                else:
                    agent_name = db_name

                # Check mutation.type. Only some types correspond to situations
                # that can be represented in INDRA; return None if we cannot
                # map to an INDRA statement (which will block processing of
                # the statement in process_relation).
                if mutation.type == 'AASite':
                    # Do not handle this
                    # Example:
                    # MedscanEntity(name='D1', urn='urn:agi-aa:D1',
                    # type='AASite', properties=None)
                    return None
                elif mutation.type == 'Mutation':
                    # Convert mutation properties to an INDRA MutCondition
                    r_old, pos, r_new = _parse_mut_string(mutation.name)
                    if r_old is None:
                        logger.warning('Could not parse mutation string: ' +
                                       mutation.name)
                        # Don't create an agent
                        return None
                    else:
                        try:
                            cond = MutCondition(pos, r_old, r_new)
                            return Agent(normalize_medscan_name(agent_name),
                                         db_refs=db_refs, mutations=[cond])
                        except BaseException:
                            logger.warning('Could not parse mutation ' +
                                           'string: ' + mutation.name)
                            return None
                elif mutation.type == 'MethSite':
                    # Convert methylation site information to an INDRA
                    # ModCondition
                    res, pos = _parse_mod_string(mutation.name)
                    if res is None:
                        return None
                    cond = ModCondition('methylation', res, pos)
                    return Agent(normalize_medscan_name(agent_name),
                                 db_refs=db_refs, mods=[cond])

                    # Example:
                    # MedscanEntity(name='R457',
                    # urn='urn:agi-s-llid:R457-2185', type='MethSite',
                    # properties=None)
                elif mutation.type == 'PhosphoSite':
                    # Convert phosphorylation site information to an INDRA
                    # ModCondition
                    res, pos = _parse_mod_string(mutation.name)
                    if res is None:
                        return None
                    cond = ModCondition('phosphorylation', res, pos)
                    return Agent(normalize_medscan_name(agent_name),
                                 db_refs=db_refs, mods=[cond])

                    # Example:
                    # MedscanEntity(name='S455',
                    # urn='urn:agi-s-llid:S455-47', type='PhosphoSite',
                    # properties=None)
                    pass
                elif mutation.type == 'Lysine':
                    # Ambiguous whether this is a methylation or
                    # demethylation; skip

                    # Example:
                    # MedscanEntity(name='K150',
                    # urn='urn:agi-s-llid:K150-5624', type='Lysine',
                    # properties=None)
                    return None
                else:
                    logger.warning('Processor currently cannot process ' +
                                   'mutations of type ' + mutation.type)
            else:
                # Handle the more common case where we just ground the entity
                # without mutation or modification information
                db_refs, db_name = _urn_to_db_refs(entity.urn)
                if db_refs is None:
                    return None
                db_refs['TEXT'] = entity.name

                if db_name is None:
                    agent_name = db_refs['TEXT']
                else:
                    agent_name = db_name

                return Agent(normalize_medscan_name(agent_name),
                             db_refs=db_refs)


class MedscanRelation(object):
    """A structure representing the information contained in a Medscan
    SVO xml element as well as associated entities and properties.

    Attributes
    ----------
    uri : str
        The URI of the current document (such as a PMID)
    sec : str
        The section of the document the relation occurs in
    entities : dict
        A dictionary mapping entity IDs from the same sentence to MedscanEntity
        objects.
    tagged_sentence : str
        The sentence from which the relation was extracted, with some tagged
        phrases and annotations.
    subj : str
        The entity ID of the subject
    verb : str
        The verb in the relationship between the subject and the object
    obj : str
        The entity ID of the object
    svo_type : str
        The type of SVO relationship (for example, CONTROL indicates
        that the verb is normalized)
    """
    def __init__(self, uri, sec, entities, tagged_sentence, subj, verb, obj,
                 svo_type):
        self.uri = uri
        self.sec = sec
        self.entities = entities
        self.tagged_sentence = tagged_sentence

        self.subj = subj
        self.verb = verb
        self.obj = obj

        self.svo_type = svo_type


def normalize_medscan_name(name):
    """Removes the "complex" and "complex complex" suffixes from a medscan
    agent name so that it better corresponds with the grounding map.

    Parameters
    ----------
    name: str
        The Medscan agent name

    Returns
    -------
    norm_name: str
        The Medscan agent name with the "complex" and "complex complex"
        suffixes removed.
    """
    suffix = ' complex'

    for i in range(2):
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name



def _parse_mod_string(s):
    """Parses a string referring to a protein modification of the form
    (residue)(position), such as T47.

    Parameters
    ----------
    s : str
        A string representation of a protein residue and position being
        modified

    Returns
    -------
    residue : str
        The residue being modified (example: T)
    position : str
        The position at which the modification is happening (example: 47)
    """
    m = re.match('([A-Za-z])+([0-9]+)', s)
    assert(m is not None)
    return (m.group(1), m.group(2))


def _parse_mut_string(s):
    """
    A string representation of a protein mutation of the form
    (old residue)(position)(new residue). Example: T34U.

    Parameters
    ----------
    s : str
        The string representation of the protein mutation

    Returns
    -------
    old_residue : str
        The old residue, or None of the mutation string cannot be parsed
    position : str
        The position at which the mutation occurs, or None if the mutation
        string cannot be parsed
    new_residue : str
        The new residue, or None if the mutation string cannot be parsed
    """
    m = re.match('([A-Za-z]+)([0-9]+)([A-Za-z]+)', s)
    if m is None:
        # Mutation string does not fit this pattern, other patterns not
        # currently supported
        return None, None, None
    else:
        return (m.group(1), m.group(2), m.group(3))


def _urn_to_db_refs(urn):
    """Converts a Medscan URN to an INDRA db_refs dictionary with grounding
    information.

    Parameters
    ----------
    url : str
        A Medscan URN

    Returns
    -------
    db_refs : dict
        A dictionary with grounding information, mapping databases to database
        identifiers. If the Medscan URN is not recognized, returns an empty
        dictionary.
    db_name : str
        The Famplex name, if available; otherwise the HGNC name if available;
        otherwise None
    """
    # Convert a urn to a db_refs dictionary
    if urn is None:
        return {}, None

    p = 'urn:([^:]+):([^:]+)'
    m = re.match(p, urn)
    if m is None:
        return None, None

    urn_type = m.group(1)
    urn_id = m.group(2)

    db_refs = {}
    db_name = None

    # TODO: support more types of URNs
    if urn_type == 'agi-cas':
        # Identifier is CAS, convert to CHEBI
        chebi_id = get_chebi_id_from_cas(urn_id)
        if chebi_id:
            db_refs['CHEBI'] = 'CHEBI:%s' % chebi_id
    elif urn_type == 'agi-llid':
        # This is an Entrez ID, convert to HGNC
        hgnc_id = get_hgnc_from_entrez(urn_id)
        if hgnc_id is not None:
            db_refs['HGNC'] = hgnc_id

            # Convert the HGNC ID to a Uniprot ID
            uniprot_id = get_uniprot_id(hgnc_id)
            db_refs['UP'] = uniprot_id

            # Try to lookup HGNC name; if it's available, set it to the
            # agent name
            db_name = get_hgnc_name(hgnc_id)
    elif urn_type == 'agi-ncimorgan':
        # Identifier is MESH
        db_refs['MESH'] = urn_id
    elif urn_type == 'agi-ncimcelltype':
        # Identifier is MESH
        db_refs['MESH'] = urn_id
    elif urn_type == 'agi-meshdis':
        # Identifier is MESH
        db_refs['MESHDIS'] = urn_id
    elif urn_type == 'agi-gocomplex':
        # Identifier is GO
        db_refs['GO'] = 'GO:%s' % urn_id
    elif urn_type == 'agi-go':
        # Identifier is GO
        db_refs['GO'] = 'GO:%s' % urn_id
    elif urn_type == 'agi-ncimtissue':
        # Identifier is MESH
        db_refs['MESH'] = urn_id

    # If we have a GO or MESH grounding, see if there is a corresponding
    # Famplex grounding
    db_sometimes_maps_to_famplex = ['GO', 'MESH']
    for db in db_sometimes_maps_to_famplex:
        if db in db_refs:
            key = (db, db_refs[db])
            if key in famplex_map:
                db_refs['FPLX'] = famplex_map[key]

    # If the urn corresponds to an eccode, groudn to famplex if that eccode
    # is in the Famplex equivalences table
    if urn.startswith('urn:agi-enz'):
        tokens = urn.split(':')
        eccode = tokens[2]
        key = ('ECCODE', eccode)
        if key in famplex_map:
            db_refs['FPLX'] = famplex_map[key]

    # If the Medscan URN itself maps to a Famplex id, add a Famplex grounding
    key = ('MEDSCAN', urn)
    if key in famplex_map:
        db_refs['FPLX'] = famplex_map[key]

    # If there is a Famplex grounding, use Famplex for entity name
    if 'FPLX' in db_refs:
        db_name = db_refs['FPLX']

    return db_refs, db_name


def _extract_id(id_string):
    """Extracts the numeric ID from the representation of the subject or
    object ID that appears as an attribute of the svo element in the Medscan
    XML document.

    Parameters
    ----------
    id_string : str
        The ID representation that appears in the svo element in the XML
        document (example: ID{123})

    Returns
    -------
    id : str
        The numeric ID, extracted from the svo element's attribute
        (example: 123)
    """
    p = 'ID\\{([0-9]+)\\}'
    matches = re.match(p, id_string)
    assert(matches is not None)
    return matches.group(1)


def _untag_sentence(s):
    """Removes all tags in the sentence, returning the original sentence
    without Medscan annotations.

    Parameters
    ----------
    s : str
        The tagged sentence

    Returns
    -------
    untagged_sentence : str
        Sentence with tags and annotations stripped out
    """
    p = 'ID{[0-9,]+=([^}]+)}'
    s = re.sub(p, '\\1', s)

    s = re.sub('CONTEXT{[^}]+}', '', s)
    s = re.sub('GLOSSARY{[^}]+}', '', s)
    return s


def _extract_sentence_tags(tagged_sentence):
    """Given a tagged sentence, extracts a dictionary mapping tags to the words
    or phrases that they tag.

    Parameters
    ----------
    tagged_sentence : str
        The sentence with Medscan annotations and tags

    Returns
    -------
    tags : dict
        A dictionary mapping tags to the words or phrases that they tag.
    """
    p = re.compile('ID{([0-9,]+)=([^}]+)}')
    tags = {}

    # Iteratively look for all matches of this pattern
    endpos = 0
    while True:
        match = p.search(tagged_sentence, pos=endpos)
        if not match:
            break
        endpos = match.end()

        tags[match.group(1)] = match.group(2)
    return tags
