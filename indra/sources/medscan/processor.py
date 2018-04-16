import re
import os
import codecs
import lxml.etree
import collections
from indra.statements import *
from indra.databases.chebi_client import get_chebi_id_from_cas
from indra.databases.hgnc_client import get_hgnc_from_entrez, get_uniprot_id


logger = logging.getLogger('medscan')

MedscanEntity = collections.namedtuple('MedscanEntity', ['name', 'urn', 'type',
                           'properties'])


MedscanProperty = collections.namedtuple('MedscanProperty',
                                         ['type', 'name', 'urn'])


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


def parse_mod_string(s):
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


def parse_mut_string(s):
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
        # Mutation string does fit this pattern, other patterns not currently
        # supported
        return None, None, None
    else:
        return (m.group(1), m.group(2), m.group(3))


def urn_to_db_refs(urn):
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
    """
    # Convert a urn to a db_refs dictionary
    if urn is None:
        return {}

    p = 'urn:([^:]+):([^:]+)'
    m = re.match(p, urn)
    assert m is not None, m

    urn_type = m.group(1)
    urn_id = m.group(2)

    db_refs = {}

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
    return db_refs


def extract_id(id_string):
    """Extracts the numeric ID from the representation of the subject or
    object ID that appears as an attribute of the svo element in the Medscna
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


def untag_sentence(s):
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


def extract_sentence_tags(tagged_sentence):
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


class MedscanProcessor(object):
    """Processes Medscan data into INDRA statements.

    Attributes
    ----------
    statements : list<str>
        A list of extracted INDRA statements
    num_entities : int
        The total number of subject or object entities the processor attempted
        to resolve
    num_entities_not_found : int
        The number of subject or object IDs which could not be resolved by
        looking in the list of entities or tagged phrases.
    unmapped_urns : set of str
        The list of URNs that are not mapped to any external database,
        optionally extracted from rnef files provided by the Medscan
        developers. Currently unused.
    """
    def __init__(self, medscan_resource_dir):
        self.statements = []
        self.num_entities_not_found = 0
        self.num_entities = 0
        self.unmapped_urns = set()
        self.relations = []

        # Read in and populate a list of unmapped urns
        if medscan_resource_dir is not None:
            fname_unmapped_complexes = os.path.join(medscan_resource_dir,
                                                    'Unmapped Complexes.rnef')
            fname_classes = os.path.join(medscan_resource_dir,
                                         'Unmapped Functional classes.rnef')
            for fname in [fname_unmapped_complexes, fname_classes]:
                with codecs.open(fname, 'rb') as f:
                    for event, elem in lxml.etree.iterparse(f,
                                                            events=('start',
                                                                    'end'),
                                                            encoding='utf-8'):
                        if event == 'start':
                            urn = elem.attrib.get('urn')
                            if urn is not None:
                                self.unmapped_urns.add(urn)

    def process_csxml_from_file_handle(self, f, num_documents):
        """Processes a filehandle to MedScan csxml input into INDRA
        statements.

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
            elif event == 'start' and elem.tag == 'match':
                match_text = elem.attrib.get('chars')
            elif event == 'start' and elem.tag == 'entity':
                if not in_prop:
                    ent_id = elem.attrib['msid']
                    ent_urn = elem.attrib.get('urn')
                    ent_type = elem.attrib['type']
                    entities[ent_id] = MedscanEntity(match_text, ent_urn,
                                                     ent_type, {})
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
                    print("Processed %d documents"
                          % doc_counter)
                if num_documents is not None and doc_counter >= num_documents:
                    break

        print("Done processing %d documents" % doc_counter)


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
        untagged_sentence = untag_sentence(relation.tagged_sentence)
        source_id = relation.uri
        m = re.match('info:pmid/([0-9]+)', source_id)
        if m is not None:
            # Extract the pmid from the URI if the URI refers to a pmid
            pmid = m.group(1)
        annotations = None
        ev = [Evidence(source_api='medscan', source_id=source_id, pmid=pmid,
                       text=untagged_sentence, annotations=None,
                       epistemics=None)]

        # These normalized verbs are mapped to IncreaseAmount statements
        increase_amount_verbs = ['ExpressionControl-positive',
                                 'MolSynthesis-positive']

        # These normalized verbs are mapped to DecreaseAmount statements
        decrease_amount_verbs = ['ExpressionControl-negative',
                                 'MolSynthesis-negative']

        if relation.verb in increase_amount_verbs:
            # If the normalized verb corresponds to an IncreaseAmount statement
            # then make one
            self.statements.append(
                                   IncreaseAmount(subj, obj, evidence=ev)
                                  )
        elif relation.verb in decrease_amount_verbs:
            # If the normalized verb corresponds to a DecreaseAmount statement
            # then make one
            self.statements.append(
                                   DecreaseAmount(subj, obj, evidence=ev)
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

            self.statements.append(statement_type(subj, obj, evidence=ev))

        elif relation.verb == 'Binding':
            # The Binding normalized verb corresponds to the INDRA Complex
            # statement.
            self.statements.append(
                                   Complex([subj, obj], evidence=ev)
                                  )

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
        tags = extract_sentence_tags(relation.tagged_sentence)

        if entity_id is None:
            return None
        self.num_entities = self.num_entities + 1

        entity_id = extract_id(entity_id)

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
            return Agent(db_refs['TEXT'], db_refs=db_refs)
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

                db_refs = urn_to_db_refs(protein.urn)
                db_refs['TEXT'] = protein.name

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
                    r_old, pos, r_new = parse_mut_string(mutation.name)
                    if r_old is None:
                        logger.warning('Could not parse mutation string: ' +
                                       mutation.name)
                        # Don't create an agent
                        return None
                    else:
                        try:
                            cond = MutCondition(pos, r_old, r_new)
                            return Agent(db_refs['TEXT'], db_refs=db_refs,
                                         mutations=[cond])
                        except BaseException:
                            logger.warning('Could not parse mutation ' +
                                           'string: ' + mutation.name)
                            return None
                elif mutation.type == 'MethSite':
                    # Convert methylation site information to an INDRA
                    # ModCondition
                    res, pos = parse_mod_string(mutation.name)
                    if res is None:
                        return None
                    cond = ModCondition('methylation', res, pos)
                    return Agent(db_refs['TEXT'], db_refs=db_refs,
                                 mods=[cond])

                    # Example:
                    # MedscanEntity(name='R457',
                    # urn='urn:agi-s-llid:R457-2185', type='MethSite',
                    # properties=None)
                elif mutation.type == 'PhosphoSite':
                    # Convert phosphorylation site information to an INDRA
                    # ModCondition
                    res, pos = parse_mod_string(mutation.name)
                    if res is None:
                        return None
                    cond = ModCondition('phosphorylation', res, pos)
                    return Agent(db_refs['TEXT'], db_refs=db_refs,
                                 mods=[cond])

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
                db_refs = urn_to_db_refs(entity.urn)
                db_refs['TEXT'] = entity.name
                return Agent(db_refs['TEXT'], db_refs=db_refs)
