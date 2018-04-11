#from indra.sources.medscan.api import MedscanRelation, MedscanEntity
#from indra.sources.medscan.api import *
import re
import sys
from indra.statements import Agent
import lxml.etree
import os
import codecs
import collections

from indra.databases.hgnc_client import get_hgnc_from_entrez
from indra.databases.chebi_client import get_chebi_id_from_cas

def urn_to_db_refs(urn):
    # Convert a urn to a db_refs dictionary
    if urn is None:
        return {}

    p = 'urn:([^:]+):([^:]+)'
    m = re.match(p, urn)
    assert m is not None, m

    urn_type = m.group(1)
    urn_id = m.group(2)

    db_refs = {}

    if urn_type == 'agi-cas':
        # Identifier is CAS, convert to CHEBI
        db_refs['CHEBI'] = get_chebi_id_from_cas(urn_id)
    elif urn_type == 'agi-llid':
        # This is an Entrez ID, convert to HGNC
        db_refs['HGNC'] = get_hgnc_from_entrez(urn_id)
    elif urn_type == 'agi-ncimorgan':
        # Identifier is MESH
        db_refs['MESH'] = urn_id
    elif urn_type == 'agi-ncimcelltype':
        # Identifier is MESH
        db_refs['MESH'] = urn_id
    elif urn_type == 'agi-meshdis':
        # Identifier is MESH
        db_refs['MESH'] = urn_id
    elif urn_type == 'agi-gocomplex':
        # Identifier is GO
        db_refs['GO'] = urn_id
    elif urn_type == 'agi-go':
        # Identifier is GO
        db_refs['GO'] = urn_id
    elif urn_type == 'agi-ncimtissue':
        # Identifier is MESH
        db_refs['MESH'] = urn_id
    return db_refs


def extract_id(id_string):
    p = 'ID\\{([0-9]+)\\}'
    matches = re.match(p, id_string)
    assert(matches is not None)
    return matches.group(1)

def untag_sentence(s):
    p = 'ID{[0-9,]+=([^}]+)}'
    s = re.sub(p, '\\1', s)
    #
    s = re.sub('CONTEXT{[^}]+}', '', s)
    s = re.sub('GLOSSARY{[^}]+}', '', s)
    return s

def extract_sentence_tags(tagged_sentence):
    p = re.compile('ID{([0-9,]+)=([^}]+)}')
    tags = {}

    endpos = 0
    while True:
        match = p.search(tagged_sentence, pos=endpos)
        if not match:
            break
        endpos = match.end()

        tags[match.group(1)] = match.group(2)
    return tags

class MedscanProcessor(object):
    def __init__(self, medscan_resource_dir):
        self.statements = []
        self.urn_examples = collections.defaultdict(set)
        self.num_entities_not_found = 0
        self.num_entities = 0
        self.unmapped_urns = set()

        # Read in and populate a list of unmapped urns
        fname_unmapped_complexes = os.path.join(medscan_resource_dir,
                                                'Unmapped Complexes.rnef')
        fname_classes = os.path.join(medscan_resource_dir,
                                     'Unmapped Functional classes.rnef')
        for fname in [fname_unmapped_complexes, fname_classes]:
            with codecs.open(fname, 'rb') as f:
                for event, elem in lxml.etree.iterparse(f, 
                                                        events=('start', 'end'),
                                                        encoding='utf-8'):
                    if event == 'start':
                        urn = elem.attrib.get('urn')
                        if urn is not None:
                            self.unmapped_urns.add(urn)
        print('# unmapped URNs:', len(self.unmapped_urns))


    def process_relation(self, relation):
        # Extract sentence tags mapping ids to the text. We refer to this
        # mapping only if the entity doesn't appear in the grounded entity
        # list
        relation.tags = extract_sentence_tags(relation.tagged_sentence)

        subj = self.agent_from_entity(relation, relation.subj)
        obj = self.agent_from_entity(relation, relation.subj)
        if subj is None or obj is None:
            # Don't extract a statement if the subject or object cannot
            # be resolved
            return
        
        if relation.verb == 'UnknownRegulation-positive':
            pass
        elif relation.verb == 'ExpressionControl-positive':
            pass
        elif relation.verb == 'Binding':
            import ipdb; ipdb.set_trace()

    def agent_from_entity(self, relation, entity_id):
        if entity_id is None:
            return None
        self.num_entities = self.num_entities + 1

        entity_id = extract_id(entity_id)

        if entity_id not in relation.entities:
            # The entity is not in the grounded entity list
            # Instead, make an ungrounded entity, with entity text that
            # appeared in the sentence

            if entity_id not in relation.tags:
                self.num_entities_not_found = self.num_entities_not_found + 1
                print('Entity not found')
                return None
            entity_text = relation.tags[entity_id]
            db_refs = {'TEXT': entity_text}
            return Agent(entity_text, db_refs=db_refs)
        else:
            entity = relation.entities[entity_id]
            db_refs = {}
            db_refs['TEXT'] = entity.match_text

            # Parse the URN
            if entity.urn is not None:
                #print(entity.urn)
                p = 'urn:([^:]+):([^:]+)'
                m = re.match(p, entity.urn)
                assert m is not None, m
                if entity.urn not in self.unmapped_urns:
                    self.urn_examples[m.group(1)].add(entity)
            else:
                # urn is None
                pass

