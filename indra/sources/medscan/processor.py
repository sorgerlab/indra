import re
import sys
from indra.statements import *
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
        hgnc_id = get_hgnc_from_entrez(urn_id)
        if hgnc_id is not None:
            db_refs['HGNC'] = hgnc_id
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
        self.modification_examples = collections.defaultdict(int)
        self.num_entities_not_found = 0
        self.num_entities = 0
        self.unmapped_urns = set()

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

    def process_relation(self, relation, last_relation):
        subj = self.agent_from_entity(relation, relation.subj)
        obj = self.agent_from_entity(relation, relation.obj)
        if subj is None or obj is None:
            # Don't extract a statement if the subject or object cannot
            # be resolved
            return

        # Make evidence
        untagged_sentence = untag_sentence(relation.tagged_sentence)
        source_id = relation.uri
        m = re.match('info:pmid/([0-9]+)', source_id)
        assert(m is not None)
        pmid = m.group(1)
        annotations = None
        ev = [Evidence(source_api='medscan', source_id=source_id, pmid=pmid,
                       text=untagged_sentence, annotations=None,
                       epistemics=None)]

        increase_amount_verbs = ['ExpressionControl-positive',
                                 'MolSynthesis-positive']
        decrease_amount_verbs = ['ExpressionControl-negative',
                                 'MolSynthesis-negative']

        if relation.verb in increase_amount_verbs:
            self.statements.append(
                                   IncreaseAmount(subj, obj, evidence=ev)
                                  )
        elif relation.verb in decrease_amount_verbs:
            self.statements.append(
                                   DecreaseAmount(subj, obj, evidence=ev)
                                  )
        elif relation.verb == 'ProtModification':
            if last_relation is not None:
                self.modification_examples[last_relation.verb] += 1

            if last_relation is None:
                # We cannot make a statement unless we have more fine-grained
                # information on the relation type from a preceding
                # unnormalized SVO
                return

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
                # This verb is not handled
                return

            self.statements.append(statement_type(subj, obj, evidence=ev))

        elif relation.verb == 'Binding':
            self.statements.append(
                                   Complex([subj, obj], evidence=ev)
                                  )

    def agent_from_entity(self, relation, entity_id):
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
        else:
            # Convert the URN grounding to an INDRA grounding
            entity = relation.entities[entity_id]
            db_refs = urn_to_db_refs(entity.urn)
            db_refs['TEXT'] = entity.match_text

        return Agent(db_refs['TEXT'], db_refs=db_refs)
