from urllib.parse import unquote

import re
import os
import glob
import time
import shutil
import tempfile
import logging
from math import floor
import lxml.etree
import collections

from indra.databases import go_client, mesh_client
from indra.statements import *
from indra.databases.chebi_client import get_chebi_id_from_cas, \
    get_chebi_name_from_id
from indra.databases.hgnc_client import get_hgnc_from_entrez, get_uniprot_id, \
        get_hgnc_name
from indra.util import read_unicode_csv
from indra.sources.reach.processor import ReachProcessor, Site

from .fix_csxml_character_encoding import fix_character_encoding

logger = logging.getLogger(__name__)


MedscanEntity = collections.namedtuple('MedscanEntity', ['name', 'urn', 'type',
                                                         'properties',
                                                         'ch_start', 'ch_end'])


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


def _fix_different_refs(a1, a2, ref_key):
    if all(ref_key in a.db_refs for a in [a1, a2]) \
           and a1.db_refs[ref_key] != a2.db_refs[ref_key]:
        a1.name = a1.db_refs[ref_key]
        a2.name = a2.db_refs[ref_key]
        return True
    return False


def _is_statement_in_list(new_stmt, old_stmt_list):
    """Return True of given statement is equivalent to on in a list

    Determines whether the statement is equivalent to any statement in the
    given list of statements, with equivalency determined by Statement's
    equals method.

    Parameters
    ----------
    new_stmt : indra.statements.Statement
        The statement to compare with
    old_stmt_list : list[indra.statements.Statement]
        The statement list whose entries we compare with statement

    Returns
    -------
    in_list : bool
        True if statement is equivalent to any statements in the list
    """
    for old_stmt in old_stmt_list:
        if old_stmt.equals(new_stmt):
            return True
        elif old_stmt.evidence_equals(new_stmt) and old_stmt.matches(new_stmt):
            # If we're comparing a complex, make sure the agents are sorted.
            if isinstance(new_stmt, Complex):
                agent_pairs = zip(old_stmt.sorted_members(),
                                  new_stmt.sorted_members())
            else:
                agent_pairs = zip(old_stmt.agent_list(), new_stmt.agent_list())

            # Compare agent-by-agent.
            for ag_old, ag_new in agent_pairs:
                s_old = set(ag_old.db_refs.items())
                s_new = set(ag_new.db_refs.items())

                # If they're equal this isn't the one we're interested in.
                if s_old == s_new:
                    continue

                # If the new statement has nothing new to offer, just ignore it
                if s_old > s_new:
                    return True

                # If the new statement does have something new, add it to the
                # existing statement. And then ignore it.
                if s_new > s_old:
                    ag_old.db_refs.update(ag_new.db_refs)
                    return True

                # If this is a case where different CHEBI ids were mapped to
                # the same entity, set the agent name to the CHEBI id.
                if _fix_different_refs(ag_old, ag_new, 'CHEBI'):
                    # Check to make sure the newly described statement does
                    # not match anything.
                    return _is_statement_in_list(new_stmt, old_stmt_list)

                # If this is a case, like above, but with UMLS IDs, do the same
                # thing as above. This will likely never be improved.
                if _fix_different_refs(ag_old, ag_new, 'UMLS'):
                    # Check to make sure the newly described statement does
                    # not match anything.
                    return _is_statement_in_list(new_stmt, old_stmt_list)

                logger.warning("Found an unexpected kind of duplicate. "
                               "Ignoring it.")
                return True

            # This means all the agents matched, which can happen if the
            # original issue was the ordering of agents in a Complex.
            return True

        elif old_stmt.get_hash(True, True) == new_stmt.get_hash(True, True):
            # Check to see if we can improve the annotation of the existing
            # statement.
            e_old = old_stmt.evidence[0]
            e_new = new_stmt.evidence[0]
            if e_old.annotations['last_verb'] is None:
                e_old.annotations['last_verb'] = e_new.annotations['last_verb']

            # If the evidence is "the same", modulo annotations, just ignore it
            if e_old.get_source_hash(True) == e_new.get_source_hash(True):
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


# These normalized verbs are mapped to IncreaseAmount statements
INCREASE_AMOUNT_VERBS = ['ExpressionControl-positive',
                         'MolSynthesis-positive',
                         'CellExpression',
                         'QuantitativeChange-positive',
                         'PromoterBinding']

# These normalized verbs are mapped to DecreaseAmount statements
DECREASE_AMOUNT_VERBS = ['ExpressionControl-negative',
                         'MolSynthesis-negative',
                         'miRNAEffect-negative',
                         'QuantitativeChange-negative']

# These normalized verbs are mapped to Activation statements (indirect)
ACTIVATION_VERBS = ['UnknownRegulation-positive',
                    'Regulation-positive']
# These normalized verbs are mapped to Activation statements (direct)
D_ACTIVATION_VERBS = ['DirectRegulation-positive',
                      'DirectRegulation-positive--direct interaction']
# All activation verbs
ALL_ACTIVATION_VERBS = ACTIVATION_VERBS + D_ACTIVATION_VERBS

# These normalized verbs are mapped to Inhibition statements (indirect)
INHIBITION_VERBS = ['UnknownRegulation-negative',
                    'Regulation-negative']
# These normalized verbs are mapped to Inhibition statements (direct)
D_INHIBITION_VERBS = ['DirectRegulation-negative',
                      'DirectRegulation-negative--direct interaction']
# All inhibition verbs
ALL_INHIBITION_VERBS = INHIBITION_VERBS + D_INHIBITION_VERBS

PMID_PATT = re.compile('info:pmid/(\d+)')


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
        self.last_site_info_in_sentence = None
        self.files_processed = 0
        self._gen = None
        self._tmp_dir = None
        self._pmids_handled = set()
        self._sentences_handled = set()
        self.__f = None
        return

    def iter_statements(self, populate=True):
        if self._gen is None and not self.statements:
            raise InputError("No generator has been initialized. Use "
                             "`process_directory` or `process_file` first.")
        if self.statements and not self._gen:
            for stmt in self.statements:
                yield stmt
        else:
            for stmt in self._gen:
                if populate:
                    self.statements.append(stmt)
                yield stmt

    def process_directory(self, directory_name, lazy=False):
        # Process each file
        glob_pattern = os.path.join(directory_name, '*.csxml')
        files = glob.glob(glob_pattern)
        self._gen = self._iter_over_files(files)
        if not lazy:
            for stmt in self._gen:
                self.statements.append(stmt)

        return

    def _iter_over_files(self, files):
        # Create temporary directory into which to put the csxml files with
        # normalized character encodings
        self.__tmp_dir = tempfile.mkdtemp('indra_medscan_processor')
        tmp_file = os.path.join(self.__tmp_dir, 'fixed_char_encoding')

        num_files = float(len(files))
        percent_done = 0
        start_time_s = time.time()

        logger.info("%d files to read" % int(num_files))
        for filename in files:
            logger.info('Processing %s' % filename)
            fix_character_encoding(filename, tmp_file)
            with open(tmp_file, 'rb') as self.__f:
                for stmt in self._iter_through_csxml_file_from_handle():
                    yield stmt

            percent_done_now = floor(100.0 * self.files_processed / num_files)
            if percent_done_now > percent_done:
                percent_done = percent_done_now
                ellapsed_s = time.time() - start_time_s
                ellapsed_min = ellapsed_s / 60.0

                msg = 'Processed %d of %d files (%f%% complete, %f minutes)' % \
                      (self.files_processed, num_files, percent_done,
                       ellapsed_min)
                logger.info(msg)

        # Delete the temporary directory
        shutil.rmtree(self.__tmp_dir)
        return

    def process_csxml_file(self, filename, interval=None, lazy=False):
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
        filename : string
            The path to a Medscan csxml file.
        interval : (start, end) or None
            Select the interval of documents to read, starting with the
            `start`th document and ending before the `end`th document. If
            either is None, the value is considered undefined. If the value
            exceeds the bounds of available documents, it will simply be
            ignored.
        lazy : bool
            If True, only create a generator which can be used by the
            `get_statements` method. If True, populate the statements list now.
        """
        if interval is None:
            interval = (None, None)

        tmp_fname = tempfile.mktemp(os.path.basename(filename))
        fix_character_encoding(filename, tmp_fname)

        self.__f = open(tmp_fname, 'rb')
        self._gen = self._iter_through_csxml_file_from_handle(*interval)
        if not lazy:
            for stmt in self._gen:
                self.statements.append(stmt)
        return

    def _iter_through_csxml_file_from_handle(self, start=None, stop=None):
        pmid = None
        sec = None
        tagged_sent = None
        doc_idx = 0
        entities = {}
        match_text = None
        in_prop = False
        last_relation = None
        property_entities = []
        property_name = None

        # Go through the document again and extract statements
        good_relations = []
        skipping_doc = False
        skipping_sent = False
        for event, elem in lxml.etree.iterparse(self.__f,
                                                events=('start', 'end'),
                                                encoding='utf-8',
                                                recover=True):
            if elem.tag in ['attr', 'toks']:
                continue
            # If opening up a new doc, set the PMID
            if event == 'start' and elem.tag == 'doc':
                if start is not None and doc_idx < start:
                    logger.info("Skipping document number %d." % doc_idx)
                    skipping_doc = True
                    continue

                if stop is not None and doc_idx >= stop:
                    logger.info("Reach the end of the allocated docs.")
                    break

                uri = elem.attrib.get('uri')
                re_pmid = PMID_PATT.match(uri)
                if re_pmid is None:
                    logger.warning("Could not extract pmid from: %s." % uri)
                    skipping_doc = True

                pmid = re_pmid.group(1)
                pmid_num = int(pmid)
                if pmid_num in self._pmids_handled:
                    logger.warning("Skipping repeated pmid: %s from %s."
                                   % (pmid, self.__f.name))
                    skipping_doc = True
            # If getting a section, set the section type
            elif event == 'start' and elem.tag == 'sec' and not skipping_doc:
                sec = elem.attrib.get('type')
            # Set the sentence context
            elif event == 'start' and elem.tag == 'sent' and not skipping_doc:
                tagged_sent = elem.attrib.get('msrc')
                h = hash(tagged_sent)
                if h in self._sentences_handled:
                    skipping_sent = True
                    continue
                skipping_sent = False

                # Reset last_relation between sentences, since we will only be
                # interested in the relation immediately preceding a CONTROL
                # statement but within the same sentence.
                last_relation = None

                entities = {}
            elif event == 'end' and elem.tag == 'sent' and not skipping_doc \
                    and not skipping_sent:
                # End of sentence; deduplicate and copy statements from this
                # sentence to the main statements list

                for s in self.sentence_statements:
                    yield s
                self.sentence_statements = []
                self._sentences_handled.add(h)
                good_relations = []

                # Reset site info
                self.last_site_info_in_sentence = None
            elif event == 'start' and elem.tag == 'match' and not skipping_doc\
                    and not skipping_sent:
                match_text = elem.attrib.get('chars')
                match_start = int(elem.attrib.get('coff'))
                match_end = int(elem.attrib.get('clen')) + match_start
            elif event == 'start' and elem.tag == 'entity' \
                    and not skipping_doc and not skipping_sent:
                if not in_prop:
                    ent_id = elem.attrib['msid']
                    ent_urn = elem.attrib.get('urn')
                    ent_type = elem.attrib['type']
                    entities[ent_id] = MedscanEntity(match_text, ent_urn,
                                                     ent_type, {},
                                                     match_start, match_end)
                else:
                    ent_type = elem.attrib['type']
                    ent_urn = elem.attrib['urn']
                    ent_name = elem.attrib['name']
                    property_entities.append(MedscanEntity(ent_name, ent_urn,
                                                           ent_type, None,
                                                           None, None))
            elif event == 'start' and elem.tag == 'svo' and not skipping_doc \
                    and not skipping_sent:
                subj = elem.attrib.get('subj')
                verb = elem.attrib.get('verb')
                obj = elem.attrib.get('obj')
                svo_type = elem.attrib.get('type')

                # Aggregate information about the relation
                relation = MedscanRelation(pmid=pmid, sec=sec, uri=uri,
                                           tagged_sentence=tagged_sent,
                                           entities=entities, subj=subj,
                                           verb=verb, obj=obj,
                                           svo_type=svo_type)
                if svo_type == 'CONTROL':
                    good_relations.append(relation)
                    self.process_relation(relation, last_relation)
                else:
                    # Sometimes a CONTROL SVO can be after an unnormalized SVO
                    # that is a more specific but less uniform version of the
                    # same extracted statement.
                    last_relation = relation
            elif event == 'start' and elem.tag == 'prop' and not skipping_doc \
                    and not skipping_sent:
                in_prop = True
                property_name = elem.attrib.get('name')
                property_entities = []
            elif event == 'end' and elem.tag == 'prop' and not skipping_doc \
                    and not skipping_sent:
                in_prop = False
                entities[ent_id].properties[property_name] = property_entities
            elif event == 'end' and elem.tag == 'doc':
                doc_idx += 1
                # Give a status update
                if doc_idx % 100 == 0:
                    logger.info("Processed %d documents" % doc_idx)
                self._pmids_handled.add(pmid_num)
                self._sentences_handled = set()

            # Solution for memory leak found here:
            # https://stackoverflow.com/questions/12160418/why-is-lxml-etree-iterparse-eating-up-all-my-memory?lq=1
            elem.clear()

        self.files_processed += 1
        self.__f.close()
        return

    def _add_statement(self, stmt):
        if not _is_statement_in_list(stmt, self.sentence_statements):
            self.sentence_statements.append(stmt)
        return

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
        subj_res = self.agent_from_entity(relation, relation.subj)
        obj_res = self.agent_from_entity(relation, relation.obj)
        if subj_res is None or obj_res is None:
            # Don't extract a statement if the subject or object cannot
            # be resolved
            return
        subj, subj_bounds = subj_res
        obj, obj_bounds = obj_res

        # Make evidence object
        untagged_sentence = _untag_sentence(relation.tagged_sentence)
        if last_relation:
            last_verb = last_relation.verb
        else:
            last_verb = None
        # Get the entity information with the character coordinates
        annotations = {'verb': relation.verb, 'last_verb': last_verb,
                       'agents': {'coords': [subj_bounds, obj_bounds]}}
        epistemics = dict()
        epistemics['direct'] = False  # Overridden later if needed
        ev = [Evidence(source_api='medscan', source_id=relation.uri,
                       pmid=relation.pmid, text=untagged_sentence,
                       annotations=annotations, epistemics=epistemics)]

        if relation.verb in INCREASE_AMOUNT_VERBS:
            # If the normalized verb corresponds to an IncreaseAmount statement
            # then make one
            self._add_statement(IncreaseAmount(subj, obj, evidence=ev))
        elif relation.verb in DECREASE_AMOUNT_VERBS:
            # If the normalized verb corresponds to a DecreaseAmount statement
            # then make one
            self._add_statement(DecreaseAmount(subj, obj, evidence=ev))
        elif relation.verb in ALL_ACTIVATION_VERBS:
            # If the normalized verb corresponds to an Activation statement,
            # then make one
            if relation.verb in D_ACTIVATION_VERBS:
                ev[0].epistemics['direction'] = True
            self._add_statement(Activation(subj, obj, evidence=ev))
        elif relation.verb in ALL_INHIBITION_VERBS:
            # If the normalized verb corresponds to an Inhibition statement,
            # then make one
            if relation.verb in D_INHIBITION_VERBS:
                ev[0].epistemics['direct'] = True
            self._add_statement(Inhibition(subj, obj, evidence=ev))

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
                    self._add_statement(s)
            else:
                self._add_statement(statement_type(subj, obj, evidence=ev))

        elif relation.verb == 'Binding':
            # The Binding normalized verb corresponds to the INDRA Complex
            # statement.
            self._add_statement(
                                   Complex([subj, obj], evidence=ev)
                                  )
        elif relation.verb == 'ProtModification-negative':
            pass  # TODO? These occur so infrequently so maybe not worth it
        elif relation.verb == 'Regulation-unknown':
            pass  # TODO? These occur so infrequently so maybe not worth it
        elif relation.verb == 'StateEffect-positive':
            pass
            # self._add_statement(
            #                       ActiveForm(subj, obj, evidence=ev)
            #                      )
            # TODO: disabling for now, since not sure whether we should set
            # the is_active flag
        elif relation.verb == 'StateEffect':
            self.last_site_info_in_sentence = \
                    ProteinSiteInfo(site_text=subj.name,
                                    object_text=obj.db_refs['TEXT'])
        return

    def agent_from_entity(self, relation, entity_id):
        """Create a (potentially grounded) INDRA Agent object from a given
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
        self.num_entities += 1

        entity_id = _extract_id(entity_id)

        if entity_id not in relation.entities and \
                entity_id not in tags:
            # Could not find the entity in either the list of grounded
            # entities of the items tagged in the sentence. Happens for
            # a very small percentage of the dataset.
            self.num_entities_not_found += 1
            return None

        if entity_id not in relation.entities:
            # The entity is not in the grounded entity list
            # Instead, make an ungrounded entity, with TEXT corresponding to
            # the words with the given entity id tagged in the sentence.
            entity_data = tags[entity_id]
            db_refs = {'TEXT': entity_data['text']}
            ag = Agent(normalize_medscan_name(db_refs['TEXT']),
                       db_refs=db_refs)
            return ag, entity_data['bounds']
        else:
            entity = relation.entities[entity_id]
            bounds = (entity.ch_start, entity.ch_end)

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
                            ag = Agent(normalize_medscan_name(agent_name),
                                       db_refs=db_refs, mutations=[cond])
                            return ag, bounds
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
                    ag = Agent(normalize_medscan_name(agent_name),
                               db_refs=db_refs, mods=[cond])
                    return ag, bounds

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
                    ag = Agent(normalize_medscan_name(agent_name),
                               db_refs=db_refs, mods=[cond])
                    return ag, bounds

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

                ag = Agent(normalize_medscan_name(agent_name),
                           db_refs=db_refs)
                return ag, bounds


class MedscanRelation(object):
    """A structure representing the information contained in a Medscan
    SVO xml element as well as associated entities and properties.

    Attributes
    ----------
    pmid : str
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
    def __init__(self, pmid, uri, sec, entities, tagged_sentence, subj, verb, obj,
                 svo_type):
        self.pmid = pmid
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


MOD_PATT = re.compile('([A-Za-z])+([0-9]+)')


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
    m = MOD_PATT.match(s)
    assert m is not None
    return m.groups()


MUT_PATT = re.compile('([A-Za-z]+)([0-9]+)([A-Za-z]+)')


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
    m = MUT_PATT.match(s)
    if m is None:
        # Mutation string does not fit this pattern, other patterns not
        # currently supported
        return None, None, None
    else:
        return m.groups()


URN_PATT = re.compile('urn:([^:]+):([^:]+)')


def _urn_to_db_refs(urn):
    """Converts a Medscan URN to an INDRA db_refs dictionary with grounding
    information.

    Parameters
    ----------
    urn : str
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

    m = URN_PATT.match(urn)
    if m is None:
        return None, None

    urn_type, urn_id = m.groups()

    db_refs = {}
    db_name = None

    # TODO: support more types of URNs
    if urn_type == 'agi-cas':
        # Identifier is CAS, convert to CHEBI
        chebi_id = get_chebi_id_from_cas(urn_id)
        if chebi_id:
            db_refs['CHEBI'] = chebi_id
            db_name = get_chebi_name_from_id(chebi_id)
    elif urn_type == 'agi-llid':
        # This is an Entrez ID, convert to HGNC
        hgnc_id = get_hgnc_from_entrez(urn_id)
        if hgnc_id is not None:
            db_refs['HGNC'] = hgnc_id

            # Convert the HGNC ID to a Uniprot ID
            uniprot_id = get_uniprot_id(hgnc_id)
            if uniprot_id is not None:
                db_refs['UP'] = uniprot_id

            # Try to lookup HGNC name; if it's available, set it to the
            # agent name
            db_name = get_hgnc_name(hgnc_id)
    elif urn_type in ['agi-meshdis', 'agi-ncimorgan', 'agi-ncimtissue',
                      'agi-ncimcelltype']:
        if urn_id.startswith('C') and urn_id[1:].isdigit():
            # Identifier is probably UMLS
            db_refs['UMLS'] = urn_id
        else:
            # Identifier is MESH
            urn_mesh_name = unquote(urn_id)
            mesh_id, mesh_name = mesh_client.get_mesh_id_name(urn_mesh_name,
                                                              offline=True)
            if mesh_id:
                db_refs['MESH'] = mesh_id
                db_name = mesh_name
            else:
                db_name = urn_mesh_name
    elif urn_type == 'agi-gocomplex':
        # Identifier is GO
        db_refs['GO'] = 'GO:%s' % urn_id
    elif urn_type == 'agi-go':
        # Identifier is GO
        db_refs['GO'] = 'GO:%s' % urn_id

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
    elif 'GO' in db_refs:
        db_name = go_client.get_go_label(db_refs['GO'])

    return db_refs, db_name


TAG_PATT = re.compile('ID{([0-9,]+)=([^}]+)}')
JUNK_PATT = re.compile('(CONTEXT|GLOSSARY){[^}]+}+')
ID_PATT = re.compile('ID\\{([0-9]+)\\}')


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
    matches = ID_PATT.match(id_string)
    assert matches is not None
    return matches.group(1)


def _untag_sentence(tagged_sentence):
    """Removes all tags in the sentence, returning the original sentence
    without Medscan annotations.

    Parameters
    ----------
    tagged_sentence : str
        The tagged sentence

    Returns
    -------
    untagged_sentence : str
        Sentence with tags and annotations stripped out
    """
    untagged_sentence = TAG_PATT.sub('\\2', tagged_sentence)
    clean_sentence = JUNK_PATT.sub('', untagged_sentence)
    return clean_sentence.strip()


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
    untagged_sentence = _untag_sentence(tagged_sentence)
    decluttered_sentence = JUNK_PATT.sub('', tagged_sentence)
    tags = {}

    # Iteratively look for all matches of this pattern
    endpos = 0
    while True:
        match = TAG_PATT.search(decluttered_sentence, pos=endpos)
        if not match:
            break
        endpos = match.end()
        text = match.group(2)
        text = text.replace('CONTEXT', '')
        text = text.replace('GLOSSARY', '')
        text = text.strip()
        start = untagged_sentence.index(text)
        stop = start + len(text)

        tag_key = match.group(1)
        if ',' in tag_key:
            for sub_key in tag_key.split(','):
                if sub_key == '0':
                    continue
                tags[sub_key] = {'text': text, 'bounds': (start, stop)}
        else:
            tags[tag_key] = {'text': text, 'bounds': (start, stop)}
    return tags
