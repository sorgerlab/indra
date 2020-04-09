import yaml
import logging
from indra.belief.wm_scorer import get_eidos_bayesian_scorer
from indra.sources.eidos import reground_texts
from indra.tools import assemble_corpus as ac
from indra.belief import BeliefEngine
from . import file_defaults, default_key_base, InvalidCorpusError
from .corpus import Corpus

logger = logging.getLogger(__name__)


class LiveCurator(object):
    """Class coordinating the real-time curation of a corpus of Statements.

    Parameters
    ----------
    scorer : indra.belief.BeliefScorer
        A scorer object to use for the curation
    corpora : dict[str, Corpus]
        A dictionary mapping corpus IDs to Corpus objects.
    """

    def __init__(self, scorer=None, corpora=None, eidos_url=None,
                 ont_manager=None):
        self.corpora = corpora if corpora else {}
        self.scorer = scorer if scorer else get_eidos_bayesian_scorer()
        self.ont_manager = ont_manager
        self.eidos_url = eidos_url

    # TODO: generalize this to other kinds of scorers
    def reset_scorer(self):
        """Reset the scorer used for curation."""
        logger.info('Resetting the scorer')
        self.scorer = get_eidos_bayesian_scorer()
        for corpus_id, corpus in self.corpora.items():
            corpus.curations = {}

    def get_corpus(self, corpus_id, check_s3=True, use_cache=True):
        """Return a corpus given an ID.

        If the corpus ID cannot be found, an InvalidCorpusError is raised.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus to return.
        check_s3 : bool
            If True, look on S3 for the corpus if it's not currently loaded.
            Default: True
        use_cache : bool
            If True, look in local cache before trying to find corpus on s3.
            If True while check_s3 if False, this option will be ignored.
            Default: False.

        Returns
        -------
        Corpus
            The corpus with the given ID.
        """
        logger.info('Getting corpus "%s"' % corpus_id)
        corpus = self.corpora.get(corpus_id)
        if corpus:
            logger.info('Found corpus loaded in memory')
        if check_s3 and corpus is None:
            logger.info('Corpus not loaded, looking on S3')
            corpus = Corpus.load_from_s3(s3key=corpus_id,
                                         force_s3_reload=not use_cache,
                                         raise_exc=True)
            logger.info('Adding corpus to loaded corpora')
            self.corpora[corpus_id] = corpus

            # Run update beliefs. The belief update needs to be inside this
            # if statement to avoid infinite recursion
            beliefs = self.update_beliefs(corpus_id)
        elif corpus is None:
            raise InvalidCorpusError

        return corpus

    def get_curations(self, corpus_id, reader):
        """Download curations for corpus id filtered to reader

        Parameters
        ----------
        corpus_id: str
            The ID of the corpus to download curations from
        reader : str
            The name of the reader to filter to. Has to be among valid
            reader names of 'all'.

        Returns
        -------
        dict
            A dict containing the requested curations
        """
        logger.info('Getting curations for corpus %s' % corpus_id)
        corpus = self.get_corpus(corpus_id, check_s3=True, use_cache=True)
        corpus_curations = corpus.get_curations(corpus_id,
                                                look_in_cache=True)
        # Get all statements that have curations
        curated_stmts = {}
        for uuid in corpus_curations:
            curated_stmts[uuid] = corpus.statements[uuid]
        if reader and reader != 'all':
            # Filter out statements and curations that don't contain material
            # from provided reader (in source api of statement)
            filtered_curations = {}
            filtered_stmts = {}
            for uuid, stmt in curated_stmts.items():
                # Check if any of the evidences are from the provided reader
                for ev in stmt.evidence:
                    if ev.source_api == reader.lower():
                        filtered_stmts[uuid] = stmt
                        filtered_curations[uuid] = corpus_curations[uuid]
                        break
            data = {'curations': filtered_curations,
                    'statements': {uuid: st.to_json() for uuid, st in
                                   filtered_stmts.items()}}
        else:
            data = {'curations': corpus_curations,
                    'statements': {uuid: st.to_json() for uuid, st in
                                   curated_stmts.items()}}
        return data

    def submit_curation(self, corpus_id, curations):
        """Submit correct/incorrect curations fo a given corpus.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus to which the curations apply.
        curations : dict
            A dict of curations with keys corresponding to Statement UUIDs and
            values corresponding to correct/incorrect feedback.
        """
        logger.info('Submitting curations for corpus "%s"' % corpus_id)
        corpus = self.get_corpus(corpus_id, check_s3=True, use_cache=True)
        # Start tabulating the curation counts
        prior_counts = {}
        subtype_counts = {}
        # Take each curation from the input
        for uuid, correct in curations.items():
            # Save the curation in the corpus
            # TODO: handle already existing curation
            stmt = corpus.statements.get(uuid)
            if stmt is None:
                logger.warning('%s is not in the corpus.' % uuid)
                continue
            corpus.curations[uuid] = correct
            # Now take all the evidences of the statement and assume that
            # they follow the correctness of the curation and contribute to
            # counts for their sources
            for ev in stmt.evidence:
                # Make the index in the curation count list
                idx = 0 if correct else 1
                extraction_rule = ev.annotations.get('found_by')
                # If there is no extraction rule then we just score the source
                if not extraction_rule:
                    try:
                        prior_counts[ev.source_api][idx] += 1
                    except KeyError:
                        prior_counts[ev.source_api] = [0, 0]
                        prior_counts[ev.source_api][idx] += 1
                # Otherwise we score the specific extraction rule
                else:
                    try:
                        subtype_counts[ev.source_api][extraction_rule][idx] \
                            += 1
                    except KeyError:
                        if ev.source_api not in subtype_counts:
                            subtype_counts[ev.source_api] = {}
                        subtype_counts[ev.source_api][extraction_rule] = [0, 0]
                        subtype_counts[ev.source_api][extraction_rule][idx] \
                            += 1
        # Finally, we update the scorer with the new curation counts
        self.scorer.update_counts(prior_counts, subtype_counts)

    def save_curation(self, corpus_id, save_to_cache=True):
        """Save the current state of curations for a corpus given its ID

        If the corpus ID cannot be found, an InvalidCorpusError is raised.

        Parameters
        ----------
        corpus_id : str
            the ID of the corpus to save the
        save_to_cache : bool
            If True, also save the current curation to the local cache.
            Default: True.
        """
        # Do NOT use cache or S3 when getting the corpus, otherwise it will
        # overwrite the current corpus
        logger.info('Saving curations for corpus "%s"' % corpus_id)
        corpus = self.get_corpus(corpus_id, check_s3=False, use_cache=False)
        corpus.upload_curations(corpus_id, save_to_cache=save_to_cache)

    def update_metadata(self, corpus_id, meta_data, save_to_cache=True):
        """Update the meta data for a given corpus

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus to update the meta data for
        meta_data : dict
            A json compatible dict containing the meta data
        save_to_cache : bool
            If True, also update the local cache of the meta data dict.
            Default: True.
        """
        logger.info('Updating meta data for corpus "%s"' % corpus_id)
        corpus = self.get_corpus(corpus_id, check_s3=True, use_cache=True)

        # Loop and add/overwrite meta data key value pairs
        for k, v in meta_data.items():
            corpus.meta_data[k] = v

        if save_to_cache:
            meta_file_key = '%s/%s/%s.json' % (default_key_base,
                                               corpus_id,
                                               file_defaults['meta'])
            corpus._save_to_cache(meta=meta_file_key)

    def update_beliefs(self, corpus_id):
        """Return updated belief scores for a given corpus.

        Parameters
        ----------
        corpus_id : str
            The ID of the corpus for which beliefs are to be updated.

        Returns
        -------
        dict
            A dictionary of belief scores with keys corresponding to Statement
            UUIDs and values to new belief scores.
        """
        logger.info('Updating beliefs for corpus "%s"' % corpus_id)
        # TODO check which options are appropriate for get_corpus
        corpus = self.get_corpus(corpus_id)
        be = BeliefEngine(self.scorer)
        stmts = list(corpus.statements.values())
        be.set_prior_probs(stmts)
        # Here we set beliefs based on actual curation
        for uuid, correct in corpus.curations.items():
            stmt = corpus.statements.get(uuid)
            if stmt is None:
                logger.warning('%s is not in the corpus.' % uuid)
                continue
            stmt.belief = correct
        belief_dict = {st.uuid: st.belief for st in stmts}
        return belief_dict

    def update_groundings(self, corpus_id):
        # TODO check which options are appropriate for get_corpus
        logger.info('Updating groundings for corpus "%s"' % corpus_id)
        corpus = self.get_corpus(corpus_id)

        # Send the latest ontology and list of concept texts to Eidos
        yaml_str = yaml.dump(self.ont_manager.yaml_root)
        concepts = []
        for stmt in corpus.raw_statements:
            for concept in stmt.agent_list():
                concept_txt = concept.db_refs.get('TEXT')
                concepts.append(concept_txt)
        groundings = reground_texts(concepts, yaml_str,
                                    webservice=self.eidos_url)
        # Update the corpus with new groundings
        idx = 0
        for stmt in corpus.raw_statements:
            for concept in stmt.agent_list():
                concept.db_refs['WM'] = groundings[idx]
                idx += 1
        assembled_statements = default_assembly(corpus.raw_statements)
        corpus.statements = {s.uuid: s for s in assembled_statements}
        return assembled_statements


def default_assembly(stmts):
    from indra.belief.wm_scorer import get_eidos_scorer
    from indra.preassembler.hierarchy_manager import get_wm_hierarchies
    hm = get_wm_hierarchies()
    scorer = get_eidos_scorer()
    stmts = ac.run_preassembly(stmts, belief_scorer=scorer,
                               return_toplevel=True,
                               flatten_evidence=True,
                               normalize_equivalences=True,
                               normalize_opposites=True,
                               normalize_ns='WM',
                               flatten_evidence_collect_from='supported_by',
                               poolsize=4,
                               hierarchies=hm)
    stmts = ac.merge_groundings(stmts)
    stmts = ac.merge_deltas(stmts)
    stmts = ac.standardize_names_groundings(stmts)
    return stmts
