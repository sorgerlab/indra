import pickle
import logging
import numpy as np
import pandas as pd
from collections import Counter
from typing import Union, Sequence, Optional, List
from sklearn.base import BaseEstimator
from indra.statements import Evidence, Statement, get_all_descendants
from indra.belief import BeliefScorer, check_extra_evidence, \
                         get_stmt_evidence, SimpleScorer


logger = logging.getLogger(__name__)


class SklearnScorer(BeliefScorer):
    """Use a pre-trained Sklearn classifier to predict belief scores.

    An implementing instance of this base class has two personalities: as a
    subclass of BeliefScorer, it implements the functions required by the
    BeliefEngine, `score_statements` and `check_prior_probs`. It also behaves
    like an sklearn model by composition, implementing methods `fit`,
    `predict`, `predict_proba`, and `predict_log_proba`, which are passed
    through to an internal sklearn model.

    A key role of this wrapper class is to implement the preprocessing of
    statement properties into a feature matrix in a standard way, so that
    a classifier trained on one corpus of statement data will still work when
    used on another corpus.

    Implementing subclasses must implement at least one of the methods for
    building the feature matrix, `stmts_to_matrix` or `df_to_matrix`.

    Parameters
    ----------
    model :
        Any instance of a classifier object supporting the methods `fit`,
        `predict_proba`, `predict`, and `predict_log_proba`.
    """
    def __init__(
        self,
        model: BaseEstimator,
    ):
        self.model = model

    def check_prior_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Empty implementation for now."""
        pass

    def score_statements(
        self,
        statements: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> Sequence[float]:
        return self.predict_proba(statements, extra_evidence)[:, 1]

    def stmts_to_matrix(
        self,
        stmts: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> np.ndarray:
        """Convert a list of Statements to a feature matrix."""
        raise NotImplementedError('Need to implement the stmts_to_matrix '
                                  'method')

    def df_to_matrix(
        self,
        df: pd.DataFrame,
    ) -> np.ndarray:
        """Convert a statement DataFrame to a feature matrix."""
        raise NotImplementedError('Need to implement the df_to_matrix '
                                   'method')

    def to_matrix(self,
        stmt_data: Union[np.ndarray, Sequence[Statement], pd.DataFrame],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> np.ndarray:
        """Get stmt feature matrix by calling appropriate method.

        If `stmt_data` is already a matrix (e.g., obtained after performing a
        train/test split on a matrix generated for a full statement corpus), it
        is returned directly; if a DataFrame of Statement metadata,
        `self.df_to_matrix` is called; if a list of Statements,
        `self.stmts_to_matrix` is called.

        Parameters
        ----------
        stmt_data :
            Statement content to be used to generate a feature matrix.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).

        Returns
        -------
        Feature matrix for the statement data.
        """
        # If we got a Numpy array, just use it!
        if isinstance(stmt_data, np.ndarray):
            stmt_arr = stmt_data
        # Otherwise check if we have a dataframe or a list of statements
        # and call the appropriate *_to_matrix method
        elif isinstance(stmt_data, pd.DataFrame):
            if extra_evidence is not None:
                raise NotImplementedError(
                   'extra_evidence cannot be used with a statement DataFrame.')
            stmt_arr = self.df_to_matrix(stmt_data)
        # Check if stmt_data is a list/tuple (i.e., of Statements):
        elif isinstance(stmt_data, (list, tuple)):
            # Check that the first entry is a Statement
            if not isinstance(stmt_data[0], Statement):
                raise ValueError('stmt_data must contain Statements.')
            stmt_arr = self.stmts_to_matrix(stmt_data, extra_evidence)
        # If it's something else, error
        else:
            raise TypeError(f'stmt_data is type {type(stmt_data)}: '
                            'must be a numpy array, DataFrame, or '
                            'list/tuple of Statements')
        return stmt_arr

    def fit(self,
        stmt_data: Union[np.ndarray, Sequence[Statement], pd.DataFrame],
        y_arr: Sequence[float],
        extra_evidence: Optional[List[List[Evidence]]] = None,
        *args,
        **kwargs,
    ):
        """Preprocess stmt data and run sklearn model `fit` method.

        Additional `args` and `kwargs` are passed to the `fit` method of the
        wrapped sklearn model.

        Parameters
        ----------
        stmt_data :
            Statement content to be used to generate a feature matrix.
        y_arr :
            Class values for the statements (e.g., a vector of 0s and 1s
            indicating correct or incorrect).
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).
        """
        # Check dimensions of stmts (x) and y_arr
        if len(stmt_data) != len(y_arr):
            raise ValueError("Number of stmts/rows must match length of y_arr.")
        # Get the data matrix based on the stmt list or stmt DataFrame
        stmt_arr = self.to_matrix(stmt_data, extra_evidence)
        # Call the fit method of the internal sklearn model
        self.model.fit(stmt_arr, y_arr, *args, **kwargs)

    def predict_proba(
        self,
        stmt_data: Union[np.ndarray, Sequence[Statement], pd.DataFrame],
        extra_evidence: Optional[List[List[Evidence]]] = None,
        *args,
        **kwargs,
    ) -> np.ndarray:
        """Preprocess stmt data and run sklearn model `predict_proba` method.

        Additional `args` and `kwargs` are passed to the `predict_proba` method
        of the wrapped sklearn model.

        Parameters
        ----------
        stmt_data :
            Statement content to be used to generate a feature matrix.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).
        """
        # Call the prediction method of the internal sklearn model
        stmt_arr = self.to_matrix(stmt_data, extra_evidence)
        return self.model.predict_proba(stmt_arr, *args, **kwargs)

    def predict(
        self,
        stmt_data: Union[np.ndarray, Sequence[Statement], pd.DataFrame],
        extra_evidence: Optional[List[List[Evidence]]] = None,
        *args,
        **kwargs,
    ) -> np.ndarray:
        """Preprocess stmt data and run sklearn model `predict` method.

        Additional `args` and `kwargs` are passed to the `predict` method of
        the wrapped sklearn model.

        Parameters
        ----------
        stmt_data :
            Statement content to be used to generate a feature matrix.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).
        """
        stmt_arr = self.to_matrix(stmt_data, extra_evidence)
        return self.model.predict(stmt_arr, *args, **kwargs)

    def predict_log_proba(
        self,
        stmt_data: Union[np.ndarray, Sequence[Statement], pd.DataFrame],
        extra_evidence: Optional[List[List[Evidence]]] = None,
        *args,
        **kwargs,
    ) -> np.ndarray:
        """Preprocess stmt data and run sklearn model `predict_log_proba`.

        Additional `args` and `kwargs` are passed to the `predict` method of
        the wrapped sklearn model.

        Parameters
        ----------
        stmt_data :
            Statement content to be used to generate a feature matrix.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).
        """
        stmt_arr = self.to_matrix(stmt_data, extra_evidence)
        return self.model.predict_log_proba(stmt_arr, *args, **kwargs)


class CountsScorer(SklearnScorer):
    """Belief model learned from evidence counts and other stmt properties.

    If using a DataFrame for Statement data, it should have the following
    columns:

    * `stmt_type`
    * `source_counts`

    Alternatively, if the DataFrame doesn't have a `source_counts` column, it
    should have columns with names matching the sources in `self.source_list`.

    Parameters
    ----------
    model :
        Any instance of a classifier object supporting the methods `fit`,
        `predict_proba`, `predict`, and `predict_log_proba`.
    source_list :
        List of strings denoting the evidence sources (evidence.source_api
        values) to be used for prediction.
    include_more_specific :
        If True, will add extra columns to the statement data matrix for the
        source counts drawn from more specific evidences; if use_num_pmids is
        True, will also add an additional column for the number of PMIDs from
        more specific evidences. If False, these columns will not be included
        even if the `extra_evidence` argument is passed to the
        `stmts_to_matrix` method. This is to ensure that the featurization of
        statements is consistent between training and prediction.
    use_stmt_type :
        Whether to include statement type as a feature.
    use_num_members :
        Whether to include a feature denoting the number of members of the
        statement. Primarily for stratifying belief predictions about Complex
        statements with more than two members. Cannot be used for statement
        data passed in as a DataFrame.
    use_num_pmids :
        Whether to include a feature for the total number of unique PMIDs
        supporting each statement. Cannot be used for statement passed in as a
        DataFrame.
    use_promoter :
        Whether to include a feature giving the fraction of evidence (0 to 1)
        containing the (case-insensitive) word "promoter". Tends to improve
        misclassification of Complex statements that actually refer to
        protein-DNA binding.
    use_avg_evidence_len :
        Whether to include a feature giving the average evidence sentence
        length (in space-separated tokens).

    Example
    -------
    .. code-block:: python

        from sklearn.linear_model import LogisticRegression
        clf = LogisticRegression()
        all_stmt_sources = CountsScorer.get_all_sources(stmts)
        scorer = CountsScorer(clf, all_stmt_sources, use_stmt_type=True,
                              use_num_pmids=True)
        scorer.fit(stmts, y_arr)
        be = BeliefEngine(scorer)
        be.set_hierarchy_probs(stmts)
    """
    def __init__(
        self,
        model: BaseEstimator,
        source_list: List[str],
        include_more_specific: bool = False,
        use_stmt_type: bool = False,
        use_num_members: bool = False,
        use_num_pmids: bool = False,
        use_promoter: bool = False,
        use_avg_evidence_len: bool = False,
    ):
        # Call superclass constructor to store the model
        super(CountsScorer, self).__init__(model)
        self.source_list = source_list
        self.include_more_specific = include_more_specific
        self.use_stmt_type = use_stmt_type
        self.use_num_members = use_num_members
        self.use_num_pmids = use_num_pmids
        self.use_promoter = use_promoter
        self.use_avg_evidence_len = use_avg_evidence_len
        # Build dictionary mapping INDRA Statement types to integers
        if use_stmt_type:
            all_stmt_types = get_all_descendants(Statement)
            self.stmt_type_map = {t.__name__: ix
                                  for ix, t in enumerate(all_stmt_types)}

    @staticmethod
    def get_all_sources(
        stmts: Sequence[Statement],
        include_more_specific: bool = True,
        include_less_specific: bool = True,
    ) -> List[str]:
        """Get a list of all the source_apis supporting the given statements.

        Useful for determining the set of sources to be used for fitting
        and prediction.

        Parameters
        ----------
        stmts :
            A list of INDRA Statements to collect source APIs for.
        include_more_specific :
            If True (default), then includes the source APIs for the more
            specific statements in the `supports` attribute of each statement.
        include_less_specific :
            If True (default), then includes the source APIs for the less
            specific statements in the `supported_by` attribute of each
            statement.

        Returns
        -------
        A list of (unique) source_apis found in the set of statements.
        """
        stmt_sources = set([ev.source_api for s in stmts for ev in s.evidence])
        if include_more_specific:
            stmt_sources.update([ev.source_api
                                 for stmt in stmts
                                 for supp_stmt in stmt.supports
                                 for ev in supp_stmt.evidence])
        if include_less_specific:
            stmt_sources.update([ev.source_api
                                 for stmt in stmts
                                 for supp_by_stmt in stmt.supported_by
                                 for ev in supp_by_stmt.evidence])
        return list(stmt_sources)

    def stmts_to_matrix(
        self,
        stmts: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> np.ndarray:
        """Convert a list of Statements to a feature matrix.

        Features are encoded as follows:

        * One column for every source listed in `self.source_list`, containing
          the number of statement evidences from that source. If
          `self.include_more_specific` is True and `extra_evidence` is
          provided, these are used in combination with the Statement's own
          evidence in determining source counts.
        * If `self.use_stmt_type` is set, statement type is included via
          one-hot encoding, with one column for each statement type.
        * If `self.use_num_members` is set, a column is added for the number
          of agents in the Statement.
        * If `self.use_num_pmids` is set, a column is added with the total
          total number of unique PMIDs supporting the Statement.  If
          `extra_evidence` is provided, these are used in combination with the
          Statement's own evidence in determining the number of PMIDs.

        Parameters
        ----------
        stmts :
            A list or tuple of INDRA Statements to be used to generate a
            feature matrix.
        extra_evidence :
            A list corresponding to the given list of statements, where
            each entry is a list of Evidence objects providing additional
            support for the corresponding statement (i.e., Evidences that
            aren't already included in the Statement's own evidence list).

        Returns
        -------
        Feature matrix for the statement data.
        """
        # Check arguments for including more specific evidences
        if self.include_more_specific and extra_evidence is None:
            logger.info("CountScorer is set to include_more_specific "
                        "evidences but no extra_evidence was included.")
            extra_evidence = [[] for stmt in stmts]
        elif not self.include_more_specific and extra_evidence is not None:
            logger.warning("extra_evidence was included but CountScorer "
                           "instance is not set to include_more_specific "
                           "evidences so extra_evidence will be ignored.")
        # Check our list of extra evidences
        check_extra_evidence(extra_evidence, len(stmts))
        # Add categorical features and collect source_apis
        cat_features = []
        stmt_sources = set()
        for ix, stmt in enumerate(stmts):
            # Collect all source_apis from stmt evidences
            dir_pmids = set()
            promoter_ct = 0
            evidence_lens = []
            for ev in stmt.evidence:
                stmt_sources.add(ev.source_api)
                dir_pmids.add(ev.pmid)
                if ev.text is not None:
                    evidence_lens.append(len(ev.text.split()))
                    if 'promoter' in ev.text.lower():
                        promoter_ct += 1

            indir_pmids = set()
            if self.include_more_specific:
                for ev in extra_evidence[ix]:
                    stmt_sources.add(ev.source_api)
                    indir_pmids.add(ev.pmid)
            # Collect non-source count features (e.g. type) from stmts
            feature_row = []
            # One-hot encoding of stmt type
            if self.use_stmt_type:
                stmt_type_ix = self.stmt_type_map[type(stmt).__name__]
                type_features = [1 if ix == stmt_type_ix else 0
                                 for ix in range(len(self.stmt_type_map))]
                feature_row.extend(type_features)
            # Add field for number of members
            if self.use_num_members:
                feature_row.append(len(stmt.agent_list()))
            # Add field with number of unique PMIDs
            if self.use_num_pmids:
                feature_row.append(len(dir_pmids))
                if self.include_more_specific:
                    feature_row.append(len(indir_pmids))
            # Add a field specifying the percentage of evidences containing
            # the word "promoter":
            if self.use_promoter:
                promoter_pct = promoter_ct / len(stmt.evidence) \
                                        if len(stmt.evidence) > 0 else 0
                feature_row.append(promoter_pct)
            # Add a field giving length of the sentence in words
            if self.use_avg_evidence_len:
                avg_evidence_len = np.mean(evidence_lens) if evidence_lens \
                                                          else 0
                feature_row.append(avg_evidence_len)
            # Only add a feature row if we're using some of the features.
            if feature_row:
                cat_features.append(feature_row)

        # Before proceeding, check whether all source_apis are in
        # source_list
        if stmt_sources.difference(set(self.source_list)):
            logger.info("source_list does not include all source_apis "
                             "in the statement data.")
        # Get source count features
        # If we have extra_evidence, we double the source count features
        if extra_evidence is None:
            num_cols = len(self.source_list)
        else:
            num_cols = len(self.source_list) * 2
        num_rows = len(stmts)
        x_arr = np.zeros((num_rows, num_cols))
        for stmt_ix, stmt in enumerate(stmts):
            # Source from the stmt itself
            direct_sources = [ev.source_api for ev in stmt.evidence]
            dsrc_ctr = Counter(direct_sources)
            for src_ix, src in enumerate(self.source_list):
                x_arr[stmt_ix, src_ix] = dsrc_ctr.get(src, 0)
            # Get indirect evidences
            if self.include_more_specific:
                indirect_sources = [ev.source_api
                                    for ev in extra_evidence[stmt_ix]]
                idsrc_ctr = Counter(indirect_sources)
                for src_ix, src in enumerate(self.source_list):
                    x_arr[stmt_ix, src_ix + len(self.source_list)] = \
                                            idsrc_ctr.get(src, 0)
        # If we have any categorical features, turn them into an array and
        # add them to matrix
        if cat_features:
            cat_arr = np.array(cat_features)
            x_arr = np.hstack((x_arr, cat_arr))
        return x_arr


    def df_to_matrix(
        self,
        df: pd.DataFrame,
    ) -> np.ndarray:
        """Convert a DataFrame of statement data to a feature matrix.

        Based on information available in a DataFrame of statement data, this
        implementation uses only source counts and statement type in building a
        feature matrix, and will raise a ValueError if either
        `self.use_num_members` or `self.use_num_pmids` is set.

        Features are encoded as follows:

        * One column for every source listed in `self.source_list`, containing
          the number of statement evidences from that source. If
          `extra_evidence` is provided, these are used in combination with the
          Statement's own evidence in determining source counts.
        * If `self.use_stmt_type` is set, statement type is included via
          one-hot encoding, with one column for each statement type.

        Parameters
        ----------
        df :
            A pandas DataFrame with statement metadata. It should have columns
            `stmt_type` and `source_counts`; alternatively, if it doesn't have
            a `source_counts` column, it should have columns with names
            matching the sources in `self.source_list`.

        Returns
        -------
        Feature matrix for the statement data.
        """
        required_cols = {'stmt_type'}
        # Currently, statement DataFrames are not expected to contain
        # number of members or num_pmids as a data column, hence we raise a
        # ValueError if either of these are set
        if self.use_num_members:
            raise ValueError('use_num_members not supported for statement '
                             'DataFrames.')
        if self.use_num_pmids:
            raise ValueError('use_num_pmids not supported for statement '
                             'DataFrames.')
        # Make sure that the dataframe contains at least all of the above
        # columns
        if not required_cols.issubset(set(df.columns)):
            raise ValueError('Statement DataFrame is missing required '
                             'columns.')
        # Check for the source_counts column. If it's there, we're good
        if 'source_counts' in df.columns:
            has_sc_col = True
        # If it's not, make sure that we have columns named for sources in
        # self.source_list:
        else:
            has_sc_col = False
            for source in self.source_list:
                if source not in df.columns:
                    raise ValueError(f'Expected column "{source}" not in the '
                                      'given statement DataFrame')

        # Add categorical features and collect source_apis
        cat_features = []
        stmt_sources = set()
        # For every statement entry in the dataframe...
        for rowtup in df.itertuples():
            # Collect statement sources
            # ...if there's a source_counts col with dicts
            if has_sc_col:
                stmt_sources |= set(rowtup.source_counts.keys())
            # Collect non-source count features (e.g. type) from stmts
            feature_row = []
            # One-hot encoding of stmt type
            if self.use_stmt_type:
                stmt_type_ix = self.stmt_type_map[rowtup.stmt_type]
                type_features = [1 if ix == stmt_type_ix else 0
                                 for ix in range(len(self.stmt_type_map))]
                feature_row.extend(type_features)
            # Only add a feature row if we're using some of the features.
            if feature_row:
                cat_features.append(feature_row)

        # Before proceeding, check whether all source_apis are in
        # source_list. If we don't have a source_counts dict, we don't look
        # for columns beyond the sources in the source list, and we are
        # guaranteed to have all of them because of the check performed above
        source_diff = stmt_sources.difference(set(self.source_list))
        if has_sc_col and source_diff:
            logger.warning("source_list does not include all source_apis "
                           f"in the statement data: {str(source_diff)}")

        # Get source count features
        num_cols = len(self.source_list)
        num_rows = len(df)
        x_arr = np.zeros((num_rows, num_cols))
        for stmt_ix, rowtup in enumerate(df.itertuples()):
            for src_ix, src in enumerate(self.source_list):
                # Get counts from the source_count dictionary
                if has_sc_col:
                    x_arr[stmt_ix, src_ix] = rowtup.source_counts.get(src, 0)
                # ...or get counts from named source column
                else:
                    x_arr[stmt_ix, src_ix] = rowtup._asdict()[src]

        # If we have any categorical features, turn them into an array and
        # add them to matrix
        if cat_features:
            cat_arr = np.array(cat_features)
            x_arr = np.hstack((x_arr, cat_arr))
        return x_arr


class HybridScorer(BeliefScorer):
    """TODO: Docstring"""
    def __init__(
        self,
        counts_scorer: CountsScorer,
        simple_scorer: SimpleScorer,
    ):
        self.counts_scorer = counts_scorer
        self.simple_scorer = simple_scorer

    def check_prior_probs(
        self,
        statements: Sequence[Statement],
    ) -> None:
        """Empty implementation for now."""
        # TODO Check that for a given list of statements, that all sources in
        # the statements are either in the CountsScorer source_list, or in
        # the prior dict from the SimpleScorer
        pass

    def score_statements(
        self,
        statements: Sequence[Statement],
        extra_evidence: Optional[List[List[Evidence]]] = None,
    ) -> Sequence[float]:
        # Get beliefs from the sklearn model, using the sources in the
        # CountScorer source_list as features
        skl_beliefs = self.counts_scorer.predict_proba(statements,
                                                       extra_evidence)[:, 1]
        skl_sources = self.counts_scorer.source_list
        hybrid_beliefs = []
        # Iterate over the statements...
        for ix, stmt in enumerate(statements):
            # ...get both the statement's own evidence and the more-specific
            # (extra) evidences
            all_evidence = get_stmt_evidence(stmt, ix, extra_evidence)
            # Next, filter out any evidences that have sources in the skl
            # model source list, leaving behind the rest. At the same time,
            # record whether we've found any sources the skl model source list.
            filt_evidence = []
            has_skl_source = False
            for ev in all_evidence:
                if ev.source_api in skl_sources:
                    has_skl_source = True
                else:
                    filt_evidence.append(ev)
            # Get the simple belief
            simple_bel = self.simple_scorer.score_evidence_list(filt_evidence)
            # Calculate hybrid belief: the probability that all sources, both
            # those evaluated by the sklearn model and the simplescorer, are
            # not jointly incorrect. If there are no sources from the skl
            # model list, we set the skl belief to 0 so the probability comes
            # only from the simple scorer
            skl_bel = skl_beliefs[ix] if has_skl_source else 0
            hybrid_bel = 1 - (1 - skl_bel) * (1 - simple_bel)
            hybrid_beliefs.append(hybrid_bel)
        return hybrid_beliefs


