from io import StringIO
import pandas
import requests
from . import SimpleScorer, BayesianScorer


def load_eidos_curation_table():
    """Return a pandas table of Eidos curation data."""
    url = 'https://raw.githubusercontent.com/clulab/eidos/master/' + \
        'src/main/resources/org/clulab/wm/eidos/english/confidence/' + \
        'rule_summary.tsv'
    # Load the table of scores from the URL above into a data frame
    res = StringIO(requests.get(url).text)
    table = pandas.read_table(res, sep='\t')
    # Drop the last "Grant total" row
    table = table.drop(table.index[len(table)-1])
    return table


def get_eidos_bayesian_scorer():
    """Return a BayesianScorer based on Eidos curation counts."""
    table = load_eidos_curation_table()
    subtype_counts = {'eidos': {r: [c, i] for r, c, i in
                              zip(table['RULE'], table['Num correct'],
                                  table['Num incorrect'])}}

    scorer = BayesianScorer(prior_counts={}, subtype_counts=subtype_counts)
    return scorer


def get_eidos_scorer():
    """Return a SimpleScorer based on Eidos curated precision estimates."""
    table = load_eidos_curation_table()

    # Get the overall precision
    total_num = table['COUNT of RULE'].sum()
    weighted_sum = table['COUNT of RULE'].dot(table['% correct'])
    precision = weighted_sum / total_num
    # We have to divide this into a random and systematic component, for now
    # in an ad-hoc manner
    syst_error = 0.05
    rand_error = 1 - precision - syst_error
    prior_probs = {'rand': {'eidos': rand_error}, 'syst': {'eidos': syst_error}}

    # Get a dict of rule-specific errors.
    subtype_probs = {'eidos':
                     {k: 1.0-min(v, 0.95)-syst_error for k, v
                      in zip(table['RULE'], table['% correct'])}}
    print(subtype_probs)
    scorer = SimpleScorer(prior_probs, subtype_probs)
    return scorer
