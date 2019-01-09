from io import StringIO
import pandas
import requests
from . import SimpleScorer


class BayesianScorer(SimpleScorer):
    def __init__(self, prior_counts, subtype_counts):
        self.prior_counts = prior_counts
        self.subtype_counts = subtype_counts
        self._update_probs()

    def update_probs(self):
        prior_probs = {source: (p / (n + p)) for source, (p, n) in
                       self.prior_counts.items()}
        # FIXME: this will need to be fixed
        subtype_probs = {source: {rule: (p / (n + p)} for rule, (p, n)
                         in st.items() for st, rule in
                         self.subtype_counts.items()}
        super()._update_probs(prior_probs, subtype_probs)

    def update_counts(self, counts):
        for source, (pos, neg) in counts.items():
            # FIXME: handle subtypes
            self.prior_counts[source][0] += pos
            self.prior_counts[source][1] += neg


def get_eidos_counts():
    url = 'https://raw.githubusercontent.com/clulab/eidos/master/' + \
        'src/main/resources/org/clulab/wm/eidos/english/confidence/' + \
        'rule_summary.tsv'
    # Load the table of scores from the URL above into a data frame
    res = StringIO(requests.get(url).text)
    table = pandas.read_table(res, sep='\t')
    # Drop the last "Grant total" row
    table = table.drop(table.index[len(table)-1])

    # Get the overall precision
    total_corr = table['Num correct'].sum()
    total_incorr = table['Num incorrect'].sum()
    precision = total_corr / (total_corr + total_incorr)
    prior_probs = {'rand': {'eidos': rand_error},
                   'syst': {'eidos': syst_error}}

    prior_counts = {'eidos': {r: [c, i] for r, c, i
                              zip(table['RULE'], table['Num correct'],
                              table['Num incorrect'])}}

    # We have to divide this into a random and systematic component, for now
    # in an ad-hoc manner
    syst_error = 0.05
    rand_error = 1 - precision - syst_error
    prior_probs = {'rand': {'eidos': rand_error}, 'syst': {'eidos': syst_error}}

    # Get a dict of rule-specific errors.
    subtype_probs = {'eidos':
                     {k: 1.0-min(v, 0.95)-syst_error for k, v
                      in zip(table['RULE'], table['% correct'])}}

def get_eidos_scorer():
    url = 'https://raw.githubusercontent.com/clulab/eidos/master/' + \
        'src/main/resources/org/clulab/wm/eidos/english/confidence/' + \
        'rule_summary.tsv'

    # Load the table of scores from the URL above into a data frame
    res = StringIO(requests.get(url).text)
    table = pandas.read_table(res, sep='\t')
    # Drop the last "Grant total" row
    table = table.drop(table.index[len(table)-1])

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

    scorer = SimpleScorer(prior_probs, subtype_probs)
    return scorer
