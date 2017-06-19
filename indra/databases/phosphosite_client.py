from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
from os.path import dirname, abspath, join
from collections import namedtuple, defaultdict
from indra.util import read_unicode_csv

PhosphoSite = namedtuple('PhosphoSite',
                         ['GENE', 'PROTEIN', 'ACC_ID', 'HU_CHR_LOC', 'MOD_RSD',
                          'SITE_GRP_ID', 'ORGANISM', 'MW_kD', 'DOMAIN',
                          'SITE_7_AA', 'LT_LIT', 'MS_LIT', 'MS_CST', 'CST_CAT'])

_data_by_up = None
_data_by_site_grp = None

def _read_phospho_site_dataset():
    global _data_by_up
    global _data_by_site_grp
    if _data_by_up is None or _data_by_site_grp is None:
        phosphosite_data_file = join(dirname(abspath(__file__)),
                                 '../resources/Phosphorylation_site_dataset.tsv')
        reader = read_unicode_csv(phosphosite_data_file, delimiter='\t',
                                  skiprows=4)
        # Build up a dict by protein
        data_by_up = {}
        data_by_site_grp = defaultdict(list)
        for row in reader:
            site = PhosphoSite(*row)
            data_by_up[site.PROTEIN] = site
            data_by_site_grp[site.SITE_GRP_ID].append(site)
        _data_by_up = data_by_up
        _data_by_site_grp = data_by_site_grp
    return (_data_by_up, _data_by_site_grp)

if __name__ == '__main__':
    (data_by_up, data_by_site_grp) = _read_phospho_site_dataset()


