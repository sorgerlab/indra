from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str
import logging
from os.path import dirname, abspath, join
from collections import namedtuple, defaultdict
from indra.util import read_unicode_csv

logger = logging.getLogger('phosphosite')

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
        data_by_up = defaultdict(dict)
        data_by_site_grp = defaultdict(list)
        for row in reader:
            site = PhosphoSite(*row)
            res_pos = site.MOD_RSD.split('-')[0]
            data_by_up[site.ACC_ID][res_pos] = site
            data_by_site_grp[site.SITE_GRP_ID].append(site)
        _data_by_up = data_by_up
        _data_by_site_grp = data_by_site_grp
    return (_data_by_up, _data_by_site_grp)

def map_to_human_site(up_id, mod_res, mod_pos):
    (data_by_up, data_by_site_grp) = _read_phospho_site_dataset()
    sites_for_up = data_by_up.get(up_id)
    # No info in Phosphosite for this Uniprot ID
    if sites_for_up is None:
        return None
    site_info = sites_for_up.get('%s%s' % (mod_res, mod_pos))
    if not site_info:
        return None
    # Lookup site group
    site_grp_list = data_by_site_grp.get(site_info.SITE_GRP_ID)
    # If an empty list, then return None (is unlikely to happen)
    if not site_grp_list:
        return None
    # Check for a human protein in the list
    human_sites = [s for s in site_grp_list if s.ORGANISM == 'human']
    if not human_sites:
        return None
    if len(human_sites) > 1:
        logger.info("More than one human site found for %s, site group %s" %
                    (up_id, site_info.SITE_GRP_ID))
    human_site = human_sites[0].MOD_RSD.split('-')[0]
    human_res = human_site[0]
    assert human_res == mod_res
    human_pos = human_site[1:]
    return human_pos

if __name__ == '__main__':
    (data_by_up, data_by_site_grp) = _read_phospho_site_dataset()


