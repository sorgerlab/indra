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
        data_by_up = defaultdict(lambda: defaultdict(list))
        data_by_site_grp = defaultdict(list)
        for row in reader:
            site = PhosphoSite(*row)
            res_pos = site.MOD_RSD.split('-')[0]
            base_acc_id = site.ACC_ID.split('-')[0]
            data_by_up[base_acc_id][res_pos].append(site)
            data_by_site_grp[site.SITE_GRP_ID].append(site)
        _data_by_up = data_by_up
        _data_by_site_grp = data_by_site_grp
    return (_data_by_up, _data_by_site_grp)

def map_to_human_site(up_id, mod_res, mod_pos):
    (data_by_up, data_by_site_grp) = _read_phospho_site_dataset()
    sites_for_up = data_by_up.get(up_id)
    # No info in Phosphosite for this Uniprot ID
    if not sites_for_up:
        return None
    site_info_list = sites_for_up.get('%s%s' % (mod_res, mod_pos))
    # If this site doesn't exist for this protein, will return an empty list
    if not site_info_list:
        return None
    # At this point site_info_list contains a list of PhosphoSite objects
    # for the given protein with an entry for the given site, however, they
    # may be different isoforms; they may also not contain the reference
    # isoform; or they may only contain a single isoform
    # If it's a single isoform, take it!
    if len(site_info_list) == 1:
        site_info = site_info_list[0]
    # If there is more than one entry for this site, take the one from the
    # reference sequence, if it's in there (for example, a site that is present
    # in a mouse isoform as well as the mouse reference sequence)
    elif len(site_info_list) > 0:
        logger.debug('More than one entry in PhosphoSite for %s, site %s%s' %
                    (up_id, mod_res, mod_pos))
        ref_site_info = None
        for si in site_info_list:
            logger.debug('\tSite info: acc_id %s, site_grp_id %s' %
                        (si.ACC_ID, si.SITE_GRP_ID))
            if si.ACC_ID == up_id:
                logger.info('\tFound entry matching reference acc id')
                ref_site_info = si
        if ref_site_info is None:
            site_info = site_info_list[0]
            logger.info('Reference sequence match not found, choosing first '
                         'entry, acc_id %s site_grp %s' %
                         (site_info.ACC_ID, site_info.SITE_GRP_ID))
        else:
            site_info = ref_site_info
    # Lookup site group
    site_grp_list = data_by_site_grp.get(site_info.SITE_GRP_ID)
    # If an empty list, then return None (is unlikely to happen)
    if not site_grp_list:
        return None
    # Check for a human protein in the list
    human_sites = [s for s in site_grp_list if s.ORGANISM == 'human']
    if not human_sites:
        return None
    # If there are multiple isoforms, choose the base one
    # (no hyphen in the accession ID)
    if len(human_sites) > 1:
        # We assume that multiple human sites only arise from multiple isoforms
        # of the same protein, which will share an accession ID, and that
        # only one of these will be the reference sequence (with no hyphen)
        base_id_sites = [site for site in human_sites
                         if site.ACC_ID.find('-') == -1]
        if base_id_sites:
            assert len(base_id_sites) == 1, 'Should only be one human ref seq'
            human_site = base_id_sites[0]
        # There is no base ID site, i.e., all mapped sites are for specific
        # isoforms only, so skip it!
        else:
            logger.info('Human isoform matches, but no human ref seq match '
                        'for %s %s %s; not mapping' % (up_id, mod_res, mod_pos))
            return None
    # If there is only one human site, take it
    else:
        human_site = human_sites[0]
    human_site_str = human_site.MOD_RSD.split('-')[0]
    human_res = human_site_str[0]
    assert human_res == mod_res
    human_pos = human_site_str[1:]
    return human_pos

if __name__ == '__main__':
    (data_by_up, data_by_site_grp) = _read_phospho_site_dataset()


