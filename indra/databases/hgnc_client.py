import os
import re
import logging
from collections import defaultdict

import requests
from typing import Set, Union
import xml.etree.ElementTree as ET
from functools import lru_cache

from indra.util import read_unicode_csv, UnicodeXMLTreeBuilder as UTB
from indra.resources import get_resource_path

logger = logging.getLogger(__name__)


hgnc_url = 'http://rest.genenames.org/fetch/'


def get_uniprot_id(hgnc_id):
    """Return the UniProt ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Note that the HGNC ID is a number that is
        passed as a string. It is not the same as the HGNC gene symbol.

    Returns
    -------
    uniprot_id : str
        The UniProt ID corresponding to the given HGNC ID.
    """
    uniprot_id = uniprot_ids.get(hgnc_id)
    # The lookup can yield an empty string. Instead return None.
    if not uniprot_id:
        return None
    return uniprot_id


def get_entrez_id(hgnc_id):
    """Return the Entrez ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Note that the HGNC ID is a number that is
        passed as a string. It is not the same as the HGNC gene symbol.

    Returns
    -------
    entrez_id : str
        The Entrez ID corresponding to the given HGNC ID.
    """
    entrez_id = entrez_ids.get(hgnc_id)
    # The lookup can yield an empty string. Instead return None.
    if not entrez_id:
        return None
    return entrez_id


def get_hgnc_from_entrez(entrez_id):
    """Return the HGNC ID corresponding to the given Entrez ID.

    Parameters
    ----------
    entrez_id : str
        The Entrez ID to be converted, a number passed as a string.

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given Entrez ID.
    """
    hgnc_id = entrez_ids_reverse.get(entrez_id)
    return hgnc_id


def get_ensembl_id(hgnc_id):
    """Return the Ensembl ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Note that the HGNC ID is a number that is
        passed as a string. It is not the same as the HGNC gene symbol.

    Returns
    -------
    ensembl_id : str
        The Ensembl ID corresponding to the given HGNC ID.
    """
    return ensembl_ids.get(hgnc_id)


def get_hgnc_from_ensembl(ensembl_id):
    """Return the HGNC ID corresponding to the given Ensembl ID.

    Parameters
    ----------
    ensembl_id : str
        The Ensembl ID to be converted, a number passed as a string.

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given Ensembl ID.
    """
    return ensembl_ids_reverse.get(ensembl_id)


def get_hgnc_name(hgnc_id):
    """Return the HGNC symbol corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted.

    Returns
    -------
    hgnc_name : str
        The HGNC symbol corresponding to the given HGNC ID.
    """
    try:
        hgnc_name = hgnc_names[hgnc_id]
    except KeyError:
        xml_tree = get_hgnc_entry(hgnc_id)
        if xml_tree is None:
            return None
        hgnc_name_tag =\
            xml_tree.find("result/doc/str[@name='symbol']")
        if hgnc_name_tag is None:
            return None
        hgnc_name = hgnc_name_tag.text.strip()
    return hgnc_name


def get_hgnc_id(hgnc_name):
    """Return the HGNC ID corresponding to the given HGNC symbol.

    Parameters
    ----------
    hgnc_name : str
        The HGNC symbol to be converted. Example: BRAF

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given HGNC symbol.
    """
    return hgnc_ids.get(hgnc_name)


def get_current_hgnc_id(hgnc_name):
    """Return HGNC ID(s) corresponding to a current or outdated HGNC symbol.

    Parameters
    ----------
    hgnc_name : str
        The HGNC symbol to be converted, possibly an outdated symbol.

    Returns
    -------
    str or list of str or None
        If there is a single HGNC ID corresponding to the given current or
        outdated HGNC symbol, that ID is returned as a string. If the symbol
        is outdated and maps to multiple current IDs, a list of these
        IDs is returned. If the given name doesn't correspond to either
        a current or an outdated HGNC symbol, None is returned.
    """
    hgnc_id = get_hgnc_id(hgnc_name)
    if hgnc_id:
        return hgnc_id
    hgnc_id = prev_sym_map.get(hgnc_name)
    return hgnc_id


def get_hgnc_from_mouse(mgi_id):
    """Return the HGNC ID corresponding to the given MGI mouse gene ID.

    Parameters
    ----------
    mgi_id : str
        The MGI ID to be converted. Example: "2444934"

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given MGI ID.
    """
    if mgi_id.startswith('MGI:'):
        mgi_id = mgi_id[4:]
    return mouse_map.get(mgi_id)


def get_hgnc_from_rat(rgd_id):
    """Return the HGNC ID corresponding to the given RGD rat gene ID.

    Parameters
    ----------
    rgd_id : str
        The RGD ID to be converted. Example: "1564928"

    Returns
    -------
    hgnc_id : str
        The HGNC ID corresponding to the given RGD ID.
    """
    if rgd_id.startswith('RGD:'):
        rgd_id = rgd_id[4:]
    return rat_map.get(rgd_id)


def get_rat_id(hgnc_id):
    """Return the RGD rat ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Example: "1097"

    Returns
    -------
    rgd_id : str
        The RGD ID corresponding to the given HGNC ID.
    """
    for k, v in rat_map.items():
        if v == hgnc_id:
            return k


def get_mouse_id(hgnc_id):
    """Return the MGI mouse ID corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted. Example: "1097"

    Returns
    -------
    mgi_id : str
        The MGI ID corresponding to the given HGNC ID.
    """
    for k, v in mouse_map.items():
        if v == hgnc_id:
            return k


@lru_cache(maxsize=1000)
def get_hgnc_entry(hgnc_id):
    """Return the HGNC entry for the given HGNC ID from the web service.

    Parameters
    ----------
    hgnc_id : str
        The HGNC ID to be converted.

    Returns
    -------
    xml_tree : ElementTree
        The XML ElementTree corresponding to the entry for the
        given HGNC ID.
    """
    url = hgnc_url + 'hgnc_id/%s' % hgnc_id
    headers = {'Accept': '*/*'}
    res = requests.get(url, headers=headers)
    if not res.status_code == 200:
        return None
    xml_tree = ET.XML(res.content, parser=UTB())
    return xml_tree


def get_gene_type(hgnc_id: str) -> Union[str, None]:
    """Return the locus type of the genve with the given HGNC ID.

    See more under Locus type at
    https://www.genenames.org/help/symbol-report/#!/#tocAnchor-1-2

    Parameters
    ----------
    hgnc_id :
        The HGNC ID of the gene to get the locus type of.

    Returns
    -------
    :
        The locus type of the given gene.
    """
    return gene_type.get(hgnc_id)


def is_kinase(gene_name):
    """Return True if the given gene name is a kinase.

    Parameters
    ----------
    gene_name : str
        The HGNC gene symbol corresponding to the protein.

    Returns
    -------
    bool
        True if the given gene name corresponds to a kinase, False otherwise.
    """
    return gene_name in kinases


def is_transcription_factor(gene_name):
    """Return True if the given gene name is a transcription factor.

    Parameters
    ----------
    gene_name : str
        The HGNC gene symbol corresponding to the protein.

    Returns
    -------
    bool
        True if the given gene name corresponds to a transcription factor,
        False otherwise.
    """
    return gene_name in tfs


def is_phosphatase(gene_name):
    """Return True if the given gene name is a phosphatase.

    Parameters
    ----------
    gene_name : str
        The HGNC gene symbol corresponding to the protein.

    Returns
    -------
    bool
        True if the given gene name corresponds to a phosphatase,
        False otherwise.
    """
    return gene_name in phosphatases


def get_enzymes(hgnc_id: str) -> Set[str]:
    """Return the EC codes corresponding to the given HGNC ID.

    Parameters
    ----------
    hgnc_id :
        The HGNC ID to be converted.

    Returns
    -------
    :
        A set of EC codes
    """
    return hgnc_to_enzymes.get(hgnc_id, set())


def get_hgncs_from_enzyme(ec_code: str) -> Set[str]:
    """Return the HGNC ids associated with a given enzyme.

    Parameters
    ----------
    ec_code :
        The EC code (e.g., 2.4.1.228)

    Returns
    -------
    :
        A set of HGNC identifiers
    """
    return enzyme_to_hgncs.get(ec_code, set())


def _read_hgnc_maps():
    hgnc_file = get_resource_path("hgnc_entries.tsv")
    csv_rows = read_unicode_csv(hgnc_file, delimiter='\t', encoding='utf-8')
    hgnc_names = {}
    hgnc_ids = {}
    hgnc_withdrawn = []
    uniprot_ids = {}
    entrez_ids = {}
    entrez_ids_reverse = {}
    mouse_map = {}
    rat_map = {}
    prev_sym_map = {}
    ensembl_ids = {}
    ensembl_ids_reverse = {}
    hgnc_withdrawn_new_ids = {}
    gene_types = {}
    hgnc_to_enzymes = defaultdict(set)
    enzyme_to_hgncs = defaultdict(set)
    # Skip the header
    next(csv_rows)
    for row in csv_rows:
        hgnc_id = row[0][5:]
        hgnc_status = row[3]
        if hgnc_status in {'Approved', 'Entry Withdrawn'}:
            hgnc_name = row[1]
            hgnc_names[hgnc_id] = hgnc_name
            # Note that withdrawn entries don't overlap with approved
            # entries at this point so it's safe to add mappings for
            # withdrawn names
            hgnc_ids[hgnc_name] = hgnc_id
        elif hgnc_status == 'Symbol Withdrawn':
            descr = row[2]
            m = re.match(r'symbol withdrawn, see \[HGNC:(?: ?)(\d+)\]', descr)
            new_id = m.groups()[0]
            hgnc_withdrawn.append(hgnc_id)
            hgnc_withdrawn_new_ids[hgnc_id] = new_id
        # Uniprot
        uniprot_id = row[6]
        if uniprot_id:
            uniprot_ids[hgnc_id] = uniprot_id
        # Entrez
        entrez_id = row[5]
        if entrez_id:
            entrez_ids[hgnc_id] = entrez_id
            entrez_ids_reverse[entrez_id] = hgnc_id
        # Mouse
        mgi_id = row[7]
        if mgi_id:
            mgi_ids = mgi_id.split(', ')
            for mgi_id in mgi_ids:
                if mgi_id.startswith('MGI:'):
                    mgi_id = mgi_id[4:]
                mouse_map[mgi_id] = hgnc_id
        # Rat
        rgd_id = row[8]
        if rgd_id:
            rgd_ids = rgd_id.split(', ')
            for rgd_id in rgd_ids:
                if rgd_id.startswith('RGD:'):
                    rgd_id = rgd_id[4:]
                rat_map[rgd_id] = hgnc_id
        # Previous symbols
        prev_sym_entry = row[9]
        if prev_sym_entry:
            prev_syms = prev_sym_entry.split(', ')
            for prev_sym in prev_syms:
                # If we already mapped this previous symbol to another ID
                if prev_sym in prev_sym_map:
                    # If we already have a list here, we just extend it
                    if isinstance(prev_sym_map[prev_sym], list):
                        prev_sym_map[prev_sym].append(hgnc_id)
                    # Otherwise we create a list and start it with the two
                    # IDs we know the symbol is mapped to
                    else:
                        prev_sym_map[prev_sym] = [prev_sym_map[prev_sym],
                                                  hgnc_id]
                # Otherwise we just make a string entry here
                else:
                    prev_sym_map[prev_sym] = hgnc_id
        ensembl_id = row[10]
        # Ensembl IDs
        if ensembl_id:
            ensembl_ids[hgnc_id] = ensembl_id
            ensembl_ids_reverse[ensembl_id] = hgnc_id
        gene_type = row[11]
        if gene_type:
            gene_types[hgnc_id] = gene_type
        enyzyme_ids = row[12]
        if enyzyme_ids:
            for enzyme_id in enyzyme_ids.split(", "):
                hgnc_to_enzymes[hgnc_id].add(enzyme_id)
                enzyme_to_hgncs[enzyme_id].add(hgnc_id)

    for old_id, new_id in hgnc_withdrawn_new_ids.items():
        hgnc_names[old_id] = hgnc_names[new_id]

    return (
        hgnc_names, hgnc_ids, hgnc_withdrawn,
        uniprot_ids, entrez_ids, entrez_ids_reverse, mouse_map, rat_map,
        prev_sym_map, ensembl_ids, ensembl_ids_reverse, gene_types,
        dict(hgnc_to_enzymes), dict(enzyme_to_hgncs),
    )


(
    hgnc_names, hgnc_ids, hgnc_withdrawn, uniprot_ids, entrez_ids,
    entrez_ids_reverse, mouse_map, rat_map, prev_sym_map, ensembl_ids,
    ensembl_ids_reverse, gene_type,
    hgnc_to_enzymes, enzyme_to_hgncs,
) = _read_hgnc_maps()


def _read_kinases():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                         'resources', 'kinases.tsv')
    kinase_table = read_unicode_csv(fname, delimiter='\t')
    gene_names = [lin[1] for lin in list(kinase_table)[1:]]
    return gene_names


def _read_phosphatases():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                         'resources', 'phosphatases.tsv')
    p_table = read_unicode_csv(fname, delimiter='\t')
    # First column is phosphatase names
    # Second column is HGNC ids
    p_names = [row[0] for row in p_table]
    return p_names


def _read_tfs():
    fname = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir,
                         'resources', 'transcription_factors.csv')
    tf_table = read_unicode_csv(fname)
    gene_names = [lin[1] for lin in list(tf_table)[1:]]
    return gene_names


kinases, phosphatases, tfs = _read_kinases(), _read_phosphatases(), _read_tfs()
