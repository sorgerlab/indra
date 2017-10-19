"""This is a tool which allows you to create a small sample of the NIH FTP
services on your local machine to test your database managers.
"""
from __future__ import absolute_import, print_function, unicode_literals
from builtins import dict, str

import xml.etree.ElementTree as ET
import tarfile
import zipfile
import gzip
import random
import os
import re
import shutil
from io import BytesIO

if __name__ == '__main__':
    # NOTE: PEP8 will complain about this, however having the args parsed up
    # here prevents a long wait just to fined out you entered a command wrong.
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument(
        dest='n',
        help=('Select the number of examples of each case to produce. '
              'A larger number will generally take longer, but will also '
              'help protect you from being attacked by special cases.'),
        type=int,
        default=2
        )
    parser.add_argument(
        dest='parent_dir',
        help=('Select the name/path of the top level directory containing '
              'your sample.')
        )
    args = parser.parse_args()

from indra.literature import pubmed_client as pub
from indra.db.manage_content import PmcOA, Medline, Manuscripts


def _get_example(case, med_pmid_list, pmc_dicts, man_dicts):
    if case == (1, 0, 0):
        pmid = random.choice(med_pmid_list)
        pmc_pmid_list = [d['PMID'] for d in pmc_dicts if d['PMID'] != '']
        # pmc has only a tiny fraction of all pmid's, so this will be very fast
        while pmid in pmc_pmid_list:
            pmid = random.choice(med_pmid_list)
        ret = (pmid, '', '')
    elif case == (1, 1, 0):
        pmc_pmid_dict_list = [d for d in pmc_dicts if d['PMID'] != '']
        man_pmid_list = [d['PMID'] for d in man_dicts]
        d = random.choice(pmc_pmid_dict_list)
        while d['PMID'] in man_pmid_list:
            d = random.choice(pmc_pmid_dict_list)
        ret = (d['PMID'], d['Accession ID'], '')
    elif case == (0, 1, 0):
        pmcid_list = [d['Accession ID'] for d in pmc_dicts if d['PMID'] == '']
        ret = ('', random.choice(pmcid_list), '')
    elif case == (1, 1, 1):
        pmcid_list = [d['Accession ID'] for d in pmc_dicts]
        d = random.choice(man_dicts)
        while d['PMCID'] not in pmcid_list:
            d = random.choice(man_dicts)
        ret = (d['PMID'], d['PMCID'], d['MID'])
    elif case == (1, 0, 1):
        pmcid_list = [d['Accession ID'] for d in pmc_dicts]
        d = random.choice(man_dicts)
        while d['PMCID'] in pmcid_list:
            d = random.choice(man_dicts)
        ret = (d['PMID'], d['PMCID'], d['MID'])
    else:
        raise Exception("Bad case: %s" % str(case))
    return ret


def build_set(n, parent_dir):
    """Create the nastiest set of content we're willing/able to handle.

    We create a small local representation of the entirety of the nih
    repositories we use, including all the nasty corner cases we can manage.
    This allows for rapid development and testing.

    Parameters
    ----------
    n : int
        The number of copies of each test case to be included. Examples are
        chosen as randomly as possible. Multiple samples could increase the
        reliability of the test.
    parent_dir : str
        The head of the tree that stands in place of the url to the nih ftp
        directory.
    """
    # Create the necessary directories.
    def get_path(sub_path):
        return os.path.join(parent_dir, sub_path)

    if os.path.exists(parent_dir):
        shutil.rmtree(parent_dir)
    os.makedirs(parent_dir)
    os.makedirs(get_path('pub/pmc'))
    os.makedirs(get_path('pubmed/baseline'))
    os.makedirs(get_path('pub/pmc/manuscript'))

    # Get the pmid data from medline (med_pmid_list)
    print("Getting medline lists...")
    med_pmid_list = []
    med = Medline()
    for i in range(1, 7):
        buf = BytesIO()
        med.ftp.ret_file("../MuId-PmId-%d.zip" % i, buf)
        zf = zipfile.ZipFile(buf)
        with zf.open(zf.namelist()[0]) as id_f:
            id_str = id_f.read().decode('utf8')
        med_pmid_list += [l.split('\t')[1] for l in id_str.splitlines()]

    # Get the data from pmc oa (pmc_dicts)
    print("Getting pmc oa lists....")
    pmc = PmcOA()
    pmc_dicts = pmc.ftp.get_csv_as_dict('oa_file_list.csv')

    # Get the data for the manuscripts (man_dicts)
    print("Getting manuscript lists...")
    man = Manuscripts()
    man_dicts = man.ftp.get_csv_as_dict('filelist.csv')

    # Get pmid, pmcid, mid tuples for the examples that we will use.
    print("Generating example sets...")
    examples = []
    for case in [(1,0,0), (1,1,0), (0,1,0), (1,1,1), (1,0,1)]:
        for _ in range(n):
            example = _get_example(case, med_pmid_list, pmc_dicts, man_dicts)
            examples.append(example)
    double_doi_info = med.get_article_info('medline17n0343.xml.gz')
    pmids_w_double_doi = [
        k for k, v in double_doi_info.items()
        if v['doi'] is not None and len(v['doi']) > 100
        ]
    examples.append((random.choice(pmids_w_double_doi), '', '',))

    # Create the test medline file.
    print("Creating medline test file...")
    pmid_list = [pmid for pmid, _, _ in examples if pmid != '']
    tree = None
    for pmid in pmid_list:
        params = {'db': 'pubmed', 'retmode': 'xml', 'id': pmid}
        if tree is None:
            tree = pub.send_request(pub.pubmed_fetch, params)
        else:
            child = pub.send_request(pub.pubmed_fetch, params).getchildren()[0]
            tree.append(child)
    if tree is not None:
        f_bts = b''
        f_bts += b"<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
        f_bts += ET.tostring(tree)
        f_path = get_path('pubmed/baseline/medline17nTEST.xml.gz')
        with open(f_path, 'wb') as gzf:
            gzf.write(gzip.compress(f_bts))

    # Create the test pmc oa article directory.
    print("Getting pmc oa xmls...")
    art_dirname = get_path('pub/pmc/articles.TEST.xml')
    if os.path.exists(art_dirname):
        shutil.rmtree(art_dirname)
    os.mkdir(art_dirname)
    pmcid_list = [pmcid for _, pmcid, _ in examples if pmcid != '']
    ex_pmc_dicts = [d for d in pmc_dicts if d['Accession ID'] in pmcid_list]
    for d in ex_pmc_dicts:
        fname = pmc.ftp.download_file(d['File'])
        with tarfile.open(fname, 'r:gz') as tar:
            mems = tar.getmembers()
            mem = [mem for mem in mems if mem.name.endswith('.nxml')][0]
            f_str = tar.extractfile(mem).read()
        fname = d['Accession ID'] + '.nxml'
        re_ret = re.findall('<journal-title>(.*?)</journal-title>',
                            f_str.decode('utf8'))
        if len(re_ret):
            sub_dir = os.path.join(
                art_dirname,
                re_ret[0].replace(' ', '_').replace('&', '')
                )
        else:
            sub_dir = os.path.join(art_dirname, 'unknown')
        if not os.path.exists(sub_dir):
            os.mkdir(sub_dir)
        path = os.path.join(sub_dir, fname)
        with open(path, 'wb') as f:
            f.write(f_str)
    with tarfile.open(art_dirname + '.tar.gz', 'w:gz') as tar:
        for dirname in os.listdir(art_dirname):
            tar.add(
                os.path.join(art_dirname, dirname),
                arcname=dirname
                )
    shutil.rmtree(art_dirname)

    # Create deleted pmids file (just make an empty file,for now.
    # TODO: Add test case to touch this.
    with open(get_path('pubmed/deleted.pmids.gz'), 'wb') as gzf:
        gzf.write(gzip.compress(b''))

    # Create the test manuscripts file.
    print('Adding manuscript directories...')
    dirfmt = get_path('pub/pmc/manuscript/%s')
    dirnames = [dirfmt % ('PMC00%dXXXXXX.xml' % i) for i in range(2, 6)]
    for dirname in dirnames:
        if os.path.exists(dirname):
            shutil.rmtree(dirname)
        os.mkdir(dirname)
    ex_man_dicts = [d for d in man_dicts if d['PMCID'] in pmcid_list]
    for d in ex_man_dicts:
        d['Tarfile'] = man.get_tarname_from_filename(d['File'])
    tar_members = dict.fromkeys(set(
        [d['Tarfile'] for d in ex_man_dicts]
        ))
    for tarname in tar_members.keys():
        if not os.path.exists(tarname):
            print("\tDownloading %s..." % tarname)
            man.ftp.download_file(tarname)
    for d in ex_man_dicts:
        parent_dir = os.path.join(
            dirfmt % tarname.replace('.tar.gz', ''),
            os.path.dirname(d['File'])
            )
        test_fname = os.path.join(
            dirfmt % tarname.replace('.tar.gz', ''),
            d['File']
            )
        with tarfile.open(d['Tarfile'], 'r:gz') as tar:
            print('\tExtracting %s from %s...' % (d['File'], d['Tarfile']))
            tar.extract(d['File'])
        if not os.path.exists(parent_dir):
            os.makedirs(parent_dir)
        os.rename(d['File'], test_fname)
    for dirname in dirnames:
        with tarfile.open(dirname + '.tar.gz', 'w:gz') as tar:
            for sub_dirname in os.listdir(dirname):
                tar.add(
                    os.path.join(dirname, sub_dirname),
                    arcname=sub_dirname
                    )
        shutil.rmtree(dirname)

    return examples


if __name__ == '__main__':
    build_set(args.n, args.parent_dir)
