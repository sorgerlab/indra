"""
Create a small sample of the NIH FTP services on your local machine to test
your database managers with representative cases involving different
combinations of IDs (PMID, PMCID, Author's Manuscript IDs) and also articles
that yield/don't yield Statements when read.

Collections of articles covering these cases are downloaded from the NIH
FTP services and repackaged locally in file structure mirroring that of the
FTP services themselves.
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
from indra.db.content_manager import PmcOA, Pubmed, Manuscripts


def _get_example(case, med_pmid_list, pmc_dicts, man_dicts):
    # PMID, no PMCID, no MS ID
    if case == (1, 0, 0):
        pmid = random.choice(med_pmid_list)
        pmc_pmid_list = [d['PMID'] for d in pmc_dicts if d['PMID'] != '']
        # pmc has only a tiny fraction of all pmid's, so this will be very fast
        while pmid in pmc_pmid_list:
            pmid = random.choice(med_pmid_list)
        ret = (pmid, '', '')
    # PMID, PMCID, no MS ID
    elif case == (1, 1, 0):
        pmc_pmid_dict_list = [d for d in pmc_dicts if d['PMID'] != '']
        man_pmid_list = [d['PMID'] for d in man_dicts]
        d = random.choice(pmc_pmid_dict_list)
        while d['PMID'] in man_pmid_list:
            d = random.choice(pmc_pmid_dict_list)
        ret = (d['PMID'], d['Accession ID'], '')
    # no PMID, PMCID, no MS ID
    elif case == (0, 1, 0):
        pmcid_list = [d['Accession ID'] for d in pmc_dicts if d['PMID'] == '']
        ret = ('', random.choice(pmcid_list), '')
    # PMID, PMCID, MS ID
    elif case == (1, 1, 1):
        pmcid_list = [d['Accession ID'] for d in pmc_dicts]
        d = random.choice(man_dicts)
        while d['PMCID'] not in pmcid_list:
            d = random.choice(man_dicts)
        ret = (d['PMID'], d['PMCID'], d['MID'])
    # PMID, no PMCID, MS ID
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

    We create a small local representation of the entirety of the NLM
    repositories we use, including all the nasty corner cases we can manage.
    This allows for rapid development and testing.

    Parameters
    ----------
    n : int
        The number of instances (distinct articles) of each test case to be
        included. Examples are chosen as randomly as possible. Multiple samples
        generally increase the reliability of the test.
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
    med = Pubmed()
    for i in range(1, 7):
        buf = BytesIO()
        med.ftp.ret_file("MuId-PmId-%d.zip" % i, buf)
        zf = zipfile.ZipFile(buf)
        with zf.open(zf.namelist()[0]) as id_f:
            id_str = id_f.read().decode('utf8')
        med_pmid_list += [l.split('\t')[1] for l in id_str.splitlines()]

    statementful_pmids = [
        '20949557', '23898069', '19801969', '21042724', '14675752', '25897078',
        '25486481', '12890751', '11251186', '20622853', '25616414', '21878640',
        '23295773', '19747910', '25778309', '25939761', '11871856', '16580132',
        '24730770', '23921085', '22018470', '19405127', '21464949', '18321309',
        '7907095', '12048232', '23751074', '18711136', '13679391', '22193543',
        '26645886', '27086966', '14570914', '20538416', '9417079', '23200589',
        '15146469', '18084123', '19265534', '19449221', '27381626', '14976202',
        '22445724', '20040392', '26039245', '17881156', '15902258', '1745350',
        '18276758', '22764095', '20652941', '25834816', '23068100', '16407218',
        '18830263', '24265318', '19752028', '8589722', '22671588', '14745431',
        '25042645', '19403642', '14707024', '23536437', '21167476', '22801439',
        '25726184', '19723643', '17409824', '28679432', '26908611', '20164468',
        '15189946', '12086229', '21900397', '12324477', '15545228', '23376846',
        '21719749', '20608972', '23583295', '23236067', '9705962', '20068183',
        '19437340', '14534726', '25731731', '15337767', '28067895', '25092803',
        '19261749', '22272295', '27121230', '23302038', '17410335', '17399955',
        '16254247', '21685363', '26598524', '25645929', '1386335', '20606534',
        '22492281', '22158902', '22022427', '24775712', '21298412', '24753544',
        '12553064', '19681600', '17912454', '17597401', '20672986', '21362231',
        '17999917', '21470928', '27334922', '16159962', '21079653', '15125833',
        '27617579', '19048115', '18687691', '27797218', '26413934', '16684954',
        '20501406', '27515963', '22784503', '25941399', '12473120', '17891137',
        '16733295', '23826126', '21427728', '8900182', '26234677', '24648515',
        '25786138', '12958678', '16998791', '19061835', '11283269', '18258923',
        '11839584', '20132317', '19158374', '23245941', '23352210', '15465819',
        '15386433', '22575647', '15966238', '23633483', '25131797', '17102080',
        '19956840', '18506362', '17961162', '1607067', '24770328', '19825990',
        '22365656', '19720761', '24435975', '26882953', '17292826', '25119113',
        '26044620', '20717925', '15316008', '16619041', '19893488', '26999786',
        '26103054', '17331464', '20022966', '24189165', '19059939', '25474223',
        '20507346', '20976540', '2810532', '15685397', '27562587', '18538673',
        '15712349', '15448517', '27467210', '7584044', '21330319', '18381962',
        '24789704', '19058873', '10523313'
        ]

    # Get the data from pmc oa (pmc_dicts)
    print("Getting pmc oa lists....")
    pmc = PmcOA()
    pmc_dicts = pmc.ftp.get_csv_as_dict('oa_file_list.csv', header=0)

    # Get the data for the manuscripts (man_dicts)
    print("Getting manuscript lists...")
    man = Manuscripts()
    man_dicts = man.ftp.get_csv_as_dict('filelist.csv', header=0)

    # Get pmid, pmcid, mid tuples for the examples that we will use.
    print("Generating example sets...")
    examples = []
    for case in [(1,0,0), (1,1,0), (0,1,0), (1,1,1), (1,0,1)]:
        for _ in range(n):
            example = _get_example(case, med_pmid_list, pmc_dicts, man_dicts)
            examples.append(example)
    # Add a few pmids that probably include some statements.
    for pmid in random.sample(statementful_pmids, n):
        examples.append((pmid, '', ''))
    # Add a special article to check article info.
    double_doi_info = med.get_article_info('baseline/pubmed18n0343.xml.gz')
    pmids_w_double_doi = [
        k for k, v in double_doi_info.items()
        if v['doi'] is not None and len(v['doi']) > 100
        ]
    assert len(pmids_w_double_doi), "No double dois found."
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
        f_path = get_path('pubmed/baseline/pubmed18nTEST.xml.gz')
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
