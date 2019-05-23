"""This tool provides a uniform method for createing a robust indra version
string, both from within python and from commandline. If possible, the version
will include the git commit hash. Otherwise, the version will be marked with
'UNHASHED'.
"""

import re
from os import devnull, chdir, curdir
from os.path import dirname, abspath
from subprocess import check_output, CalledProcessError

from indra import __version__


INDRA_GITHASH = None


def get_git_info():
    """Get a dict with useful git info."""
    start_dir = abspath(curdir)
    try:
        chdir(dirname(abspath(__file__)))
        re_patt_str = (r'commit\s+(?P<commit_hash>\w+).*?Author:\s+'
                       r'(?P<author_name>.*?)\s+<(?P<author_email>.*?)>\s+Date:\s+'
                       r'(?P<date>.*?)\n\s+(?P<commit_msg>.*?)(?:\ndiff.*?)?$')
        show_out = check_output(['git', 'show']).decode('ascii')
        revp_out = check_output(['git', 'rev-parse', '--abbrev-ref', 'HEAD'])
        revp_out = revp_out.decode('ascii').strip()
        m = re.search(re_patt_str, show_out, re.DOTALL)
        assert m is not None, \
            "Regex pattern:\n\n\"%s\"\n\n failed to match string:\n\n\"%s\"" \
            % (re_patt_str, show_out)
        ret_dict = m.groupdict()
        ret_dict['branch_name'] = revp_out
    finally:
        chdir(start_dir)
    return ret_dict


def get_version(with_git_hash=True, refresh_hash=False):
    """Get an indra version string, including a git hash."""
    version = __version__
    if with_git_hash:
        global INDRA_GITHASH
        if INDRA_GITHASH is None or refresh_hash:
            with open(devnull, 'w') as nul:
                try:
                    ret = check_output(['git', 'rev-parse', 'HEAD'],
                                       cwd=dirname(__file__), stderr=nul)
                except CalledProcessError:
                    ret = b'UNHASHED'
            INDRA_GITHASH = ret.strip()
            if isinstance(INDRA_GITHASH, bytes):
                INDRA_GITHASH = INDRA_GITHASH.decode('utf-8')
        version = '%s-%s' % (version, INDRA_GITHASH)
    return version


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('-n', '--no_hash', action='store_true',
                        help='Choose to not include the git hash.')
    args = parser.parse_args()
    version = get_version(not args.no_hash)
    print(version)
