from __future__ import absolute_import, print_function, unicode_literals, \
                       division
from builtins import dict, str

if __name__ == '__main__':
    import sys
    from glob import glob
    from os.path import join, isfile, isdir
    import subprocess

    target_dir = sys.argv[1]
    dir_files = glob(join(target_dir, '*'))
    for dir in dir_files:
        if not isdir(dir):
            continue
        # Check for the presence of the content_types dir and the output
        # directory
        content_types_path = join(dir, 'content_types.pkl')
        output_dir = join(dir, 'output')
        if not isfile(content_types_path):
            continue
        if not isdir(output_dir):
            continue
        cmd_list = ['qsub', '-b', 'y', '-V', '-cwd', 'python', '-m',
                    'indra.tools.reading.pmid_reading.read_pmids',
                    'upload_json', output_dir, content_types_path]
        print(' '.join(cmd_list))
        subprocess.call(cmd_list)
