#!/usr/bin/env python

from distutils.core import setup
import distutils.cmd
import sys, os, subprocess, traceback, re

def main():

    setup(name='indra',
          version=get_version(),
          description='Integrated Network and Dynamical Reasoning Assembler',
          long_description='INDRA is a framework ' + \
              'for assembling rule-based mathematical models of biochemical ' + \
              'systems from natural language and pathway databases.'
          author='Benjamin Gyori',
          author_email='benjamin_gyori@hms.harvard.edu',
          url='http://github.com/indra',
          packages=['indra', 'indra.bel', 'indra.biopax', 
                    'indra.reach', 'indra.trips', 'indra.databases'],
          requires=['pysb', 'objectpath', 'BeautifulSoup', 'rdflib'],
          keywords=['systems', 'biology', 'model', 'pathway', 'assembler', 'nlp', 
                    'mechanism', 'biochemistry'],
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          )

def get_version():
    rv_filename = 'RELEASE-VERSION'
    base_path = os.path.dirname(__file__)
    rv_path = os.path.join(base_path, rv_filename)
    try:
        version = get_git_version()
        version_file = open(rv_path, 'w')
        version_file.write(version)
        version_file.close() # ensure sdist build process sees new contents
    except Exception:
        try:
            version_file = open(rv_path, 'r')
            version = version_file.read()
        except IOError as e:
            if e.errno == 2 and e.filename == rv_path:
                sys.stderr.write(
"""This does not appear to be a git repository, and the file %s is not
present. In order to build or install INDRA, please either download a
distribution from http://pypi.python.org/pypi/indra or clone the git repository
at https://github.com/indra/indra.git\n""" % rv_filename)
                sys.exit(1)
    return version

class GitError(Exception):
    pass

def get_git_version():
    """Get a nice version number from git-describe"""

    # ensure that we are working in an indra git repo
    setup_path = os.path.abspath(os.path.dirname(__file__))
    if not os.path.exists(os.path.join(setup_path, '.git')):
        raise Exception("setup.py is not in the root of a git repository; "
                        "aborting")
    os.chdir(setup_path)

    # run git describe
    gitcmd = ['git', 'describe', '--always']
    try:
        gitproc = subprocess.Popen(gitcmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        retcode = gitproc.wait()
        if retcode:
            raise GitError(gitproc.stderr.read())
        version = gitproc.stdout.next().strip()
        version = re.sub(r'^v', '', version)
        return version
    except (OSError, GitError) as e:
        raise Exception("Error running 'git describe' to determine version:\n\n" +
                        "command\n=====\n" + " ".join(gitcmd) + "\n\n" +
                        "error\n====\n" + str(e) + "\n")

if __name__ == '__main__':
    main()
