from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup
import sys

def main():
    # Only install functools32 if we're in Python 2 (it's not available
    # for Python 3)
    install_list = ['pysb>=1.3.0', 'objectpath', 'rdflib==4.2.1',
                    'requests>=2.11', 'lxml', 'ipython', 'future',
                    'networkx==1.11', 'pandas']
    # We implicitly require numpy since pandas depends on numpy
    if sys.version_info[0] == 2:
        install_list.append('functools32')

    extras_require = {
                      # Inputs and outputs
                      'biopax': ['jnius-indra'],
                      'trips_offline': ['pykqml'],
                      'reach_offline': ['jnius-indra'],
                      'eidos_offline': ['pyyaml', 'jnius-indra'],
                      'geneways': ['stemming', 'nltk'],
                      'sofia': ['openpyxl'],
                      'hume': ['rdflib-jsonld'],
                      'ndex': ['ndex2'],
                      'bel': ['pybel'],
                      'sbml': ['python-libsbml'],
                      # Tools and analysis
                      'machine': ['pytz', 'tzlocal', 'tweepy', 'ndex2',
                                  'pyyaml', 'click'],
                      'explanation': ['kappy==4.0.0rc1', 'paths-graph'],
                      # AWS interface and database
                      'aws': ['boto3'],
                      'db': ['sqlalchemy', 'boto3'],
                      # Utilities
                      'graph': ['pygraphviz'],
                      'plot': ['matplotlib'],
                      'isi': ['nltk']
                      }
    extras_require['all'] = list({dep for deps in extras_require.values()
                                  for dep in deps})

    setup(name='indra',
          version='1.7.0',
          description='Integrated Network and Dynamical Reasoning Assembler',
          long_description='INDRA is a framework '
              'for assembling rule-based mathematical models and '
              'mechanistic networks of biochemical systems from natural '
              'language and pathway databases.',
          author='Benjamin Gyori',
          author_email='benjamin_gyori@hms.harvard.edu',
          url='http://github.com/sorgerlab/indra',
          packages=['indra', 'indra.assemblers', 'indra.belief',
                    'indra.benchmarks', 'indra.databases', 'indra.db',
                    'indra.explanation',
                    'indra.literature', 'indra.mechlinker',
                    'indra.preassembler', 'indra.sources',
                    'indra.sources.bel',
                    'indra.sources.biopax', 'indra.sources.cwms',
                    'indra.sources.eidos',
                    'indra.sources.geneways', 'indra.sources.hume',
                    'indra.sources.index_cards',
                    'indra.sources.indra_db_rest',
                    'indra.sources.ndex_cx', 'indra.sources.reach',
                    'indra.sources.sofia',
                    'indra.sources.sparser', 'indra.sources.tees',
                    'indra.sources.trips', 'indra.resources',
                    'indra.resources.famplex', 'indra.tests',
                    'indra.tools', 'indra.tools.reading',
                    'indra.tools.reading.db_reading',
                    'indra.tools.reading.pmid_reading',
                    'indra.tools.reading.starcluster_reading',
                    'indra.tools.reading.util',
                    'indra.tools.machine', 'indra.util'],
          install_requires=install_list,
          tests_require=['jnius-indra', 'jsonschema', 'coverage', 'matplotlib'],
          extras_require=extras_require,
          include_package_data=True,
          keywords=['systems', 'biology', 'model', 'pathway', 'assembler',
                    'nlp', 'mechanism', 'biochemistry', 'network'],
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          entry_points={'console_scripts':
                        ['indra-machine = indra.tools.machine.cli:main']}
        )


if __name__ == '__main__':
    main()
