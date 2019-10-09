import os
from setuptools import setup


readme_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                           'README.md')
with open(readme_path, 'r', encoding='utf-8') as fh:
    long_description = fh.read()


def main():
    install_list = ['pysb>=1.3.0', 'objectpath', 'rdflib',
                    'requests>=2.11', 'lxml', 'ipython', 'future',
                    'networkx>=2', 'pandas', 'ndex2==2.0.1', 'jinja2',
                    'protmapper>=0.0.14']

    extras_require = {
                      # Inputs and outputs
                      'biopax': ['cython', 'pyjnius==1.1.4'],
                      'trips_offline': ['pykqml'],
                      'reach_offline': ['cython', 'pyjnius==1.1.4'],
                      'eidos_offline': ['pyyaml', 'cython', 'pyjnius==1.1.4'],
                      'geneways': ['stemming', 'nltk'],
                      'sofia': ['openpyxl'],
                      'bel': ['pybel'],
                      'sbml': ['python-libsbml'],
                      'obo': ['obonet'],
                      # Tools and analysis
                      'machine': ['pytz', 'tzlocal', 'tweepy', 'pyyaml',
                                  'click'],
                      'explanation': ['kappy==4.0.0rc1', 'paths-graph'],
                      'adeft': ['adeft'],
                      # AWS interface and database
                      'aws': ['boto3', 'reportlab'],
                      # Utilities
                      'graph': ['pygraphviz'],
                      'plot': ['matplotlib'],
                      'isi': ['nltk'],
                      'api': ['flask']
                      }
    extras_require['all'] = list({dep for deps in extras_require.values()
                                  for dep in deps})

    setup(name='indra',
          version='1.14.1',
          description='Integrated Network and Dynamical Reasoning Assembler',
          long_description=long_description,
          long_description_content_type='text/markdown',
          author='Benjamin Gyori',
          author_email='benjamin_gyori@hms.harvard.edu',
          url='http://github.com/sorgerlab/indra',
          packages=['indra', 'indra.assemblers', 'indra.assemblers.cag',
                    'indra.assemblers.cx', 'indra.assemblers.cyjs',
                    'indra.assemblers.english', 'indra.assemblers.figaro',
                    'indra.assemblers.graph', 'indra.assemblers.html',
                    'indra.assemblers.index_card',
                    'indra.assemblers.indranet',
                    'indra.assemblers.kami', 'indra.assemblers.pybel',
                    'indra.assemblers.pysb', 'indra.assemblers.sbgn',
                    'indra.assemblers.sif', 'indra.assemblers.tsv',
                    'indra.belief',
                    'indra.benchmarks', 'indra.databases',
                    'indra.explanation', 'indra.explanation.model_checker',
                    'indra.literature', 'indra.mechlinker',
                    'indra.preassembler',
                    'indra.preassembler.grounding_mapper', 'indra.sources',
                    'indra.sources.bel',
                    'indra.sources.biopax', 'indra.sources.cwms',
                    'indra.sources.eidos',
                    'indra.sources.geneways', 'indra.sources.hprd',
                    'indra.sources.hume', 'indra.sources.index_cards',
                    'indra.sources.indra_db_rest', 'indra.sources.lincs_drug',
                    'indra.sources.ndex_cx', 'indra.sources.reach',
                    'indra.sources.rlimsp', 'indra.sources.sofia',
                    'indra.sources.sparser', 'indra.sources.tas',
                    'indra.sources.tees',
                    'indra.sources.trips', 'indra.sources.trrust',
                    'indra.resources',
                    'indra.resources.famplex', 'indra.statements',
                    'indra.tests', 'indra.tests.test_obo_clients',
                    'indra.tools', 'indra.tools.reading',
                    'indra.tools.reading.pmid_reading',
                    'indra.tools.reading.starcluster_reading',
                    'indra.tools.reading.util',
                    'indra.tools.machine', 'indra.util'],
          install_requires=install_list,
          extras_require=extras_require,
          include_package_data=True,
          keywords=['systems', 'biology', 'model', 'pathway', 'assembler',
                    'nlp', 'mechanism', 'biochemistry', 'network'],
          classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
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
