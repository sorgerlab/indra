import os
from setuptools import setup


readme_path = os.path.join(os.path.abspath(os.path.dirname(__file__)),
                           'README.md')
with open(readme_path, 'r', encoding='utf-8') as fh:
    long_description = fh.read()


def main():
    install_list = ['pysb>=1.3.0,<=1.9.1', 'objectpath', 'rdflib==4.2.2',
                    'requests>=2.11', 'lxml', 'ipython', 'future',
                    'networkx>=2', 'pandas', 'ndex2==2.0.1', 'jinja2',
                    'protmapper>=0.0.21', 'obonet', 'sympy==1.3',
                    'tqdm', 'pybiopax>=0.0.5']

    extras_require = {
                      # Inputs and outputs
                      'trips_offline': ['pykqml'],
                      'reach_offline': ['cython', 'pyjnius==1.1.4'],
                      'eidos_offline': ['pyyaml>=5.1.0', 'cython',
                                        'pyjnius==1.1.4'],
                      'hypothesis': ['gilda'],
                      'geneways': ['stemming', 'nltk<3.6'],
                      'sofia': ['openpyxl'],
                      'bel': ['pybel>=0.15.0,<0.16.0'],
                      'sbml': ['python-libsbml'],
                      # Tools and analysis
                      'machine': ['pytz', 'tzlocal', 'tweepy', 'pyyaml>=5.1.0',
                                  'click'],
                      'explanation': ['kappy==4.1.2', 'paths-graph'],
                      'grounding': ['adeft', 'gilda'],
                      # AWS interface and database
                      'aws': ['boto3', 'reportlab'],
                      # Utilities
                      'graph': ['pygraphviz'],
                      'plot': ['matplotlib'],
                      'isi': ['nltk<3.6', 'unidecode'],
                      'api': ['flask', 'flask_restx', 'flask_cors',
                              'docstring-parser']
                      }
    extras_require['all'] = list({dep for deps in extras_require.values()
                                  for dep in deps})

    setup(name='indra',
          version='1.19.0',
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
                    'indra.explanation.pathfinding',
                    'indra.literature', 'indra.mechlinker',
                    'indra.ontology', 'indra.ontology.bio',
                    'indra.ontology.world', 'indra.ontology.virtual',
                    'indra.ontology.app', 'indra.pipeline',
                    'indra.preassembler',
                    'indra.preassembler.grounding_mapper', 'indra.sources',
                    'indra.sources.bel', 'indra.sources.biofactoid',
                    'indra.sources.biopax', 'indra.sources.ctd',
                    'indra.sources.crog',
                    'indra.sources.cwms', 'indra.sources.dgi',
                    'indra.sources.drugbank', 'indra.sources.eidos',
                    'indra.sources.geneways', 'indra.sources.hprd',
                    'indra.sources.hume', 'indra.sources.hypothesis',
                    'indra.sources.index_cards',
                    'indra.sources.indra_db_rest', 'indra.sources.isi',
                    'indra.sources.minerva', 'indra.sources.ndex_cx',
                    'indra.sources.reach', 'indra.sources.omnipath',
                    'indra.sources.phosphoelm',
                    'indra.sources.rlimsp', 'indra.sources.signor',
                    'indra.sources.sofia',
                    'indra.sources.sparser', 'indra.sources.tas',
                    'indra.sources.tees',
                    'indra.sources.trips', 'indra.sources.trrust',
                    'indra.sources.virhostnet',
                    'indra.resources', 'indra.statements',
                    'indra.tests', 'indra.tests.test_obo_clients',
                    'indra.tools',
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
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            ],
          entry_points={'console_scripts':
                        ['indra-machine = indra.tools.machine.cli:main']}
        )


if __name__ == '__main__':
    main()
