from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup
import sys

def main():
    # Only install functools32 if we're in Python 2 (it's not available
    # for Python 3)
    install_list = ['pysb', 'objectpath', 'rdflib', 'requests', 'lxml',
                    'jsonpickle', 'ipython', 'future']
    if sys.version_info[0] == 2:
        install_list.append('functools32')

    setup(name='indra',
          version='1.3.0',
          description='Integrated Network and Dynamical Reasoning Assembler',
          long_description='INDRA is a framework '
              'for assembling rule-based mathematical models and '
              'mechanistic networks of biochemical systems from natural '
              'language and pathway databases.',
          author='Benjamin Gyori',
          author_email='benjamin_gyori@hms.harvard.edu',
          url='http://github.com/sorgerlab/indra',
          packages=['indra', 'indra.assemblers', 'indra.bel', 'indra.belief',
                    'indra.benchmarks',
                    'indra.biopax', 'indra.databases', 'indra.index_cards',
                    'indra.literature', 'indra.mechlinker',
                    'indra.preassembler', 'indra.reach', 'indra.resources',
                    'indra.tests', 'indra.tools', 'indra.tools.reading',
                    'indra.trips', 'indra.util'],
          install_requires=install_list,
          tests_require=['jnius-indra', 'jsonschema', 'coverage', 'matplotlib'],
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
          )
if __name__ == '__main__':
    main()
