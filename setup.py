from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup

def main():
    setup(name='indra',
          version='1.2.0',
          description='Integrated Network and Dynamical Reasoning Assembler',
          long_description='INDRA is a framework '
              'for assembling rule-based mathematical models and '
              'mechanistic networks of biochemical systems from natural '
              'language and pathway databases.',
          author='Benjamin Gyori',
          author_email='benjamin_gyori@hms.harvard.edu',
          url='http://github.com/sorgerlab/indra',
          packages=['indra', 'indra.assemblers', 'indra.bel', 'indra.biopax',
                    'indra.reach', 'indra.trips', 'indra.databases',
                    'indra.preassembler', 'indra.mechlinker',
                    'indra.tools', 'indra.tests', 'indra.resources',
                    'indra.literature'],
          install_requires=['pysb', 'objectpath', 'rdflib',
                            'functools32', 'requests', 'lxml', 'jsonpickle'],
          tests_require=['jnius-indra', 'jsonschema', 'coverage'],
          include_package_data=True,
          keywords=['systems', 'biology', 'model', 'pathway', 'assembler', 'nlp',
                    'mechanism', 'biochemistry', 'network'],
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
if __name__ == '__main__':
    main()
