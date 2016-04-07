from ez_setup import use_setuptools
use_setuptools()
from setuptools import setup

def main():
    setup(name='indra',
          version='1.0.1',
          description='Integrated Network and Dynamical Reasoning Assembler',
          long_description='INDRA is a framework '
              'for assembling rule-based mathematical models of biochemical '
              'systems from natural language and pathway databases.',
          author='Benjamin Gyori',
          author_email='benjamin_gyori@hms.harvard.edu',
          url='http://github.com/sorgerlab/indra',
          download_url='https://github.com/sorgerlab/indra/tarball/1.0.1',
          packages=['indra', 'indra.bel', 'indra.biopax', 
                    'indra.reach', 'indra.trips', 'indra.databases'],
          install_requires=['pysb', 'objectpath', 'rdflib', 'jnius-indra', 
                            'functools32', 'requests', 'lxml'],
          include_package_data=True,
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
if __name__ == '__main__':
    main()
