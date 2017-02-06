Installing INDRA
================

Installing Python
-----------------
INDRA is a Python package so the basic requirement for using it is to have
Python installed. Python is shipped with most Linux distributions and with
OSX. INDRA works with both Python 2 and 3 (tested with
2.7 and 3.5).

On Mac, the preferred way to install Python (over the built-in version) is
using `Homebrew <http://brew.sh/>`_.

.. code-block:: bash

    brew install python

On Windows, we recommend using `Anaconda <https://www.continuum.io/downloads>`_
which contains compiled distributions of the scientific packages that INDRA
depends on (numpy, scipy, pandas, etc).

Installing INDRA
----------------

Installing from Github
``````````````````````
The preferred way to install INDRA is to clone this repository into a local
folder and run setup.py from the terminal as

.. code-block:: bash

    git clone https://github.com/sorgerlab/indra.git
    cd indra
    python setup.py install

Alternatively one can use pip to install INDRA directly from Github as

.. code-block:: bash

    pip install git+https://github.com/sorgerlab/indra.git

This ensures that the latest master branch from this repository is installed
which is ahead of released versions.

Installing with pip
```````````````````
Releases of INDRA are also available via
`pip <https://pip.pypa.io/en/latest/installing/>`_. You can install the latest
released version of INDRA as

.. code-block:: bash

    pip install indra

INDRA dependencies
------------------

INDRA depends on a few standard Python packages (e.g. rdflib, requests) and
also PySB (for more information on PySB, see http://pysb.org). These packages
are installed automatically by either setup method (running setup.py install
or using pip). Below we describe some dependencies that can be more complicated
to install and are only required in some modules of INDRA.

Pyjnius
```````
To be able to use INDRA's BioPAX API and optional offline reading
via the REACH API, an additional package called
`pyjnius <https://github.com/kivy/pyjnius>`_ is needed to allow using Java/Scala
classes from Python. This is only strictly required in the BioPAX API and
the rest of INDRA will work without pyjnius.
Pyjnius needs JRE and JDK 1.8 to be installed. On Mac, you have to install both
`Legacy Java for OS X <http://support.apple.com/kb/DL1572>`_ and
`JDK and JRE from Oracle <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_.
Then set JAVA\_HOME to your JDK home directory, for instance

.. code-block:: bash

    export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Home

Then first install cython (tested with version 0.23.5) followed by jnius-indra

.. code-block:: bash

    pip install cython==0.23.5
    pip install jnius-indra

Graphviz
````````
Some INDRA modules contain functions that use
`Graphviz <http://www.graphviz.org/>`_ to visualize graphs. On most systems, doing

.. code-block:: bash

    pip install pygraphviz

works. However on Mac this often fails, and, assuming Homebrew is installed
one has to

.. code-block:: bash

    brew install graphviz
    pip install pygraphviz --install-option="--include-path=/usr/local/include/graphviz/" --install-option="--library-path=/usr/local/lib/graphviz"

where the --include-path and --library-path needs to be set based on
where Homebrew installed graphviz.

Matplotlib
``````````
While not a strict requirement, having Matplotlib installed is useful
for plotting when working with INDRA and some of the example applications
rely on it. It can be installed as

.. code-block:: bash

    pip install matplotlib

Optional additional dependencies
````````````````````````````````
Some applications built on top of INDRA (for instance The RAS Machine) have
additional dependencies. In such cases a specific `README` or
`requirements.txt` is provided in the folder to guide the set up.
