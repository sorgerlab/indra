Installation
============

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

Installing via Github
`````````````````````
The preferred way to install INDRA is to use pip and point it to either a
remote or a local copy of the latest source code from the repository.
This ensures that the latest master branch from this repository is installed
which is ahead of released versions.

To install directly from Github, do:

.. code-block:: bash

    pip install git+https://github.com/sorgerlab/indra.git

Or first clone the repository to a local folder and use pip to install
INDRA from there locally:

.. code-block:: bash

    git clone https://github.com/sorgerlab/indra.git
    cd indra
    pip install .

Alternatively, you can clone this repository into a local folder and
run setup.py from the terminal as

.. code-block:: bash

    git clone https://github.com/sorgerlab/indra.git
    cd indra
    python setup.py install

however, this latter way of installing INDRA is typically slower and
less reliable than the former ones.

Cloning the source code from Github
```````````````````````````````````
You may want to simply clone the source code without installing INDRA
as a system-wide package.

.. code-block:: bash

    git clone https://github.com/sorgerlab/indra.git

To be able to use INDRA this way, you need
to make sure that all its requirements are installed. To be able to
`import indra`, you also need the folder to be visible on your
`PYTHONPATH <https://docs.python.org/2/using/cmdline.html#envvar-PYTHONPATH>`_
environmental variable.

Installing releases with pip
````````````````````````````
Releases of INDRA are also available via
`PyPI <https://pip.pypa.io/en/latest/installing/>`_. You can install the latest
released version of INDRA as

.. code-block:: bash

    pip install indra

INDRA dependencies
------------------

INDRA depends on a few standard Python packages (e.g. rdflib, requests, pysb).
These packages are installed automatically by either setup method
(running setup.py install or using pip). Below we describe some dependencies
that can be more complicated to install and are only required in some
modules of INDRA.

PySB and BioNetGen
``````````````````
INDRA builds on the `PySB <http://pysb.org>`_ framework to assemble rule-based
models of biochemical systems. The `pysb` python package is installed by
the standard install procedure. However, to be able to generate mathematical
model equations and to export to formats such as SBML, the
`BioNetGen <http://bionetgen.org/index.php/BioNetGen_Distributions>`_
framework also needs to be installed in a way that is visible to PySB.
Detailed instructions are given in the
`PySB documentation <http://docs.pysb.org/en/latest/installation.html#option-1-install-pysb-natively-on-your-computer>`_.

.. _pyjniussetup:

Pyjnius
```````
To be able to use INDRA's BioPAX API and optional offline reading
via the REACH and Eidos APIs, an additional package called
`pyjnius <https://github.com/kivy/pyjnius>`_ is needed to allow using Java/Scala
classes from Python. This is only strictly required in these input sources and
the rest of INDRA will work without pyjnius.

1. Install `JRE and JDK 8 from Oracle <http://www.oracle.com/technetwork/java/javase/downloads/index.html>`_. Pyjnius is currently incompatible with Java 9, so
make sure to get Java 8.

2. On Mac, `install Legacy Java for OSX <http://support.apple.com/kb/DL1572>`_.
If you have trouble installing it, you can try the following as an alternative.
Edit

.. code-block:: bash

    /Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Info.plist

(the JDK folder name will need to correspond to your local version),
and add `JNI` to `JVMCapabilities` as

.. code-block:: xml

    ...
    <dict>
        <key>JVMCapabilities</key>
        <array>
            <string>CommandLine</string>
            <string>JNI</string>
        </array>
    ...

3. Set JAVA\_HOME to your JDK home directory, for instance

.. code-block:: bash

    export JAVA_HOME=/Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Home

4. Then first install cython (tested with version 0.23.5) followed by jnius-indra. These need to be
   broken up into two sequential calls to pip install.

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
additional dependencies that are encoded as
`extras in setup.py <https://setuptools.readthedocs.io/en/latest/setuptools.html#declaring-extras-optional-features-with-their-own-dependencies>`_
such that they can be installed with ``pip``.

.. code-block:: bash

    pip install indra[machine]

In other cases a specific `README` or `requirements.txt` is provided in the folder to guide the set up.
