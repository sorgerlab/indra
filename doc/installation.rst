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

INDRA depends on a few standard Python packages (e.g. rdflib, requests,
objectpath).
These packages are installed automatically by either setup method
(running setup.py install or using pip). 

Below we provide a detailed description of some extra dependencies that may
require special steps to install.

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
Some dependencies of INDRA are only needed by certain submodules or are only
used in specialized use cases. These are not installed by default but are
listed as "extra" requirements, and can be installed separately using pip.
An extra dependency list (e.g. one called extra_list) can be
installed as

.. code-block:: bash

    pip install indra[extra_list]

You can also install all extra dependencies by doing

.. code-block:: bash

   pip install indra --install-option="complete"

or 

.. code-block:: bash

   pip install indra[all]

In all of the above, you may replace `indra` with `.` (if you're in a local
copy of the `indra` folder or with the Github URL of the INDRA repo, depending
on your installation method.
See also the corresponding
`pip documentation <https://packaging.python.org/tutorials/installing-packages/#installing-setuptools-extras>`_
for more information.

The table below provides the name and the description of each "extra" list
of dependencies.

+-----------------+------------------------------------------------------+
|Extra list name  |Purpose                                               |
+=================+======================================================+
|biopax           |BioPAX input processing and Pathway Commons queries   |
+-----------------+------------------------------------------------------+
|bel              |BEL input processing and output assembly              |
+-----------------+------------------------------------------------------+
|trips_offline    |Offline reading with local instance of TRIPS system   |
+-----------------+------------------------------------------------------+
|reach_offline    |Offline reading with local instance of REACH system   |
+-----------------+------------------------------------------------------+
|eidos_offline    |Offline reading with local instance of Eidos system   |
+-----------------+------------------------------------------------------+
|geneways         |Genewayas reader input processing                     |
+-----------------+------------------------------------------------------+
|sofia            |SOFIA reader input processing                         |
+-----------------+------------------------------------------------------+
|bbn              |BBN reader input processing                           |
+-----------------+------------------------------------------------------+
|ndex             |NDEx client for network input processing and CX upload|
+-----------------+------------------------------------------------------+
|machine          |Running a local instance of a "RAS machine"           |
+-----------------+------------------------------------------------------+
|explanation      |Finding explanatory paths in rule-based models        |
+-----------------+------------------------------------------------------+
|aws              |Accessing AWS compute and storage resources           |
+-----------------+------------------------------------------------------+
|db               |Setting up and using a DB instance to store Statements|
+-----------------+------------------------------------------------------+
|graph            |Assembling into a visualizing Graphviz graphs         |
+-----------------+------------------------------------------------------+
|plot             |Create and display plots                              |
+-----------------+------------------------------------------------------+
Configuring INDRA
-----------------
Various aspects of INDRA, including API keys, dependency locations, and
Java memory limits, are parameterized by a configuration file that lives in
~/.config/indra/config.ini. The default
configuration file is provided in indra/resources/default_config.ini, and
is copied to ~/.config/indra/config.ini when INDRA starts if no configuration
already exists. Every value in the configuration can also be set as an
environment variable: for a given configuration key, INDRA will first check
for an environment variable with that name and if not present, will use
the value in the configuration file. In other words, an environment variable,
when set, takes precedence over the value set in the config file.

Configuration values include:
- REACHPATH: The location of the JAR file containing a local instance of the
REACH reading system

- EIDOSPATH: The location of the JAR file containing a local instance of the
Eidos reading system

- SPARSERPATH: The location of a local instance of the Sparser
reading system (path to a folder)

- DRUMPATH: The location of a local installation of the DRUM reading system
(path to a folder)

- NDEX_USERNAME, NDEX_PASSWORD: Credentials for accessing the NDEx web service

- ELSEVIER_API_KEY, ELSEVIER_INST_KEY: Elsevier web service API keys

- BIOGRID_API_KEY: API key for BioGRID web service (see 
  http://wiki.thebiogrid.org/doku.php/biogridrest)

- INDRA_DEFAULT_JAVA_MEM_LIMIT: Maximum memory limit for Java virtual machines
  launched by INDRA
