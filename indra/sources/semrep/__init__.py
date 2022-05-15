"""This module implements an API and processor for the SemRep system
available at https://github.com/lhncbc/SemRep (see also
https://doi.org/10.1016/j.jbi.2003.11.003). Setting up the SemRep requires
downloading several large compressed software releases and data files, and then
running installation scripts. A Dockerfile is provided in this module
which allows automating the build into a ~100GB image, or can be used as
a template for a local installation.

Currently, the XML output format of SemRep is supported which can be
produced as follows:

[SemRep base folder]/bin/semrep.v1.8 -X -L 2018 -Z 2018AA input.txt output.xml

Predicates produced by SemRep are documented in the Appendix of
https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-486#Sec26.
"""

from .api import *
